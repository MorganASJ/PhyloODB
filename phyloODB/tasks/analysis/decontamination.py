import tempfile
import time
import json
from abc import ABC, abstractmethod
from datetime import datetime
import os
import shutil
import subprocess
import gzip
import csv
import glob
import math
import re
from typing import Optional, Dict
import concurrent.futures
import threading
import hashlib
import uuid

from ..task import Task
from ..reporting import resolve_report_base_path
from ...database import DBManager
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
)

class Decontamination(Task):
    '''A class that handles the decontamination task
        Decontamination takes a set of input accessions, a tree, and a set of reference accessions.
        The most recent labeled run is recorded in Environment_Variables as
        ACTIVE_DECONT_RUN_<library_id> for downstream tasks (e.g. export).
    '''
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=16):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        # Targets (legacy accessions + new targets)
        targets_raw = []
        targets_raw.extend(self.data.get("targets", []) or [])
        targets_raw.extend(self.data.get("accessions", []) or [])
        self.accessions = list(dict.fromkeys(normalize_accessions(targets_raw)))
        # References (legacy ref_accessions + new refs)
        refs_raw = []
        refs_raw.extend(self.data.get("refs", []) or [])
        refs_raw.extend(self.data.get("ref_accessions", []) or [])
        self.ref_accessions = list(dict.fromkeys(normalize_accessions(refs_raw)))
        self.library_id = self.data.get("library_id")
        self.library_name = self.data.get("library_name")
        self.ref_taxid = self.data.get("ref_taxid")
        self.ref_clade = self.data.get("ref_clade")
        self.ref_rule_rank = self.data.get("ref_rule_rank")
        self.ref_rule_quantity = self.data.get("ref_rule_quantity")
        # New selector naming
        self.ref_select_clade = self.data.get("ref_select_clade")
        self.ref_select_rank = self.data.get("ref_select_rank")
        self.ref_select_top = self.data.get("ref_select_top")
        self.out_dir = self.data.get("out_dir")
        self.force = self.data.get("force", False)
        self.stage = checkpoint if checkpoint is not None else 0
        self.rank = (self.data.get("rank") or "phylum").lower()
        self.off_clade_fraction = float(self.data.get("off_clade_fraction", 0.1))
        self.min_buscos = int(self.data.get("min_buscos", 20))
        self.min_identity = float(self.data.get("min_identity", 0.0))
        self.min_coverage = float(self.data.get("min_coverage", 0.0))
        self.min_delta_bitscore = float(self.data.get("min_delta_bitscore", 0.0))
        self.min_bitscore = float(self.data.get("min_bitscore", 0.0))
        self.max_evalue = float(self.data.get("max_evalue", 0.0)) if self.data.get("max_evalue") is not None else None
        self.min_hits = int(self.data.get("min_hits", 1))
        self._hit_window_explicit = "hit_window" in self.data
        hit_window_raw = self.data.get("hit_window")
        if hit_window_raw in (None, ""):
            self.hit_window = 1
        else:
            try:
                self.hit_window = max(1, int(hit_window_raw))
            except (TypeError, ValueError):
                self.hit_window = 1
        self.config_path = self.data.get("config_path")
        self.run_id = self.data.get("run_id")
        self.run_label = self.data.get("run_label")
        self.config_signature = None
        self.use_paralog_filtered_refs = bool(self.data.get("use_paralog_filtered_refs", False))
        self.include_duplicated = bool(self.data.get("include_duplicated", False))
        self.ref_blastdb_path = self.data.get("_ref_blastdb_path")
        self.report_path = self.data.get("report_path")
        if not self.report_path:
            self.report_path = str(
                resolve_report_base_path(
                    self,
                    namespace="decontamination-reports",
                    default_stem="decontamination",
                    run_label=self.run_label or self.library_name or self.library_id,
                    cache_attr="_decontamination_report_dir",
                )
            )
        self.allow_same_species = bool(self.data.get("allow_same_species", False))
        self.allow_sparse_references = bool(self.data.get("allow_sparse_references", False))
        max_concurrent_raw = self.data.get("max_concurrent")
        if max_concurrent_raw in (None, ""):
            self.max_concurrent = required_threads
        else:
            try:
                self.max_concurrent = max(1, int(max_concurrent_raw))
            except (TypeError, ValueError):
                self.max_concurrent = required_threads
        self.blast_threads = 1
        self.blastp_path = None
        self.groups = []
        self.config_accessions = set()
        self.config_targets = set()
        self._taxonomy_cache_lock = threading.RLock()
        self._lineage_cache = {}
        self._rank_cache = {}
        self._descendant_cache = {}

    def _get_lineage_cached(self, taxid: int | None, db=None):
        if not taxid:
            return []
        with self._taxonomy_cache_lock:
            cached = self._lineage_cache.get(taxid)
        if cached is not None:
            return cached
        dbm = db or self.db_manager
        lineage = dbm.get_lineage_root_to_leaf(taxid) or []
        with self._taxonomy_cache_lock:
            self._lineage_cache[taxid] = lineage
        return lineage

    def _get_missing_downloads(self, accs):
        rows = self.db_manager.genomes.get_downloaded()
        present = {r[0] for r in rows} if rows else set()
        return [acc for acc in accs if acc not in present]

    def _is_descendant(self, taxid: int | None, ancestor: int | None, db=None) -> bool:
        if not taxid or not ancestor:
            return False
        if taxid == ancestor:
            return True
        key = (taxid, ancestor)
        with self._taxonomy_cache_lock:
            if key in self._descendant_cache:
                return self._descendant_cache[key]
        lineage = self._get_lineage_cached(taxid, db=db)
        out = any(t[0] == ancestor for t in lineage)
        with self._taxonomy_cache_lock:
            self._descendant_cache[key] = out
        return out

    def _busco_missing_list(self, busco_lib_id):
        present = set(self.db_manager.busco.get_processed_accessions(busco_lib_id))
        return [acc for acc in self.accessions if acc not in present]

    def _queue_busco(self, busco_lineage_name):
        missing = self._busco_missing_list(self.busco_lib_id)
        if not missing:
            return False
        for acc in missing:
            self.queue_subtask(
                job_type=4,
                status="P",
                priority=1,
                data={
                    "accession": acc,
                    "lineage": busco_lineage_name,
                    "format": "protein",
                    "force": self.force,
                },
            )
        return True

    def _get_blastdb_map(self):
        """Return accession -> blastdb path for reference accessions, preferring current library_id when available."""
        rows = self.db_manager.filtering.get_blast_dbs()
        found = {}
        for row in rows:
            _, lib_id, accession, location, _ = row
            if accession not in self.ref_accessions:
                continue
            if accession in found:
                # Prefer entries that match the target library_id when provided
                if self.library_id and found[accession][0] != self.library_id and lib_id == self.library_id:
                    found[accession] = (lib_id, location)
                continue
            found[accession] = (lib_id, location)
        return {acc: loc for acc, (_, loc) in found.items()}

    def _get_missing_blastdbs(self):
        blast_map = self._get_blastdb_map()
        return [acc for acc in self.ref_accessions if acc not in blast_map]

    def _queue_blastdb_creation(self):
        missing = self._get_missing_blastdbs()
        if missing:
            for acc in missing:
                self.queue_subtask(
                    job_type=13,
                    status="P",
                    priority=1,
                    data={
                        "accession": acc,
                        "library_id": self.library_id,
                        "force": True,
                    },
                )
            return True
        return False

    def _taxon_at_rank(self, taxid: int | None, db=None):
        if not taxid:
            return None
        rank_name = (self.rank or "").lower()
        key = (taxid, rank_name)
        with self._taxonomy_cache_lock:
            if key in self._rank_cache:
                return self._rank_cache[key]
        lineage = self._get_lineage_cached(taxid, db=db)
        for tid, _name, rank, _parent in lineage:
            if (rank or "").lower() == rank_name:
                with self._taxonomy_cache_lock:
                    self._rank_cache[key] = tid
                return tid
        with self._taxonomy_cache_lock:
            self._rank_cache[key] = None
        return None

    def _taxon_at_rank_name(self, taxid: int | None, rank_name: str, db=None):
        """Return taxid at requested rank_name within lineage; None if not found."""
        if not taxid:
            return None
        rname = (rank_name or "").lower()
        key = (taxid, rname)
        with self._taxonomy_cache_lock:
            if key in self._rank_cache:
                return self._rank_cache[key]
        lineage = self._get_lineage_cached(taxid, db=db)
        for tid, _name, rank, _parent in lineage:
            if (rank or "").lower() == rname:
                with self._taxonomy_cache_lock:
                    self._rank_cache[key] = tid
                return tid
        with self._taxonomy_cache_lock:
            self._rank_cache[key] = None
        return None

    def _resolve_taxid(self, identifier):
        """Resolve an identifier (accession/taxid/name) to a taxid, if possible."""
        if identifier is None:
            return None
        # Already int?
        try:
            return int(identifier)
        except (TypeError, ValueError):
            self.log(f"Decontam: '{identifier}' is not a numeric taxid; trying accession/clade resolution.", "DEBUG")
        ident = str(identifier).strip()
        # accession?
        genome = self.db_manager.genomes.get(ident)
        if genome and genome[1]:
            return genome[1]
        # clade name
        try:
            return resolve_clade_to_taxid(self.db_manager, ident)
        except (LookupError, ValueError):
            return None

    def _coerce_positive_int(self, value, default=1) -> int:
        try:
            return max(1, int(value))
        except (TypeError, ValueError):
            return default

    def _resolve_hit_window(self, min_hits, hit_window, *, hit_window_explicit=False, context="", warn=True):
        min_hits_val = self._coerce_positive_int(min_hits, 1)
        hit_window_val = self._coerce_positive_int(hit_window, 1)
        if hit_window_explicit:
            if min_hits_val > hit_window_val:
                if warn:
                    self.log(
                        f"Decontam: min_hits ({min_hits_val}) exceeds hit_window ({hit_window_val}){context}; "
                        "capping min_hits to hit_window.",
                        "WARNING",
                    )
                min_hits_val = hit_window_val
        else:
            if min_hits_val > hit_window_val:
                hit_window_val = min_hits_val
        return min_hits_val, hit_window_val

    @staticmethod
    def _decision_window_feasibility(possible_good, possible_bad, *, min_hits, hit_window):
        """Return whether a min_hits-of-hit_window rule can both pass and fail."""
        try:
            possible_good_val = max(0, int(possible_good))
        except (TypeError, ValueError):
            possible_good_val = 0
        try:
            possible_bad_val = max(0, int(possible_bad))
        except (TypeError, ValueError):
            possible_bad_val = 0
        try:
            min_hits_val = max(1, int(min_hits))
        except (TypeError, ValueError):
            min_hits_val = 1
        try:
            hit_window_val = max(1, int(hit_window))
        except (TypeError, ValueError):
            hit_window_val = 1
        if min_hits_val > hit_window_val:
            min_hits_val = hit_window_val

        total_possible = possible_good_val + possible_bad_val
        effective_window = min(hit_window_val, total_possible)
        win_possible = min(possible_good_val, effective_window) >= min_hits_val
        min_bad_for_loss = max(0, effective_window - min_hits_val + 1)
        lose_possible = possible_bad_val >= min_bad_for_loss

        return {
            "possible_good": possible_good_val,
            "possible_bad": possible_bad_val,
            "total_possible": total_possible,
            "min_hits": min_hits_val,
            "hit_window": hit_window_val,
            "effective_window": effective_window,
            "min_bad_for_loss": min_bad_for_loss,
            "win_possible": win_possible,
            "lose_possible": lose_possible,
        }

    def _normalize_hit_window(self):
        min_hits_val, hit_window_val = self._resolve_hit_window(
            self.min_hits,
            self.hit_window,
            hit_window_explicit=self._hit_window_explicit,
        )
        self.min_hits = min_hits_val
        self.hit_window = hit_window_val
        self.data["min_hits"] = self.min_hits
        self.data["hit_window"] = self.hit_window
        try:
            self.db_manager.tasks.update_data(self.task_id, data=self.data)
        except Exception as exc:  # boundary: persist normalized hit-window metadata only.
            self.log(f"Failed to persist normalized decontamination hit-window metadata: {exc}", "WARNING")

    def _load_config(self):
        def _accessions_from_taxid(taxid: int) -> set[str]:
            try:
                return set(
                    self.selector_candidates(
                        taxid=taxid,
                        downloaded_only=True,
                        protein_only=True,
                        status_min=1,
                    )
                )
            except Exception as exc:  # boundary: selector expansion helper failure means this config item contributes no accessions.
                self.log(f"Failed to expand config taxid {taxid} to accessions: {exc}", "WARNING")
                return set()
        if not self.config_path:
            return
        if not os.path.exists(self.config_path):
            self.handle_exception(f"Config file not found: {self.config_path}", {})
            return "ERROR"
        try:
            with open(self.config_path, "r") as handle:
                raw = handle.read()
                self.config_signature = hashlib.sha1(raw.encode()).hexdigest()
                cfg = json.loads(raw)
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            self.handle_exception(f"Failed to read config file: {exc}", {})
            return "ERROR"
        # Normalize shape
        params = {}
        groups_raw = []
        refs_override = []
        targets_override = []
        config_accessions = set()
        config_targets = set()
        if isinstance(cfg, list):
            groups_raw = cfg
        elif isinstance(cfg, dict):
            params = cfg.get("params") or {}
            groups_raw = cfg.get("groups") or cfg.get("acceptance_rules") or []
            refs_override = cfg.get("references") or []
            targets_override = cfg.get("targets") or []
        try:
            self.log(
                f"DECONTAM_TEST config: targets_override={len(targets_override)} refs_override={len(refs_override)} groups={len(groups_raw)}",
                "DEBUG",
            )
        except Exception as exc:  # boundary: diagnostic logging must not affect config loading.
            self.log(f"Failed to log decontamination config diagnostic: {exc}", "WARNING")

        # Apply param overrides
        for key, attr in {
            "rank": "rank",
            "off_clade_fraction": "off_clade_fraction",
            "min_buscos": "min_buscos",
            "min_identity": "min_identity",
            "min_coverage": "min_coverage",
            "min_delta_bitscore": "min_delta_bitscore",
            "min_bitscore": "min_bitscore",
            "max_evalue": "max_evalue",
            "min_hits": "min_hits",
            "hit_window": "hit_window",
            "window": "hit_window",
            "window_size": "hit_window",
            "ref_clade": "ref_clade",
            "ref_rule_rank": "ref_rule_rank",
            "ref_rule_quantity": "ref_rule_quantity",
            "ref_select_clade": "ref_select_clade",
            "ref_select_rank": "ref_select_rank",
            "ref_select_top": "ref_select_top",
            "use_paralog_filtered_refs": "use_paralog_filtered_refs",
            "allow_same_species": "allow_same_species",
            "allow_sparse_references": "allow_sparse_references",
        }.items():
            if key in params:
                setattr(self, attr, params[key])
                if attr == "hit_window":
                    self._hit_window_explicit = True

        if targets_override:
            config_targets.update(normalize_accessions(targets_override))

        if refs_override:
            ref_accs: list[str] = []
            for r in refs_override:
                rstr = str(r)
                normalized_ref = normalize_accessions([rstr])[0]
                genome_row = None
                genome_row = self.db_manager.genomes.get(normalized_ref)
                if genome_row is not None:
                    # Treat as explicit accession: do NOT expand to sibling assemblies
                    ref_accs.append(normalized_ref)
                    continue
                # Try resolving to taxid and expand to downloaded accessions
                taxid = self._resolve_taxid(rstr)
                if taxid:
                    accs = _accessions_from_taxid(taxid)
                    ref_accs.extend(accs)
                    continue
                # Otherwise treat as accession/name directly
                ref_accs.append(rstr)
            merged_refs = []
            if self.ref_accessions:
                merged_refs.extend(self.ref_accessions)
            merged_refs.extend(ref_accs)
            self.ref_accessions = list(dict.fromkeys(normalize_accessions(merged_refs)))
            try:
                self.log(
                    f"DECONTAM_TEST config refs merged: cli_refs={len(self.data.get('ref_accessions', []) or []) + len(self.data.get('refs', []) or [])} "
                    f"config_refs_in={len(refs_override)} final_refs={len(self.ref_accessions)}",
                    "DEBUG",
                )
            except Exception as exc:  # boundary: diagnostic logging must not affect config loading.
                self.log(f"Failed to log decontamination reference diagnostic: {exc}", "WARNING")

        # Parse groups/acceptance rules
        parsed = []
        config_accessions = set(config_accessions) if config_accessions else set()
        config_expanded_accessions = set()
        for g in groups_raw:
            if not isinstance(g, dict):
                continue
            members = g.get("members") or g.get("targets") or []
            clades = g.get("clades") or g.get("allow_clades") or []
            blacklist = g.get("blacklist") or g.get("deny_clades") or []
            min_hits = self._coerce_positive_int(g.get("min_hits", self.min_hits), self.min_hits)
            hit_window_raw = None
            for win_key in ("hit_window", "window", "window_size"):
                if win_key in g:
                    hit_window_raw = g.get(win_key)
                    break
            hit_window = None
            if hit_window_raw not in (None, ""):
                hit_window = self._coerce_positive_int(hit_window_raw, self.hit_window)
            member_accessions = []
            member_taxa = []
            for m in members if isinstance(members, list) else [members]:
                if m is None:
                    continue
                m = str(m).strip()
                # Normalize to help compare against DB accessions
                normalized_member = normalize_accessions([m])[0]
                genome_row = None
                genome_row = self.db_manager.genomes.get(normalized_member)

                is_accession = genome_row is not None
                if is_accession:
                    # Treat as explicit accession: do NOT expand to sibling assemblies
                    member_accessions.append(normalized_member)
                    config_accessions.add(normalized_member)
                else:
                    # Try resolving as taxid/name/clade and expand to downloaded accessions
                    taxid = self._resolve_taxid(m)
                    if taxid:
                        member_taxa.append(taxid)
                        accs = _accessions_from_taxid(taxid)
                        config_expanded_accessions.update(accs)
                    # Also retain the raw member string so explicit name matches still work
                    member_accessions.append(normalized_member)
                    config_accessions.add(normalized_member)
            clade_taxa = [self._resolve_taxid(c) for c in (clades if isinstance(clades, list) else [clades])]
            blacklist_taxa = [self._resolve_taxid(b) for b in (blacklist if isinstance(blacklist, list) else [blacklist])]
            parsed.append(
                {
                    "member_accessions": set(normalize_accessions(member_accessions)),
                    "member_taxa": {t for t in member_taxa if t},
                    "clades": {t for t in clade_taxa if t},
                    "blacklist": {t for t in blacklist_taxa if t},
                    "min_hits": min_hits,
                    **({"hit_window": hit_window} if hit_window is not None else {}),
                }
            )
        self.groups = parsed
        self.config_accessions = set(normalize_accessions(config_expanded_accessions or config_accessions))
        self.config_targets = set(normalize_accessions(config_targets))

    def _group_for_accession(self, accession: str, taxid: int | None, db=None):
        for g in self.groups:
            if accession in g["member_accessions"]:
                return g
            if taxid and any(self._is_descendant(taxid, mt, db=db) for mt in g["member_taxa"]):
                return g
        return None

    def _parse_blast_output(self, blast_output):
        hits = []
        for line in blast_output.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 7:
                continue
            hit = {
                "qseqid": parts[0],
                "sseqid": parts[1],
                "pident": float(parts[2]),
                "length": int(parts[3]),
                "qcovs": float(parts[4]),
                "evalue": float(parts[5]),
                "bitscore": float(parts[6]),
            }
            hits.append(hit)
        return hits

    def _run_blastp(self, query_faa, db_path):
        if not self.blastp_path:
            self.error("BLASTP_PATH not set before attempting blastp.")
            return False
        command = [
            self.blastp_path,
            "-query",
            query_faa,
            "-db",
            db_path,
            "-outfmt",
            "6 qseqid sseqid pident length qcovs evalue bitscore",
        ]
        threads = max(1, int(self.blast_threads)) if getattr(self, "blast_threads", None) else 1
        if threads > 1:
            command.extend(["-num_threads", str(threads)])
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            self.error(f"blastp failed: {result.stderr}")
            return False
        return result.stdout or ""

    def _blast_busco_against_refs(self, query_faa, ref_db_map, combined_db_path=None):
        hits = []
        if combined_db_path:
            out = self._run_blastp(query_faa, combined_db_path)
            if out is False:
                return []
            parsed = self._parse_blast_output(out)
            for hit in parsed:
                # sseqid format accession|family
                sseqid = hit["sseqid"]
                ref_acc = sseqid.split("|")[0] if "|" in sseqid else sseqid
                if self.min_identity and hit["pident"] < self.min_identity:
                    continue
                if self.min_coverage and hit["qcovs"] < self.min_coverage:
                    continue
                if self.min_bitscore and hit["bitscore"] < self.min_bitscore:
                    continue
                if self.max_evalue is not None and hit["evalue"] > self.max_evalue:
                    continue
                hit["ref_accession"] = ref_acc
                hits.append(hit)
        else:
            for ref_acc, db_path in ref_db_map.items():
                out = self._run_blastp(query_faa, db_path)
                if out is False:
                    continue
                parsed = self._parse_blast_output(out)
                for hit in parsed:
                    if self.min_identity and hit["pident"] < self.min_identity:
                        continue
                    if self.min_coverage and hit["qcovs"] < self.min_coverage:
                        continue
                    if self.min_bitscore and hit["bitscore"] < self.min_bitscore:
                        continue
                    if self.max_evalue is not None and hit["evalue"] > self.max_evalue:
                        continue
                    hit["ref_accession"] = ref_acc
                    hits.append(hit)
        hits.sort(key=lambda h: (-h["pident"], h["bitscore"]))
        # hits.sort(key=lambda h: (-h["bitscore"], h["evalue"]))
        return hits

    def _iter_query_records(self, fasta_path):
        header = None
        seq_lines = []
        with open(fasta_path) as handle:
            for raw in handle:
                line = raw.rstrip("\n")
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_lines)
                    header = line[1:].strip()
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            if header is not None:
                yield header, "".join(seq_lines)

    def _pick_duplicate_query_record(self, records, query_id, used):
        query_id = str(query_id or "")
        token = query_id.split()[0] if query_id else ""
        for idx, (header, seq) in enumerate(records):
            if idx in used:
                continue
            first = header.split()[0] if header else ""
            if query_id and (header == query_id or first == query_id or (token and token == first) or query_id in header):
                used.add(idx)
                return header, seq
        for idx, rec in enumerate(records):
            if idx in used:
                continue
            used.add(idx)
            return rec
        return None, None

    def _write_temp_query_record(self, header, sequence):
        tmp = tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False)
        try:
            tmp.write(f">{header}\n{sequence}\n")
            tmp.flush()
            return tmp.name
        finally:
            tmp.close()

    def run(self):
        """
        Execute the decontamination analysis workflow.

        This method orchestrates a multi-stage decontamination pipeline that identifies
        contaminating sequences in genomic assemblies by comparing BUSCO genes against
        reference genomes.

        Workflow Steps:
        1. **Configuration & Validation**
            - Load JSON config overrides/groups
            - Resolve accessions (from CLI or config)
            - Resolve library_id from library_name if needed
            - Determine parent BUSCO library context

        2. **Run Initialization**
            - Generate unique run_id based on reference signatures and parameters
            - Resolve reference accessions via selectors (clade, taxid, rules) or library defaults
            - Validate all accessions are downloaded

        3. **Stage 1: BUSCO Results**
            - Queue BUSCO analysis for query accessions
            - Wait for completion across all queries
            - Retry on failure if configured

        4. **Stage 2: BLAST Database Creation**
            - Queue BLAST database creation for reference genomes
            - Wait for completion
            - Retry with retry_incomplete enabled

        5. **Stage 3: Decontamination Analysis (Parallel)**
            - For each accession in parallel:
              * Retrieve taxonomy and expected taxon at target rank
              * Extract BUSCO genes and BLAST against reference databases
              * Score each BUSCO hit against reference taxonomy
              * Apply group-specific rules or default rank-based comparison
              * Calculate off-clade fraction and majority taxon
              * Assign final decision: CLEAN, CONTAMINATED, or UNCERTAIN
            - Aggregate results and report failures

        Returns:
             str or bool: "ERROR" on critical failure, False if still processing, True on completion.
        """
        # Apply JSON config overrides/groups if provided
        if self._load_config() == "ERROR":
            return "ERROR"
        self._normalize_hit_window()
        # Combine CLI targets with explicit config targets (do not infer targets from acceptance rules)
        try:
            self.log(
                f"DECONTAM_TEST targets pre-merge: cli_targets={len(self.accessions)} config_targets={len(self.config_targets or [])}",
                "DEBUG",
            )
        except Exception as exc:  # boundary: diagnostic logging must not affect task execution.
            self.log(f"Failed to log decontamination target pre-merge diagnostic: {exc}", "WARNING")
        config_targets = list(dict.fromkeys(normalize_accessions(self.config_targets or [])))
        if self.accessions:
            if config_targets:
                merged = list(self.accessions) + config_targets
                self.accessions = list(dict.fromkeys(normalize_accessions(merged)))
                self.data["accessions"] = self.accessions
                try:
                    self.db_manager.tasks.update_data(self.task_id, data=self.data)
                except Exception as exc:  # boundary: persist merged target metadata only.
                    self.log(f"Failed to persist merged decontamination targets: {exc}", "WARNING")
        elif config_targets:
            self.accessions = list(dict.fromkeys(config_targets))
            self.data["accessions"] = self.accessions
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist config target metadata only.
                self.log(f"Failed to persist config decontamination targets: {exc}", "WARNING")
        try:
            self.log(
                f"DECONTAM_TEST targets post-merge: total_targets={len(self.accessions)}",
                "DEBUG",
            )
        except Exception as exc:  # boundary: diagnostic logging must not affect task execution.
            self.log(f"Failed to log decontamination target post-merge diagnostic: {exc}", "WARNING")
        if not self.accessions:
            return self.handle_exception(
                "No accessions provided for decontamination.",
                {"accessions": self.accessions, "config_targets": list(self.config_targets or [])},
            )

        # Resolve library
        if self.library_name and not self.library_id:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
        if not self.library_id:
            return self.handle_exception("library_id (or resolvable library_name) is required.", {"library_id": self.library_id, "library_name": self.library_name})

        # Default selector placeholders (may be overridden below)
        selector_taxid = None
        selector_rank = None
        selector_qty = None

        # Resolve BUSCO context (parent library if this is derived)
        parent_library_id = self.db_manager.libraries.get_parent_id(self.library_id)
        self.busco_lib_id = parent_library_id if parent_library_id else self.library_id
        libs = self.db_manager.libraries.get(self.busco_lib_id)
        if not libs:
            return self.handle_exception("BUSCO library not found in database.", {"busco_library_id": self.busco_lib_id})
        busco_lineage_name = libs[0][1]

        # Precompute target species set (for same-species filtering)
        target_species = set()
        for acc in self.accessions:
            g = self.db_manager.genomes.get(acc)
            taxid = g[1] if g and len(g) > 1 else None
            sp = self._taxon_at_rank_name(taxid, "species")
            if sp:
                target_species.add(sp)

        # Resolve references: explicit refs + selector-based (both can be used)
        selector_taxid = self.ref_taxid
        if not selector_taxid and self.ref_clade:
            try:
                selector_taxid = resolve_clade_to_taxid(self.db_manager, self.ref_clade)
            except (LookupError, ValueError):
                selector_taxid = None
        if not selector_taxid and self.ref_select_clade:
            try:
                selector_taxid = resolve_clade_to_taxid(self.db_manager, self.ref_select_clade)
            except (LookupError, ValueError):
                selector_taxid = None
        selector_rank = self.ref_rule_rank or self.ref_select_rank
        selector_qty = self.ref_rule_quantity or self.ref_select_top

        selector_refs = []
        if selector_taxid or selector_rank or selector_qty:
            try:
                # First gather candidates
                selector_candidates = self.selector_candidates(
                    taxid=selector_taxid,
                    downloaded_only=True,
                    protein_only=True,
                    status_min=1,
                )
                if not selector_candidates:
                    raise ValueError("No accessions matched the provided selectors.")
                # Apply same-species exclusion before rule selection
                if not self.allow_same_species and target_species:
                    filtered_candidates = []
                    for c in selector_candidates:
                        g = self.db_manager.genomes.get(c)
                        taxid = g[1] if g and len(g) > 1 else None
                        sp = self._taxon_at_rank_name(taxid, "species")
                        if sp and sp in target_species:
                            continue
                        filtered_candidates.append(c)
                    selector_candidates = filtered_candidates
                if not selector_candidates and not self.allow_sparse_references:
                    raise ValueError("No selector accessions remain after same-species filtering.")
                # Apply rule selection (with BUSCO completeness weighting when available)
                selector_refs = self.apply_selector_rules(
                    selector_candidates,
                    taxid=selector_taxid,
                    rule_quantity=selector_qty,
                    rule_rank=selector_rank,
                    busco_library_id=self.busco_lib_id,
                )
                selector_refs = list(dict.fromkeys(normalize_accessions(selector_refs)))
                # Enforce minimum expected refs unless sparse allowed
                if selector_qty and selector_rank and selector_taxid:
                    # Count how many groups exist in the candidate pool at the rule rank
                    rank_token = str(selector_rank).strip().lower()
                    # Map accession -> taxid
                    acc_tax = {}
                    if selector_candidates:
                        placeholders = ",".join("?" for _ in selector_candidates)
                        self.db_manager.cursor.execute(
                            f"SELECT accession, taxid FROM Genome WHERE accession IN ({placeholders})",
                            tuple(selector_candidates),
                        )
                        rows = self.db_manager.cursor.fetchall() or []
                        acc_tax = {str(a): (int(t) if t is not None else None) for a, t in rows}
                    group_names = set()
                    for acc, tax in acc_tax.items():
                        if tax is None:
                            continue
                        lineage_rows = self.db_manager.genomes.get_lineage_root_to_leaf(tax) or []
                        for tid, name, rnk, _parent in lineage_rows:
                            if (rnk or "").lower() == rank_token:
                                group_names.add(name or str(tid))
                                break
                    required_total = int(selector_qty) * len(group_names)
                    if not self.allow_sparse_references and len(selector_refs) < required_total:
                        raise ValueError(
                            f"Selectors require at least {required_total} references at rank {selector_rank}; "
                            f"only {len(selector_refs)} available after filtering."
                        )
            except Exception as exc:  # boundary: selector/reference resolution failure becomes this task error.
                return self.handle_exception(
                    f"Failed to resolve reference selectors: {exc}",
                    {"ref_taxid": selector_taxid, "ref_rule_rank": selector_rank, "ref_rule_quantity": selector_qty},
                )

        refs_explicit_norm = list(dict.fromkeys(normalize_accessions(self.ref_accessions or [])))
        refs = list(refs_explicit_norm)
        refs.extend(selector_refs)
        try:
            self.log(
                f"DECONTAM_TEST refs before species filter: explicit_refs={len(refs_explicit_norm)} selector_refs={len(selector_refs)} merged_refs={len(refs)}",
                "DEBUG",
            )
        except Exception as exc:  # boundary: diagnostic logging must not affect task execution.
            self.log(f"Failed to log decontamination refs pre-filter diagnostic: {exc}", "WARNING")
        self._report_refs_explicit = set(refs_explicit_norm)
        self._report_refs_selected = set(selector_refs)
        self.log(f"Decontam: targets={self.accessions} refs_explicit={refs_explicit_norm} refs_selected={selector_refs}", "DEBUG")

        # Exclude same-species references unless explicitly allowed
        if not self.allow_same_species and refs:
            filtered_refs = []
            skipped = []
            for ref in refs:
                g = self.db_manager.genomes.get(ref)
                taxid = g[1] if g and len(g) > 1 else None
                sp = self._taxon_at_rank_name(taxid, "species")
                if sp and sp in target_species:
                    skipped.append(ref)
                    continue
                filtered_refs.append(ref)
            if skipped:
                self.log(
                    f"Filtered {len(skipped)} references that match target species "
                    f"(disable with --allow-same-species): {', '.join(skipped)}",
                    "DEBUG",
                )
            refs = filtered_refs
            try:
                self.log(
                    f"DECONTAM_TEST refs after species filter: kept_refs={len(refs)} skipped_refs={len(skipped)}",
                    "DEBUG",
                )
            except Exception as exc:  # boundary: diagnostic logging must not affect task execution.
                self.log(f"Failed to log decontamination refs post-filter diagnostic: {exc}", "WARNING")
            if not refs and not self.allow_sparse_references:
                return self.handle_exception(
                    "No references available after same-species filtering.",
                    {"ref_accessions": self.ref_accessions, "selector_refs": selector_refs},
                )

        if not refs:
            refs = self.db_manager.libraries.get_reference_assemblies(self.library_id) or []
        self.ref_accessions = list(dict.fromkeys(normalize_accessions(refs)))
        try:
            overlap = set(self.accessions) & set(self.ref_accessions)
            self.log(
                f"DECONTAM_TEST refs finalized: total_refs={len(self.ref_accessions)} target_ref_overlap={len(overlap)}",
                "DEBUG",
            )
        except Exception as exc:  # boundary: diagnostic logging must not affect task execution.
            self.log(f"Failed to log finalized decontamination refs diagnostic: {exc}", "WARNING")
        self.data["ref_accessions"] = self.ref_accessions
        try:
            self.db_manager.tasks.update_data(self.task_id, data=self.data)
        except Exception as exc:  # boundary: persist resolved reference metadata only.
            self.log(f"Failed to persist decontamination reference accessions: {exc}", "WARNING")
        if not self.ref_accessions:
            return self.handle_exception("No reference accessions provided or recorded for this library.", {"library_id": self.library_id})

        self.log(f"Decontam: run_id={self.run_id} run_label={self.run_label} selector_taxid={selector_taxid} selector_rank={selector_rank} selector_qty={selector_qty}", "DEBUG")

        # Generate deterministic run_id (unless provided); allow deduplication
        if not self.run_id:
            fingerprint = {
                "targets": sorted(self.accessions),
                "refs": sorted(self.ref_accessions or []),
                "rank": self.rank,
                "off_clade_fraction": self.off_clade_fraction,
                "min_buscos": self.min_buscos,
                "min_identity": self.min_identity,
                "min_coverage": self.min_coverage,
                "min_delta_bitscore": self.min_delta_bitscore,
                "min_bitscore": self.min_bitscore,
                "max_evalue": self.max_evalue,
                "min_hits": self.min_hits,
                "hit_window": self.hit_window,
                "include_duplicated": self.include_duplicated,
                "config_signature": self.config_signature,
                "selector": {
                    "taxid": selector_taxid,
                    "rank": selector_rank,
                    "qty": selector_qty,
                },
            }
            digest = hashlib.sha1(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
            self.run_id = f"dc_{digest}"
            self.data["run_id"] = self.run_id
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist generated run id metadata only.
                self.log(f"Failed to persist generated decontamination run id {self.run_id}: {exc}", "WARNING")

        # Skip duplicate runs unless force is set
        if not self.force:
            existing = self.db_manager.filtering.get_decontamination_summary(run_id=self.run_id)
            if existing:
                self.log(f"Decontamination run_id {self.run_id} already exists; skipping (use --force to rerun).", "WARNING")
                return True

        # Stage 2: build combined BUSCO BLAST DB for references (only BUSCO mode supported)
        if not self.ref_blastdb_path:
            tmpdir = tempfile.mkdtemp(prefix=f"dc_refdb_{self.run_id}_")
            self.ref_blastdb_path = os.path.join(tmpdir, "ref_buscos")
            self.data["_ref_blastdb_path"] = self.ref_blastdb_path
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist temporary reference BLAST DB path only.
                self.log(f"Failed to persist decontamination reference BLAST DB path: {exc}", "WARNING")

        def _queue_ref_busco_db():
            self.queue_subtask(
                job_type=17,
                status="P",
                priority=1,
                data={
                    "accessions": self.ref_accessions,
                    "busco_library_id": self.busco_lib_id,
                    "target_library_id": self.library_id,
                    "output_path": self.ref_blastdb_path,
                    "use_paralog_filtered": self.use_paralog_filtered_refs,
                    "force": True,
                },
            )
            return True

        def _ref_db_done():
            required_exts = ["psq", "pin", "phr"]
            return all(os.path.exists(f"{self.ref_blastdb_path}.{ext}") for ext in required_exts)

        outcome = self.manage_subtasks(
            stage=2,
            queue_fn=_queue_ref_busco_db,
            done_fn=_ref_db_done,
            wait_seconds=0,
            retry_key=None,
            max_retries=0,
            incomplete_message_fn=lambda: (f"Reference BUSCO BLAST DB not ready at {self.ref_blastdb_path}", ""),
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        # Record run metadata
        run_params = {
            "rank": self.rank,
            "off_clade_fraction": self.off_clade_fraction,
            "min_buscos": self.min_buscos,
            "min_identity": self.min_identity,
            "min_coverage": self.min_coverage,
            "min_delta_bitscore": self.min_delta_bitscore,
            "min_bitscore": self.min_bitscore,
            "max_evalue": self.max_evalue,
            "min_hits": self.min_hits,
            "hit_window": self.hit_window,
            "include_duplicated": self.include_duplicated,
            "ref_select_taxid": selector_taxid if 'selector_taxid' in locals() else None,
            "ref_select_rank": selector_rank if 'selector_rank' in locals() else None,
            "ref_select_top": selector_qty if 'selector_qty' in locals() else None,
            "config_path": self.config_path,
            "config_signature": self.config_signature,
            "run_label": self.run_label,
        }
        if not self.db_manager.filtering.add_decontamination_run(
            run_id=self.run_id,
            target_library_id=self.library_id,
            busco_library_id=self.busco_lib_id,
            targets_json=json.dumps(self.accessions),
            refs_json=json.dumps(self.ref_accessions),
            params_json=json.dumps(run_params),
            config_signature=self.config_signature,
            run_label=self.run_label,
        ):
            return self.handle_exception("Failed to record decontamination run metadata.", {"run_id": self.run_id})

        # Ensure genomes are downloaded for queries and references
        missing_downloads = self._get_missing_downloads(self.accessions + self.ref_accessions)
        if missing_downloads:
            try:
                self.log(
                    f"DECONTAM_TEST missing_downloads={len(missing_downloads)} (showing up to 10): {missing_downloads[:10]}",
                    "DEBUG",
                )
            except Exception as exc:  # boundary: diagnostic logging must not mask the primary missing-download error.
                self.log(f"Failed to log decontamination missing-download diagnostic: {exc}", "WARNING")
            return self.handle_exception(
                "Some accessions are not downloaded; download them before decontamination.",
                {"missing_accessions": missing_downloads},
            )

        # Stage 1: ensure BUSCO results for queries
        self.stage = 1
        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=lambda: self._queue_busco(busco_lineage_name),
            done_fn=lambda: len(self._busco_missing_list(self.busco_lib_id)) == 0,
            wait_seconds=int(self.data.get("busco_wait_seconds", 0)),
            retry_key="busco_retries",
            max_retries=int(self.data.get("busco_retries", 0)),
            incomplete_message_fn=lambda: (
                f"BUSCO results not yet available for all accessions ({busco_lineage_name}).",
                "",
            ),
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        # Ready to analyze
        self.stage = 3
        self.blastp_path = self.db_manager.env.get("BLASTP_PATH")
        if not self.blastp_path:
            return self.handle_exception("BLASTP_PATH is not set in environment variables.", {"variable": "BLASTP_PATH"})

        ref_db_map = {}

        # Precompute taxonomy for references
        ref_rank_taxa = {}
        ref_taxa_full = {}
        for ref in self.ref_accessions:
            genome = self.db_manager.genomes.get(ref)
            taxid = genome[1] if genome else None
            ref_taxa_full[ref] = taxid
            ref_rank_taxa[ref] = self._taxon_at_rank(taxid)

        # Filter acceptance rules to only those that match provided targets and refs
        if self.groups:
            target_taxa = {}
            for acc in self.accessions:
                g = self.db_manager.genomes.get(acc)
                taxid = g[1] if g and len(g) > 1 else None
                target_taxa[acc] = taxid
            active_groups = []
            for idx, g in enumerate(self.groups):
                target_match = False
                for acc, taxid in target_taxa.items():
                    if acc in g["member_accessions"]:
                        target_match = True
                        break
                    if taxid and any(self._is_descendant(taxid, mt, db=self.db_manager) for mt in g["member_taxa"]):
                        target_match = True
                        break
                if not target_match:
                    self.log(
                        f"Decontam: ignoring acceptance rule {idx} (no targets match members).",
                        "WARNING",
                    )
                    continue
                if g.get("clades"):
                    ref_match = False
                    for ref_taxid in ref_taxa_full.values():
                        if ref_taxid and any(self._is_descendant(ref_taxid, cl, db=self.db_manager) for cl in g["clades"]):
                            ref_match = True
                            break
                    if not ref_match:
                        self.log(
                            f"Decontam: ignoring acceptance rule {idx} (no references in allowed clades).",
                            "WARNING",
                        )
                        continue
                active_groups.append(g)
            self.groups = active_groups

        params_json = json.dumps(
            {
                "rank": self.rank,
                "off_clade_fraction": self.off_clade_fraction,
                "min_buscos": self.min_buscos,
                "min_identity": self.min_identity,
                "min_coverage": self.min_coverage,
                "min_delta_bitscore": self.min_delta_bitscore,
                "min_bitscore": self.min_bitscore,
                "max_evalue": self.max_evalue,
                "min_hits": self.min_hits,
                "hit_window": self.hit_window,
                "include_duplicated": self.include_duplicated,
                "references": self.ref_accessions,
                "config_path": self.config_path,
                "config_signature": self.config_signature,
                "run_id": self.run_id,
                "run_label": self.run_label,
            }
        )

        # Parallelization across accessions if needed
        def process_accession(acc):
            thread_db = DBManager(self.db_manager.get_path())
            try:
                thread_db.connect()
                thread_db.set_busco_run_context(
                    pipeline=self.data.get("busco_pipeline"),
                    input_mode=self.data.get("busco_input_mode"),
                    prefer_pipeline=self.data.get("prefer_busco_pipeline"),
                    prefer_input_mode=self.data.get("prefer_busco_input_mode"),
                    run_ids=self.data.get("busco_run_ids"),
                    selection=self.data.get("busco_run_selection") or "primary",
                )
            except Exception as exc:  # boundary: worker DB setup failure is reported for this accession.
                self.log(f"Failed to open DB in thread for {acc}: {exc}", "ERROR")
                return acc, "ERROR"

            genome = thread_db.genomes.get(acc)
            taxid = genome[1] if genome else None
            expected_taxon = self._taxon_at_rank(taxid, db=thread_db)
            effective_busco_run_id = thread_db.busco.get_effective_run_id_for_accession(
                acc,
                self.busco_lib_id,
                purpose="default",
            )
            group = self._group_for_accession(acc, taxid, db=thread_db)
            if not expected_taxon:
                self.log(f"No taxonomy at rank={self.rank} for {acc}; marking uncertain.", "WARNING")
                thread_db.filtering.add_decontamination_summaries(
                    [
                        {
                            "accession": acc,
                            "target_library_id": self.library_id,
                            "busco_library_id": self.busco_lib_id,
                            "run_id": self.run_id,
                            "busco_run_id": effective_busco_run_id,
                            "expected_taxid": None,
                            "majority_taxid": None,
                            "rank": self.rank,
                            "buscos_tested": 0,
                            "buscos_supporting": 0,
                            "buscos_outside": 0,
                            "off_clade_fraction": None,
                            "decision": "UNCERTAIN",
                            "params_json": params_json,
                        }
                    ]
                )
                thread_db.close()
                return acc, "UNCERTAIN"

            busco_rows = thread_db.busco.get_family_results_for_library(
                library_id=self.busco_lib_id,
                accessions=[acc],
                status=[1, 2] if self.include_duplicated else [1],
            )
            if not busco_rows:
                self.log(f"No BUSCO entries found for {acc} in library {self.busco_lib_id}.", "WARNING")
                thread_db.filtering.add_decontamination_summaries(
                    [
                        {
                            "accession": acc,
                            "target_library_id": self.library_id,
                            "busco_library_id": self.busco_lib_id,
                            "run_id": self.run_id,
                            "busco_run_id": effective_busco_run_id,
                            "expected_taxid": expected_taxon,
                            "majority_taxid": None,
                            "rank": self.rank,
                            "buscos_tested": 0,
                            "buscos_supporting": 0,
                            "buscos_outside": 0,
                            "off_clade_fraction": None,
                            "decision": "UNCERTAIN",
                            "params_json": params_json,
                        }
                    ]
                )
                thread_db.close()
                return acc, "UNCERTAIN"

            buscos_tested = 0
            buscos_supporting = 0
            buscos_outside = 0
            best_taxon_counts = {}
            vote_rows: list[dict[str, Any]] = []
            copy_vote_rows: list[dict[str, Any]] = []

            group_idx = None
            if group:
                try:
                    group_idx = self.groups.index(group)
                except ValueError:
                    group_idx = None
            group_min_hits = group.get("min_hits", self.min_hits) if group else self.min_hits
            group_hit_window_explicit = bool(group and "hit_window" in group)
            group_hit_window = group.get("hit_window") if group_hit_window_explicit else self.hit_window
            min_hits_req, hit_window = self._resolve_hit_window(
                group_min_hits,
                group_hit_window,
                hit_window_explicit=group_hit_window_explicit,
                context=f" for group {group_idx}" if group_idx is not None else "",
            )

            family_ids = [row[0] for row in busco_rows]
            ref_busco_map = {}
            if family_ids:
                ref_busco_map = thread_db.busco.get_family_presence_map(
                    library_id=self.busco_lib_id,
                    accessions=self.ref_accessions,
                    family_ids=family_ids,
                )

            def _evaluate_hits(family_id, hits):
                allowed_flags = []
                for h in hits:
                    h_taxon_full = ref_taxa_full.get(h["ref_accession"])
                    h_taxon_rank = ref_rank_taxa.get(h["ref_accession"])
                    if group and group.get("clades"):
                        if any(self._is_descendant(h_taxon_full, blk, db=thread_db) for blk in group.get("blacklist", [])):
                            allowed_flags.append(False)
                        elif group["clades"] and not any(self._is_descendant(h_taxon_full, cl, db=thread_db) for cl in group["clades"]):
                            allowed_flags.append(False)
                        else:
                            allowed_flags.append(True)
                    else:
                        allowed_flags.append(h_taxon_rank == expected_taxon)

                top_hits = hits[:hit_window]
                top_hits_report = hits[: hit_window + 1]
                top_allowed = allowed_flags[:hit_window]
                acceptable = sum(1 for ok in top_allowed if ok) >= min_hits_req

                acceptable_ref_count = 0
                refs_with_family = ref_busco_map.get(family_id, set())
                if refs_with_family:
                    if group and group.get("clades"):
                        for ref_acc in refs_with_family:
                            tax_full = ref_taxa_full.get(ref_acc)
                            if any(self._is_descendant(tax_full, blk, db=thread_db) for blk in group.get("blacklist", [])):
                                continue
                            if group["clades"] and not any(
                                self._is_descendant(tax_full, cl, db=thread_db) for cl in group["clades"]
                            ):
                                continue
                            acceptable_ref_count += 1
                    else:
                        for ref_acc in refs_with_family:
                            if ref_rank_taxa.get(ref_acc) == expected_taxon:
                                acceptable_ref_count += 1
                outside_ref_count = max(len(refs_with_family) - acceptable_ref_count, 0)
                feasibility = self._decision_window_feasibility(
                    acceptable_ref_count,
                    outside_ref_count,
                    min_hits=min_hits_req,
                    hit_window=hit_window,
                )

                best_allowed_idx = next((i for i, ok in enumerate(allowed_flags) if ok), None)
                best_outside_idx = next((i for i, ok in enumerate(allowed_flags) if not ok), None)

                decision = "outside"
                delta_bitscore = None
                best_rank_taxon = None
                best_outside_rank_taxon = None
                best_bitscore = None

                if not feasibility["win_possible"] or not feasibility["lose_possible"]:
                    if best_allowed_idx is not None:
                        best = hits[best_allowed_idx]
                        best_bitscore = best["bitscore"]
                        best_rank_taxon = ref_rank_taxa.get(best["ref_accession"])
                    if best_outside_idx is not None:
                        best_outside = hits[best_outside_idx]
                        if best_bitscore is None:
                            best_bitscore = best_outside["bitscore"]
                        best_outside_rank_taxon = ref_rank_taxa.get(best_outside["ref_accession"])
                        if best_allowed_idx is not None:
                            delta_bitscore = best_bitscore - best_outside["bitscore"]
                    decision = "unknown"
                elif best_allowed_idx is None:
                    if best_outside_idx is not None:
                        best_outside = hits[best_outside_idx]
                        best_bitscore = best_outside["bitscore"]
                        best_outside_rank_taxon = ref_rank_taxa.get(best_outside["ref_accession"])
                    decision = "outside" if acceptable_ref_count >= min_hits_req else "unknown"
                elif acceptable:
                    best = hits[best_allowed_idx]
                    best_bitscore = best["bitscore"]
                    best_rank_taxon = ref_rank_taxa.get(best["ref_accession"])
                    decision = "support"
                    if best_outside_idx is not None:
                        best_outside = hits[best_outside_idx]
                        best_outside_rank_taxon = ref_rank_taxa.get(best_outside["ref_accession"])
                        delta_bitscore = best_bitscore - best_outside["bitscore"]
                        if self.min_delta_bitscore and delta_bitscore < self.min_delta_bitscore:
                            decision = "weak"
                else:
                    best = hits[best_allowed_idx]
                    best_bitscore = best["bitscore"]
                    best_rank_taxon = ref_rank_taxa.get(best["ref_accession"])
                    decision = "outside"
                    if best_outside_idx is not None:
                        best_outside = hits[best_outside_idx]
                        best_outside_rank_taxon = ref_rank_taxa.get(best_outside["ref_accession"])
                        delta_bitscore = best_bitscore - best_outside["bitscore"]

                return (
                    decision,
                    best_rank_taxon,
                    best_outside_rank_taxon,
                    best_bitscore,
                    delta_bitscore if delta_bitscore is not None else 0,
                    top_hits_report,
                )

            duplicate_groups = {}
            if self.include_duplicated:
                for row in busco_rows:
                    if int(row[3] or 0) == 2:
                        duplicate_groups.setdefault(str(row[0]), []).append(row)

            for family_id, lib_id, accession, status, sequence, score, length in busco_rows:
                status_val = int(status or 0)
                if status_val == 2:
                    continue
                location = thread_db.busco.get_family_location(
                    family_id,
                    lib_id,
                    accession,
                    sequence_kind="prot",
                )
                if not location or not os.path.exists(location):
                    continue
                hits = self._blast_busco_against_refs(location, ref_db_map, combined_db_path=self.ref_blastdb_path)
                if not hits:
                    continue
                decision, best_rank_taxon, best_outside_rank_taxon, best_bitscore, delta_bitscore, top_hits_report = _evaluate_hits(
                    family_id,
                    hits,
                )
                buscos_tested += 1
                if decision in ("support", "weak"):
                    buscos_supporting += 1
                elif decision == "outside":
                    buscos_outside += 1

                if best_rank_taxon:
                    best_taxon_counts[best_rank_taxon] = best_taxon_counts.get(best_rank_taxon, 0) + 1

                vote_rows.append(
                    {
                        "family_id": family_id,
                        "busco_library_id": self.busco_lib_id,
                        "target_library_id": self.library_id,
                        "accession": acc,
                        "run_id": self.run_id,
                        "busco_run_id": effective_busco_run_id,
                        "expected_taxid": expected_taxon,
                        "best_taxid": best_rank_taxon,
                        "runner_taxid": best_outside_rank_taxon,
                        "rank": self.rank,
                        "best_bitscore": best_bitscore,
                        "delta_bitscore": delta_bitscore,
                        "decision": decision,
                        "top_hits_json": json.dumps(top_hits_report) if top_hits_report else None,
                    }
                )

            for family_id, dup_rows in duplicate_groups.items():
                lib_id = dup_rows[0][1]
                accession = dup_rows[0][2]
                location = thread_db.busco.get_family_location(
                    family_id,
                    lib_id,
                    accession,
                    sequence_kind="prot",
                )
                if not location or not os.path.exists(location):
                    continue
                records = list(self._iter_query_records(location))
                if not records:
                    continue
                used_records = set()
                for _family_id, _lib_id, _accession, status, sequence, _score, _length in dup_rows:
                    query_header, query_sequence = self._pick_duplicate_query_record(records, sequence, used_records)
                    if not query_header or not query_sequence:
                        continue
                    query_path = self._write_temp_query_record(query_header, query_sequence)
                    try:
                        hits = self._blast_busco_against_refs(query_path, ref_db_map, combined_db_path=self.ref_blastdb_path)
                    finally:
                        try:
                            os.remove(query_path)
                        except OSError as exc:
                            self.log(f"Failed to remove temporary decontamination query file {query_path}: {exc}", "WARNING")
                    if not hits:
                        continue
                    decision, best_rank_taxon, best_outside_rank_taxon, best_bitscore, delta_bitscore, top_hits_report = _evaluate_hits(
                        family_id,
                        hits,
                    )
                    copy_vote_rows.append(
                        {
                            "family_id": family_id,
                            "busco_library_id": self.busco_lib_id,
                            "target_library_id": self.library_id,
                            "accession": acc,
                            "run_id": self.run_id,
                            "busco_run_id": effective_busco_run_id,
                            "query_id": sequence,
                            "query_header": query_header,
                            "query_status": status,
                            "expected_taxid": expected_taxon,
                            "best_taxid": best_rank_taxon,
                            "runner_taxid": best_outside_rank_taxon,
                            "rank": self.rank,
                            "best_bitscore": best_bitscore,
                            "delta_bitscore": delta_bitscore,
                            "decision": decision,
                            "top_hits_json": json.dumps(top_hits_report) if top_hits_report else None,
                        }
                    )

            if buscos_tested == 0:
                off_frac = None
            else:
                off_frac = buscos_outside / buscos_tested
            majority_taxon = None
            if best_taxon_counts:
                majority_taxon = max(best_taxon_counts.items(), key=lambda kv: kv[1])[0]

            allowed_majority = False
            if majority_taxon:
                if group and group.get("clades"):
                    allowed_majority = any(
                        self._is_descendant(majority_taxon, cl, db=thread_db)
                        for cl in group.get("clades", [])
                    )
                else:
                    allowed_majority = majority_taxon == expected_taxon

            if buscos_tested < self.min_buscos:
                final_decision = "UNCERTAIN"
            elif off_frac is not None and off_frac <= self.off_clade_fraction and allowed_majority:
                final_decision = "CLEAN"
            elif off_frac is not None and off_frac > self.off_clade_fraction:
                final_decision = "CONTAMINATED"
            elif majority_taxon and not allowed_majority:
                final_decision = "CONTAMINATED"
            else:
                final_decision = "UNCERTAIN"

            if vote_rows and not thread_db.filtering.add_decontamination_votes(vote_rows):
                raise RuntimeError(f"Failed to persist decontamination votes for {acc}")
            if copy_vote_rows and not thread_db.filtering.add_decontamination_copy_votes(copy_vote_rows):
                raise RuntimeError(f"Failed to persist decontamination copy votes for {acc}")
            if not thread_db.filtering.add_decontamination_summaries(
                [
                    {
                        "accession": acc,
                        "target_library_id": self.library_id,
                        "busco_library_id": self.busco_lib_id,
                        "run_id": self.run_id,
                        "busco_run_id": effective_busco_run_id,
                        "expected_taxid": expected_taxon,
                        "majority_taxid": majority_taxon,
                        "rank": self.rank,
                        "buscos_tested": buscos_tested,
                        "buscos_supporting": buscos_supporting,
                        "buscos_outside": buscos_outside,
                        "off_clade_fraction": off_frac,
                        "decision": final_decision,
                        "params_json": params_json,
                    }
                ]
            ):
                raise RuntimeError(f"Failed to persist decontamination summary for {acc}")
            try:
                thread_db.close()
            except Exception as exc:  # boundary: worker cleanup failure is logged after primary work completes/fails.
                self.log(f"Failed to close worker database connection for decontamination: {exc}", "WARNING")
            return acc, final_decision

        effective_max = self.max_concurrent if self.max_concurrent and self.max_concurrent > 0 else self.REQUIRED_THREADS
        max_workers = max(1, min(len(self.accessions), effective_max, self.REQUIRED_THREADS))
        # Allow fewer workers when few targets, but allocate more threads per BLAST process
        self.blast_threads = max(1, math.ceil(self.REQUIRED_THREADS / max_workers)) if max_workers else self.REQUIRED_THREADS
        self.log(
            f"Stage 3: running decontamination with max_workers={max_workers} (blast threads per task={self.blast_threads}); targets={len(self.accessions)} refs={len(ref_db_map)}",
            "INFO",
        )
        results = []
        progress_total = len(self.accessions)
        progress_done = 0
        progress_next = 10
        progress_lock = threading.Lock()

        def _log_progress():
            nonlocal progress_done, progress_next
            if progress_total <= 0:
                return
            with progress_lock:
                progress_done += 1
                pct = int((progress_done / progress_total) * 100)
                while progress_next <= pct and progress_next <= 100:
                    self.log(
                        f"Decontam progress: {progress_next}% ({progress_done}/{progress_total})",
                        "INFO",
                    )
                    progress_next += 10
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="DC_") as executor:
            futures = {executor.submit(process_accession, acc): acc for acc in self.accessions}
            for future in concurrent.futures.as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as e:  # boundary: one accession failure is aggregated while independent targets continue.
                    acc = futures[future]
                    self.error(f"Decontamination failed for {acc}: {e}")
                    results.append((acc, "ERROR"))
                _log_progress()

        failed = [acc for acc, status in results if status in ("ERROR",)]
        if failed:
            preview = ", ".join(failed[:10])
            suffix = "" if len(failed) <= 10 else ", ..."
            self.log(
                f"Decontamination completed with issues for {len(failed)}/{len(self.accessions)} targets: {preview}{suffix}",
                "WARNING",
            )
        else:
            self.log(f"Decontamination completed for {len(self.accessions)} targets.", "INFO")

        # Optional report writing
        if self.report_path:
            try:
                base, ext = os.path.splitext(self.report_path)
                busco_report = self.report_path if ext else f"{base}_buscos.tsv"
                summary_report = f"{base}_summary.tsv" if ext else f"{self.report_path}_summary.tsv"
                refs_report = f"{base}_refs.tsv" if ext else f"{self.report_path}_refs.tsv"
                groups_report = f"{base}_groups.tsv" if ext else f"{self.report_path}_groups.tsv"
                blast_report = f"{base}_blast.tsv" if ext else f"{self.report_path}_blast.tsv"
                os.makedirs(os.path.dirname(busco_report) or ".", exist_ok=True)

                votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)
                # Build a quick map of decisions per target for unknown counts
                decisions_by_target = {}
                for v in votes:
                    acc = v[3]
                    dec = v[11]
                    decisions_by_target.setdefault(acc, []).append(dec)

                def tax_name(tid):
                    if tid is None:
                        return None
                    try:
                        self.db_manager.cursor.execute("SELECT name FROM Taxonomy WHERE taxid = ?", (tid,))
                        res = self.db_manager.cursor.fetchone()
                        return res[0] if res else None
                    except Exception as exc:  # boundary: optional taxonomy label lookup for report.
                        self.log(f"Failed to look up taxonomy name for {tid}: {exc}", "WARNING")
                        return None

                rank_cache = {}

                def acc_rank_taxon(acc, rank_override=None):
                    key = (acc, rank_override or self.rank)
                    if key in rank_cache:
                        return rank_cache[key]
                    try:
                        g = self.db_manager.genomes.get(acc)
                        taxid = g[1] if g and len(g) > 1 else None
                        tid_at_rank = self._taxon_at_rank_name(taxid, rank_override or self.rank)
                        name = tax_name(tid_at_rank)
                        rank_cache[key] = (tid_at_rank, name, taxid)
                    except Exception as exc:  # boundary: optional taxonomy label lookup for report.
                        self.log(f"Failed to resolve report taxonomy for {acc}: {exc}", "WARNING")
                        rank_cache[key] = (None, None, None)
                    return rank_cache[key]

                # Precompute which references contain each BUSCO family
                ref_busco_map = {}
                if votes and self.ref_accessions:
                    family_ids = [v[0] for v in votes]
                    ref_busco_map = self.db_manager.busco.get_family_presence_map(
                        library_id=self.busco_lib_id,
                        accessions=self.ref_accessions,
                        family_ids=family_ids,
                    )

                ref_taxa_full = {}
                ref_rank_taxa = {}
                for ref_acc in self.ref_accessions:
                    tid_at_rank, _name, raw_taxid = acc_rank_taxon(ref_acc)
                    ref_rank_taxa[ref_acc] = tid_at_rank
                    ref_taxa_full[ref_acc] = raw_taxid

                vote_keys = {(v[3], v[0]) for v in votes} if votes else set()
                busco_families_by_acc = {}
                if self.accessions:
                    def _chunk(seq, size):
                        for i in range(0, len(seq), size):
                            yield seq[i:i + size]
                    for chunk in _chunk(self.accessions, 800):
                        rows = self.db_manager.busco.get_family_results_for_library(
                            library_id=self.busco_lib_id,
                            accessions=list(chunk),
                            status=[1],
                        )
                        for fam, _lib_id, acc, *_rest in rows or []:
                            busco_families_by_acc.setdefault(str(acc), set()).add(str(fam))

                with open(busco_report, "w") as fh, open(blast_report, "w") as fh_blast:
                    fh.write(
                        "\t".join(
                            [
                                "family_id",
                                "accession",
                                "busco_status",
                                "target_rank",
                                "target_species",
                                "expected",
                                "decision",
                                "hit_window",
                                "min_hits",
                                "hits_in_window",
                                "accepted_in_window",
                                "rejected_in_window",
                                "top_hit_ranks",
                                "top_hit_species",
                            ]
                        )
                        + "\n"
                    )
                    fh_blast.write(
                        "\t".join(
                            [
                                "family_id",
                                "accession",
                                "decision",
                                "hit_window",
                                "min_hits",
                                "hit_index",
                                "ref_accession",
                                "ref_rank_taxid",
                                "ref_rank_name",
                                "ref_species_taxid",
                                "ref_species_name",
                                "allowed",
                                "pident",
                                "qcovs",
                                "length",
                                "evalue",
                                "bitscore",
                            ]
                        )
                        + "\n"
                    )
                    for row in votes:
                        (
                            family_id,
                            busco_library_id,
                            target_library_id,
                            accession,
                            run_id,
                            expected_taxid,
                            best_taxid,
                            runner_taxid,
                            rank,
                            best_bitscore,
                            delta_bitscore,
                            decision,
                            top_hits_json,
                            *_vote_extra,
                        ) = row

                        top_hits = []
                        try:
                            if top_hits_json:
                                top_hits = json.loads(top_hits_json)
                        except (TypeError, json.JSONDecodeError):
                            top_hits = []
                        top_hit_accs = []
                        top_hit_taxa = []
                        top_hit_bits = []
                        top_hit_allowed_list = []
                        group = self._group_for_accession(accession, None, db=self.db_manager)
                        group_clades = group.get("clades") if group else set()
                        group_blacklist = group.get("blacklist") if group else set()
                        group_min_hits = group.get("min_hits") if group else self.min_hits
                        group_hit_window_explicit = bool(group and "hit_window" in group)
                        group_hit_window = group.get("hit_window") if group_hit_window_explicit else self.hit_window
                        min_hits_req, hit_window = self._resolve_hit_window(
                            group_min_hits,
                            group_hit_window,
                            hit_window_explicit=group_hit_window_explicit,
                            warn=False,
                        )
                        acceptable_ref_count = 0
                        refs_with_family = ref_busco_map.get(family_id, set())
                        if refs_with_family:
                            if group_clades:
                                for ref_acc in refs_with_family:
                                    tax_full = ref_taxa_full.get(ref_acc)
                                    if any(self._is_descendant(tax_full, blk, db=self.db_manager) for blk in group_blacklist):
                                        continue
                                    if group_clades and not any(
                                        self._is_descendant(tax_full, cl, db=self.db_manager) for cl in group_clades
                                    ):
                                        continue
                                    acceptable_ref_count += 1
                            else:
                                for ref_acc in refs_with_family:
                                    if ref_rank_taxa.get(ref_acc) == expected_taxid:
                                        acceptable_ref_count += 1

                        hit_infos = []

                        for h in top_hits:
                            ref_acc = h.get("ref_accession") or h.get("sseqid")
                            if not ref_acc:
                                continue
                            top_hit_accs.append(ref_acc)
                            tid, name, _raw = acc_rank_taxon(ref_acc)
                            sp_tid, sp_name, _raw_sp = acc_rank_taxon(ref_acc, "species")
                            rank_name = name or tid or ""
                            species_name = sp_name or sp_tid or ""
                            bits_val = h.get("bitscore", "")
                            top_hit_taxa.append(rank_name)
                            top_hit_bits.append(str(bits_val))
                            # Re-evaluate allowed flag
                            ref_tax_full = self.db_manager.genomes.get_lineage_root_to_leaf(_raw) or []
                            allowed = False
                            if group_clades:
                                if any(self._is_descendant(_raw, blk, db=self.db_manager) for blk in group_blacklist):
                                    allowed = False
                                elif any(self._is_descendant(_raw, cl, db=self.db_manager) for cl in group_clades):
                                    allowed = True
                            else:
                                allowed = (tid == expected_taxid)
                            top_hit_allowed_list.append("1" if allowed else "0")
                            hit_infos.append(
                                {
                                    "acc": ref_acc,
                                    "rank_taxid": tid,
                                    "rank": rank_name,
                                    "species_taxid": sp_tid,
                                    "species": species_name,
                                    "bitscore": bits_val,
                                    "evalue": h.get("evalue", ""),
                                    "pident": h.get("pident", ""),
                                    "qcovs": h.get("qcovs", ""),
                                    "length": h.get("length", ""),
                                    "allowed": allowed,
                                }
                            )

                        window_hits = hit_infos[:hit_window]
                        hits_in_window = len(window_hits)
                        accepted_in_window = sum(1 for h in window_hits if h["allowed"])
                        rejected_in_window = hits_in_window - accepted_in_window
                        top_ranks = "|".join(str(h["rank"] or "") for h in window_hits)
                        top_species = "|".join(str(h["species"] or "") for h in window_hits)

                        target_tid, target_rank_name, _raw_tid = acc_rank_taxon(accession)
                        _sp_tid, target_species_name, _raw_sp = acc_rank_taxon(accession, "species")
                        target_species = target_species_name or _sp_tid or ""
                        decision_rule = f"rank:{self.rank}"
                        if group_clades:
                            group_idx = ""
                            try:
                                group_idx = str(self.groups.index(group))
                            except ValueError:
                                group_idx = ""
                            if group_idx:
                                decision_rule = f"rank:{self.rank};group:{group_idx}"
                            else:
                                decision_rule = f"rank:{self.rank};group:custom"
                        row_out = [
                            family_id,
                            accession,
                            "single_copy",
                            target_rank_name or target_tid or "",
                            target_species,
                            tax_name(expected_taxid) or expected_taxid or "",
                            decision,
                            hit_window,
                            min_hits_req,
                            hits_in_window,
                            accepted_in_window,
                            rejected_in_window,
                            top_ranks,
                            top_species,
                        ]
                        fh.write("\t".join("" if v is None else str(v) for v in row_out) + "\n")
                        for idx, h in enumerate(window_hits, start=1):
                            row_hit = [
                                family_id,
                                accession,
                                decision,
                                hit_window,
                                min_hits_req,
                                idx,
                                h.get("acc", ""),
                                h.get("rank_taxid") or "",
                                h.get("rank") or "",
                                h.get("species_taxid") or "",
                                h.get("species") or "",
                                "1" if h.get("allowed") else "0",
                                h.get("pident", ""),
                                h.get("qcovs", ""),
                                h.get("length", ""),
                                h.get("evalue", ""),
                                h.get("bitscore", ""),
                            ]
                            fh_blast.write("\t".join("" if v is None else str(v) for v in row_hit) + "\n")

                    # Add rows for BUSCO families with no hits
                    if busco_families_by_acc:
                        for acc, fams in busco_families_by_acc.items():
                            target_tid, target_rank_name, _raw_tid = acc_rank_taxon(acc)
                            _sp_tid, target_species_name, _raw_sp = acc_rank_taxon(acc, "species")
                            target_species = target_species_name or _sp_tid or ""
                            expected_tax = target_rank_name or target_tid or ""
                            group = self._group_for_accession(acc, None, db=self.db_manager)
                            group_min_hits = group.get("min_hits") if group else self.min_hits
                            group_hit_window_explicit = bool(group and "hit_window" in group)
                            group_hit_window = group.get("hit_window") if group_hit_window_explicit else self.hit_window
                            min_hits_req, hit_window = self._resolve_hit_window(
                                group_min_hits,
                                group_hit_window,
                                hit_window_explicit=group_hit_window_explicit,
                                warn=False,
                            )
                            for fam in sorted(fams):
                                if (acc, fam) in vote_keys:
                                    continue
                                row_out = [
                                    fam,
                                    acc,
                                    "no_hit",
                                    target_rank_name or target_tid or "",
                                    target_species,
                                    expected_tax,
                                    "no_hit",
                                    hit_window,
                                    min_hits_req,
                                    0,
                                    0,
                                    0,
                                    "",
                                    "",
                                ]
                                fh.write("\t".join("" if v is None else str(v) for v in row_out) + "\n")

                summaries = self.db_manager.filtering.get_decontamination_summary(run_id=self.run_id)
                with open(summary_report, "w") as fh:
                    fh.write(
                        "\t".join(
                            [
                                "accession",
                                "target_library_id",
                                "busco_library_id",
                                "run_id",
                                "busco_run_id",
                                "expected_taxid",
                                "majority_taxid",
                                "rank",
                                "buscos_tested",
                                "buscos_supporting",
                                "buscos_outside",
                                "buscos_unknown",
                                "single_copy_total",
                                "single_copy_missing_files",
                                "single_copy_no_hits",
                                "single_copy_after_hidden",
                                "off_clade_fraction",
                                "decision",
                                "params_json",
                            ]
                        )
                        + "\n"
                    )
                    accessions = [row[0] for row in summaries]
                    total_map = {}
                    sc_map = {}
                    have_file_map = {}
                    hidden_counts = {}
                    if accessions:
                        def _chunk(seq, size):
                            for i in range(0, len(seq), size):
                                yield seq[i:i + size]

                        for chunk in _chunk(accessions, 800):
                            counts = self.db_manager.busco.get_family_counts_by_accession(
                                library_id=self.busco_lib_id,
                                accessions=list(chunk),
                                status=[1],
                            )
                            for acc, cnt in counts.items():
                                total_map[str(acc)] = int(cnt or 0)
                                sc_map[str(acc)] = int(cnt or 0)

                        for chunk in _chunk(accessions, 300):
                            have_counts = self.db_manager.busco.count_existing_family_locations_by_accession(
                                library_id=self.busco_lib_id,
                                accessions=list(chunk),
                                status=[1],
                            )
                            for acc, cnt in have_counts.items():
                                have_file_map[str(acc)] = have_file_map.get(str(acc), 0) + int(cnt or 0)

                        hidden_counts = self.db_manager.filtering.paralog_hidden_counts(
                            target_library_id=self.library_id,
                            busco_library_id=self.busco_lib_id,
                            accessions=accessions,
                        )

                    for row in summaries:
                        unknown_ct = 0
                        (
                            accession,
                            target_library_id,
                            busco_library_id,
                            run_id,
                            expected_taxid,
                            majority_taxid,
                            rank,
                            buscos_tested,
                            buscos_supporting,
                            buscos_outside,
                            off_clade_fraction,
                            decision,
                            params_json,
                            *summary_extra,
                        ) = row
                        busco_run_id = summary_extra[0] if summary_extra else None
                        acc = accession
                        decs = decisions_by_target.get(acc, [])
                        unknown_ct = sum(1 for d in decs if str(d).lower() == "unknown")
                        total_cnt = total_map.get(acc, 0)
                        have_file_cnt = have_file_map.get(acc, 0)
                        missing_cnt = max(total_cnt - have_file_cnt, 0)
                        tested_cnt = int(buscos_tested or 0)
                        no_hits_cnt = max(total_cnt - missing_cnt - tested_cnt, 0)
                        sc_after_hidden = ""
                        sc_val = sc_map.get(acc)
                        hidden_pair = hidden_counts.get(acc) if hidden_counts else None
                        if sc_val is not None and hidden_pair:
                            hidden_cnt, total_hidden = hidden_pair
                            if int(total_hidden or 0) > 0:
                                sc_after_hidden = str(max(int(sc_val) - int(hidden_cnt or 0), 0))
                        row_out = [
                            accession,
                            target_library_id,
                            busco_library_id,
                            run_id,
                            busco_run_id if busco_run_id is not None else "",
                            expected_taxid,
                            majority_taxid,
                            rank,
                            buscos_tested,
                            buscos_supporting,
                            buscos_outside,
                            unknown_ct,
                            total_cnt,
                            missing_cnt,
                            no_hits_cnt,
                            sc_after_hidden,
                            off_clade_fraction,
                            decision,
                            params_json,
                        ]
                        fh.write("\t".join("" if v is None else str(v) for v in row_out) + "\n")

                # Reference composition report
                with open(refs_report, "w") as fh:
                    fh.write("\t".join(["accession", "taxid", f"{self.rank}_taxid", f"{self.rank}_name", "species_name", "species_taxid", "selected_by"]) + "\n")
                    for ref in self.ref_accessions:
                        tid_at_rank, rank_name, raw_taxid = acc_rank_taxon(ref)
                        species_tid = self._taxon_at_rank_name(raw_taxid, "species")
                        selected_by = "explicit" if ref in getattr(self, "_report_refs_explicit", set()) else ("selector" if ref in getattr(self, "_report_refs_selected", set()) else "")
                        fh.write(
                            "\t".join(
                                [
                                    ref,
                                    str(raw_taxid or ""),
                                    str(tid_at_rank or ""),
                                    str(rank_name or tid_at_rank or ""),
                                    tax_name(species_tid) or "",
                                    str(species_tid or ""),
                                    selected_by,
                                ]
                            )
                            + "\n"
                        )

                # Grouping report (acceptance rules)
                with open(groups_report, "w") as fh:
                    fh.write("\t".join(["group_index", "members", "member_taxa", "allowed_clades", "blacklist", "min_hits", "applied_targets"]) + "\n")
                    target_group_map = {}
                    for acc in self.accessions:
                        g = self._group_for_accession(acc, None, db=self.db_manager)
                        if g:
                            idx = self.groups.index(g)
                            target_group_map.setdefault(idx, []).append(acc)
                    for idx, g in enumerate(self.groups):
                        fh.write(
                            "\t".join(
                                [
                                    str(idx),
                                    ",".join(sorted(g.get("member_accessions") or [])),
                                    ",".join(str(t) for t in sorted(g.get("member_taxa") or [])),
                                    ",".join(str(t) for t in sorted(g.get("clades") or [])),
                                    ",".join(str(t) for t in sorted(g.get("blacklist") or [])),
                                    str(g.get("min_hits")),
                                    ",".join(sorted(target_group_map.get(idx, []))),
                                ]
                            )
                            + "\n"
                        )

                self.log(
                    f"Wrote decontamination reports under {os.path.dirname(busco_report) or '.'}.",
                    "INFO",
                )
            except Exception as e:  # boundary: optional report generation failure should not invalidate persisted analysis.
                self.log(f"Failed to write report to {self.report_path}: {e}", "WARNING")

        # Persist the active decontamination run for this library so other tasks can
        # resolve it without manual run id lookups.
        if self.run_id and self.library_id:
            env_key = f"ACTIVE_DECONT_RUN_{self.library_id}"
            payload = {
                "run_id": self.run_id,
                "run_label": self.run_label,
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
            if self.db_manager.env.set(env_key, payload):
                self.log(
                    f"Recorded active decontamination run for library {self.library_id}: {payload}",
                    "DEBUG",
                )
            else:
                self.log(
                    f"Failed to record active decontamination run env var {env_key}",
                    "WARNING",
                )

        return True
