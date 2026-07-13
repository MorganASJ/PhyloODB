import concurrent.futures
import hashlib
import json
import math
import os
import tempfile
import threading
import subprocess
from datetime import datetime

from .decontamination import Decontamination
from ..reporting import resolve_report_base_path
from ...database import DBManager


class InternalDecontaminationTask(Decontamination):
    """Internal-consistency decontamination using BUSCO sequences within the target set."""
    INTERNAL_CACHE_KEY_FILE = "__idc_cache_key.json"

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=16):
        super().__init__(db_path, task_id, checkpoint, data, required_threads=required_threads)
        if not self.data.get("report_path"):
            self.report_path = str(
                resolve_report_base_path(
                    self,
                    namespace="internal-decontamination-reports",
                    default_stem="internal_decontamination",
                    run_label=self.run_label or self.library_name or self.library_id,
                    cache_attr="_internal_decontamination_report_dir",
                )
            )
        self.p_value_threshold = float(self.data.get("p_value_threshold", 0.05))
        self.off_clade_fraction = float(self.data.get("off_clade_fraction", 0.05))
        try:
            self.min_alignment_length = max(0, int(self.data.get("min_alignment_length", 0) or 0))
        except (TypeError, ValueError):
            self.min_alignment_length = 0
        self.max_target_seqs = self.data.get("max_target_seqs")
        self.blast_program = self.data.get("blast_program")
        self.blast_db_type = self.data.get("blast_db_type")
        self.use_paralog_filtered_buscos = bool(self.data.get("use_paralog_filtered_buscos", False))
        self.internal_blastdb_path = self.data.get("_internal_blastdb_path")
        self.internal_id_map_path = self.data.get("_internal_blastdb_id_map")
        self.save_blast_output = self.data.get("save_blast_output")
        self.reuse_blast_results = self.data.get("reuse_blast_results")
        self.external_blast_db_path = self.data.get("external_blast_db_path") or self.data.get("external_blastdb_path")
        self.external_blast_db_type = self.data.get("external_blast_db_type")
        self.external_blast_program = self.data.get("external_blast_program")
        self.external_blast_output_dir = self.data.get("external_blast_output_dir")
        self.external_reuse_blast_results = self.data.get("external_reuse_blast_results")
        self.external_max_target_seqs = self.data.get("external_max_target_seqs")
        if not any(k in self.data for k in ("hit_window", "window", "window_size")):
            self.hit_window = 8
            self.data["hit_window"] = self.hit_window

    def _load_config(self):
        res = super()._load_config()
        if res == "ERROR":
            return res
        if not self.config_path or not os.path.exists(self.config_path):
            return res
        try:
            with open(self.config_path, "r") as handle:
                cfg = json.loads(handle.read())
        except (OSError, UnicodeError, json.JSONDecodeError):
            return res
        params = {}
        if isinstance(cfg, dict):
            params = cfg.get("params") or {}
        for key, attr in {
            "p_value_threshold": "p_value_threshold",
            "min_alignment_length": "min_alignment_length",
            "max_target_seqs": "max_target_seqs",
            "blast_program": "blast_program",
            "blast_db_type": "blast_db_type",
            "save_blast_output": "save_blast_output",
            "reuse_blast_results": "reuse_blast_results",
            "use_paralog_filtered_buscos": "use_paralog_filtered_buscos",
            "external_blast_db_path": "external_blast_db_path",
            "external_blast_db_type": "external_blast_db_type",
            "external_blast_program": "external_blast_program",
            "external_blast_output_dir": "external_blast_output_dir",
            "external_reuse_blast_results": "external_reuse_blast_results",
            "external_max_target_seqs": "external_max_target_seqs",
        }.items():
            if key in params:
                setattr(self, attr, params[key])
        return res

    @staticmethod
    def _log_comb(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    @classmethod
    def _hypergeom_sf(cls, x_minus_1: int, N: int, K: int, n: int) -> float:
        """Survival function P(X >= x) for Hypergeometric(N, K, n) with x = x_minus_1 + 1."""
        if N <= 0 or n <= 0 or K <= 0:
            return 1.0
        x = x_minus_1 + 1
        max_x = min(K, n)
        if x <= 0:
            return 1.0
        if x > max_x:
            return 0.0
        log_denom = cls._log_comb(N, n)
        total = 0.0
        for k in range(x, max_x + 1):
            log_p = cls._log_comb(K, k) + cls._log_comb(N - K, n - k) - log_denom
            total += math.exp(log_p)
        return min(max(total, 0.0), 1.0)

    @classmethod
    def _hypergeom_decision_feasibility(cls, Z: int, K: int, Y: int, p_value_threshold: float) -> dict:
        """Return whether the hypergeometric rule can both support and reject the expected taxon."""
        Z_val = max(0, int(Z))
        K_val = max(0, int(K))
        Y_val = max(0, int(Y))
        threshold = float(p_value_threshold)
        max_x = min(K_val, Y_val)
        min_x = max(0, Y_val - (Z_val - K_val))
        p_best = cls._hypergeom_sf(max_x - 1, Z_val, K_val, Y_val)
        p_worst = cls._hypergeom_sf(min_x - 1, Z_val, K_val, Y_val)
        return {
            "Z": Z_val,
            "K": K_val,
            "Y": Y_val,
            "min_x": min_x,
            "max_x": max_x,
            "p_best": p_best,
            "p_worst": p_worst,
            "win_possible": p_best < threshold,
            "lose_possible": p_worst >= threshold,
        }

    def _detect_busco_db_type(self) -> str:
        if self.blast_db_type in ("prot", "nucl"):
            return self.blast_db_type
        for acc in self.accessions:
            rows = self.db_manager.busco.get_family_results_for_library(
                library_id=self.busco_lib_id,
                accessions=[acc],
                status=[1],
            )
            for family_id, lib_id, accession, _status, _sequence, _score, _length in rows or []:
                loc = self.db_manager.busco.get_family_location(family_id, lib_id, accession)
                if not loc:
                    continue
                if loc.endswith(".fna") or loc.endswith(".fna.gz"):
                    return "nucl"
                return "prot"
        return "prot"

    def _resolve_blast_program(self, db_type: str) -> str | None:
        if self.blast_program:
            token = str(self.blast_program).lower()
            if token == "blastn":
                return self.db_manager.env.get("BLASTN_PATH")
            if token == "blastp":
                return self.db_manager.env.get("BLASTP_PATH")
            return self.blast_program
        if db_type == "nucl":
            return self.db_manager.env.get("BLASTN_PATH")
        return self.db_manager.env.get("BLASTP_PATH")

    def _run_blast(self, query_faa: str, db_path: str) -> str | bool:
        if not self.blastp_path:
            self.error("BLAST program not set before attempting BLAST.")
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
        command.extend(["-max_target_seqs", str(self._effective_max_target_seqs())])
        threads = max(1, int(self.blast_threads)) if getattr(self, "blast_threads", None) else 1
        if threads > 1:
            command.extend(["-num_threads", str(threads)])
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            self.error(f"BLAST failed: {result.stderr}")
            return False
        return result.stdout or ""

    def _blast_cache_path(self, base_dir: str, accession: str, family_id: str) -> str:
        safe_acc = str(accession).replace(os.sep, "_")
        safe_fam = str(family_id).replace(os.sep, "_")
        return os.path.join(base_dir, safe_acc, f"{safe_fam}.blast6")

    def _cache_key_path(self, base_dir: str) -> str:
        return os.path.join(base_dir, self.INTERNAL_CACHE_KEY_FILE)

    def _sha1_file(self, path: str) -> str | None:
        if not path or not os.path.exists(path):
            return None
        h = hashlib.sha1()
        try:
            with open(path, "rb") as fh:
                while True:
                    chunk = fh.read(1024 * 1024)
                    if not chunk:
                        break
                    h.update(chunk)
        except OSError:
            return None
        return h.hexdigest()

    def _effective_max_target_seqs(self) -> int:
        max_targets = None
        if self.max_target_seqs not in (None, ""):
            try:
                max_targets = max(1, int(self.max_target_seqs))
            except (TypeError, ValueError):
                max_targets = None
        if max_targets is None:
            return max(1, int(self.hit_window) + 1)
        return max(max_targets, int(self.hit_window) + 1)

    def _current_internal_cache_key(self) -> dict:
        id_map_sha1 = self._sha1_file(self.internal_id_map_path)
        return {
            "schema_version": 1,
            "kind": "internal_decontamination_blast_cache",
            "id_map_sha1": id_map_sha1,
            "blast_program": self.blastp_path,
            "db_type": self.blast_db_type,
            "outfmt": "6 qseqid sseqid pident length qcovs evalue bitscore",
            "effective_max_target_seqs": self._effective_max_target_seqs(),
        }

    def _load_cache_key(self, base_dir: str) -> dict | None:
        path = self._cache_key_path(base_dir)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "r") as fh:
                payload = json.load(fh)
            return payload if isinstance(payload, dict) else None
        except (OSError, UnicodeError, json.JSONDecodeError, TypeError):
            return None

    def _write_cache_key(self, base_dir: str) -> None:
        if not base_dir:
            return
        try:
            os.makedirs(base_dir, exist_ok=True)
            with open(self._cache_key_path(base_dir), "w") as fh:
                json.dump(self._current_internal_cache_key(), fh, sort_keys=True)
        except (OSError, TypeError, ValueError) as exc:
            self.log(f"Internal decontam: failed to write cache key under {base_dir}: {exc}", "WARNING")

    def _validate_or_disable_internal_cache_reuse(self) -> None:
        if not self.reuse_blast_results:
            return
        expected = self._current_internal_cache_key()
        found = self._load_cache_key(self.reuse_blast_results)
        if not found:
            self.log(
                f"Internal decontam: cache key missing at {self._cache_key_path(self.reuse_blast_results)}; "
                "disabling cache reuse for this run.",
                "WARNING",
            )
            self.reuse_blast_results = None
            return
        if found != expected:
            self.log(
                f"Internal decontam: cache key mismatch at {self._cache_key_path(self.reuse_blast_results)}; "
                "disabling cache reuse for this run.",
                "WARNING",
            )
            self.reuse_blast_results = None
            return
        self.log("Internal decontam: cache key validated.", "DEBUG")

    def _load_blast_output(self, path: str) -> str | None:
        if not path or not os.path.exists(path):
            return None
        try:
            with open(path, "r") as fh:
                return fh.read()
        except (OSError, UnicodeError):
            return None

    def _save_blast_output(self, path: str, data: str) -> None:
        if not path:
            return
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write(data or "")
        except OSError as exc:
            self.log(f"Internal decontam: failed to save BLAST cache output {path}: {exc}", "WARNING")

    def _load_id_map(self, path: str) -> dict:
        mapping = {}
        if not path or not os.path.exists(path):
            return mapping
        with open(path, "r") as fh:
            header = True
            for line in fh:
                if header:
                    header = False
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                sseqid, accession, family_id, taxid = parts[:4]
                taxid_val = None
                try:
                    if taxid not in ("", "NA", "None", None):
                        taxid_val = int(taxid)
                except (TypeError, ValueError):
                    taxid_val = None
                mapping[sseqid] = {
                    "accession": accession,
                    "family_id": family_id,
                    "taxid": taxid_val,
                    "orig_header": parts[4] if len(parts) > 4 else None,
                }
        return mapping

    def _external_check_enabled(self) -> bool:
        return bool(self.external_blast_db_path or self.external_blast_output_dir or self.external_reuse_blast_results)

    def _extract_taxid_from_staxids(self, staxids: str | None) -> int | None:
        if not staxids:
            return None
        token = str(staxids).strip()
        if not token or token.lower() in {"na", "none"}:
            return None
        for part in token.replace("|", ";").replace(",", ";").split(";"):
            part = part.strip()
            if not part or part.lower() in {"na", "n/a", "none", "0"}:
                continue
            try:
                return int(part)
            except (TypeError, ValueError):
                continue
        return None

    def _parse_blast_output_with_taxids(self, blast_output: str):
        hits = []
        for line in (blast_output or "").strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
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
            if len(parts) > 7:
                hit["hit_taxid"] = self._extract_taxid_from_staxids(parts[7])
            else:
                hit["hit_taxid"] = None
            hits.append(hit)
        return hits

    def _external_context_for_accession(self, accession: str, db=None) -> dict:
        dbm = db or self.db_manager
        genome = dbm.get_genome(accession)
        taxid = genome[1] if genome else None
        expected_taxon = self._taxon_at_rank(taxid, db=dbm)
        group = self._group_for_accession(accession, taxid, db=dbm)
        group_hit_window = self.hit_window
        if group and "hit_window" in group:
            try:
                group_hit_window = max(1, int(group.get("hit_window")))
            except (TypeError, ValueError):
                group_hit_window = self.hit_window
        min_hits_req = getattr(self, "min_hits", 1)
        if group and "min_hits" in group:
            try:
                min_hits_req = max(1, int(group.get("min_hits")))
            except (TypeError, ValueError):
                min_hits_req = getattr(self, "min_hits", 1)
        return {
            "taxid": taxid,
            "expected_taxon": expected_taxon,
            "group": group,
            "group_hit_window": group_hit_window,
            "group_clades": group.get("clades") if group else set(),
            "group_blacklist": group.get("blacklist") if group else set(),
            "min_hits_req": min_hits_req,
        }

    def _evaluate_external_hits(
        self,
        ext_hits: list[dict],
        *,
        expected_taxon: int | None,
        group_clades: set,
        group_blacklist: set,
        group_hit_window: int,
        min_hits_req: int,
        db=None,
    ) -> dict:
        # Apply thresholds
        filtered = []
        for hit in ext_hits or []:
            if self.min_identity and hit["pident"] < self.min_identity:
                continue
            if self.min_coverage and hit["qcovs"] < self.min_coverage:
                continue
            if self.min_alignment_length and hit["length"] < self.min_alignment_length:
                continue
            if self.min_bitscore and hit["bitscore"] < self.min_bitscore:
                continue
            if self.max_evalue is not None and hit["evalue"] > self.max_evalue:
                continue
            filtered.append(hit)
        filtered.sort(key=lambda h: (-h["bitscore"], h["evalue"]))
        # Deduplicate by subject id
        deduped = []
        seen_ids = set()
        for hit in filtered:
            sid = hit.get("sseqid")
            if sid in seen_ids:
                continue
            seen_ids.add(sid)
            deduped.append(hit)
        filtered = deduped

        # Remove likely self-hit before evaluating external support.
        # Priority:
        # 1) exact qseqid/sseqid match when present
        # 2) top hit with 100% identity (user-requested safety rule)
        removed_self = False
        for idx, hit in enumerate(filtered):
            if hit.get("sseqid") and hit.get("qseqid") and hit["sseqid"] == hit["qseqid"]:
                filtered.pop(idx)
                removed_self = True
                break
        if not removed_self and filtered:
            try:
                if float(filtered[0].get("pident", 0.0)) >= 100.0:
                    filtered.pop(0)
                    removed_self = True
            except (TypeError, ValueError):
                pass
        if removed_self:
            self.log("External check: removed likely self-hit from top external BLAST hits.", "DEBUG")

        if not filtered:
            return {
                "decision": "unknown",
                "hits": [],
                "allowed": [],
                "top_taxon": None,
                "top_acc": None,
                "top_bits": None,
                "top_eval": None,
                "best_rank_taxon": None,
                "runner_rank_taxon": None,
                "hit_window": 0,
            }

        window_size = min(max(1, int(group_hit_window)), len(filtered))
        window_hits = filtered[:window_size]

        def _allowed_external(hit_taxid):
            if not hit_taxid:
                return None
            hit_rank_taxid = self._taxon_at_rank(hit_taxid, db=db)
            if group_clades:
                if any(self._is_descendant(hit_taxid, blk, db=db) for blk in group_blacklist):
                    return False
                return any(self._is_descendant(hit_taxid, cl, db=db) for cl in group_clades)
            return hit_rank_taxid == expected_taxon

        all_allowed_flags = [_allowed_external(h.get("hit_taxid")) for h in filtered]
        feasibility = self._decision_window_feasibility(
            sum(1 for ok in all_allowed_flags if ok),
            sum(1 for ok in all_allowed_flags if ok is False),
            min_hits=min_hits_req,
            hit_window=group_hit_window,
        )

        allowed_flags = [_allowed_external(h.get("hit_taxid")) for h in window_hits]
        allowed_count = sum(1 for ok in allowed_flags if ok)
        outside_count = sum(1 for ok in allowed_flags if ok is False)
        unknown_count = sum(1 for ok in allowed_flags if ok is None)

        min_hits_req = max(1, int(min_hits_req)) if min_hits_req else 1
        if min_hits_req > window_size:
            min_hits_req = window_size

        if not feasibility["win_possible"] or not feasibility["lose_possible"]:
            decision = "unknown"
        elif allowed_count >= min_hits_req:
            decision = "support"
        elif allowed_count == 0 and outside_count > 0:
            decision = "outside"
        elif outside_count == 0 and unknown_count > 0:
            decision = "unknown"
        else:
            decision = "outside"

        best_allowed = next((h for h, ok in zip(window_hits, allowed_flags) if ok), None)
        runner = next((h for h, ok in zip(window_hits, allowed_flags) if ok is False), None)

        top_hit = None
        if decision == "outside":
            top_hit = runner
        if top_hit is None and decision == "support":
            top_hit = best_allowed
        if top_hit is None:
            top_hit = window_hits[0] if window_hits else None

        top_taxon = top_hit.get("hit_taxid") if top_hit else None
        top_acc = top_hit.get("sseqid") if top_hit else None
        top_bits = top_hit.get("bitscore") if top_hit else None
        top_eval = top_hit.get("evalue") if top_hit else None

        best_rank_taxon = None
        runner_rank_taxon = None
        if best_allowed and best_allowed.get("hit_taxid"):
            best_rank_taxon = self._taxon_at_rank(best_allowed.get("hit_taxid"), db=db)
        if runner and runner.get("hit_taxid"):
            runner_rank_taxon = self._taxon_at_rank(runner.get("hit_taxid"), db=db)

        return {
            "decision": decision,
            "hits": window_hits,
            "allowed": allowed_flags,
            "top_taxon": top_taxon,
            "top_acc": top_acc,
            "top_bits": top_bits,
            "top_eval": top_eval,
            "best_rank_taxon": best_rank_taxon,
            "runner_rank_taxon": runner_rank_taxon,
            "hit_window": window_size,
        }

    def _external_pending_votes(self) -> list[tuple[str, str]]:
        votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)
        pending = []
        for row in votes or []:
            family_id = row[0]
            accession = row[3]
            decision = row[11]
            payload = {}
            try:
                payload = json.loads(row[12]) if row[12] else {}
            except (TypeError, json.JSONDecodeError):
                payload = {}
            internal_decision = payload.get("internal_decision") or decision
            if str(internal_decision).lower() != "outside":
                continue
            if payload.get("external_check_run"):
                continue
            if payload.get("external_check_pending") is False:
                continue
            pending.append((accession, family_id))
        return pending

    def _external_subtasks_complete(self) -> bool:
        gen = int(self.data.get("_stage_4_gen", 0)) or 1
        tasks = self.db_manager.tasks.get_subtasks(self.task_id) or []
        phase_tasks = []
        for t in tasks:
            try:
                d = json.loads(t[6]) if t[6] else {}
            except (TypeError, json.JSONDecodeError):
                d = {}
            if d.get("__stage") == 4 and d.get("__gen") == gen:
                phase_tasks.append(t)
        if not phase_tasks:
            return True
        statuses = [t[2] for t in phase_tasks]
        return all(s == "C" for s in statuses)

    def _build_params_json(self) -> str:
        return json.dumps(
            {
                "rank": self.rank,
                "hit_window": self.hit_window,
                "p_value_threshold": self.p_value_threshold,
                "min_buscos": self.min_buscos,
                "min_identity": self.min_identity,
                "min_coverage": self.min_coverage,
                "min_alignment_length": self.min_alignment_length,
                "min_bitscore": self.min_bitscore,
                "max_evalue": self.max_evalue,
                "max_target_seqs": self.max_target_seqs,
                "save_blast_output": self.save_blast_output,
                "reuse_blast_results": self.reuse_blast_results,
                "external_blast_db_path": self.external_blast_db_path,
                "external_blast_db_type": self.external_blast_db_type,
                "external_blast_program": self.external_blast_program,
                "external_blast_output_dir": self.external_blast_output_dir,
                "external_reuse_blast_results": self.external_reuse_blast_results,
                "external_max_target_seqs": self.external_max_target_seqs,
                "config_path": self.config_path,
                "config_signature": self.config_signature,
                "run_id": self.run_id,
                "run_label": self.run_label,
            }
        )

    def _apply_external_results(self) -> bool:
        ext_base = self.external_reuse_blast_results or self.external_blast_output_dir
        if not ext_base:
            self.log("External check: no external BLAST outputs configured; skipping apply.", "WARNING")
            return False

        votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)
        if not votes:
            return False

        # Precompute context per accession
        contexts = {}
        for row in votes:
            acc = row[3]
            if acc not in contexts:
                contexts[acc] = self._external_context_for_accession(acc, db=self.db_manager)

        updated = 0
        for row in votes:
            (
                family_id,
                _busco_lib_id,
                _target_lib_id,
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
                busco_run_id,
            ) = row
            payload = {}
            try:
                payload = json.loads(top_hits_json) if top_hits_json else {}
            except (TypeError, json.JSONDecodeError):
                payload = {}
            internal_decision = payload.get("internal_decision") or decision
            if str(internal_decision).lower() != "outside":
                continue

            ctx = contexts.get(accession) or {}
            ext_path = self._blast_cache_path(ext_base, accession, family_id)
            ext_out = self._load_blast_output(ext_path)
            external_check_run = False
            external_check_pending = False
            external_decision = None
            external_hits = []
            external_allowed = []
            external_top_taxon = None
            external_top_acc = None
            external_top_bits = None
            external_top_eval = None
            external_hit_window = None
            new_decision = decision
            new_best_taxon = best_taxid
            new_runner_taxon = runner_taxid
            new_best_bits = best_bitscore

            if ext_out is None:
                external_check_pending = True
            else:
                external_check_run = True
                ext_hits = self._parse_blast_output_with_taxids(ext_out)
                result = self._evaluate_external_hits(
                    ext_hits,
                    expected_taxon=ctx.get("expected_taxon"),
                    group_clades=ctx.get("group_clades") or set(),
                    group_blacklist=ctx.get("group_blacklist") or set(),
                    group_hit_window=ctx.get("group_hit_window") or self.hit_window,
                    min_hits_req=ctx.get("min_hits_req") or 1,
                    db=self.db_manager,
                )
                external_decision = result.get("decision")
                external_hits = result.get("hits") or []
                external_allowed = result.get("allowed") or []
                external_top_taxon = result.get("top_taxon")
                external_top_acc = result.get("top_acc")
                external_top_bits = result.get("top_bits")
                external_top_eval = result.get("top_eval")
                external_hit_window = result.get("hit_window")

                if external_decision == "support":
                    new_decision = "support"
                    new_best_taxon = result.get("best_rank_taxon") or ctx.get("expected_taxon") or new_best_taxon
                    if result.get("runner_rank_taxon"):
                        new_runner_taxon = result.get("runner_rank_taxon")
                    if external_top_bits is not None:
                        new_best_bits = external_top_bits

            payload.update(
                {
                    "external_check_run": external_check_run,
                    "external_check_pending": external_check_pending,
                    "external_decision": external_decision,
                    "external_hit_window": external_hit_window,
                    "external_hits": external_hits,
                    "external_allowed": external_allowed,
                    "external_top_hit_taxon": external_top_taxon,
                    "external_top_hit_accession": external_top_acc,
                    "external_top_hit_bitscore": external_top_bits,
                    "external_top_hit_evalue": external_top_eval,
                }
            )

            self.db_manager.filtering.add_decontamination_vote(
                family_id,
                self.busco_lib_id,
                self.library_id,
                accession,
                run_id,
                expected_taxid,
                new_best_taxon,
                new_runner_taxon,
                rank,
                new_best_bits,
                delta_bitscore,
                new_decision,
                json.dumps(payload),
                busco_run_id=busco_run_id,
            )
            if new_decision != decision:
                updated += 1

        if updated:
            self.log(f"External check: updated {updated} BUSCO votes after external BLAST.", "INFO")

        # Recompute summaries with updated decisions
        votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)
        acc_votes = {}
        for row in votes or []:
            acc_votes.setdefault((row[3], row[13]), []).append(row)

        params_json = self._build_params_json()
        for (acc, busco_run_id), rows in acc_votes.items():
            ctx = contexts.get(acc) or {}
            expected_taxon = ctx.get("expected_taxon")
            best_taxon_counts = {}
            buscos_supporting = 0
            buscos_outside = 0
            tested = 0
            for row in rows:
                decision = row[11]
                if decision == "support":
                    buscos_supporting += 1
                    tested += 1
                elif decision == "outside":
                    buscos_outside += 1
                    tested += 1
                best_taxon = row[6]
                if best_taxon:
                    best_taxon_counts[best_taxon] = best_taxon_counts.get(best_taxon, 0) + 1

            majority_taxon = None
            if best_taxon_counts:
                majority_taxon = max(best_taxon_counts.items(), key=lambda kv: kv[1])[0]

            off_frac = (buscos_outside / tested) if tested else None
            if tested < self.min_buscos:
                final_decision = "UNCERTAIN"
            elif buscos_outside == 0:
                final_decision = "CLEAN"
            elif off_frac is not None and off_frac >= self.off_clade_fraction:
                final_decision = "CONTAMINATED"
            else:
                final_decision = "CLEAN"

            self.db_manager.filtering.add_decontamination_summary(
                acc,
                self.library_id,
                self.busco_lib_id,
                self.run_id,
                expected_taxon,
                majority_taxon,
                self.rank,
                tested,
                buscos_supporting,
                buscos_outside,
                off_frac,
                final_decision,
                params_json=params_json,
                busco_run_id=busco_run_id,
            )

        return True

    def _write_external_report(self, acc_taxa: dict, acc_rank_taxa: dict) -> None:
        if not self.report_path:
            return
        base, ext = os.path.splitext(self.report_path)
        external_report = f"{base}_external.tsv" if ext else f"{self.report_path}_external.tsv"
        external_hits_report = f"{base}_external_hits.tsv" if ext else f"{self.report_path}_external_hits.tsv"
        os.makedirs(os.path.dirname(external_report) or ".", exist_ok=True)

        votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)

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

        with open(external_report, "w") as fh, open(external_hits_report, "w") as fh_hits:
            fh.write(
                "\t".join(
                    [
                        "family_id",
                        "accession",
                        "query_taxon_id",
                        "query_label",
                        "internal_call",
                        "external_call",
                        "final_call",
                        "external_check_run",
                        "external_top_hit_taxon",
                        "external_top_hit_rank_taxon",
                        "external_top_hit_rank_name",
                        "external_top_hit_bitscore",
                        "external_top_hit_evalue",
                        "hits_in_window",
                        "allowed_in_window",
                        "outside_in_window",
                        "unknown_in_window",
                    ]
                )
                + "\n"
            )
            fh_hits.write(
                "\t".join(
                    [
                        "family_id",
                        "accession",
                        "decision",
                        "hit_index",
                        "ref_taxid",
                        "ref_rank_taxid",
                        "ref_rank_name",
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

            for row in votes or []:
                family_id, _busco, _tlib, accession, _rid, expected_taxid, _best, _runner, _rank, _bits, _delta, decision, top_hits_json = row
                payload = {}
                try:
                    payload = json.loads(top_hits_json) if top_hits_json else {}
                except (TypeError, json.JSONDecodeError):
                    payload = {}
                internal_decision = payload.get("internal_decision") or decision
                if str(internal_decision).lower() != "outside":
                    continue

                external_decision = payload.get("external_decision")
                external_check_run = bool(payload.get("external_check_run"))
                external_hits = payload.get("external_hits") or []
                external_allowed = payload.get("external_allowed") or []

                external_top_taxon = payload.get("external_top_hit_taxon")
                external_top_rank = self._taxon_at_rank(external_top_taxon) if external_top_taxon else None
                external_top_rank_name = tax_name(external_top_rank) or external_top_rank or ""

                hits_in_window = len(external_hits)
                allowed_in_window = sum(1 for ok in external_allowed if ok)
                outside_in_window = sum(1 for ok in external_allowed if ok is False)
                unknown_in_window = sum(1 for ok in external_allowed if ok is None)

                row_out = [
                    family_id,
                    accession,
                    expected_taxid or "",
                    tax_name(expected_taxid) or expected_taxid or "",
                    internal_decision,
                    external_decision or "",
                    decision,
                    "1" if external_check_run else "0",
                    external_top_taxon or "",
                    external_top_rank or "",
                    external_top_rank_name or "",
                    payload.get("external_top_hit_bitscore") or "",
                    payload.get("external_top_hit_evalue") or "",
                    hits_in_window,
                    allowed_in_window,
                    outside_in_window,
                    unknown_in_window,
                ]
                fh.write("\t".join("" if v is None else str(v) for v in row_out) + "\n")

                for idx, h in enumerate(external_hits, start=1):
                    allowed = external_allowed[idx - 1] if idx - 1 < len(external_allowed) else None
                    ref_taxid = h.get("hit_taxid")
                    ref_rank_taxid = self._taxon_at_rank(ref_taxid) if ref_taxid else None
                    ref_rank_name = tax_name(ref_rank_taxid) or ref_rank_taxid or ""
                    row_hit = [
                        family_id,
                        accession,
                        decision,
                        idx,
                        ref_taxid or "",
                        ref_rank_taxid or "",
                        ref_rank_name,
                        "" if allowed is None else ("1" if allowed else "0"),
                        h.get("pident", ""),
                        h.get("qcovs", ""),
                        h.get("length", ""),
                        h.get("evalue", ""),
                        h.get("bitscore", ""),
                    ]
                    fh_hits.write("\t".join("" if v is None else str(v) for v in row_hit) + "\n")

        self.log(
            f"Wrote external decontamination report: {external_report}, {external_hits_report}",
            "INFO",
        )

    def _write_external_summary(self, summaries) -> None:
        if not self.report_path or not self._external_check_enabled():
            return
        base, ext = os.path.splitext(self.report_path)
        external_summary = f"{base}_external_summary.tsv" if ext else f"{self.report_path}_external_summary.tsv"
        with open(external_summary, "w") as fh:
            fh.write(
                "\t".join(
                    [
                        "accession",
                        "buscos_tested",
                        "first_pass_support",
                        "first_pass_outside",
                        "second_pass_support",
                        "second_pass_outside",
                        "final_support",
                        "final_outside",
                        "final_outside_fraction",
                        "assembly_decision",
                    ]
                )
                + "\n"
            )
            summary_decisions = {row[0]: row[12] if len(row) > 13 else row[11] for row in summaries or []}
            votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)
            per_acc = {}
            for row in votes:
                acc = row[3]
                decision = row[11]
                payload = {}
                try:
                    payload = json.loads(row[12]) if row[12] else {}
                except (TypeError, json.JSONDecodeError):
                    payload = {}
                internal_decision = payload.get("internal_decision")
                external_decision = payload.get("external_decision")
                per_acc.setdefault(
                    acc,
                    {
                        "tested": 0,
                        "first_support": 0,
                        "first_outside": 0,
                        "second_support": 0,
                        "second_outside": 0,
                        "final_support": 0,
                        "final_outside": 0,
                    },
                )
                if decision in ("support", "outside"):
                    per_acc[acc]["tested"] += 1
                if internal_decision == "support":
                    per_acc[acc]["first_support"] += 1
                elif internal_decision == "outside":
                    per_acc[acc]["first_outside"] += 1
                if external_decision == "support":
                    per_acc[acc]["second_support"] += 1
                elif external_decision == "outside":
                    per_acc[acc]["second_outside"] += 1
                if decision == "support":
                    per_acc[acc]["final_support"] += 1
                elif decision == "outside":
                    per_acc[acc]["final_outside"] += 1

            for acc in self.accessions:
                stats = per_acc.get(
                    acc,
                    {
                        "tested": 0,
                        "first_support": 0,
                        "first_outside": 0,
                        "second_support": 0,
                        "second_outside": 0,
                        "final_support": 0,
                        "final_outside": 0,
                    },
                )
                tested = stats["tested"]
                final_outside = stats["final_outside"]
                frac = (final_outside / tested) if tested else ""
                fh.write(
                    "\t".join(
                        [
                            acc,
                            str(tested),
                            str(stats["first_support"]),
                            str(stats["first_outside"]),
                            str(stats["second_support"]),
                            str(stats["second_outside"]),
                            str(stats["final_support"]),
                            str(stats["final_outside"]),
                            str(frac) if frac != "" else "",
                            str(summary_decisions.get(acc, "")),
                        ]
                    )
                    + "\n"
                )
        self.log(
            f"Wrote external decontamination summary: {external_summary}",
            "INFO",
        )

    def _run_external_stage(self) -> str | bool:
        if not self._external_check_enabled():
            return True

        pending = self._external_pending_votes()
        if pending and not self.external_reuse_blast_results:
            if not self.external_blast_db_path:
                self.log("External check: blast_db_path not set; skipping external BLAST run.", "WARNING")
                self.data["_external_applied"] = True
                try:
                    self.db_manager.tasks.update_data(self.task_id, data=self.data)
                except Exception as exc:  # boundary: persist external-skip metadata only.
                    self.log(f"Failed to persist skipped external-check state: {exc}", "WARNING")
                try:
                    self.checkpoint(stage=5, checkpoint_data={"_external_applied": True})
                except Exception as exc:  # boundary: persist external-skip checkpoint only.
                    self.log(f"Failed to checkpoint skipped external-check state: {exc}", "WARNING")
                return True
            def _queue_external():
                available_threads = max(1, int(self.REQUIRED_THREADS))
                external_max_concurrent = 2 if available_threads >= 2 else 1
                external_threads = max(1, min(32, available_threads // external_max_concurrent))
                self.queue_subtask(
                    job_type=23,
                    status="P",
                    priority=1,
                    data={
                        "run_id": self.run_id,
                        "library_id": self.library_id,
                        "blast_db_path": self.external_blast_db_path,
                        "blast_db_type": self.external_blast_db_type,
                        "blast_program": self.external_blast_program,
                        "output_dir": self.external_blast_output_dir,
                        "reuse_blast_results": self.external_reuse_blast_results,
                        "max_target_seqs": self.external_max_target_seqs,
                        "hit_window": self.hit_window,
                        "force": bool(self.force),
                        "max_concurrent": external_max_concurrent,
                        "threads": external_threads,
                        "required_threads": external_threads,
                    },
                )
                return True

            outcome = self.manage_subtasks(
                stage=4,
                queue_fn=_queue_external,
                done_fn=self._external_subtasks_complete,
                wait_seconds=0,
                retry_key="external_retries",
                max_retries=int(self.data.get("external_retries", 0)),
                incomplete_message_fn=lambda: ("External BLAST subtask did not complete.", ""),
                retry_incomplete=False,
            )
            if outcome == "ERROR":
                return "ERROR"
            if outcome is False:
                return False

        if self.data.get("_external_applied"):
            return True

        applied = self._apply_external_results()
        if applied:
            acc_taxa = {}
            acc_rank_taxa = {}
            for acc in self.accessions:
                genome = self.db_manager.genomes.get(acc)
                taxid = genome[1] if genome else None
                acc_taxa[acc] = taxid
                acc_rank_taxa[acc] = self._taxon_at_rank(taxid)
            self._write_external_report(acc_taxa, acc_rank_taxa)
        self.data["_external_applied"] = True
        try:
            self.db_manager.tasks.update_data(self.task_id, data=self.data)
        except Exception as exc:  # boundary: persist external-apply metadata only.
            self.log(f"Failed to persist external-apply state: {exc}", "WARNING")
        try:
            self.checkpoint(stage=5, checkpoint_data={"_external_applied": True})
        except Exception as exc:  # boundary: persist external-apply checkpoint only.
            self.log(f"Failed to checkpoint external-apply state: {exc}", "WARNING")
        return True

    def _build_family_stats(self, id_map: dict, acc_rank_taxa: dict, acc_taxa: dict):
        family_counts = {}
        family_acc_counts = {}
        family_rank_counts = {}
        family_taxids = {}
        for meta in id_map.values():
            fam = meta.get("family_id")
            acc = meta.get("accession")
            taxid = meta.get("taxid") or acc_taxa.get(acc)
            if fam is None or acc is None:
                continue
            family_counts[fam] = family_counts.get(fam, 0) + 1
            family_acc_counts.setdefault(fam, {})
            family_acc_counts[fam][acc] = family_acc_counts[fam].get(acc, 0) + 1
            rank_taxid = acc_rank_taxa.get(acc)
            if rank_taxid:
                family_rank_counts.setdefault(fam, {})
                family_rank_counts[fam][rank_taxid] = family_rank_counts[fam].get(rank_taxid, 0) + 1
            if taxid:
                family_taxids.setdefault(fam, []).append(taxid)
        return family_counts, family_acc_counts, family_rank_counts, family_taxids

    def run(self):
        if self._load_config() == "ERROR":
            return "ERROR"

        # Combine CLI targets with explicit config targets
        config_targets = list(dict.fromkeys(self.config_targets or []))
        if self.accessions:
            if config_targets:
                merged = list(self.accessions) + config_targets
                self.accessions = list(dict.fromkeys(merged))
                self.data["accessions"] = self.accessions
                try:
                    self.db_manager.tasks.update_data(self.task_id, data=self.data)
                except Exception as exc:  # boundary: persist merged internal-decontam target metadata only.
                    self.log(f"Failed to persist merged internal-decontam targets: {exc}", "WARNING")
        elif config_targets:
            self.accessions = list(dict.fromkeys(config_targets))
            self.data["accessions"] = self.accessions
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist config internal-decontam target metadata only.
                self.log(f"Failed to persist config internal-decontam targets: {exc}", "WARNING")

        # Check we have accessions
        if not self.accessions:
            return self.handle_exception(
                "No accessions provided for internal decontamination.",
                {"accessions": self.accessions, "config_targets": list(self.config_targets or [])},
            )

        if self.use_paralog_filtered_buscos:
            self.log(
                "Internal decontam: use_paralog_filtered_buscos requested, but internal BLAST DB "
                "requires full 1-1 coverage. Ignoring this flag.",
                "WARNING",
            )
            self.use_paralog_filtered_buscos = False

        # Chekk we have a library or can resolve one from name
        if self.library_name and not self.library_id:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
        if not self.library_id:
            return self.handle_exception(
                "library_id (or resolvable library_name) is required.",
                {"library_id": self.library_id, "library_name": self.library_name},
            )

        parent_library_id = self.db_manager.libraries.get_parent_id(self.library_id)
        self.busco_lib_id = parent_library_id if parent_library_id else self.library_id
        libs = self.db_manager.libraries.get(self.busco_lib_id)
        if not libs:
            return self.handle_exception("BUSCO library not found in database.", {"busco_library_id": self.busco_lib_id})
        busco_lineage_name = libs[0][1]
        # Generate deterministic run_id (unless provided); allow deduplication/resume
        if not self.run_id:
            fingerprint = {
                "targets": sorted(self.accessions),
                "rank": self.rank,
                "hit_window": self.hit_window,
                "p_value_threshold": self.p_value_threshold,
                "off_clade_fraction": self.off_clade_fraction,
                "min_buscos": self.min_buscos,
                "min_identity": self.min_identity,
                "min_coverage": self.min_coverage,
                "min_alignment_length": self.min_alignment_length,
                "min_bitscore": self.min_bitscore,
                "max_evalue": self.max_evalue,
                "max_target_seqs": self.max_target_seqs,
                "config_signature": self.config_signature,
            }
            digest = hashlib.sha1(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
            self.run_id = f"idc_{digest}"
            self.data["run_id"] = self.run_id
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist generated internal-decontam run id metadata only.
                self.log(f"Failed to persist generated internal-decontam run id {self.run_id}: {exc}", "WARNING")

        external_check_enabled = self._external_check_enabled()
        if external_check_enabled and not self.external_blast_output_dir and not self.external_reuse_blast_results:
            tmpdir = tempfile.mkdtemp(prefix=f"idc_external_{self.run_id}_")
            self.external_blast_output_dir = tmpdir
            self.data["external_blast_output_dir"] = self.external_blast_output_dir
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist temporary external BLAST output dir metadata only.
                self.log(f"Failed to persist external BLAST output directory: {exc}", "WARNING")

        # Ensure genomes are downloaded for queries
        missing_downloads = self._get_missing_downloads(self.accessions)
        if missing_downloads:
            return self.handle_exception(
                "Some accessions are not downloaded; download them before internal decontamination.",
                {"missing_accessions": missing_downloads},
            )

        # Stage 1: ensure BUSCO results for queries
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

        # Check blast set up and detect DB type
        db_type = self._detect_busco_db_type()
        self.blastp_path = self._resolve_blast_program(db_type)
        if not self.blastp_path:
            missing_var = "BLASTN_PATH" if db_type == "nucl" else "BLASTP_PATH"
            return self.handle_exception(
                f"{missing_var} is not set in environment variables.",
                {"variable": missing_var},
            )

        # Stage 2: build global BUSCO BLAST DB for targets
        if not self.internal_blastdb_path:
            tmpdir = tempfile.mkdtemp(prefix=f"idc_refdb_{self.run_id}_")
            self.internal_blastdb_path = os.path.join(tmpdir, "internal_buscos")
            self.data["_internal_blastdb_path"] = self.internal_blastdb_path
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist temporary internal BLAST DB path metadata only.
                self.log(f"Failed to persist internal BLAST DB path: {exc}", "WARNING")

        if not self.internal_id_map_path:
            self.internal_id_map_path = f"{self.internal_blastdb_path}.id_map.tsv"
            self.data["_internal_blastdb_id_map"] = self.internal_id_map_path
            try:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            except Exception as exc:  # boundary: persist internal BLAST DB id-map path metadata only.
                self.log(f"Failed to persist internal BLAST DB id-map path: {exc}", "WARNING")

        def _queue_internal_busco_db():
            self.queue_subtask(
                job_type=17,
                status="P",
                priority=1,
                data={
                    "accessions": self.accessions,
                    "busco_library_id": self.busco_lib_id,
                    "target_library_id": self.library_id,
                    "output_path": self.internal_blastdb_path,
                    "use_paralog_filtered": False,
                    "force": True,
                    "id_mode": "internal",
                    "id_map_path": self.internal_id_map_path,
                    "db_type": db_type,
                },
            )
            return True

        def _internal_db_done():
            required_exts = ["psq", "pin", "phr"] if db_type == "prot" else ["nsq", "nin", "nhr"]
            return all(os.path.exists(f"{self.internal_blastdb_path}.{ext}") for ext in required_exts) and os.path.exists(
                self.internal_id_map_path
            )

        outcome = self.manage_subtasks(
            stage=2,
            queue_fn=_queue_internal_busco_db,
            done_fn=_internal_db_done,
            wait_seconds=0,
            retry_key=None,
            max_retries=0,
            incomplete_message_fn=lambda: (f"Internal BUSCO BLAST DB not ready at {self.internal_blastdb_path}", ""),
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        if getattr(self, "stage", 0) >= 3:
            return self._run_external_stage()

        # Record run metadata in db
        run_params = {
            "rank": self.rank,
            "hit_window": self.hit_window,
            "p_value_threshold": self.p_value_threshold,
            "min_buscos": self.min_buscos,
            "min_identity": self.min_identity,
            "min_coverage": self.min_coverage,
            "min_alignment_length": self.min_alignment_length,
            "min_bitscore": self.min_bitscore,
            "max_evalue": self.max_evalue,
            "max_target_seqs": self.max_target_seqs,
            "config_path": self.config_path,
            "config_signature": self.config_signature,
            "run_label": self.run_label,
            "blast_db_type": db_type,
            "external_blast_db_path": self.external_blast_db_path,
            "external_blast_db_type": self.external_blast_db_type,
            "external_blast_program": self.external_blast_program,
            "external_blast_output_dir": self.external_blast_output_dir,
            "external_reuse_blast_results": self.external_reuse_blast_results,
            "external_max_target_seqs": self.external_max_target_seqs,
        }
        if not self.db_manager.filtering.add_decontamination_run(
            run_id=self.run_id,
            target_library_id=self.library_id,
            busco_library_id=self.busco_lib_id,
            targets_json=json.dumps(self.accessions),
            refs_json=json.dumps(self.accessions),
            params_json=json.dumps(run_params),
            config_signature=self.config_signature,
            run_label=self.run_label,
        ):
            return self.handle_exception("Failed to record internal decontamination run metadata.", {"run_id": self.run_id})
        
        # Load BUSCO BLAST ID map for later use in processing; this should be a small file with metadata for each sequence in the internal BLAST DB
        id_map = self._load_id_map(self.internal_id_map_path)
        if not id_map:
            return self.handle_exception(
                "Failed to load BUSCO BLAST ID map.",
                {"id_map_path": self.internal_id_map_path},
            )

        # Precompute taxonomy for accessions and map counts per family
        acc_taxa = {}
        acc_rank_taxa = {}
        for acc in self.accessions:
            genome = self.db_manager.genomes.get(acc)
            taxid = genome[1] if genome else None
            acc_taxa[acc] = taxid
            acc_rank_taxa[acc] = self._taxon_at_rank(taxid)

        family_counts, family_acc_counts, family_rank_counts, family_taxids = self._build_family_stats(
            id_map, acc_rank_taxa, acc_taxa
        )

        params_json = self._build_params_json()

        #--------
        def process_accession(acc):

            # db connext
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
            group_hit_window = self.hit_window
            if group and "hit_window" in group:
                try:
                    group_hit_window = max(1, int(group.get("hit_window")))
                except (TypeError, ValueError):
                    group_hit_window = self.hit_window

            # If we can't determine an expected taxon at the rank, and there's no group definition that would allow us to evaluate hits, then we can't make any determination for this accession; mark as uncertain and move on
            if not expected_taxon and not (group and group.get("clades")):
                self.log(
                    f"No taxonomy at rank={self.rank} for {acc}; marking uncertain.",
                    "WARNING",
                )
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

            # get BUSCO rows for this accession; if none found, mark as uncertain and move on (this could happen if the BUSCO subtask failed or is still running, but we should have filtered those out in the waiting stage)
            busco_rows = thread_db.busco.get_family_results_for_library(
                library_id=self.busco_lib_id,
                accessions=[acc],
                status=[1],
            )
            if not busco_rows:
                self.log(
                    f"No BUSCO entries found for {acc} in library {self.busco_lib_id}.",
                    "WARNING",
                )
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

            # if we've already recorded votes for some families for this accession and run, skip those families (this can happen when re-running with force=False to only process families that haven't been done yet)
            existing_votes = thread_db.filtering.get_decontamination_votes(accession=acc, run_id=self.run_id)
            done_fams = {v[0] for v in existing_votes} if existing_votes else set()

            buscos_tested = 0
            buscos_supporting = 0
            buscos_outside = 0
            best_taxon_counts = {}
            vote_rows: list[dict[str, object]] = []

            def _queue_vote(
                family_id: str,
                *,
                best_taxid,
                runner_taxid,
                best_bitscore,
                decision: str,
                payload: dict[str, object],
            ) -> None:
                vote_rows.append(
                    {
                        "family_id": family_id,
                        "busco_library_id": self.busco_lib_id,
                        "target_library_id": self.library_id,
                        "accession": acc,
                        "run_id": self.run_id,
                        "busco_run_id": effective_busco_run_id,
                        "expected_taxid": expected_taxon,
                        "best_taxid": best_taxid,
                        "runner_taxid": runner_taxid,
                        "rank": self.rank,
                        "best_bitscore": best_bitscore,
                        "delta_bitscore": 0,
                        "decision": decision,
                        "top_hits_json": json.dumps(payload),
                    }
                )

            for family_id, lib_id, accession, status, sequence, score, length in busco_rows:
                if family_id in done_fams and not self.force:
                    continue
                location = thread_db.busco.get_family_location(
                    family_id,
                    lib_id,
                    accession,
                    sequence_kind=self.blast_db_type if self.blast_db_type in ("prot", "nucl") else None,
                )
                if not location or not os.path.exists(location):
                    continue
                out = None
                if self.reuse_blast_results:
                    cache_path = self._blast_cache_path(self.reuse_blast_results, acc, family_id)
                    out = self._load_blast_output(cache_path)
                    if out is None:
                        self.log(
                            f"Internal decontam: missing cached BLAST output for {acc} family {family_id} at {cache_path}",
                            "WARNING",
                        )
                        if self.save_blast_output:
                            out = self._run_blast(location, self.internal_blastdb_path)
                            if out is False:
                                continue
                            save_path = self._blast_cache_path(self.save_blast_output, acc, family_id)
                            self._save_blast_output(save_path, out)
                        else:
                            _queue_vote(
                                str(family_id),
                                best_taxid=None,
                                runner_taxid=None,
                                best_bitscore=None,
                                decision="unknown",
                                payload={"reason": "missing_blast_cache"},
                            )
                            continue
                else:
                    out = self._run_blast(location, self.internal_blastdb_path)
                    if out is False:
                        continue
                    if self.save_blast_output:
                        save_path = self._blast_cache_path(self.save_blast_output, acc, family_id)
                        self._save_blast_output(save_path, out)
                hits = self._parse_blast_output(out)
                if not hits:
                    # Record unknown vote for missing hits
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={"reason": "no_hits"},
                    )
                    continue

                # Filter hits to same family and apply thresholds
                filtered = []
                for hit in hits:
                    meta = id_map.get(hit["sseqid"])
                    if not meta:
                        continue
                    if meta.get("family_id") != family_id:
                        continue
                    if self.min_identity and hit["pident"] < self.min_identity:
                        continue
                    if self.min_coverage and hit["qcovs"] < self.min_coverage:
                        continue
                    if self.min_alignment_length and hit["length"] < self.min_alignment_length:
                        continue
                    if self.min_bitscore and hit["bitscore"] < self.min_bitscore:
                        continue
                    if self.max_evalue is not None and hit["evalue"] > self.max_evalue:
                        continue
                    hit["hit_accession"] = meta.get("accession")
                    hit["hit_taxid"] = meta.get("taxid")
                    hit["hit_rank_taxid"] = self._taxon_at_rank(hit["hit_taxid"], db=thread_db) if hit.get("hit_taxid") else None
                    filtered.append(hit)

                # Sort by bitscore
                filtered.sort(key=lambda h: (-h["bitscore"], h["evalue"]))
                # Deduplicate by hit accession (keep best-scoring HSP per accession)
                deduped = []
                seen_accs = set()
                for hit in filtered:
                    acc_hit = hit.get("hit_accession") or hit.get("sseqid")
                    if acc_hit in seen_accs:
                        continue
                    seen_accs.add(acc_hit)
                    deduped.append(hit)
                filtered = deduped

                # Remove self-hit (by exact query id if present; else by accession)
                removed_self = False
                for idx, hit in enumerate(filtered):
                    if hit["sseqid"] == hit.get("qseqid"):
                        filtered.pop(idx)
                        removed_self = True
                        break
                if not removed_self:
                    for idx, hit in enumerate(filtered):
                        if hit.get("hit_accession") == acc:
                            filtered.pop(idx)
                            removed_self = True
                            break
                if not removed_self:
                    self.log(
                        f"Internal decontam: self-hit not found for {acc} family {family_id}; proceeding.",
                        "WARNING",
                    )
                
                # Total number in family (excluding self if present) defines Z; if Z=0 then we have no power to make any determination, so record unknown vote and move on
                total_family = family_counts.get(family_id, 0)
                self_present = family_acc_counts.get(family_id, {}).get(acc, 0)
                Z = max(total_family - (1 if self_present else 0), 0)
                if Z <= 0:
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={"reason": "insufficient_family_size", "Z": Z},
                    )
                    continue

                available_hits = len(filtered)
                if available_hits <= 0:
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={"reason": "no_hits_after_filter"},
                    )
                    continue

                Y = min(int(group_hit_window), available_hits, Z)
                if Y <= 0:
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={"reason": "no_hits_in_window"},
                    )
                    continue

                group_clades = group.get("clades") if group else set()
                group_blacklist = group.get("blacklist") if group else set()

                # Checks if a hit is allowed by decent
                def _allowed(hit_taxid, hit_rank_taxid):
                    if not hit_taxid and not hit_rank_taxid:
                        return False
                    if group_clades:
                        if any(self._is_descendant(hit_taxid, blk, db=thread_db) for blk in group_blacklist):
                            return False
                        return any(self._is_descendant(hit_taxid, cl, db=thread_db) for cl in group_clades)
                    return hit_rank_taxid == expected_taxon

                window_hits = filtered[:Y]
                allowed_flags = [_allowed(h.get("hit_taxid"), h.get("hit_rank_taxid")) for h in window_hits]
                # X = number of observed allowed hits
                x = sum(1 for ok in allowed_flags if ok)

                # Compute K based on group definition, K is the number of possible allowed hits in the family (excluding self if present); 
                # if group_clades is defined, we count the number of family members that would be allowed by the group definition; 
                # otherwise we take the count of family members with the expected taxon at the rank; if K=0 then we have no power to make any determination, so record unknown vote and move on
                if group_clades:
                    fam_taxids = family_taxids.get(family_id, [])
                    K = 0
                    for tid in fam_taxids:
                        if any(self._is_descendant(tid, blk, db=thread_db) for blk in group_blacklist):
                            continue
                        if any(self._is_descendant(tid, cl, db=thread_db) for cl in group_clades):
                            K += 1
                    if self_present:
                        # Assume the query belongs to the group when explicitly grouped
                        K = max(K - 1, 0)
                else:
                    K = family_rank_counts.get(family_id, {}).get(expected_taxon, 0)
                    if self_present and expected_taxon:
                        K = max(K - 1, 0)

                if K <= 0:
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={"reason": "K_zero", "Z": Z, "K": K, "Y": Y, "x": x},
                    )
                    continue

                feasibility = self._hypergeom_decision_feasibility(Z, K, Y, self.p_value_threshold)
                # If even the best-case outcome cannot be significant, this family cannot support the expected taxon.
                if not feasibility["win_possible"]:
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={
                            "reason": "significance_impossible",
                            "Z": Z,
                            "K": K,
                            "Y": Y,
                            "x": x,
                            "p_best": feasibility["p_best"],
                            "p_worst": feasibility["p_worst"],
                        },
                    )
                    continue
                # Keep the complementary "must be able to lose" invariant explicit even though
                # the current hypergeometric setup will normally make this true automatically.
                if not feasibility["lose_possible"]:
                    _queue_vote(
                        str(family_id),
                        best_taxid=None,
                        runner_taxid=None,
                        best_bitscore=None,
                        decision="unknown",
                        payload={
                            "reason": "loss_impossible",
                            "Z": Z,
                            "K": K,
                            "Y": Y,
                            "x": x,
                            "p_best": feasibility["p_best"],
                            "p_worst": feasibility["p_worst"],
                        },
                    )
                    continue
                
                # calculate p-value for observed X; if significant, then we have support for the expected taxon; if not significant, then the hits are more consistent with being outside the expected taxon; record vote and metadata for this family and accession
                p_value = self._hypergeom_sf(x - 1, Z, K, Y)
                internal_decision = "support" if p_value < self.p_value_threshold else "outside"
                decision = internal_decision

                best_allowed = next((h for h, ok in zip(window_hits, allowed_flags) if ok), None)
                best_outside = next((h for h, ok in zip(window_hits, allowed_flags) if not ok), None)
                best_taxon = best_allowed.get("hit_rank_taxid") if best_allowed else None
                runner_taxon = best_outside.get("hit_rank_taxid") if best_outside else None
                best_bitscore = best_allowed.get("bitscore") if best_allowed else (best_outside.get("bitscore") if best_outside else None)

                external_check_run = False
                external_check_pending = False
                external_decision = None
                external_hits = []
                external_allowed = []
                external_top_taxon = None
                external_top_acc = None
                external_top_bits = None
                external_top_eval = None
                external_hit_window = None

                if internal_decision == "outside" and external_check_enabled:
                    ext_base = self.external_reuse_blast_results or self.external_blast_output_dir
                    if ext_base:
                        ext_path = self._blast_cache_path(ext_base, acc, family_id)
                        ext_out = self._load_blast_output(ext_path)
                        if ext_out is None:
                            external_check_pending = True
                        else:
                            external_check_run = True
                            ext_hits = self._parse_blast_output_with_taxids(ext_out)
                            min_hits_req = getattr(self, "min_hits", 1)
                            if group and "min_hits" in group:
                                try:
                                    min_hits_req = max(1, int(group.get("min_hits")))
                                except (TypeError, ValueError):
                                    min_hits_req = getattr(self, "min_hits", 1)
                            result = self._evaluate_external_hits(
                                ext_hits,
                                expected_taxon=expected_taxon,
                                group_clades=group_clades,
                                group_blacklist=group_blacklist,
                                group_hit_window=group_hit_window,
                                min_hits_req=min_hits_req,
                                db=thread_db,
                            )
                            external_decision = result.get("decision")
                            external_hits = result.get("hits") or []
                            external_allowed = result.get("allowed") or []
                            external_top_taxon = result.get("top_taxon")
                            external_top_acc = result.get("top_acc")
                            external_top_bits = result.get("top_bits")
                            external_top_eval = result.get("top_eval")
                            external_hit_window = result.get("hit_window")

                            if external_decision == "support":
                                decision = "support"
                                best_taxon = result.get("best_rank_taxon") or best_taxon
                                runner_taxon = result.get("runner_rank_taxon") or runner_taxon
                                if best_taxon is None and expected_taxon:
                                    best_taxon = expected_taxon
                    else:
                        external_check_pending = True

                if decision == "support":
                    buscos_supporting += 1
                else:
                    buscos_outside += 1

                if best_taxon:
                    best_taxon_counts[best_taxon] = best_taxon_counts.get(best_taxon, 0) + 1

                payload = {
                    "Z": Z,
                    "K": K,
                    "Y": Y,
                    "x": x,
                    "p_value": p_value,
                    "internal_decision": internal_decision,
                    "decision": decision,
                    "external_check_run": external_check_run,
                    "external_check_pending": external_check_pending,
                    "external_decision": external_decision,
                    "external_hit_window": external_hit_window,
                    "external_hits": external_hits,
                    "external_allowed": external_allowed,
                    "external_top_hit_taxon": external_top_taxon,
                    "external_top_hit_accession": external_top_acc,
                    "external_top_hit_bitscore": external_top_bits,
                    "external_top_hit_evalue": external_top_eval,
                    "hits": window_hits,
                    "allowed": allowed_flags,
                }
                _queue_vote(
                    str(family_id),
                    best_taxid=best_taxon,
                    runner_taxid=runner_taxon,
                    best_bitscore=best_bitscore,
                    decision=decision,
                    payload=payload,
                )
                buscos_tested += 1

            majority_taxon = None
            if best_taxon_counts:
                majority_taxon = max(best_taxon_counts.items(), key=lambda kv: kv[1])[0]

            off_frac = (buscos_outside / buscos_tested) if buscos_tested else None
            if buscos_tested < self.min_buscos:
                final_decision = "UNCERTAIN"
            elif buscos_outside == 0:
                final_decision = "CLEAN"
            elif off_frac is not None and off_frac >= self.off_clade_fraction:
                final_decision = "CONTAMINATED"
            else:
                final_decision = "CLEAN"

            if vote_rows and not thread_db.filtering.add_decontamination_votes(vote_rows):
                raise RuntimeError(f"Failed to persist internal decontamination votes for {acc}")
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
                raise RuntimeError(f"Failed to persist internal decontamination summary for {acc}")
            try:
                thread_db.close()
            except Exception as exc:  # boundary: worker cleanup failure is logged after primary work completes/fails.
                self.log(f"Failed to close worker database connection for internal decontamination: {exc}", "WARNING")
            return acc, final_decision

        #--------

        effective_max = self.max_concurrent if self.max_concurrent and self.max_concurrent > 0 else self.REQUIRED_THREADS
        max_workers = max(1, min(len(self.accessions), effective_max, self.REQUIRED_THREADS))
        self.blast_threads = max(1, math.ceil(self.REQUIRED_THREADS / max_workers)) if max_workers else self.REQUIRED_THREADS
        self._validate_or_disable_internal_cache_reuse()
        if self.save_blast_output:
            self._write_cache_key(self.save_blast_output)
        if self.reuse_blast_results:
            self._write_cache_key(self.reuse_blast_results)

        if self.save_blast_output:
            self.log(f"Internal decontam: saving BLAST output to {self.save_blast_output}", "INFO")
        if self.reuse_blast_results:
            self.log(f"Internal decontam: reusing BLAST output from {self.reuse_blast_results}", "INFO")

        self.log(
            f"Stage 3: running internal decontamination with max_workers={max_workers} "
            f"(blast threads per task={self.blast_threads}); targets={len(self.accessions)}",
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
                        f"Internal decontam progress: {progress_next}% ({progress_done}/{progress_total})",
                        "INFO",
                    )
                    progress_next += 10

        # Process accessions in parallel threads; each thread gets its own DB connection to avoid locking issues; results are collected in a thread-safe way and progress is logged periodically
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="IDC_") as executor:
            futures = {executor.submit(process_accession, acc): acc for acc in self.accessions}
            for future in concurrent.futures.as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as e:  # boundary: one accession failure is aggregated while independent targets continue.
                    acc = futures[future]
                    self.error(f"Internal decontamination failed for {acc}: {e}")
                    results.append((acc, "ERROR"))
                _log_progress()

        failed = [acc for acc, status in results if status in ("ERROR",)]
        if failed:
            preview = ", ".join(failed[:10])
            suffix = "" if len(failed) <= 10 else ", ..."
            self.log(
                f"Internal decontamination completed with issues for {len(failed)}/{len(self.accessions)} targets: {preview}{suffix}",
                "WARNING",
            )
        else:
            self.log(f"Internal decontamination completed for {len(self.accessions)} targets.", "INFO")

        if self.report_path:
            try:
                base, ext = os.path.splitext(self.report_path)
                busco_report = self.report_path if ext else f"{base}_buscos.tsv"
                summary_report = f"{base}_summary.tsv" if ext else f"{self.report_path}_summary.tsv"
                blast_report = f"{base}_blast.tsv" if ext else f"{self.report_path}_blast.tsv"
                refs_report = f"{base}_refs.tsv" if ext else f"{self.report_path}_refs.tsv"
                groups_report = f"{base}_groups.tsv" if ext else f"{self.report_path}_groups.tsv"
                summary_by_taxon = f"{base}_summary_by_taxon.tsv" if ext else f"{self.report_path}_summary_by_taxon.tsv"
                summary_by_family = f"{base}_summary_by_family.tsv" if ext else f"{self.report_path}_summary_by_family.tsv"
                os.makedirs(os.path.dirname(busco_report) or ".", exist_ok=True)

                votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.run_id)

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

                with open(busco_report, "w") as fh, open(blast_report, "w") as fh_blast:
                    fh.write(
                        "\t".join(
                            [
                                "family_id",
                                "accession",
                                "query_taxon_id",
                                "query_label",
                                "query_group",
                                "Z",
                                "K",
                                "Y",
                                "x",
                                "p_value",
                                "internal_call",
                                "external_check_run",
                                "external_top_hit_taxon",
                                "external_top_hit_accession",
                                "external_top_hit_bitscore",
                                "external_top_hit_evalue",
                                "final_call",
                                "reason",
                                "hits_in_window",
                                "correct_in_window",
                                "wrong_in_window",
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
                                "hit_index",
                                "ref_accession",
                                "ref_taxid",
                                "ref_rank_taxid",
                                "ref_rank_name",
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
                        payload = {}
                        try:
                            payload = json.loads(top_hits_json) if top_hits_json else {}
                        except (TypeError, json.JSONDecodeError):
                            payload = {}
                        hits = payload.get("hits") or []
                        allowed_flags = payload.get("allowed") or []
                        Z = payload.get("Z")
                        K = payload.get("K")
                        Y = payload.get("Y")
                        x = payload.get("x")
                        p_value = payload.get("p_value")
                        reason = payload.get("reason") or payload.get("decision") or ""
                        query_taxon = expected_taxid
                        query_label = tax_name(expected_taxid) or expected_taxid or ""
                        query_group = tax_name(expected_taxid) or expected_taxid or ""

                        external_check_run = bool(payload.get("external_check_run"))
                        external_top_taxon = payload.get("external_top_hit_taxon")
                        external_top_acc = payload.get("external_top_hit_accession")
                        external_top_bits = payload.get("external_top_hit_bitscore")
                        external_top_eval = payload.get("external_top_hit_evalue")

                        if not external_check_run:
                            external_top = None
                            if decision == "outside" and hits:
                                for h, ok in zip(hits, allowed_flags):
                                    if not ok:
                                        external_top = h
                                        break
                            if external_top_taxon is None and external_top is not None:
                                external_top_taxon = external_top.get("hit_taxid")
                            if external_top_acc is None and external_top is not None:
                                external_top_acc = external_top.get("hit_accession")
                            if external_top_bits is None and external_top is not None:
                                external_top_bits = external_top.get("bitscore")
                            if external_top_eval is None and external_top is not None:
                                external_top_eval = external_top.get("evalue")

                        hits_in_window = len(hits)
                        correct_in_window = sum(1 for ok in allowed_flags if ok)
                        wrong_in_window = hits_in_window - correct_in_window

                        final_call = "unknown"
                        if decision == "support":
                            final_call = "retain"
                        elif decision == "outside":
                            final_call = "contaminant"

                        row_out = [
                            family_id,
                            accession,
                            query_taxon or "",
                            query_label,
                            query_group,
                            Z if Z is not None else "",
                            K if K is not None else "",
                            Y if Y is not None else "",
                            x if x is not None else "",
                            p_value if p_value is not None else "",
                            decision,
                            "1" if external_check_run else "0",
                            external_top_taxon or "",
                            external_top_acc or "",
                            external_top_bits or "",
                            external_top_eval or "",
                            final_call,
                            reason,
                            hits_in_window,
                            correct_in_window,
                            wrong_in_window,
                        ]
                        fh.write("\t".join("" if v is None else str(v) for v in row_out) + "\n")

                        for idx, h in enumerate(hits, start=1):
                            allowed = allowed_flags[idx - 1] if idx - 1 < len(allowed_flags) else False
                            ref_taxid = h.get("hit_taxid")
                            ref_rank_taxid = self._taxon_at_rank(ref_taxid) if ref_taxid else None
                            ref_rank_name = tax_name(ref_rank_taxid) or ref_rank_taxid or ""
                            row_hit = [
                                family_id,
                                accession,
                                decision,
                                idx,
                                h.get("hit_accession") or "",
                                ref_taxid or "",
                                ref_rank_taxid or "",
                                ref_rank_name,
                                "1" if allowed else "0",
                                h.get("pident", ""),
                                h.get("qcovs", ""),
                                h.get("length", ""),
                                h.get("evalue", ""),
                                h.get("bitscore", ""),
                            ]
                            fh_blast.write("\t".join("" if v is None else str(v) for v in row_hit) + "\n")

                summaries = self.db_manager.filtering.get_decontamination_summary(run_id=self.run_id)
                implicit_rank_taxids = sorted({rt for rt in acc_rank_taxa.values() if rt})
                implicit_group_ids = {rt: idx + 1 for idx, rt in enumerate(implicit_rank_taxids)}
                custom_group_offset = len(implicit_rank_taxids)
                custom_group_ids = {idx: custom_group_offset + idx + 1 for idx in range(len(self.groups))}
                group_ids_by_acc = {}
                for acc in self.accessions:
                    acc_ids = []
                    rank_taxid = acc_rank_taxa.get(acc)
                    if rank_taxid and rank_taxid in implicit_group_ids:
                        acc_ids.append(implicit_group_ids[rank_taxid])
                    taxid = acc_taxa.get(acc)
                    for idx, g in enumerate(self.groups):
                        if acc in g.get("member_accessions", set()):
                            acc_ids.append(custom_group_ids[idx])
                            continue
                        if taxid and any(self._is_descendant(taxid, mt, db=self.db_manager) for mt in g.get("member_taxa", set())):
                            acc_ids.append(custom_group_ids[idx])
                    if acc_ids:
                        group_ids_by_acc[acc] = ",".join(str(i) for i in sorted(set(acc_ids)))
                    else:
                        group_ids_by_acc[acc] = ""
                with open(summary_report, "w") as fh:
                    fh.write(
                        "\t".join(
                            [
                                "accession",
                                "group_ids",
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
                                "off_clade_fraction",
                                "decision",
                                "params_json",
                            ]
                        )
                        + "\n"
                    )
                    decisions_by_target = {}
                    for v in votes:
                        acc = v[3]
                        dec = v[11]
                        decisions_by_target.setdefault(acc, []).append(dec)
                    for row in summaries:
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
                        unknown_ct = sum(1 for d in decisions_by_target.get(acc, []) if str(d).lower() == "unknown")
                        row_out = [
                            accession,
                            group_ids_by_acc.get(acc, ""),
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
                            off_clade_fraction,
                            decision,
                            params_json,
                        ]
                        fh.write("\t".join("" if v is None else str(v) for v in row_out) + "\n")

                with open(refs_report, "w") as fh:
                    fh.write("\t".join(["accession", "taxid", f"{self.rank}_taxid", f"{self.rank}_name"]) + "\n")
                    for acc in self.accessions:
                        taxid = acc_taxa.get(acc)
                        rank_taxid = acc_rank_taxa.get(acc)
                        fh.write(
                            "\t".join(
                                [
                                    acc,
                                    str(taxid or ""),
                                    str(rank_taxid or ""),
                                    str(tax_name(rank_taxid) or ""),
                                ]
                            )
                            + "\n"
                        )

                with open(groups_report, "w") as fh:
                    fh.write(
                        "\t".join(
                            [
                                "group_id",
                                "group_type",
                                "group_index",
                                "members",
                                "member_taxa",
                                "allowed_clades",
                                "blacklist",
                                "applied_targets",
                            ]
                        )
                        + "\n"
                    )
                    target_group_map = {}
                    implicit_group_map = {}
                    for acc in self.accessions:
                        rank_taxid = acc_rank_taxa.get(acc)
                        if rank_taxid:
                            implicit_group_map.setdefault(rank_taxid, []).append(acc)
                        for idx, g in enumerate(self.groups):
                            if acc in g.get("member_accessions", set()):
                                target_group_map.setdefault(idx, []).append(acc)
                                continue
                            taxid = acc_taxa.get(acc)
                            if taxid and any(self._is_descendant(taxid, mt, db=self.db_manager) for mt in g.get("member_taxa", set())):
                                target_group_map.setdefault(idx, []).append(acc)

                    for rank_taxid in implicit_rank_taxids:
                        group_id = implicit_group_ids[rank_taxid]
                        members = implicit_group_map.get(rank_taxid, [])
                        rank_name = tax_name(rank_taxid) or str(rank_taxid)
                        fh.write(
                            "\t".join(
                                [
                                    str(group_id),
                                    "implicit",
                                    "",
                                    rank_name,
                                    str(rank_taxid),
                                    "",
                                    "",
                                    ",".join(sorted(members)),
                                ]
                            )
                            + "\n"
                        )

                    for idx, g in enumerate(self.groups):
                        fh.write(
                            "\t".join(
                                [
                                    str(custom_group_ids[idx]),
                                    "custom",
                                    str(idx),
                                    ",".join(sorted(g.get("member_accessions") or [])),
                                    ",".join(str(t) for t in sorted(g.get("member_taxa") or [])),
                                    ",".join(str(t) for t in sorted(g.get("clades") or [])),
                                    ",".join(str(t) for t in sorted(g.get("blacklist") or [])),
                                    ",".join(sorted(target_group_map.get(idx, []))),
                                ]
                            )
                            + "\n"
                        )

                # Aggregate summaries by taxon and by family
                with open(summary_by_taxon, "w") as fh:
                    fh.write("\t".join(["taxon", "taxid", "tested", "support", "outside", "unknown", "contaminant_fraction"]) + "\n")
                    counts = {}
                    for row in votes:
                        (
                            fam_id,
                            _busco_lib,
                            _tlib,
                            acc,
                            _rid,
                            _exp,
                            _best,
                            _runner,
                            _rank,
                            _best_bits,
                            _delta,
                            decision,
                            top_hits_json,
                            *_vote_extra,
                        ) = row
                        taxid = acc_taxa.get(acc)
                        key = taxid or "unknown"
                        counts.setdefault(key, {"tested": 0, "support": 0, "outside": 0, "unknown": 0})
                        if decision == "support":
                            counts[key]["support"] += 1
                            counts[key]["tested"] += 1
                        elif decision == "outside":
                            counts[key]["outside"] += 1
                            counts[key]["tested"] += 1
                        else:
                            counts[key]["unknown"] += 1
                    for key, stats in counts.items():
                        tested = stats["tested"]
                        outside = stats["outside"]
                        frac = (outside / tested) if tested else ""
                        fh.write(
                            "\t".join(
                                [
                                    tax_name(key) if isinstance(key, int) else "unknown",
                                    str(key),
                                    str(tested),
                                    str(stats["support"]),
                                    str(outside),
                                    str(stats["unknown"]),
                                    str(frac) if frac != "" else "",
                                ]
                            )
                            + "\n"
                        )

                with open(summary_by_family, "w") as fh:
                    fh.write("\t".join(["family_id", "tested", "support", "outside", "unknown", "contaminant_fraction"]) + "\n")
                    fam_counts = {}
                    for row in votes:
                        fam_id = row[0]
                        decision = row[11]
                        fam_counts.setdefault(fam_id, {"tested": 0, "support": 0, "outside": 0, "unknown": 0})
                        if decision == "support":
                            fam_counts[fam_id]["support"] += 1
                            fam_counts[fam_id]["tested"] += 1
                        elif decision == "outside":
                            fam_counts[fam_id]["outside"] += 1
                            fam_counts[fam_id]["tested"] += 1
                        else:
                            fam_counts[fam_id]["unknown"] += 1
                    for fam_id, stats in fam_counts.items():
                        tested = stats["tested"]
                        outside = stats["outside"]
                        frac = (outside / tested) if tested else ""
                        fh.write(
                            "\t".join(
                                [
                                    str(fam_id),
                                    str(tested),
                                    str(stats["support"]),
                                    str(outside),
                                    str(stats["unknown"]),
                                    str(frac) if frac != "" else "",
                                ]
                            )
                            + "\n"
                        )

                self._write_external_summary(summaries)

                self.log(
                    f"Wrote internal decontamination reports under {os.path.dirname(busco_report) or '.'}.",
                    "INFO",
                )
            except Exception as e:  # boundary: optional report generation failure should not invalidate persisted analysis.
                self.log(f"Failed to write internal decontamination report to {self.report_path}: {e}", "WARNING")

        # Record active internal decontamination run for downstream tasks.
        if self.run_id and self.library_id:
            payload = {
                "run_id": self.run_id,
                "run_label": self.run_label,
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
            try:
                env_key = f"ACTIVE_INTERNAL_DECONT_RUN_{self.library_id}"
                if self.db_manager.env.set(env_key, payload):
                    self.log(
                        f"Recorded active internal decontamination run for library {self.library_id}: {payload}",
                        "DEBUG",
                    )
            except Exception as exc:  # boundary: active-run convenience pointer is optional; results have already been persisted.
                self.log(f"Failed to persist active internal decontamination pointer: {exc}", "WARNING")
            try:
                env_key = f"ACTIVE_DECONT_RUN_{self.library_id}"
                if self.db_manager.env.set(env_key, payload):
                    self.log(
                        f"Recorded active decontamination run for library {self.library_id}: {payload}",
                        "DEBUG",
                    )
            except Exception as exc:  # boundary: active-run convenience pointer is optional; results have already been persisted.
                self.log(f"Failed to persist active decontamination pointer from internal run: {exc}", "WARNING")

            try:
                self.checkpoint(stage=3, checkpoint_data={"_external_applied": self.data.get("_external_applied", False)})
            except Exception as exc:  # boundary: persist resumability checkpoint only.
                self.log(f"Failed to checkpoint internal decontamination completion: {exc}", "WARNING")

        return self._run_external_stage()
