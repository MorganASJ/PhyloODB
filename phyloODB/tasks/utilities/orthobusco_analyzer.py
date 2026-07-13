import os
import shutil
from Bio import SeqIO
import csv
import glob
import hashlib
import re
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Optional, List, Dict, Tuple, Set
from ete3 import Tree

from ...accession_utils import canonicalize_accession

# Inputs and setup

# BUSCO sequences:
# If busco_dir is set, it indexes those family FASTAs.
# If not set and you pass busco_results, it first collates complete BUSCOs
# (single-copy and multi-copy) into working_dir/busco_fastas and then indexes.
# Indexing uses uppercase AA sequence (md5 by default) so matches are exact by sequence.
# OrthoFinder results: reads Orthogroup_Sequences and a configured gene-tree directory
# containing Resolved_Gene_Trees.txt plus per-orthogroup *_tree.txt files.
# Mapping and basic reports

# Maps each orthogroup sequence to any BUSCO sequence with identical sequence.
# Writes:
# <id>_busco_to_orthogroup_map.tsv: row per matched sequence.
# <id>_og_to_busco_families.tsv: pairs OG -> BUSCO family.
# <id>_busco_family_to_ogs.tsv: pairs BUSCO family -> OG.
# Logs:
# OGs containing multiple BUSCO families.
# BUSCO families mapping to multiple OGs.
# 1-to-1 detection

# Relaxed 1-to-1: BUSCO family maps to exactly one OG (unique OG).
# Exact 1-to-1 (default): the OG contains BUSCO sequences from only that family and the multiset of BUSCO-matched sequences in the OG equals the BUSCO family’s multiset. Non-BUSCO extra sequences in the OG are ignored by default. (You can enforce “no extras at all” with require_no_extra_sequences=True.)
# Writes:
# <id>_busco_to_orthogroup_1to1.tsv with columns: BUSCO_family, Orthogroup, Exact_1_to_1, Relaxed_only.
# Occupancy check (applied before paralog analysis)

# For each exact 1-to-1 BUSCO↔OG, counts species present (from BUSCO-matched sequences only).
# Accepts pairs with Species_count >= min_species_in_trees; rejects the rest.
# Writes:
# <id>_busco_to_orthogroup_1to1_occupancy.tsv with Species_count and Status (accepted/rejected).
# Logs:
# Total exact 1-to-1, how many pass min_species, and the list of families rejected due to low occupancy.
# Paralog analysis (on accepted exact 1-to-1 OGs only)

# Parses per-OG gene trees. For each species with >1 copies:
# If the species’ copies are monophyletic → in-paralogs.
# If not monophyletic → out-paralogs.
# Writes:
# paralogs/OGXXXX_inparalogs.txt and paralogs/OGXXXX_outparalogs.txt (only if non-empty).
# paralogs/og_paralog_summary.tsv for the matched exact 1-to-1 set.
# Logs summaries for:
# All OGs (overall and above threshold).
# The exact 1-to-1 matched subset (overall and above threshold), including counts with in, out, both, none, and “No out-paralogs (in or none).”
# Final “good BUSCO families” list

# Defines good OGs as those exact 1-to-1/accepted by occupancy that have no out-paralogs (in-paralogs or none are fine).
# Converts those OGs back to their BUSCO families and writes:
# <id>_good_busco_families.txt (also returned as the function result).
# Meaning: each family in this list has a single, clean orthogroup that:
# Matches the BUSCO family exactly (by sequence identity of the BUSCO-matched members),
# Has sufficient species occupancy (>= min_species_in_trees),
# Shows no out-paralogs in the resolved gene tree (i.e., no evidence of ancient/lineage-duplicated copies breaking species monophyly).
# Augmented mapping output
# Writes:
# <id>_busco_to_orthogroup_map_with_paralog_class.tsv with OG paralog classification and a clean flag


class OrthoBuscoAnalyzer:
    """
    Cleaner, self-contained analyzer to compare OrthoFinder orthogroups with BUSCO families
    and summarize paralogy. Mirrors the external usage of the old analyzer whilst simplifying
    internals and outputs.

    Key behaviors:
    - Optionally collate complete BUSCO sequences (single-copy and multi-copy) into a unified directory (transform_busco_results)
    - Map Orthogroup sequences to BUSCO families by exact sequence match (hash-based by default)
    - Write mapping TSVs: <id>_busco_to_orthogroup_map.tsv and <id>_busco_to_orthogroup_1to1.tsv
    - Write augmented mapping TSV with paralog classification and a clean flag
    - Detect exact 1-1 matches (OG contains exactly the BUSCO family sequences: no duplicates, no missing, no extras)
    - Identify in- and out-paralogs per OG from OrthoFinder resolved gene trees and log summaries
    - Produce final list of BUSCO families whose OGs have no out-paralogs (in-paralogs or none) AND are exact 1-1
      (i.e., exclude BUSCOs mapping to multiple OGs or OGs containing multiple BUSCO families)
    - Return that final BUSCO list from compare_busco_orthofinder and write it to disk
    """

    def __init__(
        self,
        identifier: str,
        working_dir: str,
        orthofinder_run_folder: str,
        busco_dir: Optional[str] = None,
        *,
        gene_tree_dir: Optional[str] = None,
        append_log: bool = False,
    ):
        self.identifier = identifier
        self.working_dir = working_dir
        self.orthofinder_run_folder = orthofinder_run_folder
        self.busco_dir = busco_dir
        self.gene_tree_dir = gene_tree_dir or os.path.join(self.orthofinder_run_folder, "Resolved_Gene_Trees")
        os.makedirs(self.working_dir, exist_ok=True)
        self.log_file = os.path.join(self.working_dir, f"{identifier}.log")
        if not append_log:
            with open(self.log_file, "w", encoding="utf-8"):
                pass
        self.last_compare_results: Dict[str, object] = {}
        self._leaf_metadata_cache: Optional[Dict[str, Dict[str, str]]] = None
        self.log(f"Initialized OrthoBuscoAnalyzer with ID: {identifier}")

    # -------------------------------
    # Logging helper
    # -------------------------------
    def log(self, message: str):
        with open(self.log_file, "a") as log_f:
            log_f.write(message + "\n")
        # print(message)

    # -------------------------------
    # BUSCO collation and indexing
    # -------------------------------
    def transform_busco_results(
        self,
        busco_results: List[str],
        out_dir: Optional[str] = None,
        update: bool = False,
        force: bool = False,
        sequence_type: str = "protein",
        include_duplicated: bool = True,
    ) -> bool:
        """
        Collate complete BUSCO sequences (single-copy and multi-copy) from multiple
        BUSCO result directories
        into a single directory of family FASTAs. Each output FASTA is named <family>.fasta and
        contains records with IDs prefixed by species name (directory name above the run folder).
        """
        if out_dir is None:
            out_dir = os.path.join(self.working_dir, "busco_fastas")

        if os.path.exists(out_dir):
            if not force:
                self.log(f"Error: Output directory {out_dir} already exists. Use force=True to overwrite.")
                return False
            self.log(f"Warning: Output directory {out_dir} already exists. Overwriting.")
            shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)

        seq_kind = str(sequence_type or "protein").strip().lower()
        if seq_kind in {"nt", "nucleotide", "fna"}:
            valid_exts = (".fna", ".fa", ".fasta")
        else:
            valid_exts = (".faa", ".fa", ".fasta")
        sequence_subdirs = ["single_copy_busco_sequences"]
        if include_duplicated:
            sequence_subdirs.append("multi_copy_busco_sequences")

        for busco_result in busco_results:
            species_name = os.path.basename(os.path.dirname(busco_result))
            run_dirs = sorted(glob.glob(os.path.join(busco_result, "run*")))
            if not run_dirs:
                self.log(f"Warning: No BUSCO run dir in {busco_result}, skipping.")
                continue
            if len(run_dirs) > 1:
                self.log(f"Info: Multiple run dirs in {busco_result}, using: {run_dirs[0]}")
            run_dir = run_dirs[0]
            copied = 0
            for subdir in sequence_subdirs:
                busco_fasta_dir = os.path.join(run_dir, "busco_sequences", subdir)
                if not os.path.exists(busco_fasta_dir):
                    continue
                for fn in os.listdir(busco_fasta_dir):
                    if fn.endswith(valid_exts):
                        src_path = os.path.join(busco_fasta_dir, fn)
                        family = fn.rsplit(".", 1)[0]
                        dest_path = os.path.join(out_dir, f"{family}.fasta")
                        with open(dest_path, "a") as out_f:
                            for rec in SeqIO.parse(src_path, "fasta"):
                                rec.id = f"{species_name}_{rec.id}"
                                rec.description = ""
                                SeqIO.write(rec, out_f, "fasta")
                                copied += 1
            if copied == 0:
                self.log(f"Warning: No BUSCO sequence files found in {run_dir}, skipping.")

        if update:
            self.busco_dir = out_dir
        return True

    def _hash(self, s: str) -> str:
        return hashlib.md5(s.encode()).hexdigest()

    @dataclass(frozen=True)
    class _BuscoSequenceRecord:
        busco_file: str
        busco_id: str
        family: str
        length: int
        sequence_key: str
        protein_id: str

    @dataclass(frozen=True)
    class _OrthogroupSequenceRecord:
        orthogroup: str
        og_seq_id: str
        sequence_key: str
        sequence_length: int

    def _busco_protein_id(self, busco_id: str) -> str:
        parts = str(busco_id or "").strip().split("_", 2)
        if len(parts) == 3:
            return parts[2]
        return parts[-1] if parts else ""

    def _species_prefix_from_busco_id(self, busco_id: str) -> str:
        parts = str(busco_id).split('_')
        if len(parts) >= 2:
            return parts[0] + '_' + parts[1]
        return parts[0] if parts else ""

    def _extract_accession_token(self, token: str) -> str:
        text = str(token or "").strip()
        if not text:
            return ""
        match = re.match(r"^(GC[AF])_(\d+)[._](\d+)\b", text, flags=re.IGNORECASE)
        if match:
            prefix, digits, version = match.groups()
            return canonicalize_accession(f"{prefix}_{digits}.{version}")
        trimmed = re.sub(r"_(?:gff|cdhit|clean|raw)(?:[_-].*|\d.*)?$", "", text, flags=re.IGNORECASE)
        internal_match = re.match(r"^([A-Za-z0-9]+_\d{4})\b", trimmed)
        if internal_match:
            return canonicalize_accession(internal_match.group(1)) or internal_match.group(1)
        parts = trimmed.split("_")
        if len(parts) >= 2:
            return canonicalize_accession("_".join(parts[:2]))
        return canonicalize_accession(trimmed)

    def _load_leaf_metadata(self) -> Dict[str, Dict[str, str]]:
        if self._leaf_metadata_cache is not None:
            return self._leaf_metadata_cache

        working_dir = os.path.join(str(self.orthofinder_run_folder), "WorkingDirectory")
        species_ids_path = os.path.join(working_dir, "SpeciesIDs.txt")
        sequence_ids_path = os.path.join(working_dir, "SequenceIDs.txt")
        if not (os.path.exists(species_ids_path) and os.path.exists(sequence_ids_path)):
            self._leaf_metadata_cache = {}
            return self._leaf_metadata_cache

        species_map: Dict[str, str] = {}
        with open(species_ids_path, "r", encoding="utf-8") as handle:
            for raw in handle:
                if ":" not in raw:
                    continue
                species_idx, label = raw.split(":", 1)
                species_token = os.path.splitext(os.path.basename(label.strip().replace(".gz", "")))[0]
                accession = self._extract_accession_token(species_token)
                if accession:
                    species_map[str(species_idx).strip()] = accession

        leaf_map: Dict[str, Dict[str, str]] = {}
        with open(sequence_ids_path, "r", encoding="utf-8") as handle:
            for raw in handle:
                if ":" not in raw:
                    continue
                seq_id, leaf_name = raw.split(":", 1)
                species_idx = str(seq_id).split("_", 1)[0].strip()
                accession = species_map.get(species_idx, "")
                sequence = str(leaf_name).strip().split()[0]
                if accession and sequence:
                    metadata = {"accession": accession, "sequence": sequence}
                    leaf_map[sequence] = metadata

        self._leaf_metadata_cache = leaf_map
        return self._leaf_metadata_cache

    def _index_busco_sequences(self, use_hash: bool = True):
        """
        Build:
        - seq_index["by_hash"]: key (raw seq or hash) -> list of BUSCO records
        - seq_index["by_protein_id"]: protein ID suffix -> list of BUSCO records
        - family_to_seq_counts: family -> Counter(hash -> count)
        """
        if not self.busco_dir or not os.path.exists(self.busco_dir):
            self.log("Error: BUSCO directory not set or missing. Run transform_busco_results(update=True) or set busco_dir.")
            return None, None

        seq_index = {
            "by_hash": defaultdict(list),
            "by_protein_id": defaultdict(list),
        }
        family_to_seq_counts: Dict[str, Counter] = {}

        self.log(f"Collecting BUSCO sequences from {self.busco_dir}")
        for fn in os.listdir(self.busco_dir):
            if fn.endswith((".fa", ".faa", ".fasta")):
                # Get family name from the filename
                family = fn.rsplit(".", 1)[0]
                # Set up counter for this family
                fam_counter = family_to_seq_counts.setdefault(family, Counter())
                for rec in SeqIO.parse(os.path.join(self.busco_dir, fn), "fasta"):
                    seq = str(rec.seq).upper()
                    key = self._hash(seq) if use_hash else seq
                    busco_record = self._BuscoSequenceRecord(
                        busco_file=fn,
                        busco_id=str(rec.id),
                        family=family,
                        length=len(seq),
                        sequence_key=key,
                        protein_id=self._busco_protein_id(str(rec.id)),
                    )
                    seq_index["by_hash"][key].append(busco_record)
                    seq_index["by_protein_id"][busco_record.protein_id].append(busco_record)
                    fam_counter[key] += 1

        self.log(
            f"Collected {len(seq_index['by_hash'])} unique BUSCO sequences "
            f"({'hashed' if use_hash else 'raw sequences'})."
        )
        return seq_index, family_to_seq_counts

    def _write_exact_mapping_tsv(self, source_mapping_tsv: str, exact_mapping_tsv: str, exact_pairs: List[Tuple[str, str]]) -> None:
        exact_pair_set = {(str(fam), str(og)) for fam, og in exact_pairs}
        with open(source_mapping_tsv, "r", newline="") as in_f, open(exact_mapping_tsv, "w", newline="") as out_f:
            reader = csv.DictReader(in_f, delimiter="\t")
            fieldnames = list(reader.fieldnames) if reader.fieldnames else []
            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                pair = (str(row.get("BUSCO_family") or ""), str(row.get("Orthogroup") or ""))
                if pair in exact_pair_set:
                    writer.writerow(row)

    # -------------------------------
    # Orthogroup mapping to BUSCO
    # -------------------------------
    def _scan_orthogroups(self, seq_index: dict, use_hash: bool, mapping_tsv: str):
        """
        Map OG sequence records to BUSCO families using:
        1. BUSCO protein ID suffix -> OG sequence ID exact matches.
        2. Sequence-key fallback for unmatched records only, paired 1:1.
        Writes the per-sequence mapping TSV. Returns:
        - og_to_busco_families: og -> set of families present in og
        - busco_family_to_ogs: family -> set of ogs
        - og_family_seq_counts: og -> family -> Counter(hash -> count) for OG-BUSCO matched sequences only
    - og_total_seq_count: og -> total number of sequences in OG (for exact equality checks)
    - og_family_species: og -> family -> set(species)
        """
        og_dir = os.path.join(self.orthofinder_run_folder, "Orthogroup_Sequences")
        if not os.path.exists(og_dir):
            self.log(f"Error: Orthogroup directory {og_dir} does not exist.")
            return None, None, None, None, None

        og_to_busco_families: Dict[str, Set[str]] = defaultdict(set)
        busco_family_to_ogs: Dict[str, Set[str]] = defaultdict(set)
        og_family_seq_counts: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
        og_total_seq_count: Dict[str, int] = defaultdict(int)
        og_family_species: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))

        mapped_busco_ids: Set[str] = set()
        busco_by_hash = seq_index.get("by_hash", {})
        busco_by_protein_id = seq_index.get("by_protein_id", {})

        def write_match(
            writer: csv.writer,
            busco_record: OrthoBuscoAnalyzer._BuscoSequenceRecord,
            og_record: OrthoBuscoAnalyzer._OrthogroupSequenceRecord,
        ) -> None:
            writer.writerow(
                [
                    busco_record.busco_file,
                    busco_record.busco_id,
                    busco_record.family,
                    og_record.orthogroup,
                    og_record.og_seq_id,
                    busco_record.length,
                ]
            )
            og_to_busco_families[og_record.orthogroup].add(busco_record.family)
            busco_family_to_ogs[busco_record.family].add(og_record.orthogroup)
            og_family_seq_counts[og_record.orthogroup][busco_record.family][og_record.sequence_key] += 1
            mapped_busco_ids.add(str(busco_record.busco_id))
            sp = self._species_prefix_from_busco_id(busco_record.busco_id)
            og_family_species[og_record.orthogroup][busco_record.family].add(sp)

        with open(mapping_tsv, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["BUSCO_file", "BUSCO_id", "BUSCO_family", "Orthogroup", "OG_seq_id", "Sequence_length"])
            for fn in os.listdir(og_dir):
                if not fn.endswith((".fa", ".faa", ".fasta")):
                    continue
                og_name = fn.rsplit('.', 1)[0]
                og_path = os.path.join(og_dir, fn)
                og_records: List[OrthoBuscoAnalyzer._OrthogroupSequenceRecord] = []
                for rec in SeqIO.parse(og_path, "fasta"):
                    og_total_seq_count[og_name] += 1
                    seq = str(rec.seq).upper()
                    key = self._hash(seq) if use_hash else seq
                    og_records.append(
                        self._OrthogroupSequenceRecord(
                            orthogroup=og_name,
                            og_seq_id=str(rec.id),
                            sequence_key=key,
                            sequence_length=len(seq),
                        )
                    )

                matched_busco_ids: Set[str] = set()
                matched_og_indexes: Set[int] = set()
                matched_families: Set[str] = set()

                # Stage 1: preserve assembly-specific BUSCO pairings whenever IDs agree.
                for index, og_record in enumerate(og_records):
                    candidates = [
                        candidate
                        for candidate in busco_by_protein_id.get(og_record.og_seq_id, [])
                        if candidate.busco_id not in matched_busco_ids
                    ]
                    if not candidates:
                        continue

                    exact_seq_candidates = [
                        candidate for candidate in candidates if candidate.sequence_key == og_record.sequence_key
                    ]
                    chosen = sorted(
                        exact_seq_candidates or candidates,
                        key=lambda candidate: (candidate.family, candidate.busco_id),
                    )[0]
                    write_match(w, chosen, og_record)
                    matched_busco_ids.add(chosen.busco_id)
                    matched_og_indexes.add(index)
                    matched_families.add(chosen.family)

                # Stage 2: limited sequence fallback for unmatched records only.
                unmatched_by_hash: Dict[str, List[Tuple[int, OrthoBuscoAnalyzer._OrthogroupSequenceRecord]]] = defaultdict(list)
                for index, og_record in enumerate(og_records):
                    if index not in matched_og_indexes:
                        unmatched_by_hash[og_record.sequence_key].append((index, og_record))

                for sequence_key in sorted(unmatched_by_hash):
                    unmatched_ogs = sorted(
                        unmatched_by_hash[sequence_key],
                        key=lambda item: item[1].og_seq_id,
                    )
                    candidate_records = [
                        candidate
                        for candidate in busco_by_hash.get(sequence_key, [])
                        if candidate.busco_id not in matched_busco_ids
                    ]
                    if not candidate_records:
                        continue

                    candidates_by_family: Dict[str, List[OrthoBuscoAnalyzer._BuscoSequenceRecord]] = defaultdict(list)
                    for candidate in candidate_records:
                        candidates_by_family[candidate.family].append(candidate)

                    if matched_families:
                        eligible_families = sorted(
                            family for family in candidates_by_family if family in matched_families
                        )
                    elif len(candidates_by_family) == 1:
                        eligible_families = sorted(candidates_by_family)
                    else:
                        eligible_families = []

                    for family in eligible_families:
                        family_candidates = sorted(
                            candidates_by_family[family],
                            key=lambda candidate: (candidate.busco_id, candidate.protein_id),
                        )
                        pair_count = min(len(unmatched_ogs), len(family_candidates))
                        for pair_index in range(pair_count):
                            og_index, og_record = unmatched_ogs[pair_index]
                            busco_record = family_candidates[pair_index]
                            write_match(w, busco_record, og_record)
                            matched_busco_ids.add(busco_record.busco_id)
                            matched_og_indexes.add(og_index)
                            matched_families.add(family)
                        unmatched_ogs = unmatched_ogs[pair_count:]
                        if not unmatched_ogs:
                            break

        total_mapped = sum(sum(cnt.values()) for fam_map in og_family_seq_counts.values() for cnt in fam_map.values())
        self.log(f"Total BUSCO sequence mappings written: {total_mapped}")
        return (
            og_to_busco_families,
            busco_family_to_ogs,
            og_family_seq_counts,
            og_total_seq_count,
            og_family_species,
            mapped_busco_ids,
        )

    def _summarize_unmapped_buscos(
        self,
        mapped_busco_ids: Set[str],
        out_prefix: str,
    ) -> Dict[str, object]:
        if not self.busco_dir or not os.path.isdir(self.busco_dir):
            return {
                "total_busco_records": 0,
                "mapped_busco_records": 0,
                "unmapped_busco_records": 0,
                "per_species_unmapped_counts": {},
                "unmapped_families": {},
                "summary_tsv": "",
            }

        total_busco_records = 0
        unmapped_busco_records = 0
        per_species_unmapped_counts: Counter[str] = Counter()
        unmapped_families: Dict[str, List[str]] = defaultdict(list)

        for fn in os.listdir(self.busco_dir):
            if not fn.endswith((".fa", ".faa", ".fasta")):
                continue
            family = fn.rsplit(".", 1)[0]
            fasta_path = os.path.join(self.busco_dir, fn)
            for rec in SeqIO.parse(fasta_path, "fasta"):
                total_busco_records += 1
                busco_id = str(rec.id)
                if busco_id in mapped_busco_ids:
                    continue
                unmapped_busco_records += 1
                species = self._species_prefix_from_busco_id(busco_id)
                if species:
                    per_species_unmapped_counts[species] += 1
                unmapped_families[family].append(busco_id)

        summary_tsv = out_prefix + "_unmapped_buscos.tsv"
        with open(summary_tsv, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Scope", "Key", "Unmapped_BUSCOs"])
            writer.writerow(["total", "all_species", unmapped_busco_records])
            for species in sorted(per_species_unmapped_counts):
                writer.writerow(["species", species, per_species_unmapped_counts[species]])

        self.log(
            "BUSCOs found only in OrthoFinder unassigned genes are counted as unmapped here, "
            "because core-set mapping scans assigned Orthogroup_Sequences only."
        )
        self.log(f"Total BUSCO records examined: {total_busco_records}")
        self.log(f"Mapped BUSCO records: {len(mapped_busco_ids)}")
        self.log(f"Unmapped BUSCO records: {unmapped_busco_records}")
        self.log(f"Unmapped BUSCO summary TSV written: {summary_tsv}")
        if per_species_unmapped_counts:
            self.log("Unmapped BUSCO counts by species:")
            for species in sorted(per_species_unmapped_counts):
                self.log(f"  {species}: {per_species_unmapped_counts[species]}")

        families_with_unmapped = {
            family: sorted(busco_ids) for family, busco_ids in unmapped_families.items() if busco_ids
        }
        if families_with_unmapped:
            self.log(f"BUSCO families with at least one unmapped BUSCO record: {len(families_with_unmapped)}")

        return {
            "total_busco_records": total_busco_records,
            "mapped_busco_records": len(mapped_busco_ids),
            "unmapped_busco_records": unmapped_busco_records,
            "per_species_unmapped_counts": dict(per_species_unmapped_counts),
            "unmapped_families": families_with_unmapped,
            "summary_tsv": summary_tsv,
        }

    def _write_1to1_report(self, one_to_one_tsv: str, og_to_busco_families: dict, busco_family_to_ogs: dict,
                           og_family_seq_counts: dict, og_total_seq_count: dict, family_to_seq_counts: dict,
                           require_no_extra_sequences: bool = False):
        """
        Compute and write 1-1 statuses:
        - exact 1-1: OG contains exactly the BUSCO family sequences (multiset equality) and no extra sequences.
        - relaxed unique: BUSCO family maps to exactly one OG (regardless of other content).
        Returns a tuple:
            (exact_pairs: list[(family, og)], relaxed_map: dict[family -> og])
        """
        relaxed_one_to_one = {fam: list(ogs)[0] for fam, ogs in busco_family_to_ogs.items() if len(ogs) == 1}
        exact_pairs: List[Tuple[str, str]] = []

        for fam, og in relaxed_one_to_one.items():
            # Exact BUSCO-only equality: OG contains BUSCOs from one family and the BUSCO-matched sequences match the family multiset
            if len(og_to_busco_families.get(og, set())) != 1:
                continue
            og_busco_counts = og_family_seq_counts.get(og, {}).get(fam, Counter())
            fam_counts = family_to_seq_counts.get(fam, Counter())
            if not fam_counts:
                continue
            if require_no_extra_sequences:
                # Ultra-strict: the OG must contain no non-BUSCO sequences at all
                if og_total_seq_count.get(og, 0) != sum(og_busco_counts.values()):
                    continue
            # Multiset equality on BUSCO-matched sequences
            if og_busco_counts == fam_counts:
                exact_pairs.append((fam, og))

        with open(one_to_one_tsv, "w", newline="") as out2:
            w2 = csv.writer(out2, delimiter="\t")
            w2.writerow(["BUSCO_family", "Orthogroup", "Exact_1_to_1", "Relaxed_only"])  # mirrors older style
            exact_set = {f for f, _ in exact_pairs}
            for fam, og in relaxed_one_to_one.items():
                exact = "YES" if fam in exact_set else "NO"
                relaxed_only = "YES" if fam not in exact_set else "NO"
                w2.writerow([fam, og, exact, relaxed_only])

        self.log(f"Strict 1-1 families: {len(exact_pairs)}; Relaxed (unique OG) families: {len(relaxed_one_to_one)}")
        return exact_pairs, relaxed_one_to_one

    # -------------------------------
    # Paralog analysis
    # -------------------------------
    def _parse_tree(self, newick: str):
        try:
            return Tree(newick, format=1)
        except Exception as exc:  # boundary: try alternate ETE tree parser format.
            self.log(f"parse_tree: format=1 parse failed, trying default parser: {exc}")
            try:
                return Tree(newick)
            except Exception as e:  # boundary: malformed tree input is reported and skipped by caller.
                self.log(f"parse_tree: Failed to parse tree: {e}")
                return None

    def _header_accession(self, header: str) -> str:
        token = str(header or "").strip()
        if not token:
            return ""
        match = re.match(r"^(GC[AF])_(\d+)_(\d+)(?:_|$)", token, flags=re.IGNORECASE)
        if match:
            prefix, digits, version = match.groups()
            return canonicalize_accession(f"{prefix.upper()}_{digits}.{version}")
        internal_match = re.match(r"^([A-Za-z0-9]+_\d{4})(?:_|$)", token)
        if internal_match:
            return canonicalize_accession(internal_match.group(1)) or internal_match.group(1)
        parts = token.split('_')
        if len(parts) >= 2:
            return canonicalize_accession('_'.join(parts[:2]))
        return canonicalize_accession(parts[0])

    def _species_key(self, leaf_name: str) -> str:
        metadata = self._load_leaf_metadata().get(str(leaf_name).strip())
        if metadata and metadata.get("accession"):
            return str(metadata["accession"])
        return self._header_accession(leaf_name)

    def _species_tuple_accession(self, leaves: List[str]) -> str:
        if not leaves:
            return ""
        return self._header_accession(leaves[0])

    def _paralogs_for_tree(self, og_name: str, tree: Tree, out_dir: str):
        species_to_leaves: Dict[str, List[str]] = defaultdict(list)
        for leaf in tree.get_leaf_names():
            species_to_leaves[self._species_key(leaf)].append(leaf)

        inparalogs: List[str] = []
        outparalogs: List[str] = []
        in_species: Set[str] = set()
        out_species: Set[str] = set()
        for sp, leaves in species_to_leaves.items():
            if len(leaves) > 1:
                try:
                    is_mono, _, _ = tree.check_monophyly(values=leaves, target_attr="name", unrooted=True)
                except Exception as e:  # boundary: one monophyly check failure marks that OG non-monophyletic.
                    is_mono = False
                    self.log(f"identify_paralogs: Failed to check monophyly for {og_name}: {e}")
                if is_mono:
                    inparalogs.append(str(tuple(leaves)))
                    if sp:
                        in_species.add(sp)
                else:
                    outparalogs.append(str(tuple(leaves)))
                    if sp:
                        out_species.add(sp)

        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            in_file = os.path.join(out_dir, f"{og_name}_inparalogs.txt")
            out_file = os.path.join(out_dir, f"{og_name}_outparalogs.txt")
            if inparalogs:
                with open(in_file, 'w') as f:
                    f.write("\n".join(inparalogs) + "\n")
            if outparalogs:
                with open(out_file, 'w') as f:
                    f.write("\n".join(outparalogs) + "\n")

        return inparalogs, outparalogs, len(species_to_leaves), in_species, out_species

    def _ensure_resolved_gene_trees_file(self) -> Optional[str]:
        """
        Resolve the consolidated Resolved_Gene_Trees.txt path expected by
        paralog parsing code. If missing, synthesize it from per-OG *_tree.txt
        files in the configured gene tree directory.
        """
        tree_dir = str(self.gene_tree_dir)
        trees_file = os.path.join(tree_dir, "Resolved_Gene_Trees.txt")
        if os.path.exists(trees_file):
            return trees_file

        if not os.path.isdir(tree_dir):
            self.log(f"identify_paralogs: Gene tree directory not found: {tree_dir}")
            return None

        per_og_files = sorted(glob.glob(os.path.join(tree_dir, "*_tree.txt")))
        if not per_og_files:
            self.log(
                "identify_paralogs: No consolidated or per-OG resolved tree files found "
                f"in {tree_dir}"
            )
            return None

        try:
            with open(trees_file, "w") as out:
                written = 0
                for tree_path in per_og_files:
                    base = os.path.basename(tree_path)
                    og_name = base[:-9] if base.endswith("_tree.txt") else os.path.splitext(base)[0]
                    with open(tree_path, "r") as fh:
                        tree_newick = fh.read().strip()
                    if not tree_newick:
                        continue
                    out.write(f"{og_name}:{tree_newick}\n")
                    written += 1
            self.log(
                "identify_paralogs: Built fallback Resolved_Gene_Trees.txt "
                f"from {written} per-OG tree files."
            )
            return trees_file
        except (OSError, UnicodeError) as exc:
            self.log(f"identify_paralogs: Failed building fallback resolved tree file: {exc}")
            return None

    def _identify_paralogs(self, matched_exact_ogs: Set[str], min_species: int, write_tsv: bool):
        """
        Parse Resolved_Gene_Trees.txt and compute paralog summaries:
        - Overall across all OGs (for the logs)
        - For matched_exact_ogs subset separately (for the logs)
        Returns:
            dict with keys:
                above_species_threshold (dict og->n)
                below_species_threshold (dict og->n)
                og_has_inparalogs (list og)
                og_has_outparalogs (list og)
                og_has_no_paralogs (list og)
                matched* variants filtered to matched_exact_ogs
        Also writes og_paralog_summary.tsv when write_tsv is True.
        """
        trees_file = self._ensure_resolved_gene_trees_file()
        if not trees_file:
            return None

        paralog_dir = os.path.join(self.working_dir, "paralogs")
        os.makedirs(paralog_dir, exist_ok=True)

        og_has_inparalogs: List[str] = []
        og_has_outparalogs: List[str] = []
        og_has_no_paralogs: List[str] = []
        above_species_threshold: Dict[str, int] = {}
        below_species_threshold: Dict[str, int] = {}
        matched_ogs: Set[str] = set()
        og_species_inparalogs: Dict[str, Set[str]] = {}
        og_species_outparalogs: Dict[str, Set[str]] = {}

        processed_count = 0

        with open(trees_file, 'r') as fh:
            for raw in fh:
                line = raw.strip()
                if not line or ':' not in line:
                    continue
                og_name, tree_newick = line.split(':', 1)
                og_name = og_name.strip()
                tree_newick = tree_newick.strip()
                if og_name in matched_exact_ogs:
                    matched_ogs.add(og_name)

                tree = self._parse_tree(tree_newick)
                if not tree:
                    self.log(f"identify_paralogs: Skipping {og_name} due to tree parse failure.")
                    continue

                inpars, outpars, sp_count, in_species, out_species = self._paralogs_for_tree(og_name, tree, out_dir=paralog_dir)
                og_species_inparalogs[og_name] = set(in_species)
                og_species_outparalogs[og_name] = set(out_species)

                if inpars:
                    og_has_inparalogs.append(og_name)
                if outpars:
                    og_has_outparalogs.append(og_name)
                if not inpars and not outpars:
                    og_has_no_paralogs.append(og_name)

                if sp_count >= int(min_species):
                    above_species_threshold[og_name] = sp_count
                else:
                    below_species_threshold[og_name] = sp_count

                processed_count += 1

        if processed_count == 0:
            self.log("identify_paralogs: No orthogroup trees processed. Check your mapping and tree file.")
        else:
            self.log(
                f"identify_paralogs: processed={processed_count}, matched_exact={len(matched_ogs)}, "
                f"above_min_species={len(above_species_threshold)}, below_min_species={len(below_species_threshold)}, "
                f"with_inparalogs={len(og_has_inparalogs)}, with_outparalogs={len(og_has_outparalogs)}, "
                f"no_paralogs={len(og_has_no_paralogs)}"
            )

        og_paralog_classification: Dict[str, str] = {}
        og_paralog_counts: Dict[str, Tuple[int, int]] = {}
        if write_tsv:
            tsv_path = os.path.join(paralog_dir, "og_paralog_summary.tsv")
            with open(tsv_path, 'w', newline='') as tsv_f:
                writer = csv.writer(tsv_f, delimiter='\t')
                writer.writerow(["Orthogroup", "In_Paralog_Sets", "Out_Paralog_Sets", "Classification"])
                # Only summarize matched_exact_ogs for compactness like older implementation
                for og in matched_exact_ogs:
                    in_file = os.path.join(paralog_dir, f"{og}_inparalogs.txt")
                    out_file = os.path.join(paralog_dir, f"{og}_outparalogs.txt")
                    in_count = 0
                    out_count = 0
                    if os.path.exists(in_file):
                        with open(in_file) as f:
                            in_count = sum(1 for _ in f if _.strip())
                    if os.path.exists(out_file):
                        with open(out_file) as f:
                            out_count = sum(1 for _ in f if _.strip())
                    if out_count == 0 and in_count > 0:
                        classification = "Only In-Paralogs"
                    elif out_count > 0 and in_count == 0:
                        classification = "Only Out-Paralogs"
                    elif out_count > 0 and in_count > 0:
                        classification = "Both In- and Out-Paralogs"
                    else:
                        classification = "No Paralogs"
                    og_paralog_classification[og] = classification
                    og_paralog_counts[og] = (in_count, out_count)
                    writer.writerow([og, in_count, out_count, classification])
            self.log(f"Paralog summary written to {tsv_path}")

        # Return a compact dict of essentials
        return {
            "above_species_threshold": above_species_threshold,
            "below_species_threshold": below_species_threshold,
            "og_has_inparalogs": og_has_inparalogs,
            "og_has_outparalogs": og_has_outparalogs,
            "og_has_no_paralogs": og_has_no_paralogs,
            "og_paralog_classification": og_paralog_classification,
            "og_paralog_counts": og_paralog_counts,
            "og_species_inparalogs": og_species_inparalogs,
            "og_species_outparalogs": og_species_outparalogs,
        }

    # -------------------------------
    # Orchestration
    # -------------------------------
    def compare_busco_orthofinder(self, strict: bool = False, min_species_in_trees: int = 1, out_prefix: Optional[str] = None, use_hash: bool = True, busco_results: Optional[List[str]] = None, require_no_extra_sequences: bool = False):
        """
        Orchestrates BUSCO/OG mapping, 1-1 detection, paralog analysis and writes reports.

        Returns the list of BUSCO families that map to 'good' OGs defined as:
          - Exact 1-1 BUSCO<->OG (OG contains only those sequences)
          - OG has no out-paralogs (it may have in-paralogs or none)
          - Optionally above min_species_in_trees threshold in the summary/logging
        """
        self.log(f"compare_busco_orthofinder: Using OrthoFinder results at {self.orthofinder_run_folder}")
        if out_prefix is None:
            out_prefix = os.path.join(self.working_dir, f"{self.identifier}_busco_to_orthogroup")
        mapping_tsv = out_prefix + "_map.tsv"
        one_to_one_tsv = out_prefix + "_1to1.tsv"

        # Ensure BUSCO directory is available (collate if needed and locations provided)
        if not self.busco_dir:
            if busco_results:
                self.log("BUSCO directory not set; collating BUSCO results from provided locations...")
                ok = self.transform_busco_results(busco_results, out_dir=None, update=True, force=False)
                if not ok:
                    return []
            else:
                self.log("Error: BUSCO directory not set. Provide busco_results to collate or set busco_dir.")
                return []

        # Index BUSCO sequences
        seq_index, family_to_seq_counts = self._index_busco_sequences(use_hash=use_hash)
        if seq_index is None:
            return []

        # Map OG sequences to BUSCOs
        self.log(f"Mapping Orthogroups in {os.path.join(self.orthofinder_run_folder, 'Orthogroup_Sequences')} to BUSCO sequences.")
        (
            og_to_busco_families,
            busco_family_to_ogs,
            og_family_seq_counts,
            og_total_seq_count,
            og_family_species,
            mapped_busco_ids,
        ) = self._scan_orthogroups(
            seq_index=seq_index,
            use_hash=use_hash,
            mapping_tsv=mapping_tsv,
        )
        if og_to_busco_families is None:
            return []
        self.log(f"Mapping completed. Outputs: {mapping_tsv}, {one_to_one_tsv}")
        unmapped_summary = self._summarize_unmapped_buscos(mapped_busco_ids, out_prefix)

        # Write explicit OG<->BUSCO mapping TSVs
        og_to_busco_tsv = out_prefix + "_og_to_busco_families.tsv"
        busco_to_og_tsv = out_prefix + "_busco_family_to_ogs.tsv"
        with open(og_to_busco_tsv, 'w', newline='') as f1:
            w1 = csv.writer(f1, delimiter='\t')
            w1.writerow(["Orthogroup", "BUSCO_family"])  # one pair per row
            for og in sorted(og_to_busco_families.keys()):
                for fam in sorted(og_to_busco_families[og]):
                    w1.writerow([og, fam])
        with open(busco_to_og_tsv, 'w', newline='') as f2:
            w2 = csv.writer(f2, delimiter='\t')
            w2.writerow(["BUSCO_family", "Orthogroup"])  # one pair per row
            for fam in sorted(busco_family_to_ogs.keys()):
                for og in sorted(busco_family_to_ogs[fam]):
                    w2.writerow([fam, og])
        self.log(f"Mapping pair TSVs written: {og_to_busco_tsv}, {busco_to_og_tsv}")

        # 1-1 computation and report
        exact_pairs, relaxed_map = self._write_1to1_report(
            one_to_one_tsv=one_to_one_tsv,
            og_to_busco_families=og_to_busco_families,
            busco_family_to_ogs=busco_family_to_ogs,
            og_family_seq_counts=og_family_seq_counts,
            og_total_seq_count=og_total_seq_count,
            family_to_seq_counts=family_to_seq_counts,
            require_no_extra_sequences=require_no_extra_sequences,
        )

        # Detailed mapping summaries
        multi_busco_ogs = [og for og, fams in og_to_busco_families.items() if len(fams) > 1]
        self.log(f"Orthogroups containing multiple BUSCO families: {len(multi_busco_ogs)}")

        multi_og_buscos = [fam for fam, ogs in busco_family_to_ogs.items() if len(ogs) > 1]
        self.log(f"BUSCO families mapped to multiple Orthogroups: {len(multi_og_buscos)}")

        # Compute BUSCO occupancy (number of species with that BUSCO present in the OG) for exact 1-1s
        occupancy_records = []
        exact_pairs_filtered = []
        below_threshold_families = []
        for fam, og in exact_pairs:
            species_set = og_family_species.get(og, {}).get(fam, set())
            occ = len(species_set)
            occupancy_records.append((fam, og, occ))
            if occ >= min_species_in_trees:
                exact_pairs_filtered.append((fam, og))
            else:
                below_threshold_families.append((fam, og, occ))

        # Log and write occupancy summary
        self.log(f"Exact 1-1 BUSCO families total: {len(exact_pairs)}")
        self.log(f"Exact 1-1 families with >= min_species ({min_species_in_trees}): {len(exact_pairs_filtered)}")
        self.log(f"Exact 1-1 families rejected due to low occupancy: {len(below_threshold_families)}")
        occupancy_tsv = out_prefix + "_1to1_occupancy.tsv"
        with open(occupancy_tsv, 'w', newline='') as occ_f:
            w = csv.writer(occ_f, delimiter='\t')
            w.writerow(["BUSCO_family", "Orthogroup", "Species_count", "Status"])
            accepted = { (fam, og) for fam, og in exact_pairs_filtered }
            for fam, og, occ in occupancy_records:
                status = "accepted" if (fam, og) in accepted else "rejected"
                w.writerow([fam, og, occ, status])
        self.log(f"1-1 occupancy TSV written: {occupancy_tsv}")

        if below_threshold_families:
            self.log(f"BUSCO families rejected due to low occupancy: {len(below_threshold_families)}")

        # Paralog reporting should see all exact 1-1 OGs so the "regardless of min_species"
        # summary can differ from the thresholded summary. Final family selection still uses
        # only the occupancy-passing subset below.
        matched_exact_ogs = {og for _, og in exact_pairs}
        matched_exact_ogs_filtered = {og for _, og in exact_pairs_filtered}
        paralog_results = self._identify_paralogs(matched_exact_ogs=matched_exact_ogs, min_species=min_species_in_trees, write_tsv=True)
        if paralog_results is None:
            return []

        # Determine 'good' OGs: no out-paralogs (in-paralogs allowed) among matched exact 1-1 OGs
        og_has_out = set(paralog_results["og_has_outparalogs"])
        og_has_in = set(paralog_results["og_has_inparalogs"])  # not directly used but informative
        og_has_none = set(paralog_results["og_has_no_paralogs"])  # across all OGs
        good_ogs = {og for og in matched_exact_ogs_filtered if og not in og_has_out}

        # missing = matched_exact_ogs - good_ogs
        # self.log(f"Final 'good' Orthogroups (exact 1-1 and no out-paralogs): {len(good_ogs)}")
        # self.log(f"Orthogroups excluded due to out-paralogs: {len(missing)}")
        # for og in sorted(missing):
        #     self.log(f"  EXCLUDED OG {og}: BUSCO family {next((fam for fam, o in exact_pairs if o == og), 'N/A')}")

        # Final BUSCO family list mapped to these 'good' OGs (must be exact 1-1 by construction)
        fam_to_og_exact = {fam: og for fam, og in exact_pairs_filtered}
        good_families = sorted([fam for fam, og in fam_to_og_exact.items() if og in good_ogs])

        # Write final list
        final_list_path = os.path.join(self.working_dir, f"{self.identifier}_good_busco_families.txt")
        with open(final_list_path, 'w') as fh:
            for fam in good_families:
                fh.write(fam + "\n")
        self.log(f"Wrote final list of good BUSCO families ({len(good_families)}): {final_list_path}")

        exact_mapping_tsv = out_prefix + "_exact_map.tsv"
        self._write_exact_mapping_tsv(mapping_tsv, exact_mapping_tsv, exact_pairs)
        self.log(f"Exact BUSCO mapping written: {exact_mapping_tsv}")

        # Write augmented mapping with paralog classification + clean (exact 1-1, occupancy-pass, no paralogs)
        paralog_classification = paralog_results.get("og_paralog_classification", {})
        clean_pairs = {
            (fam, og)
            for fam, og in exact_pairs_filtered
            if paralog_classification.get(og) in ("No Paralogs", "Only In-Paralogs")
        }
        augmented_mapping_tsv = out_prefix + "_map_with_paralog_class.tsv"
        with open(mapping_tsv, "r", newline="") as in_f, open(augmented_mapping_tsv, "w", newline="") as out_f:
            reader = csv.DictReader(in_f, delimiter="\t")
            base_fields = list(reader.fieldnames) if reader.fieldnames else []
            fieldnames = base_fields + ["OG_Paralog_Classification", "Clean_1to1_NoParalogs"]
            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                og = row.get("Orthogroup", "")
                fam = row.get("BUSCO_family", "")
                row["OG_Paralog_Classification"] = paralog_classification.get(og, "NA")
                row["Clean_1to1_NoParalogs"] = "YES" if (fam, og) in clean_pairs else "NO"
                writer.writerow(row)
        self.log(f"Augmented mapping written: {augmented_mapping_tsv}")

        family_species_paralog_status: Dict[str, Dict[str, Dict[str, bool]]] = {}
        og_species_inparalogs = paralog_results.get("og_species_inparalogs", {}) or {}
        og_species_outparalogs = paralog_results.get("og_species_outparalogs", {}) or {}
        for fam, og in exact_pairs:
            species_rows: Dict[str, Dict[str, bool]] = {}
            for accession in sorted(set(og_species_inparalogs.get(og, set())) | set(og_species_outparalogs.get(og, set()))):
                species_rows[str(accession)] = {
                    "has_inparalogs": str(accession) in set(og_species_inparalogs.get(og, set())),
                    "has_outparalogs": str(accession) in set(og_species_outparalogs.get(og, set())),
                }
            family_species_paralog_status[str(fam)] = species_rows

        species_tsv = out_prefix + "_species_paralog_status.tsv"
        with open(species_tsv, "w", newline="") as species_f:
            writer = csv.writer(species_f, delimiter="\t")
            writer.writerow(["BUSCO_family", "Orthogroup", "Accession", "Has_In_Paralogs", "Has_Out_Paralogs"])
            for fam, og in sorted(exact_pairs):
                species_rows = family_species_paralog_status.get(str(fam), {})
                for accession in sorted(species_rows):
                    row = species_rows[accession]
                    writer.writerow([
                        fam,
                        og,
                        accession,
                        "YES" if row.get("has_inparalogs") else "NO",
                        "YES" if row.get("has_outparalogs") else "NO",
                    ])
        self.log(f"Species-specific paralog mapping written: {species_tsv}")

        self.last_compare_results = {
            "good_families": list(good_families),
            "exact_pairs": list(exact_pairs),
            "exact_pairs_filtered": list(exact_pairs_filtered),
            "family_to_orthogroup_exact": {str(fam): str(og) for fam, og in exact_pairs},
            "family_to_orthogroups_all": {
                str(fam): sorted(str(og) for og in ogs)
                for fam, ogs in busco_family_to_ogs.items()
            },
            "family_species_paralog_status": family_species_paralog_status,
            "unmapped_busco_summary": unmapped_summary,
            "paralog_results": paralog_results,
            "files": {
                "mapping_tsv": exact_mapping_tsv,
                "all_mapping_tsv": mapping_tsv,
                "exact_mapping_tsv": exact_mapping_tsv,
                "one_to_one_tsv": one_to_one_tsv,
                "augmented_mapping_tsv": augmented_mapping_tsv,
                "species_paralog_tsv": species_tsv,
                "final_list_path": final_list_path,
                "unmapped_busco_summary_tsv": unmapped_summary.get("summary_tsv", ""),
            },
        }

        # Return BUSCO families list as requested
        return good_families
