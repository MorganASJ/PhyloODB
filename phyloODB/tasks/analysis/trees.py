from __future__ import annotations

import csv
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from Bio import Phylo, SeqIO
from ete3 import Tree

from ..misc.export_library import ExportLibraryTask
from ..task import Task
from ...accession_utils import canonicalize_accession

DEFAULT_MAFFT_TASK_THREADS = 2
DEFAULT_IQTREE_TASK_THREADS = 4


def _resolve_executable(raw_path: object, fallback_names: Iterable[str]) -> Optional[str]:
    candidate = str(raw_path or "").strip()
    if candidate:
        expanded = os.path.expanduser(candidate)
        if os.path.isfile(expanded) and os.access(expanded, os.X_OK):
            return expanded
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    for name in fallback_names:
        resolved = shutil.which(name)
        if resolved:
            return resolved
    return None


def _flags_from_env_or_task(*, env_value: object, task_value: object, default_flags: list[str]) -> list[str]:
    raw = str(task_value or "").strip() or str(env_value or "").strip()
    if not raw:
        return list(default_flags)
    return shlex.split(raw)


def _ensure_thread_flag(
    flags: list[str],
    *,
    thread_tokens: set[str],
    inject_tokens: list[str],
    thread_count: int,
) -> list[str]:
    if any(token in thread_tokens for token in flags):
        return list(flags)
    return [*inject_tokens, str(max(1, int(thread_count))), *flags]


def valid_mafft_alignment(path: str, input_fasta: Optional[str] = None) -> bool:
    """Return whether *path* is a complete FASTA alignment.

    When the source FASTA is available, also require the aligned output to contain
    the same record IDs. This prevents a truncated-but-parseable file from being
    adopted after a machine or filesystem interruption.
    """
    if not path or not os.path.isfile(path) or os.path.getsize(path) <= 0:
        return False
    try:
        records = list(SeqIO.parse(path, "fasta"))
        if not records:
            return False
        lengths = {len(record.seq) for record in records}
        if len(lengths) != 1 or not next(iter(lengths)):
            return False
        output_ids = [record.id for record in records]
        if len(output_ids) != len(set(output_ids)):
            return False
        if input_fasta:
            if not os.path.isfile(input_fasta):
                return False
            input_ids = [record.id for record in SeqIO.parse(input_fasta, "fasta")]
            if not input_ids or sorted(input_ids) != sorted(output_ids):
                return False
    except (OSError, ValueError):
        return False
    return True


def valid_iqtree_tree(path: str) -> bool:
    """Return whether *path* contains a complete, parseable tree."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) <= 0:
        return False
    if Path(path).suffix.lower() in {".nex", ".nexus"}:
        try:
            tree = Phylo.read(path, "nexus")
        except Exception:
            return False
        return len(tree.get_terminals()) >= 2
    try:
        with open(path, "r", encoding="utf-8") as handle:
            newick = handle.read().strip()
    except (OSError, UnicodeError):
        return False
    if not newick or not newick.endswith(";"):
        return False
    for tree_format in (0, 1):
        try:
            tree = Tree(newick, format=tree_format)
        except Exception:
            continue
        if len(tree.get_leaf_names()) >= 2:
            return True
    return False


def _read_best_tree_path(tree_dir: str, prefix: str) -> Optional[str]:
    if not os.path.isdir(tree_dir):
        return None
    candidates = [
        os.path.join(tree_dir, f"{prefix}.treefile"),
        os.path.join(tree_dir, f"{prefix}.contree"),
    ]
    for path in candidates:
        if valid_iqtree_tree(path):
            return path
    for filename in sorted(os.listdir(tree_dir)):
        if filename.endswith(".treefile") or filename.endswith(".contree"):
            path = os.path.join(tree_dir, filename)
            if valid_iqtree_tree(path):
                return path
    return None


def expected_mafft_output_path(*, input_fasta: str, out_dir: str, output_name: Optional[str] = None) -> str:
    base_name = output_name or (Path(str(input_fasta)).stem + ".aln.fasta")
    return os.path.join(str(out_dir), base_name)


def expected_iqtree_tree_dir(*, input_alignment: str, out_dir: str, prefix: Optional[str] = None) -> tuple[str, str]:
    token = prefix or Path(str(input_alignment)).stem
    return os.path.join(str(out_dir), token), token


def run_mafft_alignment(
    task: Task,
    *,
    input_fasta: str,
    out_dir: str,
    output_name: Optional[str] = None,
    mafft_flags: object = None,
    thread_count: Optional[int] = None,
) -> tuple[str, list[str]]:
    mafft_path = _resolve_executable(
        task.db_manager.env.get("MAFFT_PATH"),
        fallback_names=("mafft",),
    )
    if not mafft_path:
        raise FileNotFoundError("MAFFT executable not found. Set MAFFT_PATH or install 'mafft' on PATH.")
    os.makedirs(out_dir, exist_ok=True)
    output_path = expected_mafft_output_path(
        input_fasta=str(input_fasta),
        out_dir=str(out_dir),
        output_name=output_name,
    )
    flags = _flags_from_env_or_task(
        env_value=task.db_manager.env.get("MAFFT_FLAGS"),
        task_value=mafft_flags,
        default_flags=["--localpair", "--maxiterate", "1000"],
    )
    flags = _ensure_thread_flag(
        flags,
        thread_tokens={"--thread"},
        inject_tokens=["--thread"],
        thread_count=thread_count if thread_count is not None else task.REQUIRED_THREADS,
    )
    command = [mafft_path, *flags, str(input_fasta)]
    task.log(f"Running MAFFT: {' '.join(command)}", "DEBUG")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {result.stderr.strip()}")
    temporary_path = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=out_dir,
            prefix=f".{os.path.basename(output_path)}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = handle.name
            handle.write(result.stdout)
        if not valid_mafft_alignment(temporary_path, str(input_fasta)):
            raise RuntimeError("MAFFT produced an invalid or incomplete alignment.")
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
    return output_path, command


def run_iqtree_analysis(
    task: Task,
    *,
    input_alignment: str,
    out_dir: str,
    prefix: Optional[str] = None,
    iqtree_flags: object = None,
    thread_count: Optional[int] = None,
) -> tuple[str, str, list[str]]:
    iqtree_path = _resolve_executable(
        task.db_manager.env.get("IQTREE_PATH"),
        fallback_names=("iqtree2", "iqtree"),
    )
    if not iqtree_path:
        raise FileNotFoundError("IQ-TREE executable not found. Set IQTREE_PATH or install 'iqtree2'/'iqtree' on PATH.")
    os.makedirs(out_dir, exist_ok=True)
    tree_dir, token = expected_iqtree_tree_dir(
        input_alignment=str(input_alignment),
        out_dir=str(out_dir),
        prefix=prefix,
    )
    os.makedirs(tree_dir, exist_ok=True)
    flags = _flags_from_env_or_task(
        env_value=task.db_manager.env.get("IQTREE_FLAGS"),
        task_value=iqtree_flags,
        default_flags=["-m", "MFP", "-B", "1000"],
    )
    flags = _ensure_thread_flag(
        flags,
        thread_tokens={"-nt", "--threads", "--threads-max"},
        inject_tokens=["-nt"],
        thread_count=thread_count if thread_count is not None else task.REQUIRED_THREADS,
    )
    if task.payload_bool("force_restart", False) and "-redo" not in flags and "--redo" not in flags:
        flags.append("-redo")
    command = [iqtree_path, "-s", str(input_alignment), "--prefix", os.path.join(tree_dir, token), *flags]
    task.log(f"Running IQ-TREE: {' '.join(command)}", "DEBUG")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed: {result.stderr.strip()}")
    best_tree = _read_best_tree_path(tree_dir, token)
    if not best_tree:
        raise FileNotFoundError(f"IQ-TREE completed but no tree file was found in {tree_dir}.")
    return tree_dir, best_tree, command


class MafftTask(Task):
    def run(self):
        input_fasta = str(self.data.get("input_fasta") or "")
        out_dir = str(self.data.get("out_dir") or "")
        output_name = self.data.get("output_name")
        if not input_fasta or not os.path.exists(input_fasta):
            return self.handle_exception("Input FASTA does not exist.", {"input_fasta": input_fasta})
        if not out_dir:
            return self.handle_exception("Output directory is required.", {})
        output_path = expected_mafft_output_path(
            input_fasta=input_fasta,
            out_dir=out_dir,
            output_name=output_name,
        )
        command: list[str] = []
        if valid_mafft_alignment(output_path, input_fasta):
            self.log(f"Using completed MAFFT alignment at {output_path}.", "INFO")
        else:
            if os.path.exists(output_path):
                self.log(f"Existing MAFFT alignment is invalid; rebuilding {output_path}.", "WARNING")
            self.log(f"Running MAFFT on {os.path.basename(input_fasta)}.", "INFO")
            try:
                output_path, command = run_mafft_alignment(
                    self,
                    input_fasta=input_fasta,
                    out_dir=out_dir,
                    output_name=output_name,
                    mafft_flags=self.data.get("mafft_flags"),
                )
            except Exception as exc:  # boundary: external MAFFT/filesystem failure becomes this task error
                return self.handle_exception("Failed to run MAFFT.", {"error": str(exc), "input_fasta": input_fasta})
        if not valid_mafft_alignment(output_path, input_fasta):
            return self.handle_exception(
                "MAFFT output is invalid or incomplete.",
                {"input_fasta": input_fasta, "output_path": output_path},
            )
        try:
            self.db_manager.artifacts.register(
                owner_type="task",
                owner_id=self.task_id,
                artifact_type="mafft_alignment",
                path=output_path,
                format="fasta",
                sequence_kind="prot",
                metadata={"command": command},
            )
        except Exception as exc:  # boundary: optional artifact catalog metadata; alignment file already exists.
            self.log(f"Failed to register MAFFT artifact: {exc}", "WARNING")
        self.log(f"MAFFT alignment written to {output_path}.", "INFO")
        return True


class IQTreeTask(Task):
    def run(self):
        input_alignment = str(self.data.get("input_alignment") or "")
        out_dir = str(self.data.get("out_dir") or "")
        prefix = self.data.get("prefix")
        if not input_alignment or not os.path.exists(input_alignment):
            return self.handle_exception("Input alignment does not exist.", {"input_alignment": input_alignment})
        if not valid_mafft_alignment(input_alignment):
            return self.handle_exception("Input alignment is invalid or incomplete.", {"input_alignment": input_alignment})
        if not out_dir:
            return self.handle_exception("Output directory is required.", {})
        tree_dir, token = expected_iqtree_tree_dir(
            input_alignment=input_alignment,
            out_dir=out_dir,
            prefix=prefix,
        )
        force_restart = self.payload_bool("force_restart", False)
        best_tree = _read_best_tree_path(tree_dir, token)
        command: list[str] = []
        if best_tree and not force_restart:
            self.log(f"Using completed IQ-TREE tree at {best_tree}.", "INFO")
        else:
            if force_restart:
                self.log(f"Rebuilding IQ-TREE output at {tree_dir}.", "WARNING")
                # Consume the clean-restart request before launching IQ-TREE. If the
                # machine dies during this attempt, startup recovery can resume the
                # new IQ-TREE checkpoint instead of repeatedly applying -redo.
                persisted_data = dict(self.data)
                persisted_data["force_restart"] = False
                self.db_manager.tasks.update_data(self.task_id, data=persisted_data)
            self.log(f"Running IQ-TREE on {os.path.basename(input_alignment)}.", "INFO")
            try:
                tree_dir, best_tree, command = run_iqtree_analysis(
                    self,
                    input_alignment=input_alignment,
                    out_dir=out_dir,
                    prefix=prefix,
                    iqtree_flags=self.data.get("iqtree_flags"),
                )
            except Exception as exc:  # boundary: external IQ-TREE/filesystem failure becomes this task error
                return self.handle_exception("Failed to run IQ-TREE.", {"error": str(exc), "input_alignment": input_alignment})
        if not valid_iqtree_tree(best_tree):
            return self.handle_exception(
                "IQ-TREE output is invalid or incomplete.",
                {"input_alignment": input_alignment, "tree_path": best_tree},
            )
        try:
            self.db_manager.artifacts.register(
                owner_type="task",
                owner_id=self.task_id,
                artifact_type="iqtree_results_dir",
                path=tree_dir,
                is_dir=True,
                format="directory",
                metadata={"command": command},
            )
            self.db_manager.artifacts.register(
                owner_type="task",
                owner_id=self.task_id,
                artifact_type="iqtree_best_tree",
                path=best_tree,
                format="newick",
                metadata={"command": command},
            )
        except Exception as exc:  # boundary: optional artifact catalog metadata; tree output already exists.
            self.log(f"Failed to register IQ-TREE artifacts: {exc}", "WARNING")
        self.log(f"IQ-TREE best tree written to {best_tree}.", "INFO")
        return True


class BuildBuscoTreesTask(ExportLibraryTask):
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=4):
        super().__init__(db_path, task_id, checkpoint, data, required_threads=required_threads)
        self.alignments_dir: Optional[str] = self.data.get("alignments_dir")
        self.trees_dir: Optional[str] = self.data.get("trees_dir")
        self.manifest_path: Optional[str] = self.data.get("manifest_path")
        self.busco_families_dir: Optional[str] = self.data.get("busco_families_dir")
        self.mafft_threads: int = int(self.data.get("mafft_threads") or DEFAULT_MAFFT_TASK_THREADS)
        self.iqtree_threads: int = int(self.data.get("iqtree_threads") or DEFAULT_IQTREE_TASK_THREADS)

    def _default_export_dir_name(self) -> str:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        suffix = re.sub(r"[^A-Za-z0-9._-]+", "_", str(self.library_name or self.library_id or "busco_trees")).strip("._-")
        return f"task_{self.task_id}_{stamp}_{suffix or 'busco_trees'}"

    def _fasta_files(self) -> list[str]:
        if not self.busco_families_dir or not os.path.isdir(self.busco_families_dir):
            return []
        return sorted(
            filename
            for filename in os.listdir(self.busco_families_dir)
            if filename.endswith(".fasta")
        )

    def _family_rows(self) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        for filename in self._fasta_files():
            family_id = os.path.splitext(filename)[0]
            raw_fasta = os.path.join(str(self.busco_families_dir), filename)
            alignment_path = expected_mafft_output_path(
                input_fasta=raw_fasta,
                out_dir=str(self.alignments_dir),
                output_name=f"{family_id}.aln.fasta",
            )
            tree_dir, prefix = expected_iqtree_tree_dir(
                input_alignment=alignment_path,
                out_dir=str(self.trees_dir),
                prefix=family_id,
            )
            rows.append(
                {
                    "family_id": family_id,
                    "raw_fasta": raw_fasta,
                    "alignment_path": alignment_path,
                    "tree_dir": tree_dir,
                    "prefix": prefix,
                }
            )
        return rows

    def _queue_mafft_subtasks(self) -> bool:
        queued = False
        for row in self._family_rows():
            if valid_mafft_alignment(row["alignment_path"], row["raw_fasta"]):
                continue
            self.queue_subtask(
                job_type=32,
                status="P",
                priority=1,
                data={
                    "input_fasta": row["raw_fasta"],
                    "out_dir": self.alignments_dir,
                    "output_name": os.path.basename(row["alignment_path"]),
                    "mafft_flags": self.data.get("mafft_flags"),
                    "required_threads": int(self.mafft_threads),
                },
            )
            queued = True
        return queued

    def _mafft_done(self) -> bool:
        rows = self._family_rows()
        return bool(rows) and all(
            valid_mafft_alignment(row["alignment_path"], row["raw_fasta"])
            for row in rows
        )

    def _queue_iqtree_subtasks(self) -> bool:
        queued = False
        for row in self._family_rows():
            best_tree = _read_best_tree_path(row["tree_dir"], row["prefix"])
            if best_tree:
                continue
            self.queue_subtask(
                job_type=33,
                status="P",
                priority=1,
                data={
                    "input_alignment": row["alignment_path"],
                    "out_dir": self.trees_dir,
                    "prefix": row["family_id"],
                    "iqtree_flags": self.data.get("iqtree_flags"),
                    "force_restart": bool(getattr(self, "_phase_meta", {}).get("gen", 1) > 1),
                    "required_threads": int(self.iqtree_threads),
                },
            )
            queued = True
        return queued

    def _iqtree_done(self) -> bool:
        rows = self._family_rows()
        if not rows:
            return False
        for row in rows:
            best_tree = _read_best_tree_path(row["tree_dir"], row["prefix"])
            if not best_tree:
                return False
        return True

    def _write_manifest(self) -> None:
        with open(self.manifest_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["family_id", "raw_fasta", "alignment_path", "tree_dir", "tree_path"])
            for row in self._family_rows():
                best_tree = _read_best_tree_path(row["tree_dir"], row["prefix"])
                if not best_tree:
                    raise FileNotFoundError(f"IQ-TREE output missing for family {row['family_id']}.")
                writer.writerow([
                    row["family_id"],
                    row["raw_fasta"],
                    row["alignment_path"],
                    row["tree_dir"],
                    best_tree,
                ])

    def run(self):
        if self.stage < 1:
            exported = super().run()
            if exported is not True:
                return exported
            self.alignments_dir = os.path.join(self.out_dir, "alignments")
            self.trees_dir = os.path.join(self.out_dir, "trees")
            self.manifest_path = os.path.join(self.out_dir, "manifest.tsv")
            self.busco_families_dir = os.path.join(self.out_dir, "busco_families")
            try:
                os.makedirs(self.alignments_dir, exist_ok=True)
                os.makedirs(self.trees_dir, exist_ok=True)
            except OSError as exc:
                return self.handle_exception("Failed to create BUSCO tree output directories.", {"error": str(exc)})
            self.checkpoint(
                1,
                {
                    "alignments_dir": self.alignments_dir,
                    "trees_dir": self.trees_dir,
                    "manifest_path": self.manifest_path,
                    "busco_families_dir": self.busco_families_dir,
                },
            )

        if not self._fasta_files():
            return self.handle_exception("No BUSCO family FASTAs were available for tree building.", {})
        self.log(
            f"Building BUSCO family trees for {len(self._fasta_files())} families into {self.trees_dir or self.out_dir}.",
            "INFO",
        )

        outcome = self.manage_subtasks(
            stage=2,
            queue_fn=self._queue_mafft_subtasks,
            done_fn=self._mafft_done,
            wait_seconds=0,
            max_retries=int(self.data.get("mafft_retries", 1)),
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        outcome = self.manage_subtasks(
            stage=3,
            queue_fn=self._queue_iqtree_subtasks,
            done_fn=self._iqtree_done,
            wait_seconds=0,
            max_retries=int(self.data.get("iqtree_retries", 1)),
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        try:
            self._write_manifest()
        except (OSError, csv.Error, FileNotFoundError) as exc:
            return self.handle_exception("Failed while building BUSCO family trees.", {"error": str(exc)})
        self.db_manager.artifacts.register(
            owner_type="export_run",
            owner_id=self.task_id,
            artifact_type="busco_tree_alignments_dir",
            path=self.alignments_dir,
            is_dir=True,
            format="directory",
        )
        self.db_manager.artifacts.register(
            owner_type="export_run",
            owner_id=self.task_id,
            artifact_type="busco_tree_results_dir",
            path=self.trees_dir,
            is_dir=True,
            format="directory",
        )
        self.db_manager.artifacts.register(
            owner_type="export_run",
            owner_id=self.task_id,
            artifact_type="busco_tree_manifest",
            path=self.manifest_path,
            format="tsv",
        )
        self.log(f"Built BUSCO family trees and wrote manifest to {self.manifest_path}.", "INFO")
        return True


class _OrthogroupTreeAnnotationMixin:
    _TREE_COLOR_BUSCO = "#2e7d32"
    _TREE_COLOR_INPARALOG = "#1565c0"
    _TREE_COLOR_OUTPARALOG = "#c62828"
    _TREE_COLOR_HIDDEN_PARALOG = "#c62828"
    _TREE_COLOR_NON_BUSCO = "#616161"
    _BUSCO_STATUS_MARKERS = {
        1: "BUSCO_SC",
        2: "BUSCO_DUP",
        3: "BUSCO_FRAG",
        4: "BUSCO_MISSING",
    }

    def _read_mapping_manifest(self) -> dict[str, dict[str, str]]:
        manifest_path = str(self.data.get("manifest_tsv") or "")
        if not manifest_path or not os.path.exists(manifest_path):
            return {}
        rows: dict[str, dict[str, str]] = {}
        with open(manifest_path, "r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                orthogroup = str(row.get("orthogroup") or row.get("family_id") or "").strip()
                if orthogroup:
                    rows[orthogroup] = {str(k): str(v or "") for k, v in row.items()}
        return rows

    def _extract_accession_token(self, token: str) -> str:
        text = str(token or "").strip()
        if not text:
            return ""
        match = re.match(r"^(GC[AF])_(\d+)[._](\d+)\b", text, flags=re.IGNORECASE)
        if match:
            prefix, digits, version = match.groups()
            return canonicalize_accession(f"{prefix}_{digits}.{version}")
        trimmed = re.sub(r"_(?:gff|cdhit|clean|raw)(?:[_-].*|\d.*)?$", "", text, flags=re.IGNORECASE)
        return canonicalize_accession(trimmed)

    def _load_leaf_metadata(self, orthofinder_location: str) -> dict[str, dict[str, str]]:
        working_dir = os.path.join(str(orthofinder_location), "WorkingDirectory")
        species_ids_path = os.path.join(working_dir, "SpeciesIDs.txt")
        sequence_ids_path = os.path.join(working_dir, "SequenceIDs.txt")
        if not (os.path.exists(species_ids_path) and os.path.exists(sequence_ids_path)):
            return {}
        species_map: dict[str, str] = {}
        with open(species_ids_path, "r", encoding="utf-8") as handle:
            for raw in handle:
                if ":" not in raw:
                    continue
                species_idx, label = raw.split(":", 1)
                species_token = os.path.splitext(os.path.basename(label.strip().replace(".gz", "")))[0]
                accession = self._extract_accession_token(species_token)
                if accession:
                    species_map[str(species_idx).strip()] = species_token
        leaf_map: dict[str, dict[str, str]] = {}
        with open(sequence_ids_path, "r", encoding="utf-8") as handle:
            for raw in handle:
                if ":" not in raw:
                    continue
                seq_id, leaf_name = raw.split(":", 1)
                species_idx = str(seq_id).split("_", 1)[0].strip()
                species_token = species_map.get(species_idx, "")
                accession = self._extract_accession_token(species_token)
                sequence = str(leaf_name).strip().split()[0]
                if accession:
                    metadata = {"accession": accession, "sequence": sequence}
                    leaf_map[sequence] = metadata
                    if species_token:
                        normalized_species = species_token.replace(".", "_")
                        leaf_map[f"{normalized_species}_{sequence}"] = metadata
        return leaf_map

    def _parse_paralog_leaf_names(self, path: str) -> set[str]:
        leaves: set[str] = set()
        if not path or not os.path.exists(path):
            return leaves
        with open(path, "r", encoding="utf-8") as handle:
            for raw in handle:
                tokens = re.findall(r"[A-Za-z0-9._-]+", str(raw or ""))
                for token in tokens:
                    leaves.add(token)
        return leaves

    def _taxonomy_labels_for_accession(self, accession: str) -> dict[str, str]:
        row = self.db_manager.genomes.get(str(accession))
        if not row:
            return {"taxon": "NA", "phylum": "NA"}
        taxid = row[1]
        lineage_rows = self.db_manager.genomes.get_lineage_root_to_leaf(int(taxid)) if taxid is not None else []
        lineage_map = {str(rank).lower(): str(name) for (_tid, name, rank, _parent) in lineage_rows if rank and name}
        taxon_name = next((str(name) for (tid, name, _rank, _parent) in lineage_rows if tid == taxid and name), "")
        return {
            "taxon": re.sub(r"[^A-Za-z0-9.]+", "_", taxon_name).strip("_") or "NA",
            "phylum": re.sub(r"[^A-Za-z0-9.]+", "_", lineage_map.get("phylum", "")).strip("_") or "NA",
        }

    def _resolve_leaf_metadata(self, leaf_name: str, leaf_metadata: dict[str, dict[str, str]]) -> dict[str, str]:
        metadata = leaf_metadata.get(leaf_name)
        if metadata:
            return metadata
        accession = self._extract_accession_token(leaf_name)
        sequence = leaf_name
        if accession and leaf_name.startswith(accession + "_"):
            sequence = leaf_name[len(accession) + 1:]
        elif accession and leaf_name.startswith(accession.replace(".", "_") + "_"):
            sequence = leaf_name[len(accession.replace(".", "_")) + 1:]
        else:
            match = re.match(r"^(GC[AF]_\d+_\d+)(?:_.+)?$", leaf_name, flags=re.IGNORECASE)
            if match and leaf_name.startswith(match.group(1) + "_"):
                sequence = leaf_name[len(match.group(1)) + 1:]
        return {
            "accession": accession,
            "sequence": sequence,
        }

    def _load_busco_markers(
        self,
        *,
        family_id: str,
        leaf_metadata: dict[str, dict[str, str]],
        source_run_ids: list[int],
    ) -> dict[tuple[str, str], str]:
        family_token = str(family_id or "").strip()
        if not family_token:
            return {}
        if source_run_ids:
            placeholders = ",".join("?" for _ in source_run_ids)
            rows = self.db_manager.cursor.execute(
                f"""
                SELECT run_id, accession, status, sequence
                FROM BUSCO_Run_Family_Data
                WHERE family_id = ? AND run_id IN ({placeholders})
                ORDER BY accession ASC, run_id DESC, sequence ASC
                """,
                (family_token, *source_run_ids),
            ).fetchall()
        else:
            rows = self.db_manager.cursor.execute(
                """
                SELECT run_id, accession, status, sequence
                FROM BUSCO_Run_Family_Data
                WHERE family_id = ?
                ORDER BY accession ASC, run_id DESC, sequence ASC
                """,
                (family_token,),
            ).fetchall()
        latest_by_accession: dict[str, tuple[int, list[tuple[int, str]]]] = {}
        for run_id, accession, status, sequence in rows or []:
            acc_token = canonicalize_accession(accession)
            if not acc_token:
                continue
            if sequence is None and int(status or 0) == 4:
                continue
            current = latest_by_accession.get(acc_token)
            if current is None or int(run_id) > current[0]:
                latest_by_accession[acc_token] = (int(run_id), [(int(status or 0), str(sequence or ""))])
            elif int(run_id) == current[0]:
                current[1].append((int(status or 0), str(sequence or "")))
        markers: dict[tuple[str, str], str] = {}
        for accession, (_run_id, entries) in latest_by_accession.items():
            for status, sequence in entries:
                if not sequence:
                    continue
                markers[(accession, sequence)] = self._BUSCO_STATUS_MARKERS.get(status, "BUSCO_OTHER")
        return markers

    def _load_hidden_paralog_leaf_names(
        self,
        *,
        family_id: str,
        source_run_ids: list[int],
    ) -> set[str]:
        target_library_id = self.data.get("target_library_id")
        busco_library_id = self.data.get("busco_library_id")
        try:
            target_library_id = int(target_library_id)
            busco_library_id = int(busco_library_id)
        except (TypeError, ValueError):
            return set()
        return self.db_manager.filtering.get_hidden_paralog_sequence_ids(
            target_library_id=target_library_id,
            busco_library_id=busco_library_id,
            family_id=family_id,
            source_run_ids=source_run_ids,
        )

    def _render_nexus(
        self,
        *,
        tree_path: str,
        family_id: str,
        source_run_ids: list[int],
        in_leaves: set[str],
        out_leaves: set[str],
        leaf_metadata: dict[str, dict[str, str]],
    ) -> str:
        tree = Tree(open(tree_path, "r", encoding="utf-8").read().strip(), format=1)
        tree.ladderize()
        busco_markers = self._load_busco_markers(
            family_id=family_id,
            leaf_metadata=leaf_metadata,
            source_run_ids=source_run_ids,
        )
        hidden_paralog_leaves = self._load_hidden_paralog_leaf_names(
            family_id=family_id,
            source_run_ids=source_run_ids,
        )

        def _leaf_marker(leaf_name: str) -> tuple[str, bool]:
            metadata = self._resolve_leaf_metadata(leaf_name, leaf_metadata)
            accession = metadata.get("accession") or "NA"
            sequence = metadata.get("sequence") or leaf_name
            marker = busco_markers.get((accession, sequence), "NON_BUSCO")
            return marker, marker != "NON_BUSCO"

        def _display_name(leaf_name: str) -> str:
            metadata = self._resolve_leaf_metadata(leaf_name, leaf_metadata)
            accession = metadata.get("accession") or "NA"
            taxonomy = self._taxonomy_labels_for_accession(accession)
            marker, _is_busco = _leaf_marker(leaf_name)
            hidden_suffix = "[*HP]" if leaf_name in hidden_paralog_leaves else ""
            return f"{leaf_name}_{taxonomy['taxon']}_{taxonomy['phylum']}_[{marker}]{hidden_suffix}"

        def _format_node(node) -> str:
            children = ""
            if not node.is_leaf():
                children = "(" + ",".join(_format_node(child) for child in node.children) + ")"
            label = ""
            if node.name:
                display_name = _display_name(str(node.name))
                label = "'" + display_name.replace("'", "''") + "'"
            color = ""
            if node.is_leaf():
                leaf_name = str(node.name)
                _marker, is_busco = _leaf_marker(leaf_name)
                if leaf_name in hidden_paralog_leaves:
                    color = self._TREE_COLOR_HIDDEN_PARALOG
                elif leaf_name in out_leaves:
                    color = self._TREE_COLOR_OUTPARALOG
                elif leaf_name in in_leaves:
                    color = self._TREE_COLOR_INPARALOG
                elif not is_busco:
                    color = self._TREE_COLOR_NON_BUSCO
                else:
                    color = self._TREE_COLOR_BUSCO
            annotation = f"[&!color={color}]" if color else ""
            branch = f":{node.dist}" if node.dist is not None else ""
            return f"{children}{label}{annotation}{branch}"

        leaves = [str(leaf.name or "") for leaf in tree.iter_leaves()]
        lines = [
            "#NEXUS",
            "",
            "begin taxa;",
            f"    dimensions ntax={len(leaves)};",
            "    taxlabels",
        ]
        for leaf in leaves:
            display_name = _display_name(leaf).replace("'", "''")
            lines.append(f"        '{display_name}'")
        lines.extend(
            [
                "    ;",
                "end;",
                "",
                "begin trees;",
                f"    tree annotated = [&R] {_format_node(tree)};",
                "end;",
                "",
            ]
        )
        return "\n".join(lines)


class AnnotateOrthogroupTreeTask(Task, _OrthogroupTreeAnnotationMixin):
    def _tree_paths_from_dir(self, input_dir: str) -> list[str]:
        paths: list[str] = []
        for filename in sorted(os.listdir(input_dir)):
            if filename == "Resolved_Gene_Trees.txt":
                continue
            if not (filename.endswith("_tree.txt") or filename.endswith(".treefile")):
                continue
            paths.append(os.path.join(input_dir, filename))
        return paths

    def run(self):
        input_tree = str(self.data.get("input_tree") or "")
        input_dir = str(self.data.get("input_dir") or "")
        output_dir = str(self.data.get("out_dir") or "")
        if bool(input_tree) == bool(input_dir):
            return self.handle_exception("Provide exactly one of input_tree or input_dir.", {})
        if not output_dir:
            source_dir = os.path.dirname(input_tree) if input_tree else input_dir
            if source_dir and os.access(source_dir, os.W_OK):
                output_dir = os.path.join(source_dir, "annotated-og-trees")
            else:
                return self.handle_exception("Output directory is required when the input location is not writable.", {})
        os.makedirs(output_dir, exist_ok=True)
        orthofinder_location = str(self.data.get("orthofinder_location") or "")
        leaf_metadata = self._load_leaf_metadata(orthofinder_location) if orthofinder_location else {}
        manifest = self._read_mapping_manifest()
        tree_paths = [input_tree] if input_tree else self._tree_paths_from_dir(input_dir)
        for tree_path in tree_paths:
            base = os.path.basename(tree_path)
            orthogroup = base[:-9] if base.endswith("_tree.txt") else os.path.splitext(base)[0]
            row = manifest.get(orthogroup, {})
            family_id = str(row.get("family_id") or row.get("busco_family_id") or "").strip()
            source_run_ids = [int(token) for token in str(row.get("source_run_ids") or "").split(",") if token.strip().isdigit()]
            in_file = row.get("paralog_in_file") or row.get("inparalog_file") or ""
            out_file = row.get("paralog_out_file") or row.get("outparalog_file") or ""
            in_leaves = self._parse_paralog_leaf_names(in_file)
            out_leaves = self._parse_paralog_leaf_names(out_file)
            try:
                nexus = self._render_nexus(
                    tree_path=tree_path,
                    family_id=family_id,
                    source_run_ids=source_run_ids,
                    in_leaves=in_leaves,
                    out_leaves=out_leaves,
                    leaf_metadata=leaf_metadata,
                )
                with open(os.path.join(output_dir, f"{orthogroup}.nex"), "w", encoding="utf-8") as handle:
                    handle.write(nexus)
            except Exception as exc:  # boundary: one tree annotation failure becomes this task error
                return self.handle_exception("Failed to annotate orthogroup tree.", {"tree_path": tree_path, "error": str(exc)})
        return True
