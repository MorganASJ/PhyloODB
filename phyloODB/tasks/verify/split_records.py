import os
import shutil
import gzip
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

from ..task import Task


class SplitRecordsTask(Task):
    """Split multiple fasta files in a genome folder into separate accessions."""

    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.accession = self.data.get("accession")
        self.folder = self.data.get("folder")
        self.split_isolated_proteomes = bool(self.data.get("split_isolated_proteomes", False))
        self.check_only = bool(self.data.get("check_only", False))
        self.report = self.data.get("report")

    def _append_comment(self, accession: str, note: str) -> None:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        try:
            self.db_manager.cursor.execute("SELECT comments FROM Genome WHERE accession = ?", (accession,))
            row = self.db_manager.cursor.fetchone()
            existing = row[0] if row else ""
            new_comment = f"[{ts}] {note}"
            combined = f"{existing}\n{new_comment}" if existing else new_comment
            self.db_manager.cursor.execute("UPDATE Genome SET comments = ? WHERE accession = ?", (combined, accession))
            self.db_manager.commit()
        except Exception as exc:  # boundary: comment enrichment failure should not block split operation.
            self.log(f"Failed to append split comment for {accession}: {exc}", "WARNING")

    def _list_pairs(self) -> Tuple[List[str], List[str]]:
        nuc = []
        prot = []
        for fname in os.listdir(self.folder):
            if fname.endswith(".fna") or fname.endswith(".fna.gz"):
                nuc.append(fname)
            elif fname.endswith(".faa") or fname.endswith(".faa.gz"):
                prot.append(fname)
        return sorted(nuc), sorted(prot)

    def _next_accession(self, base: str, idx: int) -> str:
        return f"{base}.{idx}"

    def run(self):
        if not self.accession or not self.folder:
            return self.handle_exception("accession and folder are required", {})
        if not os.path.isdir(self.folder):
            return self.handle_exception("folder does not exist", {"folder": self.folder})

        nuc_files, prot_files = self._list_pairs()
        if len(nuc_files) + len(prot_files) <= 1:
            return True  # nothing to split

        report_lines: List[str] = []

        # Pair proteins by basename
        prot_by_base = {}
        for p in prot_files:
            base = p.replace(".gz", "").replace(".faa", "")
            prot_by_base.setdefault(base, []).append(p)

        # Decide primary (keep with original accession): oldest/first nuc file
        primary_nuc = nuc_files[0] if nuc_files else None
        primary_base = primary_nuc.replace(".gz", "").replace(".fna", "") if primary_nuc else None

        splits: List[Tuple[str, List[str], List[str]]] = []  # (new_acc, nuc_files, prot_files)

        # Additional nuc files become new accessions
        for idx, nuc in enumerate(nuc_files[1:], start=1):
            new_acc = self._next_accession(self.accession, idx)
            base = nuc.replace(".gz", "").replace(".fna", "")
            related_prot = prot_by_base.pop(base, [])
            splits.append((new_acc, [nuc], related_prot))

        # Protein-only entries if requested
        if self.split_isolated_proteomes:
            for base, files in list(prot_by_base.items()):
                new_acc = self._next_accession(self.accession, len(splits) + 1)
                splits.append((new_acc, [], files))
                prot_by_base.pop(base, None)

        if self.check_only:
            self.log(f"Would split {self.accession} into {len(splits)+1} records", "DEBUG")
            return True

        # Perform splits: move files to new folders and clone Genome/Assembly rows
        for idx, (new_acc, nucs, prots) in enumerate(splits, start=1):
            new_folder = os.path.join(os.path.dirname(self.folder), new_acc)
            os.makedirs(new_folder, exist_ok=True)
            for f in nucs + prots:
                src = os.path.join(self.folder, f)
                if os.path.exists(src):
                    shutil.move(src, os.path.join(new_folder, f))
            # clone Genome row
            try:
                self.db_manager.cursor.execute("SELECT * FROM Genome WHERE accession = ?", (self.accession,))
                row = self.db_manager.cursor.fetchone()
                if row:
                    cols = [c[1] for c in self.db_manager.cursor.execute("PRAGMA table_info(Genome)")]  # type: ignore
                    data = dict(zip(cols, row))
                    data["accession"] = new_acc
                    data["location"] = new_folder
                    data["status"] = 1 if nucs else 0
                    data["protein"] = 1 if prots else 0
                    self.db_manager.genomes.upsert(data)
                    self._append_comment(new_acc, f"Split from {self.accession}")
            except Exception as exc:  # boundary: one cloned Genome row failure is logged; remaining split records continue.
                self.log(f"Failed to clone Genome for {new_acc}: {exc}", "ERROR")
            # clone Assembly row
            try:
                self.db_manager.cursor.execute("SELECT * FROM Assembly WHERE accession = ?", (self.accession,))
                row = self.db_manager.cursor.fetchone()
                if row:
                    cols = [c[1] for c in self.db_manager.cursor.execute("PRAGMA table_info(Assembly)")]  # type: ignore
                    data = dict(zip(cols, row))
                    data["accession"] = new_acc
                    self.db_manager.genomes.insert_assembly(data)
                    self._append_comment(new_acc, f"Split assembly from {self.accession}")
            except Exception as exc:  # boundary: one cloned Assembly row failure is logged; remaining split records continue.
                self.log(f"Failed to clone Assembly for {new_acc}: {exc}", "ERROR")

            self._append_comment(self.accession, f"Split out {new_acc} ({', '.join(nucs+prots)})")
            report_lines.append(f"Split {self.accession} -> {new_acc}: {', '.join(nucs+prots)}")

        # Remove moved files from original folder record status flags
        try:
            protein_left = any(f.endswith(".faa") or f.endswith(".faa.gz") for f in os.listdir(self.folder))
            self.db_manager.genomes.update_status(self.accession, 1, (datetime.now(), self.folder, protein_left))
        except Exception as exc:  # boundary: status refresh failure is logged after file split work.
            self.log(f"Failed to refresh original split record status for {self.accession}: {exc}", "WARNING")

        if self.report and report_lines:
            try:
                Path(self.report).parent.mkdir(parents=True, exist_ok=True)
                Path(self.report).write_text("\n".join(report_lines), encoding="utf-8")
            except OSError as exc:
                self.log(f"Failed to write split-record report {self.report}: {exc}", "WARNING")

        return True
