from __future__ import annotations

import os
from .base import BaseRepository
from ...accession_utils import canonicalize_accessions


class OrthoFinderRepository(BaseRepository):
    @staticmethod
    def _normalize_mcl_inflation(value):
        if value in (None, ""):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _parse_logged_mcl_inflation(command_line):
        import re

        text = str(command_line or "").strip()
        if not text:
            return None
        match = re.search(r"(?:^|\s)-I\s+([0-9]*\.?[0-9]+)(?:\s|$)", text)
        if not match:
            return None
        try:
            return float(match.group(1))
        except ValueError:
            return None

    def _read_logged_command_line(self, location):
        if not location:
            return None
        log_path = os.path.join(str(location), "Log.txt")
        if not os.path.exists(log_path):
            return None
        try:
            with open(log_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if line.startswith("Command Line:"):
                        return line.split("Command Line:", 1)[1].strip()
        except OSError:
            return None
        return None

    @staticmethod
    def _has_core_orthogroup_outputs(location):
        if not location:
            return False
        orthogroup_sequences_dir = os.path.join(str(location), "Orthogroup_Sequences")
        return os.path.isdir(orthogroup_sequences_dir)

    @staticmethod
    def _canonical_accessions(accessions):
        cleaned = []
        seen = set()
        for accession in canonicalize_accessions(accessions or []):
            if accession in seen:
                continue
            seen.add(accession)
            cleaned.append(accession)
        return cleaned

    def add_results(self, library_id, datetime=None, location=None, status: str = "ready", mcl_inflation=None, command_line=None):
        with self.core.transaction(operation=f"add OrthoFinder results for library {library_id}"):
            columns = ["library_id"]
            values = [library_id]
            if location is not None:
                columns.append("location")
                values.append(location)
            if datetime is not None:
                columns.append("date")
                values.append(datetime)
            columns.append("mcl_inflation")
            values.append(self._normalize_mcl_inflation(mcl_inflation))
            columns.append("command_line")
            values.append(str(command_line) if command_line is not None else None)
            columns.append("status")
            values.append(status)
            col_str = ", ".join(columns)
            placeholders = ", ".join(["?"] * len(values))
            self.core.execute(f"INSERT INTO OrthoFinder_Results ({col_str}) VALUES ({placeholders})", tuple(values))
            return self.cursor.lastrowid

    def add_accessions(self, orthofinder_id, accessions):
        canonical_accessions = self._canonical_accessions(accessions)
        with self.core.transaction(operation=f"add OrthoFinder accessions {orthofinder_id}"):
            if canonical_accessions:
                self.core.executemany(
                    "INSERT OR IGNORE INTO OrthoFinder_Accessions (orthofinder_id, accession) VALUES (?, ?)",
                    [(orthofinder_id, acc) for acc in canonical_accessions],
                )
            return True

    def update_location(self, orthofinder_id, location):
        with self.core.transaction(operation=f"update OrthoFinder location {orthofinder_id}"):
            self.core.execute(
                "UPDATE OrthoFinder_Results SET location = ? WHERE orthofinder_id = ?",
                (location, orthofinder_id),
            )
            return True

    def update_run_metadata(self, orthofinder_id, *, mcl_inflation=None, command_line=None):
        with self.core.transaction(operation=f"update OrthoFinder metadata {orthofinder_id}"):
            self.core.execute(
                "UPDATE OrthoFinder_Results SET mcl_inflation = ?, command_line = ? WHERE orthofinder_id = ?",
                (self._normalize_mcl_inflation(mcl_inflation), str(command_line) if command_line is not None else None, int(orthofinder_id)),
            )
            return True

    def set_status(self, orthofinder_id, status: str):
        with self.core.transaction(operation=f"set OrthoFinder status {orthofinder_id}"):
            self.core.execute(
                "UPDATE OrthoFinder_Results SET status = ? WHERE orthofinder_id = ?",
                (str(status), int(orthofinder_id)),
            )
            return True

    def get_status(self, orthofinder_id):
        row = self.core.fetchone(
            "SELECT COALESCE(status, 'ready') FROM OrthoFinder_Results WHERE orthofinder_id = ?",
            (int(orthofinder_id),),
        )
        return str(row[0]) if row and row[0] is not None else None

    def get(self, orthofinder_id: int):
        return self.core.fetchone(
            "SELECT orthofinder_id, library_id, location, COALESCE(status, 'ready'), date, mcl_inflation, command_line FROM OrthoFinder_Results WHERE orthofinder_id = ?",
            (int(orthofinder_id),),
        )

    def get_many(self, *, library_id=None, include_inactive: bool = True):
        clauses = []
        params = []
        if library_id is not None:
            clauses.append("library_id = ?")
            params.append(int(library_id))
        if not include_inactive:
            clauses.append("COALESCE(status, 'ready') = 'ready'")
        sql = "SELECT orthofinder_id, library_id, location, COALESCE(status, 'ready'), date, mcl_inflation, command_line FROM OrthoFinder_Results"
        if clauses:
            sql += " WHERE " + " AND ".join(clauses)
        sql += " ORDER BY orthofinder_id"
        return self.core.fetchall(sql, tuple(params))

    def get_accessions(self, orthofinder_id: int):
        rows = self.core.fetchall(
            "SELECT accession FROM OrthoFinder_Accessions WHERE orthofinder_id = ? ORDER BY accession",
            (int(orthofinder_id),),
        )
        return [row[0] for row in rows]

    def replace_accessions(self, orthofinder_id: int, accessions):
        accessions = self._canonical_accessions(accessions)
        with self.core.transaction(operation=f"replace OrthoFinder accessions {orthofinder_id}"):
            self.core.execute("DELETE FROM OrthoFinder_Accessions WHERE orthofinder_id = ?", (int(orthofinder_id),))
            if accessions:
                self.core.executemany(
                    "INSERT OR IGNORE INTO OrthoFinder_Accessions (orthofinder_id, accession) VALUES (?, ?)",
                    [(int(orthofinder_id), acc) for acc in accessions],
                )
            return True

    def delete_results(self, orthofinder_id):
        with self.core.transaction(operation=f"delete OrthoFinder results {orthofinder_id}"):
            self.core.execute("DELETE FROM OrthoFinder_Accessions WHERE orthofinder_id = ?", (orthofinder_id,))
            self.core.execute("DELETE FROM OrthoFinder_Results WHERE orthofinder_id = ?", (orthofinder_id,))
            return True

    def _effective_mcl_inflation(self, row):
        stored = self._normalize_mcl_inflation(row[5] if len(row) > 5 else None)
        if stored is not None:
            return stored
        command_line = row[6] if len(row) > 6 else None
        if not command_line:
            command_line = self._read_logged_command_line(row[2] if len(row) > 2 else None)
        parsed = self._parse_logged_mcl_inflation(command_line)
        if parsed is not None:
            self.update_run_metadata(int(row[0]), mcl_inflation=parsed, command_line=command_line)
        return parsed

    def _matching_ready_runs(self, accessions, mcl_inflation=None):
        requested = set(self._canonical_accessions(accessions))
        requested_inflation = self._normalize_mcl_inflation(mcl_inflation)
        if not requested:
            return []
        rows = self.core.fetchall(
            """
            SELECT orthofinder_id, library_id, location, COALESCE(status, 'ready'), date, mcl_inflation, command_line
            FROM OrthoFinder_Results
            WHERE COALESCE(status, 'ready') = 'ready'
            ORDER BY orthofinder_id DESC
            """
        )
        matches = []
        for row in rows:
            orthofinder_id = int(row[0])
            if not self._has_core_orthogroup_outputs(row[2]):
                continue
            stored = set(self._canonical_accessions(self.get_accessions(orthofinder_id)))
            if stored == requested:
                if self._effective_mcl_inflation(row) != requested_inflation:
                    continue
                matches.append(row)
        return matches

    def assert_results_exist(self, accessions, mcl_inflation=None):
        matches = self._matching_ready_runs(accessions, mcl_inflation=mcl_inflation)
        if not matches:
            return None
        return matches[0][0], matches[0][2]

    def get_by_reference_accessions(self, accessions, mcl_inflation=None):
        return self._matching_ready_runs(accessions, mcl_inflation=mcl_inflation)
