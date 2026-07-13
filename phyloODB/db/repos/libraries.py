from __future__ import annotations

from .base import BaseRepository


class LibraryRepository(BaseRepository):
    def get(self, library_id=None, *, include_inactive: bool = True):
        if library_id:
            if include_inactive:
                return self.core.fetchall("SELECT * FROM Libraries WHERE library_id = ?", (library_id,))
            return self.core.fetchall(
                "SELECT * FROM Libraries WHERE library_id = ? AND COALESCE(status, 'ready') = 'ready'",
                (library_id,),
            )
        if include_inactive:
            return self.core.fetchall("SELECT * FROM Libraries")
        return self.core.fetchall("SELECT * FROM Libraries WHERE COALESCE(status, 'ready') = 'ready'")

    def get_id(self, name, *, include_inactive: bool = True):
        sql = "SELECT library_id FROM Libraries WHERE lower(library_name) = lower(?)"
        if not include_inactive:
            sql += " AND COALESCE(status, 'ready') = 'ready'"
        row = self.core.fetchone(sql, (name,))
        return row[0] if row else None

    def get_name(self, library_id):
        row = self.core.fetchone(
            "SELECT library_name FROM Libraries WHERE library_id = ?",
            (library_id,),
        )
        return row[0] if row else None

    def get_parent_id(self, library_id):
        row = self.core.fetchone(
            "SELECT parent_id FROM Libraries WHERE library_id = ?",
            (library_id,),
        )
        return row[0] if row and row[0] else None

    def get_status(self, library_id):
        row = self.core.fetchone("SELECT COALESCE(status, 'ready') FROM Libraries WHERE library_id = ?", (library_id,))
        return str(row[0]) if row and row[0] is not None else None

    def set_status(self, library_id, status: str):
        with self.core.transaction(operation=f"set library status {library_id}"):
            self.core.execute(
                "UPDATE Libraries SET status = ? WHERE library_id = ?",
                (str(status), int(library_id)),
            )
        return True

    def assert_has_parent(self, library_id):
        return self.get_parent_id(library_id)

    def add(self, library_name, taxid, size, location, parent_id=None, ref_accessions=None, odb_version=None, status: str = "ready"):
        with self.core.transaction(operation=f"add library {library_name}"):
            normalized_status = str(status or "ready")
            existing = self.core.fetchone(
                "SELECT library_id, library_name FROM Libraries WHERE lower(library_name) = lower(?)",
                (library_name,),
            )
            if existing:
                library_id, existing_name = existing
                self.core.execute(
                    """
                    UPDATE Libraries
                    SET taxid = ?, size = ?, location = ?, parent_id = ?, odb_version = COALESCE(?, odb_version), status = ?
                    WHERE library_id = ?
                    """,
                    (taxid, size, location, parent_id, odb_version, normalized_status, library_id),
                )
                library_name = existing_name
            else:
                self.core.execute(
                    """
                    INSERT INTO Libraries (library_name, taxid, size, location, parent_id, odb_version, status)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (library_name, taxid, size, location, parent_id, odb_version, normalized_status),
                )
                library_id = self.cursor.lastrowid

            if not library_id:
                raise RuntimeError("insert did not produce a library id")

            if ref_accessions is not None:
                self.core.execute("DELETE FROM Reference_Assemblies WHERE library_id = ?", (library_id,))
                if ref_accessions:
                    self.core.executemany(
                        "INSERT OR IGNORE INTO Reference_Assemblies (accession, library_id) VALUES (?, ?)",
                        [(acc, library_id) for acc in ref_accessions],
                    )

            if location and self.manager.storage.list_roots(kind="libraries"):
                self.manager.storage.bind_library_location(int(library_id), location, kind="libraries")
            self.manager.busco.invalidate_adjusted_results_for_library(
                int(library_id),
                reason="library_definition_updated",
            )
            return library_id

    def purge(self, library_id, *, preserve_orthofinder: bool = False):
        with self.core.transaction(operation=f"purge library data {library_id}"):
            self.core.execute("DELETE FROM Reference_Assemblies WHERE library_id = ?", (library_id,))
            self.core.execute("DELETE FROM Proteome_BlastDBs WHERE library_id = ?", (library_id,))
            self.core.execute("DELETE FROM BUSCO_Results WHERE library_id = ?", (library_id,))
            self.core.execute("DELETE FROM BUSCO_Family_Data WHERE library_id = ?", (library_id,))
            self.core.execute("DELETE FROM BUSCO_descriptions WHERE library_id = ?", (library_id,))

            if not preserve_orthofinder:
                rows = self.core.fetchall(
                    "SELECT orthofinder_id FROM OrthoFinder_Results WHERE library_id = ?",
                    (library_id,),
                )
                for (orthofinder_id,) in rows:
                    self.core.execute(
                        "DELETE FROM OrthoFinder_Accessions WHERE orthofinder_id = ?",
                        (orthofinder_id,),
                    )
                self.core.execute("DELETE FROM OrthoFinder_Results WHERE library_id = ?", (library_id,))
            self.core.execute("DELETE FROM BUSCO_Adjusted_Results WHERE library_id = ?", (library_id,))
            return True

    def set_binding(self, library_id, path, *, kind="libraries"):
        return self.manager.storage.bind_library_location(library_id, path, kind=kind)

    def resolve_path(self, library_id):
        return self.manager.storage.resolve_library_location(library_id)

    def get_reference_assemblies(self, library_id):
        rows = self.core.fetchall(
            "SELECT accession FROM Reference_Assemblies WHERE library_id = ? ORDER BY accession",
            (library_id,),
        )
        return [row[0] for row in rows]

    def get_by_reference_accessions(self, accessions):
        accs = list(dict.fromkeys(str(acc) for acc in (accessions or []) if acc))
        if not accs:
            return []
        placeholders = ",".join("?" for _ in accs)
        sql = f"""
            SELECT l.library_id, l.library_name, l.parent_id, COALESCE(l.status, 'ready') AS status
            FROM Libraries l
            JOIN Reference_Assemblies ra ON ra.library_id = l.library_id
            WHERE ra.accession IN ({placeholders})
            GROUP BY l.library_id
            HAVING COUNT(DISTINCT ra.accession) = ?
        """
        return self.core.fetchall(sql, tuple(accs + [len(accs)]))

    def add_reference_assemblies(self, library_id, accessions):
        accessions = list(accessions or [])
        if not accessions:
            return True
        with self.core.transaction(operation=f"add reference assemblies to library {library_id}"):
            placeholders = ",".join(["?"] * len(accessions))
            rows = self.core.fetchall(
                f"SELECT accession FROM Genome WHERE accession IN ({placeholders})",
                tuple(accessions),
            )
            present = [row[0] for row in rows]
            if not present:
                return True
            self.core.executemany(
                "INSERT OR IGNORE INTO Reference_Assemblies (accession, library_id) VALUES (?, ?)",
                [(acc, library_id) for acc in present],
            )
            return True

    def get_busco_descriptions(self, library_id, family_ids=None):
        params = [library_id]
        query = "SELECT family_id, library_id, description, link FROM BUSCO_descriptions WHERE library_id = ?"
        if family_ids:
            family_ids = list(family_ids)
            placeholders = ",".join("?" for _ in family_ids)
            query += f" AND family_id IN ({placeholders})"
            params.extend(family_ids)
        return self.core.fetchall(query, tuple(params))

    def add_busco_descriptions(self, rows):
        rows = list(rows)
        with self.core.transaction(operation="add BUSCO descriptions"):
            self.core.executemany(
                "INSERT OR REPLACE INTO BUSCO_descriptions (family_id, library_id, description, link) VALUES (?, ?, ?, ?)",
                rows,
            )
            touched = sorted({int(row[1]) for row in rows if len(row) >= 2 and row[1] is not None})
            for library_id in touched:
                self.manager.busco.invalidate_adjusted_results_for_library(
                    int(library_id),
                    reason="library_busco_descriptions_updated",
                )
            return True

    def update_size(self, library_id, size):
        with self.core.transaction(operation=f"update library size {library_id}"):
            self.core.execute(
                "UPDATE Libraries SET size = ? WHERE library_id = ?",
                (size, library_id),
            )
            self.manager.busco.invalidate_adjusted_results_for_library(
                int(library_id),
                reason="library_size_updated",
            )
            return True

    def get_view(self):
        return self.core.fetchall("SELECT * FROM Libraries_View WHERE COALESCE(status, 'ready') = 'ready'")
