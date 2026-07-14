from __future__ import annotations

import os
import shutil
import sqlite3
import tarfile
from contextlib import suppress
from datetime import datetime

from ..core import sqlite_busy_timeout_ms
from .base import BaseRepository

SP_TOKENS = {"sp", "sp.", "spp", "spp."}


class GenomeRepository(BaseRepository):
    @staticmethod
    def _normalize_dt(value):
        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d %H:%M:%S")
        if value is None:
            return None
        return str(value).split(".")[0]

    def get(self, accession):
        return self.core.fetchone("SELECT * FROM Genome WHERE accession = ?", (accession,))

    def get_many(self, genome_ids=None, status=None):
        if genome_ids is None:
            if status is None:
                return self.core.fetchall("SELECT * FROM Genome ORDER BY accession")
            return self.core.fetchall("SELECT * FROM Genome WHERE status >= ? ORDER BY accession", (status,))
        genome_ids = list(genome_ids or [])
        if not genome_ids:
            return []
        placeholders = ", ".join("?" for _ in genome_ids)
        if status is not None:
            return self.core.fetchall(
                f"SELECT * FROM Genome WHERE accession IN ({placeholders}) AND status >= ?",
                (*genome_ids, status),
            )
        return self.core.fetchall(
            f"SELECT * FROM Genome WHERE accession IN ({placeholders})",
            tuple(genome_ids),
        )

    def get_accessions_by_taxid(self, taxid: int, include_descendants: bool = True, status_min: int | None = 1, protein_only: bool = False):
        if include_descendants:
            base_cte = """
                WITH RECURSIVE desc(taxid) AS (
                    SELECT taxid FROM Taxonomy WHERE taxid = ?
                    UNION ALL
                    SELECT t.taxid FROM Taxonomy t JOIN desc d ON t.parent_taxid = d.taxid
                )
            """
            where_terms = ["g.taxid IN (SELECT taxid FROM desc)"]
            params = [taxid]
            sql_prefix = base_cte + "\nSELECT g.accession, g.taxid, g.protein FROM Genome g LEFT JOIN Hidden_Genomes h ON h.accession = g.accession"
        else:
            where_terms = ["g.taxid = ?"]
            params = [taxid]
            sql_prefix = "SELECT g.accession, g.taxid, g.protein FROM Genome g LEFT JOIN Hidden_Genomes h ON h.accession = g.accession"

        where_terms.append("h.accession IS NULL")
        if status_min is not None:
            where_terms.append("g.status >= ?")
            params.append(status_min)
        if protein_only:
            where_terms.append("g.protein = 1")

        sql = sql_prefix + "\nWHERE " + " AND ".join(where_terms)
        return self.core.fetchall(sql, tuple(params)) or []

    def get_downloaded(self):
        return self.core.fetchall(
            """
            SELECT * FROM Genome_quick_view
            WHERE accession NOT IN (SELECT accession FROM Hidden_Genomes)
            """
        )

    def get_path(self, accession):
        row = self.core.fetchone(
            "SELECT storage_root_id, relative_path, location FROM Genome WHERE accession = ?",
            (accession,),
        )
        if not row:
            return None
        return self.manager.storage.resolve_path(
            storage_root_id=row[0],
            relative_path=row[1],
            fallback_path=row[2],
        )

    def resolve_path(self, accession):
        return self.manager.storage.resolve_genome_location(accession)

    def set_binding(self, accession, path, *, kind="genomes"):
        return self.manager.storage.bind_genome_location(accession, path, kind=kind)

    def insert(self, data):
        with self.core.transaction(operation=f"insert genome {data.get('accession')}"):
            gd = dict(data)
            accession = gd.get("accession")
            if accession:
                # Genome has always been a convenient registration entrypoint,
                # so preserve that contract now that foreign keys are enforced.
                self.core.execute(
                    "INSERT OR IGNORE INTO Assembly (accession) VALUES (?)",
                    (accession,),
                )
            location = gd.get("location")
            if gd.get("dl_date") is not None:
                gd["dl_date"] = self._normalize_dt(gd.get("dl_date"))
            keys = list(gd.keys())
            columns = ", ".join(keys)
            placeholders = ", ".join(["?"] * len(keys))
            self.core.execute(
                f"INSERT OR IGNORE INTO Genome ({columns}) VALUES ({placeholders})",
                tuple(gd[k] for k in keys),
            )
            if location and gd.get("accession") and self.manager.storage.list_roots(kind="genomes"):
                self.manager.storage.bind_genome_location(gd["accession"], location, kind="genomes")
            return True

    def upsert(self, data):
        with self.core.transaction(operation=f"upsert genome {data.get('accession')}"):
            gd = dict(data)
            accession = gd.get("accession")
            if accession:
                self.core.execute(
                    "INSERT OR IGNORE INTO Assembly (accession) VALUES (?)",
                    (accession,),
                )
            if gd.get("dl_date") is not None:
                gd["dl_date"] = self._normalize_dt(gd.get("dl_date"))
            keys = list(gd.keys())
            columns = ", ".join(keys)
            placeholders = ", ".join(["?"] * len(keys))
            protected = {"dl_date", "location", "status", "protein"}
            update_clause = ", ".join(
                f"{k}=excluded.{k}" for k in keys if k != "accession" and k not in protected
            )
            self.core.execute(
                f"INSERT INTO Genome ({columns}) VALUES ({placeholders}) ON CONFLICT(accession) DO UPDATE SET {update_clause}",
                tuple(gd[k] for k in keys),
            )
            if gd.get("location") and gd.get("accession") and self.manager.storage.list_roots(kind="genomes"):
                self.manager.storage.bind_genome_location(gd["accession"], gd["location"], kind="genomes")
            return True

    def insert_assembly(self, data):
        with self.core.transaction(operation=f"insert assembly {data.get('accession')}"):
            values = dict(data)
            keys = list(values.keys())
            columns = ", ".join(keys)
            placeholders = ", ".join(["?"] * len(keys))
            self.core.execute(
                f"INSERT OR IGNORE INTO Assembly ({columns}) VALUES ({placeholders})",
                tuple(values[k] for k in keys),
            )
            return True

    def upsert_assembly(self, data):
        with self.core.transaction(operation=f"upsert assembly {data.get('accession')}"):
            values = dict(data)
            keys = list(values.keys())
            columns = ", ".join(keys)
            placeholders = ", ".join(["?"] * len(keys))
            update_clause = ", ".join(f"{k}=excluded.{k}" for k in keys if k != "accession")
            self.core.execute(
                f"INSERT INTO Assembly ({columns}) VALUES ({placeholders}) ON CONFLICT(accession) DO UPDATE SET {update_clause}",
                tuple(values[k] for k in keys),
            )
            return True

    def update_status(self, accession, status, details=None):
        with self.core.transaction(operation=f"update genome status {accession}"):
            if status == 1:
                if not details or len(details) != 3:
                    raise ValueError("Data must contain [dl_date, location, protein]")
                dl_date, location, protein = details
                dl_date_str = self._normalize_dt(dl_date)
                self.core.execute(
                    """
                    UPDATE Genome
                    SET status = ?, dl_date = ?, location = ?, protein = ?
                    WHERE accession = ?
                    """,
                    (status, dl_date_str, location, protein, accession),
                )
                if location and self.manager.storage.list_roots(kind="genomes"):
                    self.manager.storage.bind_genome_location(accession, location, kind="genomes")
            else:
                self.core.execute(
                    "UPDATE Genome SET status = ? WHERE accession = ?",
                    (status, accession),
                )
            return True

    def set_isoforms_cleaned(self, accession, value):
        with self.core.transaction(operation=f"set isoform-cleaned flag {accession}"):
            self.core.execute(
                "UPDATE Genome SET isoforms_cleaned = ? WHERE accession = ?",
                (1 if bool(value) else 0, accession),
            )
            return True

    def set_protein(self, accession, value):
        with self.core.transaction(operation=f"set protein flag {accession}"):
            self.core.execute(
                "UPDATE Genome SET protein = ? WHERE accession = ?",
                (1 if bool(value) else 0, accession),
            )
            return True

    def get_status(self, accession):
        row = self.core.fetchone("SELECT status FROM Genome WHERE accession = ?", (accession,))
        return int(row[0]) if row and row[0] is not None else None

    def hide(self, accession: str, status: str | None = None, reason: str | None = None):
        with self.core.transaction(operation=f"hide genome {accession}"):
            self.core.execute(
                "INSERT OR REPLACE INTO Hidden_Genomes (accession, status, reason) VALUES (?, ?, ?)",
                (accession, status, reason),
            )
            return True

    def get_column_map(self):
        cached = getattr(self.manager, "_genome_column_map", None)
        if cached:
            return cached
        rows = self.core.fetchall("PRAGMA table_info(Genome)")
        col_map = {row[1]: idx for idx, row in enumerate(rows)}
        self.manager._genome_column_map = col_map
        return col_map

    def get_taxid_for_genus(self, genus):
        token = str(genus or "").strip()
        if not token:
            return None
        row = self.core.fetchone(
            "SELECT taxid FROM Taxonomy WHERE name = ? AND rank = 'genus'",
            (token,),
        )
        return int(row[0]) if row and row[0] is not None else None

    def get_taxid_for_name(self, name):
        token = str(name or "").strip()
        if not token:
            return None
        row = self.core.fetchone(
            "SELECT taxid FROM Taxonomy WHERE lower(name) = lower(?) ORDER BY CASE WHEN rank = 'genus' THEN 0 ELSE 1 END, taxid LIMIT 1",
            (token,),
        )
        return int(row[0]) if row and row[0] is not None else None

    def get_taxid_for_genus_species(self, name):
        token = str(name or "").strip()
        if not token:
            return None
        parts = token.split()
        genus = parts[0]
        species = " ".join(parts[1:]).strip() if len(parts) > 1 else None
        if species and species.lower() in SP_TOKENS:
            species = None

        if species:
            row = self.core.fetchone(
                "SELECT taxid FROM Taxonomy WHERE name = ? AND rank = 'species'",
                (f"{genus} {species}".strip(),),
            )
            if row:
                return {"taxid": int(row[0])}

        genus_taxid = self.get_taxid_for_genus(genus)
        if genus_taxid is None and not species:
            fallback_taxid = self.get_taxid_for_name(genus)
            if fallback_taxid is not None:
                return {"taxid": fallback_taxid}
        if genus_taxid is None:
            return None
        if not species:
            return {"taxid": genus_taxid}

        row = self.core.fetchone(
            "SELECT taxid FROM Taxonomy WHERE name = ? AND parent_taxid = ? AND rank = 'species'",
            (species, genus_taxid),
        )
        if row:
            return {"taxid": int(row[0])}
        return None

    def get_taxid_for_species(self, genus, species=None):
        token = str(genus or "").strip()
        if not token:
            return None
        species_token = str(species or "").strip() or None
        if species_token and species_token.lower() in SP_TOKENS:
            species_token = None
        result = self.get_taxid_for_genus_species(
            f"{token} {species_token}".strip() if species_token else token
        )
        return int(result["taxid"]) if result and result.get("taxid") is not None else None

    def get_lineage_root_to_leaf(self, taxid):
        rows = self.core.fetchall(
            """
            WITH RECURSIVE lineage(taxid, name, rank, parent_taxid, depth, path) AS (
                SELECT taxid, name, rank, parent_taxid, 0 AS depth, printf('/%d/', taxid) AS path
                FROM Taxonomy
                WHERE taxid = ?
                UNION ALL
                SELECT t.taxid, t.name, t.rank, t.parent_taxid, l.depth + 1,
                       l.path || printf('%d/', t.taxid)
                FROM Taxonomy t
                JOIN lineage l ON t.taxid = l.parent_taxid
                WHERE instr(l.path, printf('/%d/', t.taxid)) = 0
            )
            SELECT taxid, name, rank, parent_taxid
            FROM lineage
            ORDER BY depth DESC
            """,
            (taxid,),
        )
        return [(r[0], r[1], r[2], r[3]) for r in rows]

    def insert_taxonomy_information(self, taxonomy):
        with self.core.transaction(operation="insert taxonomy information"):
            for taxid_val, data in taxonomy.items():
                if isinstance(data, dict):
                    keys = list(data.keys())
                    values = [data[k] for k in keys]
                else:
                    keys = ["name", "rank", "parent_taxid"]
                    values = list(data)
                columns = ", ".join(["taxid"] + keys)
                placeholders = ", ".join(["?"] * (len(keys) + 1))
                self.core.execute(
                    f"INSERT OR REPLACE INTO Taxonomy ({columns}) VALUES ({placeholders})",
                    [taxid_val] + values,
                )
            return True

    def insert_taxdump(self, path):
        self.core.execute("PRAGMA synchronous = NORMAL")
        self.core.execute("PRAGMA temp_store = MEMORY")

        if path.endswith(".tar.gz"):
            extract_dir = os.path.join(os.path.dirname(path), "taxdump_extracted")
            os.makedirs(extract_dir, exist_ok=True)
            with tarfile.open(path, "r:gz") as tar:
                tar.extractall(extract_dir)
            taxdump_dir = extract_dir
        else:
            taxdump_dir = path

        names_path = os.path.join(taxdump_dir, "names.dmp")
        nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
        if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
            raise FileNotFoundError("names.dmp or nodes.dmp not found in taxdump directory")

        scientific_names = {}
        with open(names_path, encoding="utf-8") as handle:
            for line in handle:
                parts = [x.strip() for x in line.split("|")]
                taxid, name_txt, _, name_class = parts[:4]
                if name_class == "scientific name":
                    scientific_names[int(taxid)] = name_txt

        insert_query = (
            "INSERT INTO Taxonomy (taxid, name, rank, parent_taxid) "
            "VALUES (?, ?, ?, ?) "
            "ON CONFLICT(taxid) DO UPDATE SET "
            "name=excluded.name, rank=excluded.rank, parent_taxid=excluded.parent_taxid"
        )

        self.core.execute("BEGIN")
        try:
            batch = []
            batch_size = 50000
            with open(nodes_path, encoding="utf-8") as handle:
                for line in handle:
                    parts = [x.strip() for x in line.split("|")]
                    try:
                        taxid = int(parts[0])
                        parent_taxid = int(parts[1])
                    except (ValueError, IndexError):
                        continue
                    rank = parts[2]
                    name = scientific_names.get(taxid)
                    if not name:
                        continue
                    parent_val = None if parent_taxid == taxid else parent_taxid
                    batch.append((taxid, name, rank, parent_val))
                    if len(batch) >= batch_size:
                        self.core.executemany(insert_query, batch)
                        batch.clear()
                if batch:
                    self.core.executemany(insert_query, batch)
            self.core.commit()
        except Exception:  # boundary: taxonomy bulk import must rollback and re-raise original failure.
            self.core.rollback()
            raise
        finally:
            if path.endswith(".tar.gz"):
                shutil.rmtree(extract_dir, ignore_errors=True)
            with suppress(sqlite3.Error):
                self.core.execute(f"PRAGMA busy_timeout = {sqlite_busy_timeout_ms()}")
                self.core.execute("PRAGMA read_uncommitted = true")
