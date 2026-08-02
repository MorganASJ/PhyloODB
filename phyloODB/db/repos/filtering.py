from __future__ import annotations

import datetime as datetime_module
import sqlite3
from typing import Any, Optional, Sequence

from .base import BaseRepository, transactional
from ..errors import RepositoryWriteError


class FilteringRepository(BaseRepository):
    def _invalidate_adjusted_results_for_accessions(
        self,
        *,
        target_library_id: Optional[int],
        accessions: Sequence[str],
        reason: str,
    ) -> None:
        if target_library_id is None:
            return
        normalized = sorted({str(accession) for accession in accessions if accession is not None})
        if not normalized:
            return
        self.manager.busco.invalidate_adjusted_results_for_library(
            int(target_library_id),
            accessions=normalized,
            reason=reason,
        )

    def get_blast_dbs(self, accession=None, library_id=None):
        if accession and library_id:
            return self.core.fetchall("SELECT * FROM Proteome_BlastDBs WHERE accession = ? AND library_id = ?", (accession, library_id))
        if accession:
            return self.core.fetchall("SELECT * FROM Proteome_BlastDBs WHERE accession = ?", (accession,))
        if library_id:
            return self.core.fetchall("SELECT * FROM Proteome_BlastDBs WHERE library_id = ?", (library_id,))
        return self.core.fetchall("SELECT * FROM Proteome_BlastDBs")

    @transactional("add proteome BLAST database")
    def add_proteome_blastdb(self, accession, location, library_id=None):
        try:
            self.core.execute(
                "INSERT OR REPLACE INTO Proteome_BlastDBs (accession, library_id, location) VALUES (?, ?, ?)",
                (accession, library_id, location),
            )
            self.conn.commit()
            return self.cursor.lastrowid
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding proteome BLASTDB: {exc}") from exc

    @transactional("delete proteome BLAST databases")
    def delete_proteome_blastdb_ids(self, blastdb_ids: Sequence[int]) -> int:
        ids = [int(value) for value in blastdb_ids if value is not None]
        if not ids:
            return 0
        placeholders = ",".join("?" for _ in ids)
        self.core.execute(f"DELETE FROM Proteome_BlastDBs WHERE blastdb_id IN ({placeholders})", tuple(ids))
        deleted = int(self.cursor.rowcount or 0)
        self.conn.commit()
        return deleted

    def supports_target_library(self):
        cached = getattr(self.manager, "_paralog_filtering_has_target", None)
        if cached is not None:
            return cached
        rows = self.core.fetchall("PRAGMA table_info(Paralog_Filtering)")
        columns = {row[1] for row in rows}
        cached = "target_library_id" in columns
        self.manager._paralog_filtering_has_target = cached
        return cached

    def supports_paralog_runs(self):
        cached = getattr(self.manager, "_paralog_filtering_has_runs", None)
        if cached is not None:
            return cached
        rows = self.core.fetchall("PRAGMA table_info(Paralog_Filtering)")
        columns = {row[1] for row in rows}
        cached = "run_id" in columns
        self.manager._paralog_filtering_has_runs = cached
        return cached

    def _active_paralog_env_key(self, target_library_id: int) -> str:
        return f"ACTIVE_PARALOG_RUN_{int(target_library_id)}"

    def _latest_paralog_run_id(self, target_library_id: int) -> Optional[str]:
        row = self.core.fetchone(
            """
            SELECT run_id
            FROM Paralog_Filtering_Runs
            WHERE target_library_id = ?
            ORDER BY COALESCE(date, '') DESC, rowid DESC
            LIMIT 1
            """,
            (int(target_library_id),),
        )
        return str(row[0]) if row and row[0] is not None else None

    def resolve_paralog_run_id(self, *, target_library_id: int, run_id: Optional[str] = None) -> Optional[str]:
        if run_id is not None:
            return str(run_id)
        if not self.supports_paralog_runs():
            return None
        payload = self.manager.env.get(self._active_paralog_env_key(int(target_library_id)))
        active_run_id = payload.get("run_id") if isinstance(payload, dict) else None
        if active_run_id:
            return str(active_run_id)
        return self._latest_paralog_run_id(int(target_library_id))

    @transactional("add paralog filtering run")
    def add_paralog_filtering_run(
        self,
        run_id,
        target_library_id,
        busco_library_id,
        targets_json,
        accessions_json,
        ref_accessions_json,
        selection_mode,
        selection_params_json,
        config_signature=None,
        run_label=None,
        report_dir=None,
        datetime=None,
    ):
        try:
            self.core.execute(
                """
                INSERT OR REPLACE INTO Paralog_Filtering_Runs (
                    run_id, target_library_id, busco_library_id, targets_json, accessions_json,
                    ref_accessions_json, selection_mode, selection_params_json, config_signature,
                    run_label, report_dir, date
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    run_id,
                    target_library_id,
                    busco_library_id,
                    targets_json,
                    accessions_json,
                    ref_accessions_json,
                    selection_mode,
                    selection_params_json,
                    config_signature,
                    run_label,
                    report_dir,
                    datetime,
                ),
            )
            self.conn.commit()
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding paralog filtering run: {exc}") from exc

    def get_paralog_filtering_run(self, run_id):
        return self.core.fetchone(
            """
            SELECT run_id, target_library_id, busco_library_id, targets_json, accessions_json,
                   ref_accessions_json, selection_mode, selection_params_json, config_signature,
                   run_label, report_dir, date
            FROM Paralog_Filtering_Runs
            WHERE run_id = ?
            """,
            (str(run_id),),
        )

    def list_paralog_filtering_runs(self, *, target_library_id: Optional[int] = None):
        sql = """
            SELECT run_id, target_library_id, busco_library_id, targets_json, accessions_json,
                   ref_accessions_json, selection_mode, selection_params_json, config_signature,
                   run_label, report_dir, date
            FROM Paralog_Filtering_Runs
        """
        params = []
        if target_library_id is not None:
            sql += " WHERE target_library_id = ?"
            params.append(int(target_library_id))
        sql += " ORDER BY COALESCE(date, '') DESC, rowid DESC"
        return self.core.fetchall(sql, tuple(params))

    def _chunked(self, values, size: int = 900):
        if not values:
            return
        for idx in range(0, len(values), size):
            yield values[idx:idx + size]

    def _get_library_row(self, library_id: int):
        return self.core.fetchone(
            "SELECT library_id, library_name, taxid, location, size, odb_version, parent_id FROM Libraries WHERE library_id = ?",
            (int(library_id),),
        )

    def _get_library_size(self, library_id: int) -> int:
        size = 0
        row = self._get_library_row(int(library_id))
        if row and row[4] is not None:
            try:
                size = int(row[4])
            except (TypeError, ValueError):
                size = 0
        if size < 0:
            size = 0
        actual = self.core.fetchone(
            "SELECT COUNT(*) FROM BUSCO_descriptions WHERE library_id = ?",
            (int(library_id),),
        )
        actual_size = int(actual[0]) if actual and actual[0] is not None else 0
        if actual_size > size:
            size = actual_size
        return size

    def latest_decont_summary(self, *, target_library_id: int, accessions=None, run_id: Optional[str] = None):
        sql = "SELECT accession, run_id, decision, date FROM Decontamination_Summary WHERE target_library_id = ?"
        params = [target_library_id]
        if run_id:
            sql += " AND run_id = ?"
            params.append(run_id)
        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            sql += f" AND accession IN ({placeholders})"
            params.extend(accessions)
        rows = self.core.fetchall(sql, tuple(params))
        if run_id:
            return {str(acc): (str(run_id), decision, date) for acc, _rid, decision, date in rows}
        latest = {}
        for acc, rid, decision, date in rows:
            key = str(acc)
            prev = latest.get(key)
            if prev is None or (date or "") > (prev[2] or ""):
                latest[key] = (str(rid), decision, date)
        return latest

    def paralog_hidden_counts(
        self,
        *,
        target_library_id: int,
        busco_library_id: int,
        family_library_id: Optional[int] = None,
        accessions=None,
        run_id: Optional[str] = None,
    ):
        supports_target = self.supports_target_library()
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(target_library_id), run_id=run_id)
        sql = """
            SELECT pf.accession,
                   COUNT(DISTINCT CASE WHEN pf.clean = 0 THEN pf.family_id END) AS hidden_cnt,
                   COUNT(DISTINCT pf.family_id) AS total_cnt
            FROM Paralog_Filtering pf
            JOIN BUSCO_descriptions bd
              ON bd.family_id = pf.family_id AND bd.library_id = ?
            WHERE pf.library_id = ?
        """
        params = [int(family_library_id or busco_library_id), busco_library_id]
        if supports_target:
            sql += " AND pf.target_library_id = ?"
            params.append(target_library_id)
        if self.supports_paralog_runs() and paralog_run_id:
            sql += " AND pf.run_id = ?"
            params.append(paralog_run_id)
        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            sql += f" AND pf.accession IN ({placeholders})"
            params.extend(accessions)
        sql += " GROUP BY pf.accession"
        rows = self.core.fetchall(sql, tuple(params))
        return {str(acc): (int(hidden or 0), int(total or 0)) for acc, hidden, total in rows}

    @transactional("delete paralog filtering records")
    def delete_paralog_records(self, accession: str, *, busco_library_id: Optional[int] = None, target_library_id: Optional[int] = None, run_id: Optional[str] = None):
        supports_target = self.supports_target_library()
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(target_library_id), run_id=run_id) if target_library_id is not None else run_id
        try:
            query = "DELETE FROM Paralog_Filtering WHERE accession = ?"
            params = [accession]
            if busco_library_id is not None:
                query += " AND library_id = ?"
                params.append(busco_library_id)
            if supports_target and target_library_id is not None:
                query += " AND target_library_id = ?"
                params.append(target_library_id)
            if self.supports_paralog_runs() and paralog_run_id is not None:
                query += " AND run_id = ?"
                params.append(paralog_run_id)
            self.core.execute(query, tuple(params))
            self.conn.commit()
            if target_library_id is not None:
                self.manager.busco.invalidate_adjusted_results_for_library(
                    int(target_library_id),
                    accessions=[str(accession)],
                    reason="paralog_filtering_deleted",
                )
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(
                f"Error deleting Paralog_Filtering records for {accession} "
                f"(lib={busco_library_id}, target={target_library_id}): {exc}"
            ) from exc

    def get_paralog_results(self, target_library_id=None, accession=None, run_id=None, include_metadata: bool = False):
        supports_target = self.supports_target_library()
        supports_runs = self.supports_paralog_runs()
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(target_library_id), run_id=run_id) if target_library_id is not None else (str(run_id) if run_id is not None else None)
        params = []
        conditions = []
        if supports_target:
            if include_metadata and supports_runs:
                query = """
                    SELECT family_id, library_id, target_library_id, accession, clean, run_id,
                           selected_ref_count, selection_threshold, reused, reason_code, selection_signature
                    FROM Paralog_Filtering
                """
            else:
                query = "SELECT family_id, library_id, target_library_id, accession, clean FROM Paralog_Filtering"
            if target_library_id is not None:
                conditions.append("target_library_id = ?")
                params.append(target_library_id)
        else:
            query = "SELECT family_id, library_id, accession, clean FROM Paralog_Filtering"
            if target_library_id is not None:
                conditions.append("library_id = ?")
                params.append(target_library_id)
        if supports_runs and paralog_run_id is not None:
            conditions.append("run_id = ?")
            params.append(paralog_run_id)
        if accession is not None:
            conditions.append("accession = ?")
            params.append(accession)
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        rows = self.core.fetchall(query, tuple(params))
        if supports_target:
            return rows
        return [(family_id, library_id, library_id, accession, clean) for family_id, library_id, accession, clean in rows]

    def get_paralog_results_history(self, *, target_library_id: int, accession: Optional[str] = None):
        if not self.supports_paralog_runs():
            return self.get_paralog_results(target_library_id=target_library_id, accession=accession, include_metadata=True)
        sql = """
            SELECT family_id, library_id, target_library_id, accession, clean, run_id,
                   selected_ref_count, selection_threshold, reused, reason_code, selection_signature
            FROM Paralog_Filtering
            WHERE target_library_id = ?
        """
        params = [int(target_library_id)]
        if accession is not None:
            sql += " AND accession = ?"
            params.append(str(accession))
        sql += " ORDER BY COALESCE(date, '') DESC, rowid DESC"
        return self.core.fetchall(sql, tuple(params))

    def get_paralog_unclean_families(self, *, target_library_id: int, busco_library_id: int, accessions=None, run_id: Optional[str] = None):
        supports_target = self.supports_target_library()
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(target_library_id), run_id=run_id)
        sql = """
            SELECT pf.accession, pf.family_id
            FROM Paralog_Filtering pf
            JOIN BUSCO_descriptions bd
              ON bd.family_id = pf.family_id AND bd.library_id = ?
            WHERE pf.library_id = ?
              AND pf.clean = 0
        """
        params = [busco_library_id, busco_library_id]
        if supports_target:
            sql += " AND pf.target_library_id = ?"
            params.append(target_library_id)
        if self.supports_paralog_runs() and paralog_run_id:
            sql += " AND pf.run_id = ?"
            params.append(paralog_run_id)
        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            sql += f" AND pf.accession IN ({placeholders})"
            params.extend(accessions)
        rows = self.core.fetchall(sql, tuple(params))
        return {(str(acc), str(fam)) for acc, fam in rows}

    def get_hidden_paralog_sequence_ids(
        self,
        *,
        target_library_id: int,
        busco_library_id: int,
        family_id: str,
        source_run_ids: Optional[Sequence[int]] = None,
        run_id: Optional[str] = None,
    ) -> set[str]:
        family_token = str(family_id or "").strip()
        if not family_token:
            return set()
        supports_target = self.supports_target_library()

        def _collect(candidate_target_library_id: int) -> set[str]:
            paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(candidate_target_library_id), run_id=run_id)
            def _collect_once(filter_by_source_runs: bool) -> set[str]:
                hidden_ids: set[str] = set()

                copy_sql = """
                    SELECT DISTINCT query_id, query_header
                    FROM Paralog_Filtering_Copy
                    WHERE family_id = ?
                      AND library_id = ?
                      AND clean = 0
                """
                copy_params: list[object] = [family_token, int(busco_library_id)]
                if supports_target:
                    copy_sql += " AND target_library_id = ?"
                    copy_params.append(int(candidate_target_library_id))
                if self.supports_paralog_runs() and paralog_run_id:
                    copy_sql += " AND run_id = ?"
                    copy_params.append(paralog_run_id)
                if filter_by_source_runs and source_run_ids:
                    placeholders = ",".join("?" for _ in source_run_ids)
                    copy_sql += f" AND busco_run_id IN ({placeholders})"
                    copy_params.extend(int(value) for value in source_run_ids)
                for query_id, query_header in self.core.fetchall(copy_sql, tuple(copy_params)) or []:
                    for token in (query_id, query_header):
                        text = str(token or "").strip()
                        if not text:
                            continue
                        hidden_ids.add(text)
                        hidden_ids.add(text.split()[0])

                single_sql = """
                    SELECT DISTINCT d.sequence
                    FROM Paralog_Filtering pf
                    JOIN BUSCO_Run_Family_Data d
                      ON d.family_id = pf.family_id
                     AND d.accession = pf.accession
                     AND d.library_id = pf.library_id
                     AND d.run_id = pf.busco_run_id
                    WHERE pf.family_id = ?
                      AND pf.library_id = ?
                      AND pf.clean = 0
                """
                single_params: list[object] = [family_token, int(busco_library_id)]
                if supports_target:
                    single_sql += " AND pf.target_library_id = ?"
                    single_params.append(int(candidate_target_library_id))
                if self.supports_paralog_runs() and paralog_run_id:
                    single_sql += " AND pf.run_id = ?"
                    single_params.append(paralog_run_id)
                if filter_by_source_runs and source_run_ids:
                    placeholders = ",".join("?" for _ in source_run_ids)
                    single_sql += f" AND pf.busco_run_id IN ({placeholders})"
                    single_params.extend(int(value) for value in source_run_ids)
                for (sequence,) in self.core.fetchall(single_sql, tuple(single_params)) or []:
                    text = str(sequence or "").strip()
                    if text:
                        hidden_ids.add(text)
                return hidden_ids

            hidden_ids = _collect_once(filter_by_source_runs=True)
            if not source_run_ids:
                return hidden_ids
            return hidden_ids | _collect_once(filter_by_source_runs=False)

        primary_hidden_ids = _collect(int(target_library_id))
        if int(target_library_id) == int(busco_library_id):
            return primary_hidden_ids
        parent_hidden_ids = _collect(int(busco_library_id))
        return primary_hidden_ids | parent_hidden_ids

    @transactional("add decontamination run")
    def add_decontamination_run(self, run_id, target_library_id, busco_library_id, targets_json, refs_json, params_json, config_signature=None, run_label=None, datetime=None):
        try:
            self.core.execute(
                """
                INSERT OR REPLACE INTO Decontamination_Runs (
                    run_id, target_library_id, busco_library_id, targets_json, refs_json, params_json,
                    config_signature, run_label, date
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (run_id, target_library_id, busco_library_id, targets_json, refs_json, params_json, config_signature, run_label, datetime),
            )
            self.conn.commit()
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding decontamination run: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding decontamination run: {exc}") from exc

    def get_decontamination_run(self, run_id):
        return self.core.fetchone(
            """
            SELECT run_id, target_library_id, busco_library_id, targets_json, refs_json, params_json, config_signature, run_label, date
            FROM Decontamination_Runs
            WHERE run_id = ?
            """,
            (run_id,),
        )

    def get_decontamination_votes(self, target_library_id=None, accession=None, run_id=None, busco_run_id=None):
        params = []
        conditions = []
        query = """
            SELECT family_id, busco_library_id, target_library_id, accession, run_id, expected_taxid,
                   best_taxid, runner_taxid, rank, best_bitscore, delta_bitscore, decision, top_hits_json,
                   busco_run_id
            FROM Decontamination_Busco_Votes
        """
        if run_id is not None:
            conditions.append("run_id = ?")
            params.append(run_id)
        if busco_run_id is not None:
            conditions.append("busco_run_id = ?")
            params.append(int(busco_run_id))
        if target_library_id is not None:
            conditions.append("target_library_id = ?")
            params.append(target_library_id)
        if accession is not None:
            conditions.append("accession = ?")
            params.append(accession)
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        return self.core.fetchall(query, tuple(params))

    def add_decontamination_vote(
        self,
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
        top_hits_json=None,
        busco_run_id=None,
        datetime=None,
    ):
        return self.add_decontamination_votes(
            [
                {
                    "family_id": family_id,
                    "busco_library_id": busco_library_id,
                    "target_library_id": target_library_id,
                    "accession": accession,
                    "run_id": run_id,
                    "busco_run_id": busco_run_id,
                    "expected_taxid": expected_taxid,
                    "best_taxid": best_taxid,
                    "runner_taxid": runner_taxid,
                    "rank": rank,
                    "best_bitscore": best_bitscore,
                    "delta_bitscore": delta_bitscore,
                    "decision": decision,
                    "top_hits_json": top_hits_json,
                    "date": datetime,
                }
            ]
        )

    @transactional("add decontamination votes")
    def add_decontamination_votes(self, rows: Sequence[dict[str, Any]]) -> bool:
        if not rows:
            return True
        payload = []
        affected_accessions: set[str] = set()
        target_library_id = None
        for row in rows:
            target_library_id = row.get("target_library_id", target_library_id)
            accession = str(row["accession"])
            affected_accessions.add(accession)
            timestamp = row.get("date")
            if timestamp is None:
                timestamp = datetime_module.datetime.utcnow()
            payload.append(
                (
                    row["family_id"],
                    row["busco_library_id"],
                    row["target_library_id"],
                    accession,
                    row["run_id"],
                    row.get("busco_run_id"),
                    row.get("expected_taxid"),
                    row.get("best_taxid"),
                    row.get("runner_taxid"),
                    row["rank"],
                    row.get("best_bitscore"),
                    row.get("delta_bitscore"),
                    row["decision"],
                    row.get("top_hits_json"),
                    timestamp,
                )
            )
        try:
            self.core.executemany(
                """
                INSERT OR REPLACE INTO Decontamination_Busco_Votes (
                    family_id, busco_library_id, target_library_id, accession, run_id, busco_run_id,
                    expected_taxid, best_taxid, runner_taxid, rank, best_bitscore,
                    delta_bitscore, decision, top_hits_json, date
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                payload,
            )
            self.conn.commit()
            self._invalidate_adjusted_results_for_accessions(
                target_library_id=target_library_id,
                accessions=sorted(affected_accessions),
                reason="decontamination_votes_updated",
            )
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding decontamination votes: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding decontamination votes: {exc}") from exc

    def get_decontamination_summary(self, target_library_id=None, accession=None, run_id=None, busco_run_id=None):
        params = []
        conditions = []
        query = """
            SELECT accession, target_library_id, busco_library_id, run_id, expected_taxid, majority_taxid,
                   rank, buscos_tested, buscos_supporting, buscos_outside, off_clade_fraction, decision, params_json,
                   busco_run_id
            FROM Decontamination_Summary
        """
        if run_id is not None:
            conditions.append("run_id = ?")
            params.append(run_id)
        if busco_run_id is not None:
            conditions.append("busco_run_id = ?")
            params.append(int(busco_run_id))
        if target_library_id is not None:
            conditions.append("target_library_id = ?")
            params.append(target_library_id)
        if accession is not None:
            conditions.append("accession = ?")
            params.append(accession)
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        return self.core.fetchall(query, tuple(params))

    def add_decontamination_summary(
        self,
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
        params_json=None,
        busco_run_id=None,
        datetime=None,
    ):
        return self.add_decontamination_summaries(
            [
                {
                    "accession": accession,
                    "target_library_id": target_library_id,
                    "busco_library_id": busco_library_id,
                    "run_id": run_id,
                    "busco_run_id": busco_run_id,
                    "expected_taxid": expected_taxid,
                    "majority_taxid": majority_taxid,
                    "rank": rank,
                    "buscos_tested": buscos_tested,
                    "buscos_supporting": buscos_supporting,
                    "buscos_outside": buscos_outside,
                    "off_clade_fraction": off_clade_fraction,
                    "decision": decision,
                    "params_json": params_json,
                    "date": datetime,
                }
            ]
        )

    @transactional("add decontamination summaries")
    def add_decontamination_summaries(self, rows: Sequence[dict[str, Any]]) -> bool:
        if not rows:
            return True
        payload = []
        affected_accessions: set[str] = set()
        target_library_id = None
        for row in rows:
            target_library_id = row.get("target_library_id", target_library_id)
            accession = str(row["accession"])
            affected_accessions.add(accession)
            timestamp = row.get("date")
            if timestamp is None:
                timestamp = datetime_module.datetime.utcnow()
            payload.append(
                (
                    accession,
                    row["target_library_id"],
                    row["busco_library_id"],
                    row["run_id"],
                    row.get("busco_run_id"),
                    row.get("expected_taxid"),
                    row.get("majority_taxid"),
                    row["rank"],
                    row["buscos_tested"],
                    row["buscos_supporting"],
                    row["buscos_outside"],
                    row.get("off_clade_fraction"),
                    row["decision"],
                    row.get("params_json"),
                    timestamp,
                )
            )
        try:
            self.core.executemany(
                """
                INSERT OR REPLACE INTO Decontamination_Summary (
                    accession, target_library_id, busco_library_id, run_id, busco_run_id, expected_taxid, majority_taxid,
                    rank, buscos_tested, buscos_supporting, buscos_outside, off_clade_fraction,
                    decision, params_json, date
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                payload,
            )
            self.conn.commit()
            self._invalidate_adjusted_results_for_accessions(
                target_library_id=target_library_id,
                accessions=sorted(affected_accessions),
                reason="decontamination_summary_updated",
            )
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding decontamination summaries: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding decontamination summaries: {exc}") from exc

    def add_paralog_filtering_copy_result(
        self,
        family_id,
        library_id,
        target_library_id,
        accession,
        run_id,
        query_id,
        query_header,
        query_status,
        clean,
        selected_ref_count=None,
        reused=False,
        reason_code=None,
        selection_signature=None,
        busco_run_id=None,
        datetime=None,
    ):
        return self.add_paralog_filtering_copy_results(
            [
                {
                    "family_id": family_id,
                    "library_id": library_id,
                    "target_library_id": target_library_id,
                    "accession": accession,
                    "run_id": run_id,
                    "busco_run_id": busco_run_id,
                    "query_id": query_id,
                    "query_header": query_header,
                    "query_status": query_status,
                    "clean": clean,
                    "selected_ref_count": selected_ref_count,
                    "reused": reused,
                    "reason_code": reason_code,
                    "selection_signature": selection_signature,
                    "date": datetime,
                }
            ]
        )

    @transactional("add paralog filtering copy results")
    def add_paralog_filtering_copy_results(self, rows: Sequence[dict[str, Any]]) -> bool:
        if not rows:
            return True
        payload = []
        affected_accessions: set[str] = set()
        target_library_id = None
        for row in rows:
            target_library_id = row.get("target_library_id", target_library_id)
            accession = str(row["accession"])
            affected_accessions.add(accession)
            timestamp = row.get("date")
            if timestamp is None:
                timestamp = datetime_module.datetime.utcnow()
            payload.append(
                (
                    row["family_id"],
                    row["library_id"],
                    row["target_library_id"],
                    accession,
                    row["run_id"],
                    row.get("busco_run_id"),
                    row["query_id"],
                    row.get("query_header"),
                    row.get("query_status"),
                    row["clean"],
                    row.get("selected_ref_count"),
                    int(bool(row.get("reused"))),
                    row.get("reason_code"),
                    row.get("selection_signature"),
                    timestamp,
                )
            )
        try:
            self.core.executemany(
                """
                INSERT OR REPLACE INTO Paralog_Filtering_Copy (
                    family_id, library_id, target_library_id, accession, run_id, busco_run_id, query_id,
                    query_header, query_status, clean, selected_ref_count, reused,
                    reason_code, selection_signature, date
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                payload,
            )
            self.conn.commit()
            self._invalidate_adjusted_results_for_accessions(
                target_library_id=target_library_id,
                accessions=sorted(affected_accessions),
                reason="paralog_filtering_copy_updated",
            )
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding Paralog_Filtering_Copy results: {exc}") from exc

    def add_paralog_filtering_result(
        self,
        family_id,
        busco_library_id,
        target_library_id,
        accession,
        run_id,
        clean,
        selected_ref_count=None,
        selection_threshold=None,
        reused=False,
        reason_code=None,
        selection_signature=None,
        busco_run_id=None,
        datetime=None,
    ):
        return self.add_paralog_filtering_results(
            [
                {
                    "family_id": family_id,
                    "busco_library_id": busco_library_id,
                    "target_library_id": target_library_id,
                    "accession": accession,
                    "run_id": run_id,
                    "busco_run_id": busco_run_id,
                    "clean": clean,
                    "selected_ref_count": selected_ref_count,
                    "selection_threshold": selection_threshold,
                    "reused": reused,
                    "reason_code": reason_code,
                    "selection_signature": selection_signature,
                    "date": datetime,
                }
            ]
        )

    @transactional("add paralog filtering results")
    def add_paralog_filtering_results(self, rows: Sequence[dict[str, Any]]) -> bool:
        if not rows:
            return True
        supports_target = self.supports_target_library()
        payload = []
        affected_accessions: set[str] = set()
        target_library_id = None
        for row in rows:
            target_library_id = row.get("target_library_id", target_library_id)
            accession = str(row["accession"])
            affected_accessions.add(accession)
            timestamp = row.get("date")
            if timestamp is None:
                timestamp = datetime_module.datetime.utcnow()
            if supports_target:
                payload.append(
                    (
                        row["family_id"],
                        row["busco_library_id"],
                        row["target_library_id"],
                        accession,
                        row["run_id"],
                        row.get("busco_run_id"),
                        row["clean"],
                        row.get("selected_ref_count"),
                        row.get("selection_threshold"),
                        int(bool(row.get("reused"))),
                        row.get("reason_code"),
                        row.get("selection_signature"),
                        timestamp,
                    )
                )
            else:
                payload.append(
                    (
                        row["family_id"],
                        row["busco_library_id"],
                        accession,
                        row["clean"],
                        timestamp,
                    )
                )
        try:
            if supports_target:
                self.core.executemany(
                    """
                    INSERT OR REPLACE INTO Paralog_Filtering
                        (family_id, library_id, target_library_id, accession, run_id, busco_run_id, clean,
                         selected_ref_count, selection_threshold, reused, reason_code,
                         selection_signature, date)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    payload,
                )
            else:
                self.core.executemany(
                    """
                    INSERT OR REPLACE INTO Paralog_Filtering
                        (family_id, library_id, accession, clean, date)
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    payload,
                )
            self.conn.commit()
            self._invalidate_adjusted_results_for_accessions(
                target_library_id=target_library_id,
                accessions=sorted(affected_accessions),
                reason="paralog_filtering_updated",
            )
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding Paralog Filtering results: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding Paralog Filtering results: {exc}") from exc

    def add_decontamination_copy_vote(
        self,
        family_id,
        busco_library_id,
        target_library_id,
        accession,
        run_id,
        query_id,
        query_header,
        query_status,
        expected_taxid,
        best_taxid,
        runner_taxid,
        rank,
        best_bitscore,
        delta_bitscore,
        decision,
        top_hits_json=None,
        busco_run_id=None,
        datetime=None,
    ):
        return self.add_decontamination_copy_votes(
            [
                {
                    "family_id": family_id,
                    "busco_library_id": busco_library_id,
                    "target_library_id": target_library_id,
                    "accession": accession,
                    "run_id": run_id,
                    "busco_run_id": busco_run_id,
                    "query_id": query_id,
                    "query_header": query_header,
                    "query_status": query_status,
                    "expected_taxid": expected_taxid,
                    "best_taxid": best_taxid,
                    "runner_taxid": runner_taxid,
                    "rank": rank,
                    "best_bitscore": best_bitscore,
                    "delta_bitscore": delta_bitscore,
                    "decision": decision,
                    "top_hits_json": top_hits_json,
                    "date": datetime,
                }
            ]
        )

    @transactional("add decontamination copy votes")
    def add_decontamination_copy_votes(self, rows: Sequence[dict[str, Any]]) -> bool:
        if not rows:
            return True
        payload = []
        affected_accessions: set[str] = set()
        target_library_id = None
        for row in rows:
            target_library_id = row.get("target_library_id", target_library_id)
            accession = str(row["accession"])
            affected_accessions.add(accession)
            timestamp = row.get("date")
            if timestamp is None:
                timestamp = datetime_module.datetime.utcnow()
            payload.append(
                (
                    row["family_id"],
                    row["busco_library_id"],
                    row["target_library_id"],
                    accession,
                    row["run_id"],
                    row.get("busco_run_id"),
                    row["query_id"],
                    row.get("query_header"),
                    row.get("query_status"),
                    row.get("expected_taxid"),
                    row.get("best_taxid"),
                    row.get("runner_taxid"),
                    row["rank"],
                    row.get("best_bitscore"),
                    row.get("delta_bitscore"),
                    row["decision"],
                    row.get("top_hits_json"),
                    timestamp,
                )
            )
        try:
            self.core.executemany(
                """
                INSERT OR REPLACE INTO Decontamination_Busco_Copy_Votes (
                    family_id, busco_library_id, target_library_id, accession, run_id, busco_run_id,
                    query_id, query_header, query_status, expected_taxid, best_taxid,
                    runner_taxid, rank, best_bitscore, delta_bitscore, decision, top_hits_json, date
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                payload,
            )
            self.conn.commit()
            self._invalidate_adjusted_results_for_accessions(
                target_library_id=target_library_id,
                accessions=sorted(affected_accessions),
                reason="decontamination_copy_votes_updated",
            )
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding decontamination copy votes: {exc}") from exc

    def get_latest_decontamination_summary(self, *, target_library_id: int, accessions=None, run_id: Optional[str] = None):
        return self.latest_decont_summary(target_library_id=target_library_id, accessions=accessions, run_id=run_id)

    def get_latest_decontamination_summary_with_fallback(self, *, target_library_id: int, parent_library_id: Optional[int] = None, accessions=None, run_id: Optional[str] = None):
        primary = self.latest_decont_summary(target_library_id=target_library_id, accessions=accessions, run_id=run_id)
        if not parent_library_id:
            return primary
        fallback = self.latest_decont_summary(target_library_id=parent_library_id, accessions=accessions, run_id=run_id)
        merged = dict(primary)
        for acc, entry in fallback.items():
            if acc not in merged:
                merged[acc] = entry
        return merged

    def get_decontamination_unclean_families(self, *, target_library_id: int, busco_library_id: int, accessions_by_run: dict[str, list[str]], supported_decisions):
        unclean = set()
        supported = tuple(str(dec) for dec in supported_decisions if dec is not None) or ("support", "weak")
        supported_placeholders = ",".join("?" for _ in supported)
        for run_id, acc_list in accessions_by_run.items():
            if not acc_list:
                continue
            placeholders = ",".join("?" for _ in acc_list)
            sql = f"""
                SELECT v.accession, v.family_id
                FROM Decontamination_Busco_Votes v
                WHERE v.target_library_id = ?
                  AND v.busco_library_id = ?
                  AND v.run_id = ?
                  AND v.accession IN ({placeholders})
                  AND COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders})
            """
            params = [target_library_id, busco_library_id, str(run_id), *acc_list, *supported]
            rows = self.core.fetchall(sql, tuple(params))
            for acc, fam in rows:
                unclean.add((str(acc), str(fam)))
        return unclean

    def get_decontamination_accessions(self, *, target_library_id: int, run_id: Optional[str] = None, accessions=None):
        sql = "SELECT DISTINCT accession FROM Decontamination_Summary WHERE target_library_id = ?"
        params = [target_library_id]
        if run_id:
            sql += " AND run_id = ?"
            params.append(run_id)
        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            sql += f" AND accession IN ({placeholders})"
            params.extend(accessions)
        rows = self.core.fetchall(sql, tuple(params))
        return {str(row[0]) for row in rows}

    def get_paralog_filtering_accessions(self, *, target_library_id: int, busco_library_id: int, accessions=None, run_id: Optional[str] = None):
        supports_target = self.supports_target_library()
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(target_library_id), run_id=run_id)
        sql = """
            SELECT DISTINCT pf.accession
            FROM Paralog_Filtering pf
            JOIN BUSCO_descriptions bd
              ON bd.family_id = pf.family_id AND bd.library_id = ?
            WHERE pf.library_id = ?
        """
        params = [busco_library_id, busco_library_id]
        if supports_target:
            sql += " AND pf.target_library_id = ?"
            params.append(target_library_id)
        if self.supports_paralog_runs() and paralog_run_id:
            sql += " AND pf.run_id = ?"
            params.append(paralog_run_id)
        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            sql += f" AND pf.accession IN ({placeholders})"
            params.extend(accessions)
        rows = self.core.fetchall(sql, tuple(params))
        return {str(row[0]) for row in rows}

    def _duplicate_family_copy_map(self, *, busco_library_id: int, accessions: Sequence[str]) -> dict[tuple[str, str], list[str]]:
        if not accessions:
            return {}
        out: dict[tuple[str, str], list[str]] = {}
        for chunk in self._chunked(list(accessions), 800):
            placeholders = ",".join("?" for _ in chunk)
            rows = self.core.fetchall(
                f"""
                SELECT accession, family_id, sequence
                FROM BUSCO_Family_Data
                WHERE library_id = ?
                  AND status = 2
                  AND accession IN ({placeholders})
                ORDER BY accession, family_id, sequence
                """,
                tuple([busco_library_id] + list(chunk)),
            )
            for acc, fam, seq in rows:
                if seq is None:
                    continue
                out.setdefault((str(acc), str(fam)), []).append(str(seq))
        return out

    def _duplicate_paralog_copy_lookup(
        self,
        *,
        target_library_id: int,
        busco_library_id: int,
        accessions: Sequence[str],
        run_id: Optional[str] = None,
    ) -> dict[tuple[str, str, str], tuple[bool, Optional[str]]]:
        if not accessions:
            return {}
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(target_library_id), run_id=run_id)
        out: dict[tuple[str, str, str], tuple[bool, Optional[str]]] = {}
        for chunk in self._chunked(list(accessions), 800):
            placeholders = ",".join("?" for _ in chunk)
            sql = f"""
                SELECT accession, family_id, query_id, clean, query_header
                FROM Paralog_Filtering_Copy
                WHERE target_library_id = ?
                  AND library_id = ?
            """
            params = [target_library_id, busco_library_id]
            if self.supports_paralog_runs() and paralog_run_id:
                sql += " AND run_id = ?"
                params.append(paralog_run_id)
            sql += f" AND accession IN ({placeholders})"
            params.extend(list(chunk))
            rows = self.core.fetchall(sql, tuple(params))
            for acc, fam, query_id, clean, query_header in rows:
                out[(str(acc), str(fam), str(query_id))] = (bool(clean), str(query_header) if query_header else None)
        return out

    def _duplicate_decont_copy_lookup(
        self,
        *,
        busco_library_id: int,
        accessions_by_run: dict[tuple[str, int], list[str]],
    ) -> dict[tuple[str, str, str], tuple[str, Optional[str]]]:
        out: dict[tuple[str, str, str], tuple[str, Optional[str]]] = {}
        for (run_id, source_lib), accs in (accessions_by_run or {}).items():
            for chunk in self._chunked(list(accs), 800):
                placeholders = ",".join("?" for _ in chunk)
                rows = self.core.fetchall(
                    f"""
                    SELECT accession, family_id, query_id, decision, query_header
                    FROM Decontamination_Busco_Copy_Votes
                    WHERE target_library_id = ?
                      AND busco_library_id = ?
                      AND run_id = ?
                      AND accession IN ({placeholders})
                    """,
                    tuple([source_lib, busco_library_id, run_id] + list(chunk)),
                )
                for acc, fam, query_id, decision, query_header in rows:
                    out[(str(acc), str(fam), str(query_id))] = (str(decision), str(query_header) if query_header else None)
        return out

    def get_rescued_duplicate_copies(
        self,
        *,
        target_library_id: int,
        busco_library_id: int,
        accessions: Sequence[str],
        include_paralog: bool = False,
        include_decontam: bool = False,
        paralog_run_id: Optional[str] = None,
        accessions_by_run: Optional[dict[tuple[str, int], list[str]]] = None,
        supported_decisions: Optional[Sequence[str]] = None,
    ) -> dict[tuple[str, str], dict[str, Any]]:
        if not accessions or (not include_paralog and not include_decontam):
            return {}

        copy_map = self._duplicate_family_copy_map(
            busco_library_id=busco_library_id,
            accessions=accessions,
        )
        if not copy_map:
            return {}

        supported = {str(token) for token in (supported_decisions or ("support", "weak"))}
        paralog_lookup = (
            self._duplicate_paralog_copy_lookup(
                target_library_id=target_library_id,
                busco_library_id=busco_library_id,
                accessions=accessions,
                run_id=paralog_run_id,
            )
            if include_paralog
            else {}
        )
        decont_lookup = (
            self._duplicate_decont_copy_lookup(
                busco_library_id=busco_library_id,
                accessions_by_run=accessions_by_run or {},
            )
            if include_decontam
            else {}
        )

        rescued: dict[tuple[str, str], dict[str, Any]] = {}
        for (acc, fam), query_ids in copy_map.items():
            if len(query_ids) < 2:
                continue
            clean_ids: list[tuple[str, Optional[str]]] = []
            dirty_ids: list[str] = []
            unresolved = False
            for query_id in query_ids:
                rejected = False
                has_unknown = False
                chosen_header = None
                if include_paralog:
                    pentry = paralog_lookup.get((acc, fam, query_id))
                    if pentry is None:
                        has_unknown = True
                    else:
                        pflag, pheader = pentry
                        if chosen_header is None and pheader:
                            chosen_header = pheader
                        if not pflag:
                            rejected = True
                if include_decontam:
                    dentry = decont_lookup.get((acc, fam, query_id))
                    if dentry is None:
                        has_unknown = True
                    else:
                        decision, dheader = dentry
                        if chosen_header is None and dheader:
                            chosen_header = dheader
                        if decision not in supported:
                            rejected = True
                if rejected:
                    dirty_ids.append(query_id)
                elif has_unknown:
                    unresolved = True
                else:
                    clean_ids.append((query_id, chosen_header))
            if unresolved:
                continue
            if len(clean_ids) == 1 and len(dirty_ids) == len(query_ids) - 1 and dirty_ids:
                rescued[(acc, fam)] = {
                    "query_id": clean_ids[0][0],
                    "query_header": clean_ids[0][1],
                    "dirty_query_ids": tuple(dirty_ids),
                }
        return rescued

    def get_rescued_duplicate_counts(
        self,
        *,
        target_library_id: int,
        busco_library_id: int,
        accessions: Sequence[str],
        include_paralog: bool = False,
        include_decontam: bool = False,
        paralog_run_id: Optional[str] = None,
        accessions_by_run: Optional[dict[tuple[str, int], list[str]]] = None,
        supported_decisions: Optional[Sequence[str]] = None,
    ) -> dict[str, int]:
        rescued = self.get_rescued_duplicate_copies(
            target_library_id=target_library_id,
            busco_library_id=busco_library_id,
            accessions=accessions,
            include_paralog=include_paralog,
            include_decontam=include_decontam,
            paralog_run_id=paralog_run_id,
            accessions_by_run=accessions_by_run,
            supported_decisions=supported_decisions,
        )
        counts: dict[str, int] = {}
        for acc, _fam in rescued.keys():
            counts[acc] = counts.get(acc, 0) + 1
        return counts

    def _decontam_contaminated_counts(
        self,
        *,
        family_library_id: int,
        target_library_id: int,
        busco_library_id: int,
        accessions_by_run: dict[tuple[str, int], list[str]],
        supported_decisions: Sequence[str],
    ):
        contaminated: dict[str, int] = {}
        supported = tuple(str(dec) for dec in supported_decisions if dec is not None) or ("support", "weak")
        for (run_id, target_lib), acc_list in accessions_by_run.items():
            if not acc_list:
                continue
            for chunk in self._chunked(acc_list):
                placeholders = ",".join("?" for _ in chunk)
                supported_placeholders = ",".join("?" for _ in supported)
                sql = f"""
                    SELECT v.accession,
                           COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders}) THEN v.family_id END) AS contam_cnt
                    FROM Decontamination_Busco_Votes v
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = v.family_id AND bd.library_id = ?
                    WHERE v.target_library_id = ?
                      AND v.busco_library_id = ?
                      AND v.run_id = ?
                      AND v.accession IN ({placeholders})
                    GROUP BY v.accession
                """
                params = [*supported, family_library_id, target_lib, busco_library_id, run_id, *chunk]
                rows = self.core.fetchall(sql, tuple(params))
                for acc, cnt in rows:
                    contaminated[str(acc)] = int(cnt or 0)
        return contaminated

    def _decontam_overlap_counts(
        self,
        *,
        family_library_id: int,
        target_library_id: int,
        busco_library_id: int,
        accessions_by_run: dict[tuple[str, int], list[str]],
        supported_decisions: Sequence[str],
    ):
        supports_target = self.supports_target_library()
        paralog_run_id = self.resolve_paralog_run_id(target_library_id=int(family_library_id), run_id=None)
        overlap: dict[str, int] = {}
        supported = tuple(str(dec) for dec in supported_decisions if dec is not None) or ("support", "weak")
        for (run_id, target_lib), acc_list in accessions_by_run.items():
            if not acc_list:
                continue
            for chunk in self._chunked(acc_list):
                placeholders = ",".join("?" for _ in chunk)
                supported_placeholders = ",".join("?" for _ in supported)
                sql = f"""
                    SELECT v.accession,
                           COUNT(DISTINCT v.family_id) AS overlap_cnt
                    FROM Decontamination_Busco_Votes v
                    JOIN Paralog_Filtering pf
                      ON pf.family_id = v.family_id AND pf.accession = v.accession
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = v.family_id AND bd.library_id = ?
                    WHERE v.target_library_id = ?
                      AND v.busco_library_id = ?
                      AND v.run_id = ?
                      AND COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders})
                      AND pf.clean = 0
                """
                params = [*supported, family_library_id, target_lib, busco_library_id, run_id]
                if supports_target:
                    sql += " AND pf.target_library_id = ?"
                    params.append(family_library_id)
                if self.supports_paralog_runs() and paralog_run_id:
                    sql += " AND pf.run_id = ?"
                    params.append(paralog_run_id)
                sql += " AND pf.library_id = ?"
                params.append(busco_library_id)
                sql += f" AND v.accession IN ({placeholders}) GROUP BY v.accession"
                params.extend(chunk)
                rows = self.core.fetchall(sql, tuple(params))
                for acc, cnt in rows:
                    overlap[str(acc)] = int(cnt or 0)
        return overlap

    def _decontam_decision_counts(
        self,
        *,
        family_library_id: int,
        busco_library_id: int,
        accessions_by_run: dict[tuple[str, int], list[str]],
    ):
        counts: dict[str, tuple[int, int, int]] = {}
        for (run_id, target_lib), acc_list in accessions_by_run.items():
            if not acc_list:
                continue
            for chunk in self._chunked(acc_list):
                placeholders = ",".join("?" for _ in chunk)
                sql = f"""
                    SELECT v.accession,
                           COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'support' THEN v.family_id END) AS support_cnt,
                           COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'weak' THEN v.family_id END) AS weak_cnt,
                           COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'unknown' THEN v.family_id END) AS unknown_cnt
                    FROM Decontamination_Busco_Votes v
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = v.family_id AND bd.library_id = ?
                    WHERE v.target_library_id = ?
                      AND v.busco_library_id = ?
                      AND v.run_id = ?
                      AND v.accession IN ({placeholders})
                    GROUP BY v.accession
                """
                params = [family_library_id, target_lib, busco_library_id, run_id, *chunk]
                rows = self.core.fetchall(sql, tuple(params))
                for acc, support_cnt, weak_cnt, unknown_cnt in rows:
                    counts[str(acc)] = (
                        int(support_cnt or 0),
                        int(weak_cnt or 0),
                        int(unknown_cnt or 0),
                    )
        return counts

    def get_decontamination_decision_percentages(
        self,
        *,
        library_id: int,
        accessions=None,
        run_id: Optional[str] = None,
    ):
        lib_row = self._get_library_row(int(library_id))
        if not lib_row:
            return {}
        parent_id = lib_row[6]
        busco_library_id = int(parent_id) if parent_id else int(library_id)
        size = self._get_library_size(int(library_id))
        if size <= 0:
            return {}

        accessions_list = list(dict.fromkeys([str(a) for a in (accessions or []) if a is not None])) if accessions else None
        if accessions_list is not None and not accessions_list:
            return {}

        primary = self.latest_decont_summary(
            target_library_id=int(library_id),
            accessions=accessions_list,
            run_id=run_id,
        )
        latest: dict[str, tuple[str, Optional[str], Optional[str], int]] = {
            acc: (rid, decision, date, int(library_id))
            for acc, (rid, decision, date) in primary.items()
        }
        if parent_id:
            fallback = self.latest_decont_summary(
                target_library_id=int(parent_id),
                accessions=accessions_list,
                run_id=run_id,
            )
            for acc, (rid, decision, date) in fallback.items():
                if acc not in latest:
                    latest[acc] = (rid, decision, date, int(parent_id))

        if not latest:
            return {}

        accessions_by_run: dict[tuple[str, int], list[str]] = {}
        for acc, (rid, _decision, _date, source_lib) in latest.items():
            accessions_by_run.setdefault((str(rid), int(source_lib)), []).append(acc)

        raw_counts = self._decontam_decision_counts(
            family_library_id=int(library_id),
            busco_library_id=busco_library_id,
            accessions_by_run=accessions_by_run,
        )

        results: dict[str, tuple[float, float, float]] = {}
        for acc, (support_cnt, weak_cnt, unknown_cnt) in raw_counts.items():
            results[acc] = (
                round(100.0 * support_cnt / size, 2),
                round(100.0 * weak_cnt / size, 2),
                round(100.0 * unknown_cnt / size, 2),
            )
        return results
