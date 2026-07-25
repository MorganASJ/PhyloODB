from __future__ import annotations

import json
import os
import sqlite3
from typing import Any, Optional, Sequence

from .base import BaseRepository, transactional
from ..errors import RepositoryWriteError
from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE


class BuscoRepository(BaseRepository):
    MANUAL_PRIMARY_POLICIES = {"manual_override"}

    def get_processed_accessions(self, library=None):
        if library is not None:
            rows = self.core.fetchall(
                """
                SELECT DISTINCT accession
                FROM (
                    SELECT br.accession
                    FROM BUSCO_Results br
                    WHERE br.library_id = ?
                      AND NOT EXISTS (
                          SELECT 1
                          FROM BUSCO_Runs r
                          WHERE r.accession = br.accession
                            AND r.library_id = br.library_id
                      )
                    UNION
                    SELECT accession FROM BUSCO_Runs WHERE library_id = ? AND status = 'completed'
                )
                """,
                (library, library),
            )
        else:
            rows = self.core.fetchall(
                """
                SELECT DISTINCT accession
                FROM (
                    SELECT br.accession
                    FROM BUSCO_Results br
                    WHERE br.library_id IS NULL
                      AND NOT EXISTS (
                          SELECT 1
                          FROM BUSCO_Runs r
                          WHERE r.accession = br.accession
                      )
                    UNION
                    SELECT accession FROM BUSCO_Runs WHERE status = 'completed'
                )
                """
            )
        return [row[0] for row in rows]

    def get_processed_accessions_any(self):
        rows = self.core.fetchall(
            """
            SELECT DISTINCT accession
            FROM (
                SELECT br.accession
                FROM BUSCO_Results br
                WHERE NOT EXISTS (
                    SELECT 1
                    FROM BUSCO_Runs r
                    WHERE r.accession = br.accession
                      AND (
                          r.library_id = br.library_id
                          OR (r.library_id IS NULL AND br.library_id IS NULL)
                      )
                )
                UNION
                SELECT accession FROM BUSCO_Runs WHERE status = 'completed'
            )
            """
        )
        return [row[0] for row in rows]

    def _get_library_row(self, library_id: int):
        return self.core.fetchone(
            "SELECT library_id, library_name, taxid, location, size, odb_version, parent_id, storage_root_id, relative_path FROM Libraries WHERE library_id = ?",
            (int(library_id),),
        )

    def _chunked(self, values, size: int = 900):
        if not values:
            return
        for idx in range(0, len(values), size):
            yield values[idx:idx + size]

    def _record_compat_event(
        self,
        event: str,
        *,
        count: int = 1,
        accession: Optional[str] = None,
        library_id: Optional[int] = None,
        run_id: Optional[int] = None,
        detail: Optional[str] = None,
    ) -> None:
        recorder = getattr(self.manager, "record_busco_compat_event", None)
        if callable(recorder):
            recorder(
                event,
                count=count,
                accession=accession,
                library_id=library_id,
                run_id=run_id,
                detail=detail,
            )

    def _normalize_accessions(self, accessions: Optional[Sequence[str]]) -> Optional[list[str]]:
        if accessions is None:
            return None
        normalized = [str(a) for a in accessions if a is not None]
        return list(dict.fromkeys(normalized))

    def _normalize_statuses(self, status: Optional[Sequence[int]]) -> Optional[list[int]]:
        if status is None:
            return None
        normalized: list[int] = []
        for item in status:
            try:
                normalized.append(int(item))
            except (TypeError, ValueError):
                continue
        return list(dict.fromkeys(normalized))

    def _fetch_family_rows(
        self,
        *,
        table: str,
        library_id: int,
        accessions: Optional[Sequence[str]] = None,
        family_ids: Optional[Sequence[str]] = None,
        status: Optional[Sequence[int]] = None,
        run_ids: Optional[Sequence[int]] = None,
    ) -> list[tuple]:
        clauses = ["library_id = ?"]
        params: list[Any] = [int(library_id)]
        if run_ids is not None:
            run_vals = self._normalize_run_ids(run_ids)
            if not run_vals:
                return []
            placeholders = ",".join("?" for _ in run_vals)
            clauses.append(f"run_id IN ({placeholders})")
            params.extend(run_vals)
        acc_vals = self._normalize_accessions(accessions)
        if acc_vals:
            placeholders = ",".join("?" for _ in acc_vals)
            clauses.append(f"accession IN ({placeholders})")
            params.extend(acc_vals)
        family_vals = [str(fam) for fam in family_ids or [] if fam is not None]
        if family_vals:
            placeholders = ",".join("?" for _ in family_vals)
            clauses.append(f"family_id IN ({placeholders})")
            params.extend(family_vals)
        status_vals = self._normalize_statuses(status)
        if status_vals:
            placeholders = ",".join("?" for _ in status_vals)
            clauses.append(f"status IN ({placeholders})")
            params.extend(status_vals)
        sql = (
            "SELECT family_id, library_id, accession, status, sequence, score, length "
            f"FROM {table} WHERE "
            + " AND ".join(clauses)
        )
        return self.core.fetchall(sql, tuple(params)) or []

    def _run_accessions_with_family_rows(self, run_map: dict[str, int]) -> set[str]:
        if not run_map:
            return set()
        seen: set[str] = set()
        for chunk in self._chunked(list(dict.fromkeys(int(run_id) for run_id in run_map.values()))):
            placeholders = ",".join("?" for _ in chunk)
            rows = self.core.fetchall(
                f"SELECT DISTINCT accession, run_id FROM BUSCO_Run_Family_Data WHERE run_id IN ({placeholders})",
                tuple(chunk),
            ) or []
            for accession, run_id in rows:
                if accession is None or run_id is None:
                    continue
                if run_map.get(str(accession)) == int(run_id):
                    seen.add(str(accession))
        return seen

    def _legacy_result_dir_guess(self, accession: str, library_id: int) -> Optional[str]:
        genome_path = self.manager.genomes.resolve_path(str(accession))
        library_name = self.manager.libraries.get_name(int(library_id))
        if not genome_path or not library_name:
            return None
        candidate = os.path.join(str(genome_path), f"{library_name}_results")
        if os.path.isdir(candidate):
            self._record_compat_event(
                "legacy_result_dir_filesystem_fallback",
                accession=str(accession),
                library_id=int(library_id),
                detail=candidate,
            )
            return candidate
        return None

    def _strict_result_dir_for_run_row(self, run_row: Optional[tuple]) -> Optional[str]:
        if not run_row:
            return None
        if len(run_row) >= 20 and run_row[7]:
            return str(run_row[7])
        if len(run_row) >= 7 and run_row[6]:
            return str(run_row[6])
        if run_row[0] is None:
            return None
        artifacts = self.manager.artifacts.find(
            owner_type="busco_run",
            owner_id=int(run_row[0]),
            artifact_type="busco_result_root",
        )
        if artifacts:
            resolved = self.manager.artifacts.resolve_path(artifacts[0])
            if resolved:
                return resolved
        return None

    def get_accessions_for_run_ids(self, run_ids: Sequence[int], *, library_id: Optional[int] = None) -> list[str]:
        run_vals = self._normalize_run_ids(run_ids)
        if not run_vals:
            return []
        rows: list[Any] = []
        for chunk in self._chunked(run_vals):
            placeholders = ",".join("?" for _ in chunk)
            sql = f"SELECT DISTINCT accession FROM BUSCO_Runs WHERE run_id IN ({placeholders})"
            params: list[Any] = list(chunk)
            if library_id is not None:
                sql += " AND library_id = ?"
                params.append(int(library_id))
            rows.extend(self.core.fetchall(sql, tuple(params)) or [])
        return list(dict.fromkeys(str(row[0]) for row in rows if row and row[0] is not None))

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

    def _custom_busco_base_counts(self, *, library_id: int, parent_id: int, accessions=None):
        sql = """
            SELECT bfd.accession,
                   t.name AS species,
                   COUNT(DISTINCT CASE WHEN bfd.status = 1 THEN bfd.family_id END) AS no_sc_complete,
                   COUNT(DISTINCT CASE WHEN bfd.status = 2 THEN bfd.family_id END) AS no_duplicated_complete,
                   COUNT(DISTINCT CASE WHEN bfd.status = 3 THEN bfd.family_id END) AS no_fragmented,
                   COUNT(DISTINCT CASE WHEN bfd.status = 4 THEN bfd.family_id END) AS no_missing
            FROM BUSCO_Family_Data bfd
            JOIN BUSCO_descriptions bd
              ON bd.family_id = bfd.family_id AND bd.library_id = ?
            JOIN Genome g ON g.accession = bfd.accession
            JOIN Taxonomy t ON g.taxid = t.taxid
            WHERE bfd.library_id = ?
        """
        params = [int(library_id), int(parent_id)]
        if accessions:
            acc_vals = [str(a) for a in accessions]
            placeholders = ",".join("?" for _ in acc_vals)
            sql += f" AND bfd.accession IN ({placeholders})"
            params.extend(acc_vals)
        sql += " GROUP BY bfd.accession, t.name ORDER BY bfd.accession"
        return self.core.fetchall(sql, tuple(params))

    def _legacy_busco_summary_counts(self, *, library_id: int, accessions=None):
        sql = """
            SELECT br.accession,
                   t.name AS species,
                   br.no_sc_complete,
                   br.no_duplicated_complete,
                   br.no_fragmented,
                   br.no_missing
            FROM BUSCO_Results br
            JOIN Genome g ON g.accession = br.accession
            JOIN Taxonomy t ON g.taxid = t.taxid
            WHERE br.library_id = ?
        """
        params = [int(library_id)]
        if accessions:
            acc_vals = [str(a) for a in accessions]
            placeholders = ",".join("?" for _ in acc_vals)
            sql += f" AND br.accession IN ({placeholders})"
            params.extend(acc_vals)
        sql += " ORDER BY br.accession"
        return self.core.fetchall(sql, tuple(params))

    def _run_summary_counts(self, *, library_id: int, accessions: Sequence[str], purpose: str = "default"):
        acc_vals = [str(a) for a in accessions if a is not None]
        if not acc_vals:
            return []
        run_map = self._resolve_busco_runs_for_query(
            int(library_id),
            accessions=acc_vals,
            purpose=purpose,
        )
        if not run_map:
            return []
        run_ids = list(dict.fromkeys(int(run_id) for run_id in run_map.values()))
        rows: list[tuple[Any, ...]] = []
        for chunk in self._chunked(run_ids):
            placeholders = ",".join("?" for _ in chunk)
            rows.extend(
                self.core.fetchall(
                    f"""
                    SELECT r.accession,
                           r.run_id,
                           t.name AS species,
                           r.no_sc_complete,
                           r.no_duplicated_complete,
                           r.no_fragmented,
                           r.no_missing
                    FROM BUSCO_Runs r
                    JOIN Genome g ON g.accession = r.accession
                    JOIN Taxonomy t ON g.taxid = t.taxid
                    WHERE r.library_id = ?
                      AND r.run_id IN ({placeholders})
                    """,
                    tuple([int(library_id), *chunk]),
                )
                or []
            )
        results = []
        for accession, run_id, species, sc_complete, dup_complete, fragmented, missing in rows:
            if run_map.get(str(accession)) != int(run_id):
                continue
            results.append((accession, species, sc_complete, dup_complete, fragmented, missing))
        return results

    def _custom_busco_base_counts_for_runs(
        self,
        *,
        library_id: int,
        parent_id: int,
        run_refs: Sequence[tuple[str, int]],
    ) -> dict[tuple[str, int], tuple[str | None, int, int, int, int]]:
        normalized: list[tuple[str, int]] = []
        for accession, run_id in run_refs or []:
            if accession is None or run_id is None:
                continue
            try:
                normalized.append((str(accession), int(run_id)))
            except (TypeError, ValueError):
                continue
        if not normalized:
            return {}
        run_ids = list(dict.fromkeys(int(run_id) for _acc, run_id in normalized))
        placeholders = ",".join("?" for _ in run_ids)
        sql = f"""
            SELECT d.accession,
                   d.run_id,
                   t.name AS species,
                   COUNT(DISTINCT CASE WHEN d.status = 1 THEN d.family_id END) AS no_sc_complete,
                   COUNT(DISTINCT CASE WHEN d.status = 2 THEN d.family_id END) AS no_duplicated_complete,
                   COUNT(DISTINCT CASE WHEN d.status = 3 THEN d.family_id END) AS no_fragmented,
                   COUNT(DISTINCT CASE WHEN d.status = 4 THEN d.family_id END) AS no_missing
            FROM BUSCO_Run_Family_Data d
            JOIN BUSCO_descriptions bd
              ON bd.family_id = d.family_id AND bd.library_id = ?
            JOIN Genome g ON g.accession = d.accession
            JOIN Taxonomy t ON g.taxid = t.taxid
            WHERE d.library_id = ?
              AND d.run_id IN ({placeholders})
            GROUP BY d.accession, d.run_id, t.name
            ORDER BY d.accession, d.run_id
        """
        rows = self.core.fetchall(sql, tuple([int(library_id), int(parent_id), *run_ids]))
        results: dict[tuple[str, int], tuple[str | None, int, int, int, int]] = {}
        for accession, run_id, species, sc_complete, dup_complete, fragmented, missing in rows or []:
            results[(str(accession), int(run_id))] = (
                species,
                int(sc_complete or 0),
                int(dup_complete or 0),
                int(fragmented or 0),
                int(missing or 0),
            )
        return results

    def get_results_percentages(
        self,
        accession=None,
        library_id=None,
        library_name=None,
        accessions=None,
    ):
        if accession is not None:
            accessions = [accession]
        if accessions:
            accessions = [str(a) for a in accessions if a is not None]

        if library_id is None and library_name:
            library_id = self.manager.libraries.get_id(library_name)

        if library_id is not None:
            lib_row = self._get_library_row(int(library_id))
            if not lib_row:
                return []
            lib_name = lib_row[1]
            lib_size = lib_row[4]
            parent_id = lib_row[6]
            if parent_id:
                try:
                    size = int(lib_size) if lib_size is not None else 0
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
                selected_accessions = self._normalize_accessions(accessions)
                run_map = self._resolve_busco_runs_for_query(
                    int(library_id),
                    accessions=selected_accessions,
                    purpose="default",
                )
                strict_rows: list[tuple[Any, ...]] = []
                strict_accessions: set[str] = set()
                if run_map:
                    strict_counts = self._custom_busco_base_counts_for_runs(
                        library_id=int(library_id),
                        parent_id=int(parent_id),
                        run_refs=[(acc, run_id) for acc, run_id in run_map.items()],
                    )
                    strict_rows = [
                        (acc, species, sc_complete, dup_complete, fragmented, missing)
                        for (acc, _run_id), (species, sc_complete, dup_complete, fragmented, missing) in strict_counts.items()
                    ]
                    strict_accessions = {str(acc) for acc, *_rest in strict_rows}
                    if strict_accessions:
                        self._record_compat_event(
                            "strict_results_percentages_accession",
                            count=len(strict_accessions),
                            library_id=int(library_id),
                        )
                    missing_run_family = set(run_map) - strict_accessions
                    for acc in sorted(missing_run_family):
                        self._record_compat_event(
                            "missing_run_family_data_for_percentages",
                            accession=acc,
                            library_id=int(library_id),
                            run_id=run_map.get(acc),
                        )
                fallback_scope = None
                if selected_accessions is not None:
                    fallback_scope = [acc for acc in selected_accessions if acc not in strict_accessions]
                legacy_rows = self._custom_busco_base_counts(
                    library_id=int(library_id),
                    parent_id=int(parent_id),
                    accessions=fallback_scope,
                )
                if selected_accessions is None and strict_accessions:
                    legacy_rows = [row for row in legacy_rows if str(row[0]) not in strict_accessions]
                if legacy_rows:
                    self._record_compat_event(
                        "legacy_family_data_fallback_for_percentages",
                        count=len({str(row[0]) for row in legacy_rows if row and row[0] is not None}),
                        library_id=int(library_id),
                    )
                rows = strict_rows + list(legacy_rows)
                results = []
                for acc, species, sc_complete, dup_complete, fragmented, missing in rows:
                    if size > 0:
                        sc_val = int(sc_complete or 0)
                        dup_val = int(dup_complete or 0)
                        frag_val = int(fragmented or 0)
                        miss_val = max(size - (sc_val + dup_val + frag_val), 0)
                        complete = round(100.0 * (sc_val + dup_val) / size, 2)
                        single_copy_complete = round(100.0 * sc_val / size, 2)
                        duplicated = round(100.0 * dup_val / size, 2)
                        frag_pct = round(100.0 * frag_val / size, 2)
                        miss_pct = round(100.0 * miss_val / size, 2)
                    else:
                        complete = single_copy_complete = duplicated = frag_pct = miss_pct = None
                    results.append(
                        (
                            acc,
                            species,
                            lib_name,
                            complete,
                            single_copy_complete,
                            duplicated,
                            frag_pct,
                            miss_pct,
                        )
                    )
                return results

            sql = """
                SELECT brp.accession, brp.species, brp.library_name,
                       brp.complete, brp.single_copy_complete, brp.duplicated,
                       brp.fragmented, brp.missing
                FROM BUSCO_Results_Percentages brp
                JOIN Libraries l ON l.library_name = brp.library_name
                WHERE l.library_id = ?
            """
            params = [int(library_id)]
            if accessions:
                placeholders = ",".join("?" for _ in accessions)
                sql += f" AND brp.accession IN ({placeholders})"
                params.extend(accessions)
            return self.core.fetchall(sql, tuple(params))

        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            return self.core.fetchall(
                f"SELECT * FROM BUSCO_Results_Percentages WHERE accession IN ({placeholders})",
                tuple(accessions),
            )
        return self.core.fetchall("SELECT * FROM BUSCO_Results_Percentages")

    def _normalize_adjusted_flags(
        self,
        *,
        include_paralog: Optional[bool] = None,
        include_decontam: Optional[bool] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
    ) -> tuple[bool, bool, bool, bool, bool]:
        return (
            True if include_paralog is None else bool(include_paralog),
            True if include_decontam is None else bool(include_decontam),
            bool(allow_ambiguous_contaminants) if allow_ambiguous_contaminants is not None else False,
            bool(strict_decontamination) if strict_decontamination is not None else False,
            bool(rescue_duplicates) if rescue_duplicates is not None else False,
        )

    def _resolve_adjusted_accessions(self, *, library_id: int, parent_id: Optional[int], accessions=None) -> list[str]:
        if accessions:
            return list(dict.fromkeys([str(a) for a in accessions if a is not None]))
        busco_library_id = int(parent_id) if parent_id else int(library_id)
        if parent_id:
            rows = self._custom_busco_base_counts(
                library_id=int(library_id),
                parent_id=busco_library_id,
                accessions=None,
            )
        else:
            return self.get_processed_accessions(int(library_id))
        return list(dict.fromkeys([str(row[0]) for row in rows if row and row[0] is not None]))

    def _resolve_effective_decont_context(
        self,
        *,
        target_library_id: int,
        parent_library_id: Optional[int],
        accessions: Sequence[str],
        decont_run_id: Optional[str],
    ) -> dict[str, tuple[str, int, Optional[str]]]:
        if not accessions:
            return {}
        primary = self.manager.filtering.latest_decont_summary(
            target_library_id=int(target_library_id),
            accessions=accessions,
            run_id=decont_run_id,
        )
        resolved: dict[str, tuple[str, int, Optional[str]]] = {
            str(acc): (str(run_id), int(target_library_id), decision)
            for acc, (run_id, decision, _date) in primary.items()
        }
        if parent_library_id:
            fallback = self.manager.filtering.latest_decont_summary(
                target_library_id=int(parent_library_id),
                accessions=accessions,
                run_id=decont_run_id,
            )
            for acc, (run_id, decision, _date) in fallback.items():
                key = str(acc)
                if key not in resolved:
                    resolved[key] = (str(run_id), int(parent_library_id), decision)
        return resolved

    def _make_adjusted_cache_key(
        self,
        *,
        library_id: int,
        accession: str,
        effective_busco_run_id: Optional[int],
        effective_paralog_run_id: Optional[str],
        effective_decont_run_id: Optional[str],
        effective_decont_library_id: Optional[int],
        include_paralog: bool,
        include_decontam: bool,
        allow_ambiguous_contaminants: bool,
        strict_decontamination: bool,
        rescue_duplicates: bool,
    ) -> str:
        payload = {
            "library_id": int(library_id),
            "accession": str(accession),
            "effective_busco_run_id": int(effective_busco_run_id) if effective_busco_run_id is not None else None,
            "effective_paralog_run_id": str(effective_paralog_run_id) if effective_paralog_run_id is not None else None,
            "effective_decont_run_id": str(effective_decont_run_id) if effective_decont_run_id is not None else None,
            "effective_decont_library_id": int(effective_decont_library_id) if effective_decont_library_id is not None else None,
            "include_paralog": bool(include_paralog),
            "include_decontam": bool(include_decontam),
            "allow_ambiguous_contaminants": bool(allow_ambiguous_contaminants),
            "strict_decontamination": bool(strict_decontamination),
            "rescue_duplicates": bool(rescue_duplicates),
        }
        return json.dumps(payload, sort_keys=True, separators=(",", ":"))

    def _read_adjusted_cache(self, cache_keys: Sequence[str]) -> dict[str, tuple]:
        rows_by_key: dict[str, tuple] = {}
        for chunk in self._chunked(list(cache_keys), 800):
            placeholders = ",".join("?" for _ in chunk)
            rows = self.core.fetchall(
                f"""
                SELECT cache_key, accession, species, library_name,
                       complete, single_copy_complete, duplicated, fragmented, missing,
                       hidden_paralog, contaminated,
                       effective_decont_run_id, effective_decont_decision,
                       has_paralog, has_decont, status
                FROM BUSCO_Adjusted_Results
                WHERE cache_key IN ({placeholders})
                """,
                tuple(chunk),
            )
            for row in rows or []:
                rows_by_key[str(row[0])] = row
        return rows_by_key

    def _cache_row_to_adjusted_result(self, row: tuple):
        return (
            str(row[1]),
            row[2],
            row[3],
            row[4],
            row[5],
            row[6],
            row[7],
            row[8],
            row[9],
            row[10],
            row[11],
            row[12],
            bool(row[13]),
            bool(row[14]),
        )

    @transactional("update adjusted BUSCO cache")
    def _upsert_adjusted_cache_rows(
        self,
        *,
        library_id: int,
        cache_keys: dict[str, str],
        effective_busco_run_ids: dict[str, int],
        effective_decont_context: dict[str, tuple[str, int, Optional[str]]],
        include_paralog: bool,
        include_decontam: bool,
        allow_ambiguous_contaminants: bool,
        strict_decontamination: bool,
        rescue_duplicates: bool,
        rows,
    ) -> None:
        payload = []
        for row in rows or []:
            acc = str(row[0])
            decont_info = effective_decont_context.get(acc)
            payload.append(
                (
                    cache_keys[acc],
                    int(library_id),
                    acc,
                    row[1],
                    row[2],
                    int(effective_busco_run_ids[acc]) if acc in effective_busco_run_ids else None,
                    str(decont_info[0]) if decont_info is not None else None,
                    int(decont_info[1]) if decont_info is not None else None,
                    str(decont_info[2]) if decont_info is not None and decont_info[2] is not None else None,
                    1 if include_paralog else 0,
                    1 if include_decontam else 0,
                    1 if allow_ambiguous_contaminants else 0,
                    1 if strict_decontamination else 0,
                    1 if rescue_duplicates else 0,
                    row[3],
                    row[4],
                    row[5],
                    row[6],
                    row[7],
                    row[8],
                    row[9],
                    1 if row[12] else 0,
                    1 if row[13] else 0,
                )
            )
        if not payload:
            return
        self.core.executemany(
            """
            INSERT INTO BUSCO_Adjusted_Results (
                cache_key, library_id, accession, species, library_name,
                effective_busco_run_id, effective_decont_run_id, effective_decont_library_id,
                effective_decont_decision,
                include_paralog, include_decontam, allow_ambiguous_contaminants,
                strict_decontamination, rescue_duplicates,
                complete, single_copy_complete, duplicated, fragmented, missing,
                hidden_paralog, contaminated, has_paralog, has_decont,
                status, updated_at, invalidated_at, invalidation_reason
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ready', datetime('now'), NULL, NULL)
            ON CONFLICT(cache_key) DO UPDATE SET
                species = excluded.species,
                library_name = excluded.library_name,
                effective_busco_run_id = excluded.effective_busco_run_id,
                effective_decont_run_id = excluded.effective_decont_run_id,
                effective_decont_library_id = excluded.effective_decont_library_id,
                effective_decont_decision = excluded.effective_decont_decision,
                include_paralog = excluded.include_paralog,
                include_decontam = excluded.include_decontam,
                allow_ambiguous_contaminants = excluded.allow_ambiguous_contaminants,
                strict_decontamination = excluded.strict_decontamination,
                rescue_duplicates = excluded.rescue_duplicates,
                complete = excluded.complete,
                single_copy_complete = excluded.single_copy_complete,
                duplicated = excluded.duplicated,
                fragmented = excluded.fragmented,
                missing = excluded.missing,
                hidden_paralog = excluded.hidden_paralog,
                contaminated = excluded.contaminated,
                has_paralog = excluded.has_paralog,
                has_decont = excluded.has_decont,
                status = 'ready',
                updated_at = datetime('now'),
                invalidated_at = NULL,
                invalidation_reason = NULL
            """,
            payload,
        )

    def get_results_adjusted_cached(
        self,
        *,
        library_id: int,
        accessions=None,
        include_paralog: Optional[bool] = None,
        paralog_run_id: Optional[str] = None,
        include_decontam: Optional[bool] = None,
        decont_run_id: Optional[str] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
    ):
        lib_row = self._get_library_row(int(library_id))
        if not lib_row:
            return []
        parent_id = lib_row[6]
        accessions_list = self._resolve_adjusted_accessions(
            library_id=int(library_id),
            parent_id=int(parent_id) if parent_id else None,
            accessions=accessions,
        )
        if not accessions_list:
            return []

        (
            include_paralog_val,
            include_decontam_val,
            allow_ambiguous_val,
            strict_decontam_val,
            rescue_dup_val,
        ) = self._normalize_adjusted_flags(
            include_paralog=include_paralog,
            include_decontam=include_decontam,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )

        busco_scope_id = int(parent_id) if parent_id else int(library_id)
        effective_busco_run_ids = self._resolve_busco_runs_for_query(
            busco_scope_id,
            accessions=accessions_list,
            purpose="default",
        )
        effective_paralog_run_id = self.manager.filtering.resolve_paralog_run_id(
            target_library_id=int(library_id),
            run_id=paralog_run_id,
        )
        effective_decont_context = self._resolve_effective_decont_context(
            target_library_id=int(library_id),
            parent_library_id=int(parent_id) if parent_id else None,
            accessions=accessions_list,
            decont_run_id=decont_run_id,
        )

        cache_keys = {
            acc: self._make_adjusted_cache_key(
                library_id=int(library_id),
                accession=acc,
                effective_busco_run_id=effective_busco_run_ids.get(acc),
                effective_paralog_run_id=effective_paralog_run_id,
                effective_decont_run_id=effective_decont_context.get(acc, (None, None, None))[0],
                effective_decont_library_id=effective_decont_context.get(acc, (None, None, None))[1],
                include_paralog=include_paralog_val,
                include_decontam=include_decontam_val,
                allow_ambiguous_contaminants=allow_ambiguous_val,
                strict_decontamination=strict_decontam_val,
                rescue_duplicates=rescue_dup_val,
            )
            for acc in accessions_list
        }
        cached_rows = self._read_adjusted_cache(list(cache_keys.values()))
        results_by_acc = {
            str(row[1]): self._cache_row_to_adjusted_result(row)
            for row in cached_rows.values()
            if str(row[15] or "ready") == "ready"
        }
        missing = [acc for acc in accessions_list if acc not in results_by_acc]
        if missing:
            live_rows = self._compute_results_adjusted_live(
                library_id=int(library_id),
                accessions=missing,
                include_paralog=include_paralog_val,
                paralog_run_id=effective_paralog_run_id,
                include_decontam=include_decontam_val,
                decont_run_id=decont_run_id,
                allow_ambiguous_contaminants=allow_ambiguous_val,
                strict_decontamination=strict_decontam_val,
                rescue_duplicates=rescue_dup_val,
            )
            if not self.manager.read_only:
                self._upsert_adjusted_cache_rows(
                    library_id=int(library_id),
                    cache_keys=cache_keys,
                    effective_busco_run_ids=effective_busco_run_ids,
                    effective_decont_context=effective_decont_context,
                    include_paralog=include_paralog_val,
                    include_decontam=include_decontam_val,
                    allow_ambiguous_contaminants=allow_ambiguous_val,
                    strict_decontamination=strict_decontam_val,
                    rescue_duplicates=rescue_dup_val,
                    rows=live_rows,
                )
                self.conn.commit()
            for row in live_rows or []:
                results_by_acc[str(row[0])] = row
        return [results_by_acc[acc] for acc in accessions_list if acc in results_by_acc]

    def get_results_adjusted(
        self,
        *,
        library_id: int,
        accessions=None,
        include_paralog: Optional[bool] = None,
        paralog_run_id: Optional[str] = None,
        include_decontam: Optional[bool] = None,
        decont_run_id: Optional[str] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
    ):
        return self.get_results_adjusted_cached(
            library_id=library_id,
            accessions=accessions,
            include_paralog=include_paralog,
            paralog_run_id=paralog_run_id,
            include_decontam=include_decontam,
            decont_run_id=decont_run_id,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )

    def get_display_results_for_runs(
        self,
        *,
        library_id: int,
        run_refs: Sequence[tuple[str, int]],
        include_paralog: Optional[bool] = None,
        paralog_run_id: Optional[str] = None,
        include_decontam: Optional[bool] = None,
        decont_run_id: Optional[str] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
    ) -> dict[tuple[str, int], tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], bool, bool, Optional[float], Optional[float], Optional[float], Optional[str], Optional[str], Optional[str]]]:
        normalized_refs: list[tuple[str, int]] = []
        for accession, run_id in run_refs or []:
            if accession is None or run_id is None:
                continue
            try:
                normalized_refs.append((str(accession), int(run_id)))
            except (TypeError, ValueError):
                continue
        if not normalized_refs:
            return {}

        lib_row = self._get_library_row(int(library_id))
        if not lib_row:
            return {}
        parent_id = int(lib_row[6]) if lib_row[6] else None
        busco_library_id = parent_id if parent_id else int(library_id)
        displayed_runs_by_accession: dict[str, list[int]] = {}
        for accession, run_id in normalized_refs:
            displayed_runs_by_accession.setdefault(str(accession), []).append(int(run_id))
        size = self._get_library_size(int(library_id))
        (
            include_paralog_val,
            include_decontam_val,
            allow_ambiguous_val,
            strict_decontam_val,
            rescue_dup_val,
        ) = self._normalize_adjusted_flags(
            include_paralog=include_paralog,
            include_decontam=include_decontam,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )

        base_by_ref: dict[tuple[str, int], tuple[Any, int, int, int, int]] = {}
        if parent_id:
            base_by_ref = self._custom_busco_base_counts_for_runs(
                library_id=int(library_id),
                parent_id=int(parent_id),
                run_refs=normalized_refs,
            )
        else:
            run_ids = [run_id for _acc, run_id in normalized_refs]
            run_rows: dict[int, tuple] = {}
            for chunk in self._chunked(list(dict.fromkeys(run_ids))):
                placeholders = ",".join("?" for _ in chunk)
                rows = self.core.fetchall(
                    f"""
                    SELECT run_id, accession, no_sc_complete, no_duplicated_complete, no_fragmented, no_missing
                    FROM BUSCO_Runs
                    WHERE library_id = ? AND run_id IN ({placeholders})
                    """,
                    tuple([int(library_id), *chunk]),
                )
                for row in rows or []:
                    run_rows[int(row[0])] = row
            for accession, run_id in normalized_refs:
                row = run_rows.get(int(run_id))
                if not row or str(row[1]) != accession:
                    continue
                sc_val = int(row[2] or 0)
                dup_val = int(row[3] or 0)
                frag_val = int(row[4] or 0)
                miss_val = int(row[5] or 0)
                if miss_val == 0 and size > 0:
                    miss_val = max(size - (sc_val + dup_val + frag_val), 0)
                base_by_ref[(accession, int(run_id))] = (None, sc_val, dup_val, frag_val, miss_val)

        if not base_by_ref:
            return {}

        supported_decisions = ["support"]
        if not strict_decontam_val:
            supported_decisions.append("weak")
        if allow_ambiguous_val:
            supported_decisions.append("unknown")
        supported = tuple(str(dec) for dec in supported_decisions if dec is not None) or ("support", "weak")
        supported_placeholders = ",".join("?" for _ in supported)

        def _pct(count: Optional[int]) -> Optional[float]:
            if count is None:
                return None
            if size <= 0:
                return None
            return round((100.0 * int(count)) / float(size), 2)

        def _latest_paralog_run_for_ref(accession: str, run_id: int) -> Optional[str]:
            has_busco_link = self.manager._column_exists("Paralog_Filtering", "busco_run_id")

            def _has_legacy_unlinked_rows() -> bool:
                if not has_busco_link:
                    return True
                params: list[Any] = [int(library_id), str(accession), int(busco_library_id)]
                sql = """
                    SELECT 1
                    FROM Paralog_Filtering pf
                    LEFT JOIN BUSCO_Runs br ON br.run_id = pf.busco_run_id
                    WHERE pf.target_library_id = ?
                      AND pf.accession = ?
                      AND pf.library_id = ?
                      AND (pf.busco_run_id IS NULL OR br.run_id IS NULL)
                """
                if paralog_run_id:
                    sql += " AND pf.run_id = ?"
                    params.append(str(paralog_run_id))
                sql += " LIMIT 1"
                row = self.core.fetchone(sql, tuple(params))
                return bool(row)

            params: list[Any] = [int(library_id), str(accession), int(busco_library_id)]
            sql = """
                SELECT pf.run_id
                FROM Paralog_Filtering pf
                WHERE pf.target_library_id = ?
                  AND pf.accession = ?
                  AND pf.library_id = ?
            """
            if has_busco_link:
                sql += " AND pf.busco_run_id = ?"
                params.append(int(run_id))
            else:
                sql += """
                  AND EXISTS (
                        SELECT 1
                        FROM BUSCO_Run_Family_Data d
                        JOIN BUSCO_descriptions bd
                          ON bd.family_id = d.family_id AND bd.library_id = ?
                        WHERE d.run_id = ?
                          AND d.accession = pf.accession
                          AND d.family_id = pf.family_id
                    )
                """
                params.extend([int(library_id), int(run_id)])
            if paralog_run_id:
                sql += " AND pf.run_id = ?"
                params.append(str(paralog_run_id))
            sql += " ORDER BY COALESCE(pf.date, '') DESC, pf.rowid DESC LIMIT 1"
            row = self.core.fetchone(sql, tuple(params))
            if (not row or not row[0]) and has_busco_link and _has_legacy_unlinked_rows():
                params = [int(library_id), str(accession), int(busco_library_id)]
                sql = """
                    SELECT pf.run_id
                    FROM Paralog_Filtering pf
                    JOIN BUSCO_Run_Family_Data d
                      ON d.family_id = pf.family_id
                     AND d.accession = pf.accession
                     AND d.run_id = ?
                     AND d.library_id = ?
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = pf.family_id AND bd.library_id = ?
                    WHERE pf.target_library_id = ?
                      AND pf.accession = ?
                      AND pf.library_id = ?
                """
                params = [int(run_id), int(busco_library_id), int(library_id), int(library_id), str(accession), int(busco_library_id)]
                if paralog_run_id:
                    sql += " AND pf.run_id = ?"
                    params.append(str(paralog_run_id))
                sql += " ORDER BY COALESCE(pf.date, '') DESC, pf.rowid DESC LIMIT 1"
                row = self.core.fetchone(sql, tuple(params))
            return str(row[0]) if row and row[0] else None

        def _paralog_counts(accession: str, run_id: int) -> tuple[Optional[int], bool, Optional[str]]:
            selected_run = _latest_paralog_run_for_ref(accession, run_id)
            if selected_run is None:
                return None, False, None
            has_busco_link = self.manager._column_exists("Paralog_Filtering", "busco_run_id")
            used_legacy_fallback = False
            if has_busco_link:
                row = self.core.fetchone(
                    """
                    SELECT
                        COUNT(DISTINCT CASE WHEN pf.clean = 0 THEN pf.family_id END) AS hidden_cnt,
                        COUNT(DISTINCT pf.family_id) AS total_cnt
                    FROM Paralog_Filtering pf
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = pf.family_id AND bd.library_id = ?
                    WHERE pf.target_library_id = ?
                      AND pf.library_id = ?
                      AND pf.accession = ?
                      AND pf.run_id = ?
                      AND pf.busco_run_id = ?
                    """,
                    (int(library_id), int(library_id), int(busco_library_id), str(accession), str(selected_run), int(run_id)),
                )
                if (not row or int(row[1] or 0) <= 0):
                    legacy_row = self.core.fetchone(
                        """
                        SELECT 1
                        FROM Paralog_Filtering pf
                        LEFT JOIN BUSCO_Runs br ON br.run_id = pf.busco_run_id
                        WHERE pf.target_library_id = ?
                          AND pf.accession = ?
                          AND pf.library_id = ?
                          AND pf.run_id = ?
                          AND (pf.busco_run_id IS NULL OR br.run_id IS NULL)
                        LIMIT 1
                        """,
                        (int(library_id), str(accession), int(busco_library_id), str(selected_run)),
                    )
                    if not legacy_row:
                        total_cnt = int(row[1] or 0) if row else 0
                        if total_cnt <= 0:
                            return None, False, selected_run
                    used_legacy_fallback = True
                    row = self.core.fetchone(
                        """
                        SELECT
                            COUNT(DISTINCT CASE WHEN pf.clean = 0 THEN pf.family_id END) AS hidden_cnt,
                            COUNT(DISTINCT pf.family_id) AS total_cnt
                        FROM Paralog_Filtering pf
                        JOIN BUSCO_Run_Family_Data d
                          ON d.family_id = pf.family_id
                         AND d.accession = pf.accession
                         AND d.run_id = ?
                         AND d.library_id = ?
                        JOIN BUSCO_descriptions bd
                          ON bd.family_id = pf.family_id AND bd.library_id = ?
                        WHERE pf.target_library_id = ?
                          AND pf.library_id = ?
                          AND pf.accession = ?
                          AND pf.run_id = ?
                        """,
                        (int(run_id), int(busco_library_id), int(library_id), int(library_id), int(busco_library_id), str(accession), str(selected_run)),
                    )
            else:
                row = self.core.fetchone(
                    """
                    SELECT
                        COUNT(DISTINCT CASE WHEN pf.clean = 0 THEN pf.family_id END) AS hidden_cnt,
                        COUNT(DISTINCT pf.family_id) AS total_cnt
                    FROM Paralog_Filtering pf
                    JOIN BUSCO_Run_Family_Data d
                      ON d.family_id = pf.family_id
                     AND d.accession = pf.accession
                     AND d.run_id = ?
                     AND d.library_id = ?
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = pf.family_id AND bd.library_id = ?
                    WHERE pf.target_library_id = ?
                      AND pf.library_id = ?
                      AND pf.accession = ?
                      AND pf.run_id = ?
                    """,
                    (int(run_id), int(busco_library_id), int(library_id), int(library_id), int(busco_library_id), str(accession), str(selected_run)),
                )
            total_cnt = int(row[1] or 0) if row else 0
            if total_cnt <= 0:
                return None, False, selected_run
            return int(row[0] or 0), True, selected_run

        def _latest_decont_context_for_ref(accession: str, run_id: int) -> tuple[Optional[str], Optional[int], Optional[str]]:
            target_candidates = [int(library_id)]
            if parent_id:
                target_candidates.append(int(parent_id))
            summary_has_link = self.manager._column_exists("Decontamination_Summary", "busco_run_id")
            votes_has_link = self.manager._column_exists("Decontamination_Busco_Votes", "busco_run_id")

            for target_library_id in target_candidates:
                def _has_legacy_unlinked_decont_rows() -> bool:
                    if not summary_has_link and not votes_has_link:
                        return True
                    if summary_has_link:
                        params: list[Any] = [str(accession), int(target_library_id), int(busco_library_id)]
                        sql = """
                            SELECT 1
                            FROM Decontamination_Summary s
                            LEFT JOIN BUSCO_Runs br ON br.run_id = s.busco_run_id
                            WHERE s.accession = ?
                              AND s.target_library_id = ?
                              AND s.busco_library_id = ?
                              AND (s.busco_run_id IS NULL OR br.run_id IS NULL)
                        """
                        if decont_run_id:
                            sql += " AND s.run_id = ?"
                            params.append(str(decont_run_id))
                        sql += " LIMIT 1"
                        if self.core.fetchone(sql, tuple(params)):
                            return True
                    if votes_has_link:
                        params = [str(accession), int(target_library_id), int(busco_library_id)]
                        sql = """
                            SELECT 1
                            FROM Decontamination_Busco_Votes v
                            LEFT JOIN BUSCO_Runs br ON br.run_id = v.busco_run_id
                            WHERE v.accession = ?
                              AND v.target_library_id = ?
                              AND v.busco_library_id = ?
                              AND (v.busco_run_id IS NULL OR br.run_id IS NULL)
                        """
                        if decont_run_id:
                            sql += " AND v.run_id = ?"
                            params.append(str(decont_run_id))
                        sql += " LIMIT 1"
                        if self.core.fetchone(sql, tuple(params)):
                            return True
                    return False

                params: list[Any] = [str(accession), int(target_library_id), int(busco_library_id)]
                sql = """
                    SELECT run_id, decision
                    FROM Decontamination_Summary
                    WHERE accession = ?
                      AND target_library_id = ?
                      AND busco_library_id = ?
                """
                if summary_has_link:
                    sql += " AND busco_run_id = ?"
                    params.append(int(run_id))
                if decont_run_id:
                    sql += " AND run_id = ?"
                    params.append(str(decont_run_id))
                sql += " ORDER BY COALESCE(date, '') DESC, rowid DESC LIMIT 1"
                row = self.core.fetchone(sql, tuple(params))
                if row and row[0]:
                    return str(row[0]), int(target_library_id), str(row[1]) if row[1] is not None else None
                if not _has_legacy_unlinked_decont_rows():
                    continue
                # Fallback to vote linkage if summary rows are missing.
                params = [str(accession), int(target_library_id), int(busco_library_id)]
                sql = """
                    SELECT run_id
                    FROM Decontamination_Busco_Votes
                    WHERE accession = ?
                      AND target_library_id = ?
                      AND busco_library_id = ?
                """
                if self.manager._column_exists("Decontamination_Busco_Votes", "busco_run_id"):
                    sql += " AND busco_run_id = ?"
                    params.append(int(run_id))
                if decont_run_id:
                    sql += " AND run_id = ?"
                    params.append(str(decont_run_id))
                sql += " ORDER BY COALESCE(date, '') DESC, rowid DESC LIMIT 1"
                row = self.core.fetchone(sql, tuple(params))
                if row and row[0]:
                    return str(row[0]), int(target_library_id), None
                if votes_has_link:
                    params = [int(run_id), int(busco_library_id), str(accession), int(target_library_id), int(busco_library_id)]
                    sql = """
                        SELECT v.run_id
                        FROM Decontamination_Busco_Votes v
                        JOIN BUSCO_Run_Family_Data d
                          ON d.family_id = v.family_id
                         AND d.accession = v.accession
                         AND d.run_id = ?
                         AND d.library_id = ?
                        WHERE v.accession = ?
                          AND v.target_library_id = ?
                          AND v.busco_library_id = ?
                    """
                    if decont_run_id:
                        sql += " AND v.run_id = ?"
                        params.append(str(decont_run_id))
                    sql += " ORDER BY COALESCE(v.date, '') DESC, v.rowid DESC LIMIT 1"
                    row = self.core.fetchone(sql, tuple(params))
                    if row and row[0]:
                        return str(row[0]), int(target_library_id), None
                if len(displayed_runs_by_accession.get(str(accession), [])) == 1 and _has_legacy_unlinked_decont_rows():
                    params = [str(accession), int(target_library_id), int(busco_library_id)]
                    sql = """
                        SELECT run_id, decision
                        FROM Decontamination_Summary
                        WHERE accession = ?
                          AND target_library_id = ?
                          AND busco_library_id = ?
                          AND busco_run_id IS NULL
                    """
                    if decont_run_id:
                        sql += " AND run_id = ?"
                        params.append(str(decont_run_id))
                    sql += " ORDER BY COALESCE(date, '') DESC, rowid DESC LIMIT 1"
                    row = self.core.fetchone(sql, tuple(params))
                    if row and row[0]:
                        return str(row[0]), int(target_library_id), str(row[1]) if row[1] is not None else None
                    params = [str(accession), int(target_library_id), int(busco_library_id)]
                    sql = """
                        SELECT run_id
                        FROM Decontamination_Busco_Votes
                        WHERE accession = ?
                          AND target_library_id = ?
                          AND busco_library_id = ?
                          AND busco_run_id IS NULL
                    """
                    if decont_run_id:
                        sql += " AND run_id = ?"
                        params.append(str(decont_run_id))
                    sql += " ORDER BY COALESCE(date, '') DESC, rowid DESC LIMIT 1"
                    row = self.core.fetchone(sql, tuple(params))
                    if row and row[0]:
                        return str(row[0]), int(target_library_id), None
            return None, None, None

        def _decont_counts(accession: str, run_id: int, selected_run: str, source_library_id: int) -> tuple[Optional[int], Optional[int], Optional[int], Optional[int], bool]:
            if self.manager._column_exists("Decontamination_Busco_Votes", "busco_run_id"):
                row = self.core.fetchone(
                    f"""
                    SELECT
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'support' THEN v.family_id END),
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'weak' THEN v.family_id END),
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'unknown' THEN v.family_id END),
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders}) THEN v.family_id END),
                        COUNT(DISTINCT v.family_id)
                    FROM Decontamination_Busco_Votes v
                    {"JOIN BUSCO_descriptions bd ON bd.family_id = v.family_id AND bd.library_id = ?" if parent_id else ""}
                    WHERE v.accession = ?
                      AND v.target_library_id = ?
                      AND v.busco_library_id = ?
                      AND v.run_id = ?
                      AND v.busco_run_id = ?
                    """,
                    tuple(
                        [
                            *supported,
                            *((int(library_id),) if parent_id else ()),
                            str(accession),
                            int(source_library_id),
                            int(busco_library_id),
                            str(selected_run),
                            int(run_id),
                        ]
                    ),
                )
                total_cnt = int(row[4] or 0) if row else 0
                if total_cnt <= 0:
                    row = self.core.fetchone(
                        f"""
                        SELECT
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'support' THEN v.family_id END),
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'weak' THEN v.family_id END),
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'unknown' THEN v.family_id END),
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders}) THEN v.family_id END),
                            COUNT(DISTINCT v.family_id)
                        FROM Decontamination_Busco_Votes v
                        JOIN BUSCO_Run_Family_Data d
                          ON d.family_id = v.family_id
                         AND d.accession = v.accession
                         AND d.run_id = ?
                         AND d.library_id = ?
                        {"JOIN BUSCO_descriptions bd ON bd.family_id = v.family_id AND bd.library_id = ?" if parent_id else ""}
                        WHERE v.accession = ?
                          AND v.target_library_id = ?
                          AND v.busco_library_id = ?
                          AND v.run_id = ?
                        """,
                        tuple(
                            [
                                *supported,
                                int(run_id),
                                int(busco_library_id),
                                *((int(library_id),) if parent_id else ()),
                                str(accession),
                                int(source_library_id),
                                int(busco_library_id),
                                str(selected_run),
                            ]
                        ),
                    )
                    total_cnt = int(row[4] or 0) if row else 0
                if total_cnt <= 0 and len(displayed_runs_by_accession.get(str(accession), [])) == 1:
                    row = self.core.fetchone(
                        f"""
                        SELECT
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'support' THEN v.family_id END),
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'weak' THEN v.family_id END),
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'unknown' THEN v.family_id END),
                            COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders}) THEN v.family_id END),
                            COUNT(DISTINCT v.family_id)
                        FROM Decontamination_Busco_Votes v
                        {"JOIN BUSCO_descriptions bd ON bd.family_id = v.family_id AND bd.library_id = ?" if parent_id else ""}
                        WHERE v.accession = ?
                          AND v.target_library_id = ?
                          AND v.busco_library_id = ?
                          AND v.run_id = ?
                          AND v.busco_run_id IS NULL
                        """,
                        tuple(
                            [
                                *supported,
                                *((int(library_id),) if parent_id else ()),
                                str(accession),
                                int(source_library_id),
                                int(busco_library_id),
                                str(selected_run),
                            ]
                        ),
                    )
            else:
                row = self.core.fetchone(
                    f"""
                    SELECT
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'support' THEN v.family_id END),
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'weak' THEN v.family_id END),
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') = 'unknown' THEN v.family_id END),
                        COUNT(DISTINCT CASE WHEN COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders}) THEN v.family_id END),
                        COUNT(DISTINCT v.family_id)
                    FROM Decontamination_Busco_Votes v
                    JOIN BUSCO_Run_Family_Data d
                      ON d.family_id = v.family_id
                     AND d.accession = v.accession
                     AND d.run_id = ?
                     AND d.library_id = ?
                    {"JOIN BUSCO_descriptions bd ON bd.family_id = v.family_id AND bd.library_id = ?" if parent_id else ""}
                    WHERE v.accession = ?
                      AND v.target_library_id = ?
                      AND v.busco_library_id = ?
                      AND v.run_id = ?
                    """,
                    tuple(
                        [
                            *supported,
                            int(run_id),
                            int(busco_library_id),
                            *((int(library_id),) if parent_id else ()),
                            str(accession),
                            int(source_library_id),
                            int(busco_library_id),
                            str(selected_run),
                        ]
                    ),
                )
            total_cnt = int(row[4] or 0) if row else 0
            if total_cnt <= 0:
                return None, None, None, None, False
            return int(row[0] or 0), int(row[1] or 0), int(row[2] or 0), int(row[3] or 0), True

        results: dict[tuple[str, int], tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], bool, bool, Optional[float], Optional[float], Optional[float], Optional[str], Optional[str], Optional[str]]] = {}
        for accession, run_id in normalized_refs:
            base_entry = base_by_ref.get((accession, int(run_id)))
            if not base_entry:
                continue

            _species, sc_val, dup_val, frag_val, miss_val = base_entry
            if miss_val == 0 and size > 0:
                miss_val = max(size - (sc_val + dup_val + frag_val), 0)

            hidden_val: Optional[int] = None
            has_paralog = False
            selected_paralog_run: Optional[str] = None
            if include_paralog_val or paralog_run_id:
                hidden_val, has_paralog, selected_paralog_run = _paralog_counts(accession, int(run_id))

            support_val: Optional[int] = None
            weak_val: Optional[int] = None
            unknown_val: Optional[int] = None
            contam_val: Optional[int] = None
            has_decont = False
            decont_decision: Optional[str] = None
            selected_decont_run: Optional[str] = None
            if include_decontam_val or decont_run_id:
                selected_decont_run, source_library_id, decont_decision = _latest_decont_context_for_ref(accession, int(run_id))
                if selected_decont_run and source_library_id is not None:
                    support_val, weak_val, unknown_val, contam_val, has_decont = _decont_counts(
                        accession,
                        int(run_id),
                        str(selected_decont_run),
                        int(source_library_id),
                    )

            overlap_val = 0
            if include_paralog_val and include_decontam_val and hidden_val is not None and contam_val is not None and selected_paralog_run and selected_decont_run:
                if self.manager._column_exists("Paralog_Filtering", "busco_run_id") and self.manager._column_exists("Decontamination_Busco_Votes", "busco_run_id"):
                    overlap_row = self.core.fetchone(
                        f"""
                        SELECT COUNT(DISTINCT v.family_id)
                        FROM Decontamination_Busco_Votes v
                        JOIN Paralog_Filtering pf
                          ON pf.family_id = v.family_id
                         AND pf.accession = v.accession
                         AND pf.target_library_id = ?
                         AND pf.run_id = ?
                         AND pf.busco_run_id = v.busco_run_id
                        {"JOIN BUSCO_descriptions bd ON bd.family_id = v.family_id AND bd.library_id = ?" if parent_id else ""}
                        WHERE v.accession = ?
                          AND v.target_library_id = ?
                          AND v.busco_library_id = ?
                          AND v.run_id = ?
                          AND v.busco_run_id = ?
                          AND COALESCE(v.decision, 'unknown') NOT IN ({supported_placeholders})
                          AND pf.clean = 0
                        """,
                        tuple(
                            [
                                int(library_id),
                                str(selected_paralog_run),
                                *((int(library_id),) if parent_id else ()),
                                str(accession),
                                int(parent_id if parent_id else library_id),
                                int(busco_library_id),
                                str(selected_decont_run),
                                int(run_id),
                                *supported,
                            ]
                        ),
                    )
                else:
                    overlap_row = None
                overlap_val = int(overlap_row[0] or 0) if overlap_row and overlap_row[0] is not None else 0

            rescued_val = 0
            if rescue_dup_val:
                accessions_by_run = {}
                if selected_decont_run:
                    accessions_by_run[(str(selected_decont_run), int(parent_id if parent_id else library_id))] = [str(accession)]
                rescue_map = self.manager.filtering.get_rescued_duplicate_counts(
                    target_library_id=int(library_id),
                    busco_library_id=int(busco_library_id),
                    accessions=[str(accession)],
                    include_paralog=include_paralog_val,
                    include_decontam=include_decontam_val,
                    paralog_run_id=selected_paralog_run,
                    accessions_by_run=accessions_by_run,
                    supported_decisions=supported_decisions,
                )
                rescued_val = int(rescue_map.get(accession, 0) or 0)

            contam_effective = contam_val
            if include_decontam_val and contam_val is not None and include_paralog_val and hidden_val is not None:
                contam_effective = max(contam_val - overlap_val, 0)

            sc_adj = int(sc_val or 0)
            dup_adj = int(dup_val or 0)
            if include_paralog_val and hidden_val is not None:
                sc_adj -= hidden_val
            if include_decontam_val and contam_effective is not None:
                sc_adj -= contam_effective
            if rescued_val:
                sc_adj += rescued_val
                dup_adj = max(dup_adj - rescued_val, 0)
            if sc_adj < 0:
                sc_adj = 0

            hidden_display = hidden_val if include_paralog_val else (0 if hidden_val is not None else None)
            contam_display = contam_effective if include_decontam_val else (0 if contam_effective is not None else None)

            complete_pct = _pct(sc_adj + dup_adj)
            single_copy_pct = _pct(sc_adj)
            duplicated_pct = _pct(dup_adj)
            frag_pct = _pct(frag_val)
            miss_final = max(size - (sc_adj + dup_adj + int(frag_val or 0) + (hidden_display or 0) + (contam_display or 0)), 0) if size > 0 else None
            miss_pct = _pct(miss_final) if miss_final is not None else None
            hidden_pct = _pct(hidden_display) if hidden_display is not None else None
            contam_pct = _pct(contam_display) if contam_display is not None else None
            support_pct = _pct(support_val) if support_val is not None else None
            weak_pct = _pct(weak_val) if weak_val is not None else None
            unknown_pct = _pct(unknown_val) if unknown_val is not None else None

            results[(accession, int(run_id))] = (
                complete_pct,
                single_copy_pct,
                duplicated_pct,
                frag_pct,
                miss_pct,
                hidden_pct,
                contam_pct,
                bool(has_paralog),
                bool(has_decont),
                support_pct,
                weak_pct,
                unknown_pct,
                selected_paralog_run,
                selected_decont_run,
                decont_decision,
            )
        return results

    def _compute_results_adjusted_live(
        self,
        *,
        library_id: int,
        accessions=None,
        include_paralog: Optional[bool] = None,
        paralog_run_id: Optional[str] = None,
        include_decontam: Optional[bool] = None,
        decont_run_id: Optional[str] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
    ):
        lib_row = self._get_library_row(int(library_id))
        if not lib_row:
            return []
        lib_name = lib_row[1]
        parent_id = lib_row[6]
        busco_library_id = int(parent_id) if parent_id else int(library_id)

        size = self._get_library_size(int(library_id))
        if parent_id:
            rows = self._custom_busco_base_counts(
                library_id=int(library_id),
                parent_id=busco_library_id,
                accessions=accessions,
            )
        else:
            rows = self._run_summary_counts(
                library_id=int(library_id),
                accessions=accessions,
            )
        base: dict[str, tuple[str | None, int, int, int, int]] = {}
        for acc, species, sc_complete, dup_complete, fragmented, missing in rows:
            sc_val = int(sc_complete or 0)
            dup_val = int(dup_complete or 0)
            frag_val = int(fragmented or 0)
            miss_val = int(missing or 0)
            if miss_val == 0 and size > 0:
                miss_val = max(size - (sc_val + dup_val + frag_val), 0)
            base[str(acc)] = (species, sc_val, dup_val, frag_val, miss_val)

        accessions_list = list(dict.fromkeys([str(a) for a in (accessions or []) if a is not None])) if accessions else list(base.keys())
        if not accessions_list:
            return []

        include_paralog, include_decontam, allow_ambiguous, strict_decontam, rescue_dup = self._normalize_adjusted_flags(
            include_paralog=include_paralog,
            include_decontam=include_decontam,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )

        supported_decisions = ["support"]
        if not strict_decontam:
            supported_decisions.append("weak")
        if allow_ambiguous:
            supported_decisions.append("unknown")

        hidden_counts: dict[str, tuple[int, int]] = {}
        has_paralog: set[str] = set()
        hidden_counts = self.manager.filtering.paralog_hidden_counts(
            target_library_id=int(library_id),
            busco_library_id=busco_library_id,
            accessions=accessions_list,
            run_id=paralog_run_id,
        )
        has_paralog = {acc for acc, (_hidden, total) in hidden_counts.items() if total > 0}

        decont_primary = self.manager.filtering.latest_decont_summary(
            target_library_id=int(library_id),
            accessions=accessions_list,
            run_id=decont_run_id,
        )
        decont_latest: dict[str, tuple[str, Optional[str], Optional[str], int]] = {
            acc: (run_id, decision, date, int(library_id))
            for acc, (run_id, decision, date) in decont_primary.items()
        }
        if parent_id:
            decont_fallback = self.manager.filtering.latest_decont_summary(
                target_library_id=int(parent_id),
                accessions=accessions_list,
                run_id=decont_run_id,
            )
            for acc, (run_id, decision, date) in decont_fallback.items():
                if acc not in decont_latest:
                    decont_latest[acc] = (run_id, decision, date, int(parent_id))

        has_decont = set(decont_latest.keys())
        accessions_by_run: dict[tuple[str, int], list[str]] = {}
        for acc, (run_id, _decision, _date, source_lib) in decont_latest.items():
            accessions_by_run.setdefault((str(run_id), int(source_lib)), []).append(acc)

        contaminated_raw = (
            self.manager.filtering._decontam_contaminated_counts(
                family_library_id=int(library_id),
                target_library_id=int(library_id),
                busco_library_id=busco_library_id,
                accessions_by_run=accessions_by_run,
                supported_decisions=supported_decisions,
            )
            if accessions_by_run
            else {}
        )

        overlap = {}
        if include_paralog and include_decontam and accessions_by_run:
            overlap = self.manager.filtering._decontam_overlap_counts(
                family_library_id=int(library_id),
                target_library_id=int(library_id),
                busco_library_id=busco_library_id,
                accessions_by_run=accessions_by_run,
                supported_decisions=supported_decisions,
            )

        rescued_counts = {}
        if rescue_dup:
            rescued_counts = self.manager.filtering.get_rescued_duplicate_counts(
                target_library_id=int(library_id),
                busco_library_id=busco_library_id,
                accessions=accessions_list,
                include_paralog=include_paralog,
                include_decontam=include_decontam,
                paralog_run_id=paralog_run_id,
                accessions_by_run=accessions_by_run,
                supported_decisions=supported_decisions,
            )

        results = []
        for acc in accessions_list:
            if acc not in base:
                continue
            species, sc_val, dup_val, frag_val, _miss_val = base.get(acc, (None, 0, 0, 0, size))
            hidden_val = hidden_counts.get(acc, (0, 0))[0] if acc in has_paralog else None
            contam_val = contaminated_raw.get(acc) if acc in has_decont else None
            rescued_val = rescued_counts.get(acc, 0) if rescue_dup else 0

            if include_decontam and contam_val is not None and include_paralog and hidden_val is not None:
                contam_effective = max(contam_val - overlap.get(acc, 0), 0)
            else:
                contam_effective = contam_val if contam_val is not None else None

            sc_adj = sc_val
            dup_adj = dup_val
            if include_paralog and hidden_val is not None:
                sc_adj -= hidden_val
            if include_decontam and contam_effective is not None:
                sc_adj -= contam_effective
            if rescued_val:
                sc_adj += rescued_val
                dup_adj = max(dup_adj - rescued_val, 0)
            if sc_adj < 0:
                sc_adj = 0

            hidden_display = hidden_val if include_paralog else (0 if hidden_val is not None else None)
            contam_display = contam_effective if include_decontam else (0 if contam_effective is not None else None)

            if size > 0:
                complete_pct = round(100.0 * (sc_adj + dup_adj) / size, 2)
                single_copy_pct = round(100.0 * sc_adj / size, 2)
                duplicated_pct = round(100.0 * dup_adj / size, 2)
                frag_pct = round(100.0 * frag_val / size, 2)
                miss_final = max(size - (sc_adj + dup_adj + frag_val + (hidden_display or 0) + (contam_display or 0)), 0)
                miss_pct = round(100.0 * miss_final / size, 2)
                hidden_pct = round(100.0 * hidden_display / size, 2) if hidden_display is not None else None
                contam_pct = round(100.0 * contam_display / size, 2) if contam_display is not None else None
            else:
                complete_pct = single_copy_pct = duplicated_pct = frag_pct = miss_pct = None
                hidden_pct = contam_pct = None

            run_id = decont_latest.get(acc, (None, None, None, None))[0] if acc in has_decont else None
            decision = decont_latest.get(acc, (None, None, None, None))[1] if acc in has_decont else None
            results.append(
                (
                    acc,
                    species,
                    lib_name,
                    complete_pct,
                    single_copy_pct,
                    duplicated_pct,
                    frag_pct,
                    miss_pct,
                    hidden_pct,
                    contam_pct,
                    run_id,
                    decision,
                    acc in has_paralog,
                    acc in has_decont,
                )
            )
        return results

    @transactional("invalidate adjusted BUSCO results")
    def _mark_adjusted_results_stale(self, where_sql: str, params: Sequence[Any], *, reason: Optional[str]) -> bool:
        try:
            self.core.execute(
                f"""
                UPDATE BUSCO_Adjusted_Results
                SET status = 'stale',
                    invalidated_at = datetime('now'),
                    invalidation_reason = ?,
                    updated_at = datetime('now')
                WHERE {where_sql}
                """,
                tuple([str(reason or "dependency_changed")] + list(params)),
            )
            self.conn.commit()
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Error invalidating adjusted BUSCO results cache: {exc}") from exc

    def _dependent_library_ids_for_busco_library(self, busco_library_id: int) -> list[int]:
        rows = self.core.fetchall(
            """
            SELECT library_id
            FROM Libraries
            WHERE library_id = ? OR parent_id = ?
            """,
            (int(busco_library_id), int(busco_library_id)),
        )
        return list(dict.fromkeys(int(row[0]) for row in rows or [] if row and row[0] is not None))

    def invalidate_adjusted_results_for_library(
        self,
        library_id: int,
        *,
        accessions: Optional[Sequence[str]] = None,
        reason: Optional[str] = None,
    ) -> bool:
        clauses = ["library_id = ?"]
        params: list[Any] = [int(library_id)]
        if accessions:
            acc_vals = [str(a) for a in accessions if a is not None]
            if acc_vals:
                placeholders = ",".join("?" for _ in acc_vals)
                clauses.append(f"accession IN ({placeholders})")
                params.extend(acc_vals)
        return self._mark_adjusted_results_stale(" AND ".join(clauses), params, reason=reason)

    def invalidate_adjusted_results_for_busco_scope(
        self,
        busco_library_id: int,
        *,
        accessions: Optional[Sequence[str]] = None,
        reason: Optional[str] = None,
    ) -> bool:
        library_ids = self._dependent_library_ids_for_busco_library(int(busco_library_id))
        if not library_ids:
            return True
        placeholders = ",".join("?" for _ in library_ids)
        clauses = [f"library_id IN ({placeholders})"]
        params: list[Any] = list(library_ids)
        if accessions:
            acc_vals = [str(a) for a in accessions if a is not None]
            if acc_vals:
                acc_placeholders = ",".join("?" for _ in acc_vals)
                clauses.append(f"accession IN ({acc_placeholders})")
                params.extend(acc_vals)
        return self._mark_adjusted_results_stale(" AND ".join(clauses), params, reason=reason)

    def invalidate_adjusted_results_for_busco_run(self, run_id: int, *, reason: Optional[str] = None) -> bool:
        row = self.core.fetchone(
            "SELECT accession, library_id FROM BUSCO_Runs WHERE run_id = ?",
            (int(run_id),),
        )
        clauses = ["effective_busco_run_id = ?"]
        params: list[Any] = [int(run_id)]
        if row and row[0] is not None and row[1] is not None:
            dependent = self._dependent_library_ids_for_busco_library(int(row[1]))
            if dependent:
                lib_placeholders = ",".join("?" for _ in dependent)
                clauses.append(f"(accession = ? AND library_id IN ({lib_placeholders}))")
                params.append(str(row[0]))
                params.extend(dependent)
        return self._mark_adjusted_results_stale(" OR ".join(clauses), params, reason=reason)

    def invalidate_adjusted_results_for_decont_scope(
        self,
        library_id: int,
        *,
        accessions: Optional[Sequence[str]] = None,
        run_id: Optional[str] = None,
        reason: Optional[str] = None,
    ) -> bool:
        clauses = ["library_id = ?"]
        params: list[Any] = [int(library_id)]
        if accessions:
            acc_vals = [str(a) for a in accessions if a is not None]
            if acc_vals:
                placeholders = ",".join("?" for _ in acc_vals)
                clauses.append(f"accession IN ({placeholders})")
                params.extend(acc_vals)
        if run_id is not None:
            clauses.append("effective_decont_run_id = ?")
            params.append(str(run_id))
        return self._mark_adjusted_results_stale(" AND ".join(clauses), params, reason=reason)

    @transactional("delete adjusted BUSCO results")
    def delete_adjusted_results_for_library(self, library_id: int) -> bool:
        try:
            self.core.execute("DELETE FROM BUSCO_Adjusted_Results WHERE library_id = ?", (int(library_id),))
            self.conn.commit()
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(
                f"Error deleting adjusted BUSCO results cache for library {library_id}: {exc}"
            ) from exc

    def is_lineage_downloaded(self, lineage):
        library_id = self.manager.libraries.get_id(lineage)
        if library_id is None:
            return False
        location = self.manager.libraries.resolve_path(int(library_id))
        return bool(location and os.path.isdir(location) and os.path.exists(os.path.join(location, "dataset.cfg")))

    @transactional("add BUSCO results")
    def add_results(self, accession, library_id, busco_result, datetime=None):
        try:
            self.core.execute(
                """
                INSERT INTO BUSCO_Results (accession, library_id, no_sc_complete, no_duplicated_complete, no_fragmented, no_missing, date)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(accession, library_id) DO UPDATE SET
                    no_sc_complete = excluded.no_sc_complete,
                    no_duplicated_complete = excluded.no_duplicated_complete,
                    no_fragmented = excluded.no_fragmented,
                    no_missing = excluded.no_missing,
                    date = excluded.date
                """,
                (
                    accession,
                    int(library_id),
                    busco_result["Single copy BUSCOs"],
                    busco_result["Multi copy BUSCOs"],
                    busco_result["Fragmented BUSCOs"],
                    busco_result["Missing BUSCOs"],
                    datetime,
                ),
            )
            self.conn.commit()
            self.invalidate_adjusted_results_for_busco_scope(
                int(library_id),
                accessions=[str(accession)],
                reason="busco_results_updated",
            )
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding BUSCO results: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding BUSCO results: {exc}") from exc

    @transactional("delete BUSCO records")
    def delete_records(self, accession, library_id):
        try:
            self.core.execute(
                "DELETE FROM BUSCO_Family_Locations WHERE accession = ? AND library_id = ?",
                (accession, int(library_id)),
            )
            self.core.execute(
                "DELETE FROM BUSCO_Family_Data WHERE accession = ? AND library_id = ?",
                (accession, int(library_id)),
            )
            self.core.execute(
                "DELETE FROM BUSCO_Results WHERE accession = ? AND library_id = ?",
                (accession, int(library_id)),
            )
            self.conn.commit()
            self.invalidate_adjusted_results_for_busco_scope(
                int(library_id),
                accessions=[str(accession)],
                reason="busco_records_deleted",
            )
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(
                f"Error deleting BUSCO records for {accession} (lib={library_id}): {exc}"
            ) from exc

    @transactional("create BUSCO run")
    def create_run(
        self,
        *,
        accession: str,
        library_id: int,
        lineage_name: str,
        input_mode: str,
        pipeline: str,
        pipeline_params_effective: Optional[dict] = None,
        pipeline_params_source: Optional[dict] = None,
        busco_cli_args: Optional[list] = None,
        busco_version: Optional[str] = None,
        result_dir: Optional[str] = None,
        proteome_profile_id: Optional[int] = None,
        status: str = "running",
    ):
        try:
            self.core.execute(
                """
                INSERT INTO BUSCO_Runs (
                    accession, library_id, lineage_name, input_mode, pipeline,
                    pipeline_params_effective_json, pipeline_params_source_json,
                    busco_cli_args_json, busco_version, result_dir, proteome_profile_id, status, started_at, updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'), datetime('now'))
                """,
                (
                    accession,
                    int(library_id),
                    lineage_name,
                    input_mode,
                    pipeline,
                    json.dumps(pipeline_params_effective or {}),
                    json.dumps(pipeline_params_source or {}),
                    json.dumps(busco_cli_args or []),
                    busco_version,
                    result_dir,
                    int(proteome_profile_id) if proteome_profile_id is not None else None,
                    status,
                ),
            )
            self.conn.commit()
            return int(self.cursor.lastrowid)
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Error creating BUSCO run: {exc}") from exc

    def get_run(self, run_id: int):
        return self.core.fetchone(
            """
            SELECT r.run_id, r.accession, r.library_id, r.lineage_name, r.input_mode, r.pipeline,
                   r.result_dir, r.status, r.no_sc_complete, r.no_duplicated_complete,
                   r.no_fragmented, r.no_missing, r.started_at, r.completed_at,
                   r.proteome_profile_id, pp.profile_name, COALESCE(pp.is_default, 0) AS is_default_profile
            FROM BUSCO_Runs r
            LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
            WHERE r.run_id = ?
            """,
            (int(run_id),),
        )

    @transactional("update BUSCO run")
    def update_run(
        self,
        run_id: int,
        *,
        status: Optional[str] = None,
        result_dir: Optional[str] = None,
        lineage_name: Optional[str] = None,
        input_mode: Optional[str] = None,
        pipeline: Optional[str] = None,
        counts: Optional[dict] = None,
        completed: bool = False,
        busco_version: Optional[str] = None,
        proteome_profile_id: Optional[int] = None,
    ):
        try:
            assignments = []
            params: list[Any] = []
            if status is not None:
                assignments.append("status = ?")
                params.append(status)
            if result_dir is not None:
                assignments.append("result_dir = ?")
                params.append(result_dir)
            if lineage_name is not None:
                assignments.append("lineage_name = ?")
                params.append(lineage_name)
            if input_mode is not None:
                assignments.append("input_mode = ?")
                params.append(input_mode)
            if pipeline is not None:
                assignments.append("pipeline = ?")
                params.append(pipeline)
            if busco_version is not None:
                assignments.append("busco_version = ?")
                params.append(busco_version)
            if proteome_profile_id is not None:
                assignments.append("proteome_profile_id = ?")
                params.append(int(proteome_profile_id))
            if counts is not None:
                assignments.extend(
                    [
                        "no_sc_complete = ?",
                        "no_duplicated_complete = ?",
                        "no_fragmented = ?",
                        "no_missing = ?",
                    ]
                )
                params.extend(
                    [
                        int(counts.get("Single copy BUSCOs", 0) or 0),
                        int(counts.get("Multi copy BUSCOs", 0) or 0),
                        int(counts.get("Fragmented BUSCOs", 0) or 0),
                        int(counts.get("Missing BUSCOs", 0) or 0),
                    ]
                )
            assignments.append("updated_at = datetime('now')")
            if completed:
                assignments.append("completed_at = datetime('now')")
            params.append(int(run_id))
            self.core.execute(
                f"UPDATE BUSCO_Runs SET {', '.join(assignments)} WHERE run_id = ?",
                tuple(params),
            )
            self.conn.commit()
            if any(
                value is not None
                for value in (status, result_dir, lineage_name, input_mode, pipeline, counts, busco_version, proteome_profile_id)
            ):
                self.invalidate_adjusted_results_for_busco_run(int(run_id), reason="busco_run_updated")
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Error updating BUSCO run {run_id}: {exc}") from exc

    @staticmethod
    def _normalize_input_mode(value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        token = str(value).strip().lower()
        if token in {"nucl", "nucleotide"}:
            return "genome"
        if token == "prot":
            return "protein"
        return token or None

    @staticmethod
    def _normalize_pipeline(value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        return str(value).strip().lower() or None

    @classmethod
    def _matches_pipeline_alias(cls, requested: Optional[str], *, row_pipeline: Optional[str], row_input_mode: Optional[str]) -> bool:
        requested_norm = cls._normalize_pipeline(requested)
        if requested_norm is None:
            return True
        row_pipeline_norm = cls._normalize_pipeline(row_pipeline)
        row_mode_norm = cls._normalize_input_mode(row_input_mode)
        if requested_norm == "busco":
            # "busco" denotes a run produced by BUSCO itself, irrespective of
            # which BUSCO predictor was recorded.  Keep derived OrthoFinder
            # runs out of this group: add-library uses this alias when checking
            # that its source BUSCO evidence exists.
            return row_pipeline_norm in {"miniprot", "metaeuk", "augustus"}
        if requested_norm == "proteome":
            return row_mode_norm == "protein"
        return row_pipeline_norm == requested_norm

    @staticmethod
    def _normalize_proteome_profile(value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        return str(value).strip() or None

    def _normalize_run_ids(self, values: Optional[Sequence[Any]]) -> list[int]:
        resolved: list[int] = []
        for value in values or []:
            text = str(value).strip()
            if not text:
                continue
            if text.startswith("@") and len(text) > 1:
                env_value = self.manager.get_environment_variable(text[1:])
                if isinstance(env_value, str):
                    tokens = [part.strip() for part in env_value.split(",") if part.strip()]
                elif isinstance(env_value, (list, tuple, set)):
                    tokens = [str(part).strip() for part in env_value if str(part).strip()]
                else:
                    tokens = []
                for token in tokens:
                    if token.isdigit():
                        resolved.append(int(token))
                continue
            if text.isdigit():
                resolved.append(int(text))
        return list(dict.fromkeys(resolved))

    @staticmethod
    def _primary_sort_key(row: Sequence[Any]) -> tuple[Any, ...]:
        completed_at = str(row[13] or "")
        single_copy = int(row[9] or 0)
        total_complete = int(row[9] or 0) + int(row[10] or 0)
        duplicated = int(row[10] or 0)
        return (single_copy, total_complete, -duplicated, completed_at, int(row[0]))

    @transactional("delete BUSCO run")
    def delete_run(self, run_id: int) -> bool:
        try:
            self.invalidate_adjusted_results_for_busco_run(int(run_id), reason="busco_run_deleted")
            self.core.execute("DELETE FROM BUSCO_Primary WHERE run_id = ?", (int(run_id),))
            self.core.execute("DELETE FROM BUSCO_Run_Family_Artifacts WHERE run_id = ?", (int(run_id),))
            self.core.execute("DELETE FROM BUSCO_Run_Family_Locations WHERE run_id = ?", (int(run_id),))
            self.core.execute("DELETE FROM BUSCO_Run_Family_Data WHERE run_id = ?", (int(run_id),))
            self.core.execute(
                "DELETE FROM Artifacts WHERE owner_type = 'busco_run' AND owner_id = ?",
                (str(int(run_id)),),
            )
            self.core.execute("DELETE FROM BUSCO_Runs WHERE run_id = ?", (int(run_id),))
            self.conn.commit()
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Error deleting BUSCO run {run_id}: {exc}") from exc

    def get_run_status(self, run_id: int) -> Optional[str]:
        row = self.core.fetchone("SELECT status FROM BUSCO_Runs WHERE run_id = ?", (int(run_id),))
        return str(row[0]) if row and row[0] is not None else None

    def set_run_status(self, run_id: int, status: str) -> bool:
        return self.update_run(int(run_id), status=str(status))

    def _sync_legacy_from_run(self, *, accession: str, library_id: int, run_id: int) -> None:
        row = self.core.fetchone(
            """
            SELECT no_sc_complete, no_duplicated_complete, no_fragmented, no_missing, completed_at
            FROM BUSCO_Runs
            WHERE run_id = ? AND accession = ? AND library_id = ?
            """,
            (int(run_id), accession, int(library_id)),
        )
        if row is None:
            return
        sc, dup, frag, miss, completed_at = row
        self.core.execute(
            "DELETE FROM BUSCO_Family_Locations WHERE accession = ? AND library_id = ?",
            (accession, int(library_id)),
        )
        self.core.execute(
            "DELETE FROM BUSCO_Family_Data WHERE accession = ? AND library_id = ?",
            (accession, int(library_id)),
        )
        self.core.execute(
            "DELETE FROM BUSCO_Results WHERE accession = ? AND library_id = ?",
            (accession, int(library_id)),
        )
        self.core.execute(
            """
            INSERT INTO BUSCO_Results (accession, library_id, no_sc_complete, no_duplicated_complete, no_fragmented, no_missing, date)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                accession,
                int(library_id),
                int(sc or 0),
                int(dup or 0),
                int(frag or 0),
                int(miss or 0),
                completed_at,
            ),
        )
        self.core.execute(
            """
            INSERT INTO BUSCO_Family_Data (family_id, library_id, accession, status, sequence, score, length)
            SELECT family_id, library_id, accession, status, sequence, score, length
            FROM BUSCO_Run_Family_Data
            WHERE run_id = ? AND accession = ? AND library_id = ?
            """,
            (int(run_id), accession, int(library_id)),
        )
        self.core.execute(
            """
            INSERT INTO BUSCO_Family_Locations (family_id, library_id, accession, location)
            SELECT family_id, library_id, accession, location
            FROM BUSCO_Run_Family_Locations
            WHERE run_id = ? AND accession = ? AND library_id = ?
            """,
            (int(run_id), accession, int(library_id)),
        )

    @transactional("set BUSCO primary run")
    def set_primary_run(
        self,
        *,
        accession: str,
        library_id: int,
        run_id: int,
        purpose: str = "default",
        policy: Optional[str] = None,
        updated_by: Optional[str] = None,
    ):
        try:
            self.core.execute(
                """
                INSERT INTO BUSCO_Primary (accession, library_id, purpose, run_id, policy, updated_at, updated_by)
                VALUES (?, ?, ?, ?, ?, datetime('now'), ?)
                ON CONFLICT(accession, library_id, purpose) DO UPDATE SET
                    run_id = excluded.run_id,
                    policy = excluded.policy,
                    updated_at = excluded.updated_at,
                    updated_by = excluded.updated_by
                """,
                (accession, int(library_id), purpose, int(run_id), policy, updated_by),
            )
            if str(purpose).strip().lower() == "default":
                self._sync_legacy_from_run(accession=str(accession), library_id=int(library_id), run_id=int(run_id))
            self.conn.commit()
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(
                f"Error setting BUSCO primary run for {accession}/{library_id}: {exc}"
            ) from exc

    @transactional("clear BUSCO primary run")
    def clear_primary_run(self, accession: str, library_id: int, purpose: str = "default") -> bool:
        try:
            self.core.execute(
                """
                DELETE FROM BUSCO_Primary
                WHERE accession = ? AND library_id = ? AND purpose = ?
                """,
                (str(accession), int(library_id), str(purpose)),
            )
            self.conn.commit()
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(
                f"Error clearing BUSCO primary run for {accession}/{library_id}/{purpose}: {exc}"
            ) from exc

    def get_primary_run(self, accession, library_id, purpose: str = "default"):
        return self.core.fetchone(
            """
            SELECT r.run_id, r.accession, r.library_id, r.lineage_name, r.input_mode, r.pipeline,
                   r.result_dir, r.status, r.no_sc_complete, r.no_duplicated_complete,
                   r.no_fragmented, r.no_missing, r.completed_at,
                   r.proteome_profile_id, pp.profile_name, COALESCE(pp.is_default, 0) AS is_default_profile
            FROM BUSCO_Primary p
            JOIN BUSCO_Runs r ON r.run_id = p.run_id
            LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
            WHERE p.accession = ? AND p.library_id = ? AND p.purpose = ?
              AND COALESCE(r.status, 'completed') = 'completed'
            """,
            (accession, int(library_id), purpose),
        )

    def get_primary_assignment(self, accession: str, library_id: int, purpose: str = "default"):
        return self.core.fetchone(
            """
            SELECT p.run_id, p.policy, p.updated_by, p.updated_at
            FROM BUSCO_Primary p
            WHERE p.accession = ? AND p.library_id = ? AND p.purpose = ?
            """,
            (str(accession), int(library_id), str(purpose)),
        )

    def is_manual_primary_override(self, accession: str, library_id: int, purpose: str = "default") -> bool:
        row = self.get_primary_assignment(str(accession), int(library_id), purpose=purpose)
        if not row or row[1] is None:
            return False
        return str(row[1]).strip().lower() in self.MANUAL_PRIMARY_POLICIES

    def _run_supported_sequence_kinds(self, run_id: int) -> set[str]:
        kinds: set[str] = set()
        rows = self.core.fetchall(
            """
            SELECT DISTINCT COALESCE(sequence_kind, '')
            FROM BUSCO_Run_Family_Artifacts
            WHERE run_id = ?
            """,
            (int(run_id),),
        ) or []
        for (kind,) in rows:
            token = str(kind or "").strip().lower()
            if token in {"prot", "nucl"}:
                kinds.add(token)

        if not kinds:
            for _fam, _lib, _acc, location in self.get_run_family_locations(int(run_id)) or []:
                token = str(location or "").strip().lower()
                if not token:
                    continue
                if any(token.endswith(suffix) for suffix in (".faa", ".faa.gz", ".pep", ".pep.gz", ".aa", ".aa.gz")):
                    kinds.add("prot")
                elif any(token.endswith(suffix) for suffix in (".fna", ".fna.gz", ".ffn", ".ffn.gz", ".cds", ".cds.gz")):
                    kinds.add("nucl")

        run_row = self.get_run(int(run_id))
        pipeline = self._normalize_pipeline(run_row[5] if run_row and len(run_row) > 5 else None)
        input_mode = self._normalize_input_mode(run_row[4] if run_row and len(run_row) > 4 else None)
        if input_mode == "protein":
            kinds.add("prot")
        elif input_mode == "genome":
            kinds.add("prot")
            if pipeline in {"augustus", "metaeuk"}:
                kinds.add("nucl")
        return kinds

    def run_supports_purpose(self, run_id: int, purpose: str = "default") -> bool:
        run_row = self.get_run(int(run_id))
        if not run_row:
            return False
        if str(run_row[7] or "").strip().lower() != "completed":
            return False
        purpose_token = str(purpose or "default").strip().lower()
        if purpose_token == "default":
            return True
        kinds = self._run_supported_sequence_kinds(int(run_id))
        if purpose_token == "export_protein":
            return "prot" in kinds
        if purpose_token == "export_nucleotide":
            return "nucl" in kinds
        return False

    def get_runs_for_primary_choice(
        self,
        accession: str,
        library_id: int,
        *,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
    ):
        rows = self.get_runs_for_accessions([str(accession)], library_id=int(library_id), purpose="default") or []
        required_pipeline = self._normalize_pipeline(pipeline)
        required_mode = self._normalize_input_mode(input_mode)
        required_profile = self._normalize_proteome_profile(proteome_profile)
        allowed_run_ids = set(self._normalize_run_ids(run_ids))
        filtered = []
        for row in rows:
            run_id = int(row[0])
            row_pipeline = self._normalize_pipeline(row[6])
            row_mode = self._normalize_input_mode(row[5])
            row_status = str(row[8] or "").strip().lower()
            row_profile = self._normalize_proteome_profile(row[20] if len(row) > 20 else None)
            if row_status != "completed":
                continue
            if allowed_run_ids and run_id not in allowed_run_ids:
                continue
            if not self._matches_pipeline_alias(required_pipeline, row_pipeline=row_pipeline, row_input_mode=row_mode):
                continue
            if required_mode and row_mode != required_mode:
                continue
            if required_profile and row_profile != required_profile:
                continue
            filtered.append(row)
        return filtered

    def choose_best_run(
        self,
        accession: str,
        library_id: int,
        *,
        purpose: str = "default",
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
    ) -> Optional[tuple]:
        preferred_pipeline = self._normalize_pipeline(preferred_pipeline)
        preferred_mode = self._normalize_input_mode(preferred_input_mode)
        preferred_profile = self._normalize_proteome_profile(preferred_proteome_profile)
        usable = [
            row
            for row in self.get_runs_for_primary_choice(
                accession=str(accession),
                library_id=int(library_id),
                run_ids=run_ids,
                pipeline=pipeline,
                input_mode=input_mode,
                proteome_profile=proteome_profile,
            )
            if self.run_supports_purpose(int(row[0]), purpose=purpose)
        ]
        if not usable:
            return None

        def _sort_key(row: Sequence[Any]) -> tuple[Any, ...]:
            row_pipeline = self._normalize_pipeline(row[6])
            row_mode = self._normalize_input_mode(row[5])
            row_profile = self._normalize_proteome_profile(row[20] if len(row) > 20 else None)
            prefer_hits = (
                1 if self._matches_pipeline_alias(preferred_pipeline, row_pipeline=row_pipeline, row_input_mode=row_mode) else 0,
                1 if preferred_mode and row_mode == preferred_mode else 0,
                1 if preferred_profile and row_profile == preferred_profile else 0,
            )
            return prefer_hits + self._primary_sort_key(row)

        return sorted(usable, key=_sort_key, reverse=True)[0]

    def refresh_auto_primary_runs_for_accession(
        self,
        accession: str,
        library_id: int,
        *,
        updated_by: Optional[str] = None,
        policy: str = "auto_best",
    ) -> dict[str, Optional[int]]:
        results: dict[str, Optional[int]] = {}
        default_profile = self.manager.proteomes.get_default_profile_name(str(accession))
        for purpose in ("default", "export_protein", "export_nucleotide"):
            if self.is_manual_primary_override(str(accession), int(library_id), purpose=purpose):
                current = self.get_primary_assignment(str(accession), int(library_id), purpose=purpose)
                results[purpose] = int(current[0]) if current and current[0] is not None else None
                continue
            best = None
            if default_profile:
                profile_runs = [
                    row
                    for row in self.get_runs_for_primary_choice(
                        str(accession),
                        int(library_id),
                        proteome_profile=default_profile,
                    )
                    if self.run_supports_purpose(int(row[0]), purpose=purpose)
                ]
                if profile_runs:
                    best = sorted(
                        profile_runs,
                        key=lambda row: (
                            str(row[13] or ""),
                            int(row[0]),
                        ),
                        reverse=True,
                    )[0]
            if best is None:
                best = self.choose_best_run(
                    str(accession),
                    int(library_id),
                    purpose=purpose,
                    preferred_proteome_profile=default_profile,
                )
            if best is None:
                self.clear_primary_run(str(accession), int(library_id), purpose=purpose)
                results[purpose] = None
                continue
            run_id = int(best[0])
            self.set_primary_run(
                accession=str(accession),
                library_id=int(library_id),
                run_id=run_id,
                purpose=purpose,
                policy=policy,
                updated_by=updated_by,
            )
            results[purpose] = run_id
        return results

    def get_runs_for_accessions(
        self,
        accessions: Sequence[str],
        *,
        library_id: Optional[int] = None,
        purpose: str = "default",
    ):
        if not accessions:
            return []
        params: list[Any] = [str(a) for a in accessions]
        placeholders = ",".join("?" for _ in params)
        where = [f"r.accession IN ({placeholders})", "COALESCE(r.status, 'completed') = 'completed'"]
        if library_id is not None:
            where.append("r.library_id = ?")
            params.append(int(library_id))
        sql = f"""
            SELECT
                r.run_id,
                r.accession,
                r.library_id,
                l.library_name,
                r.lineage_name,
                r.input_mode,
                r.pipeline,
                r.result_dir,
                r.status,
                r.no_sc_complete,
                r.no_duplicated_complete,
                r.no_fragmented,
                r.no_missing,
                r.completed_at,
                CASE WHEN p.run_id IS NULL THEN 0 ELSE 1 END AS is_primary,
                CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * (COALESCE(r.no_sc_complete,0) + COALESCE(r.no_duplicated_complete,0)) / l.size, 2) END AS complete,
                CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * COALESCE(r.no_sc_complete,0) / l.size, 2) END AS single_copy_complete,
                CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * COALESCE(r.no_duplicated_complete,0) / l.size, 2) END AS duplicated,
                CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * COALESCE(r.no_fragmented,0) / l.size, 2) END AS fragmented,
                CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * COALESCE(r.no_missing,0) / l.size, 2) END AS missing,
                pp.profile_name AS proteome_profile,
                COALESCE(pp.is_default, 0) AS is_default_profile
            FROM BUSCO_Runs r
            JOIN Libraries l ON l.library_id = r.library_id
            LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
            LEFT JOIN BUSCO_Primary p
                ON p.accession = r.accession
               AND p.library_id = r.library_id
               AND p.purpose = ?
               AND p.run_id = r.run_id
            WHERE {" AND ".join(where)}
            ORDER BY r.accession, is_primary DESC, r.completed_at DESC, r.run_id DESC
        """
        return self.core.fetchall(sql, tuple([purpose] + params)) or []

    def count_run_family_rows(self, run_id: int) -> int:
        row = self.core.fetchone(
            "SELECT COUNT(*) FROM BUSCO_Run_Family_Data WHERE run_id = ?",
            (int(run_id),),
        )
        return int(row[0]) if row and row[0] is not None else 0

    def get_run_family_locations(self, run_id: int):
        return self.core.fetchall(
            """
            SELECT family_id, library_id, accession, location
            FROM BUSCO_Run_Family_Locations
            WHERE run_id = ?
            ORDER BY family_id, accession
            """,
            (int(run_id),),
        )

    def get_run_family_data(self, run_id: int):
        return self.core.fetchall(
            """
            SELECT family_id, library_id, accession, status, sequence, score, length
            FROM BUSCO_Run_Family_Data
            WHERE run_id = ?
            ORDER BY family_id, accession
            """,
            (int(run_id),),
        ) or []

    def get_primary_result_dir(
        self,
        accession,
        library_id,
        purpose: str = "default",
        *,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
    ):
        row = self.get_primary_run(accession, library_id, purpose=purpose)
        if proteome_profile is not None:
            row_profile = self._normalize_proteome_profile(row[14] if row and len(row) > 14 else None)
            if row_profile != self._normalize_proteome_profile(proteome_profile):
                row = self.choose_best_run(
                    str(accession),
                    int(library_id),
                    purpose=purpose,
                    proteome_profile=proteome_profile,
                    preferred_proteome_profile=preferred_proteome_profile,
                )
        rdir = self._strict_result_dir_for_run_row(row)
        if rdir:
            return rdir
        fallback = self.core.fetchone(
            """
            SELECT result_dir
            FROM BUSCO_Runs
            WHERE accession = ? AND library_id = ? AND status = 'completed'
            ORDER BY completed_at DESC, run_id DESC
            LIMIT 1
            """,
            (accession, int(library_id)),
        )
        if fallback and fallback[0]:
            return str(fallback[0])
        return None

    def get_result_dir_for_query(
        self,
        accession: str,
        library_id: int,
        *,
        purpose: str = "default",
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
        selection: Optional[str] = None,
        allow_legacy_fallback: bool = True,
    ) -> Optional[str]:
        run_map = self._resolve_busco_runs_for_query(
            int(library_id),
            accessions=[str(accession)],
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            proteome_profile=proteome_profile,
            preferred_proteome_profile=preferred_proteome_profile,
            selection=selection,
            purpose=purpose,
        )
        run_id = run_map.get(str(accession))
        if run_id is None:
            return None
        row = self.get_run(int(run_id))
        rdir = self._strict_result_dir_for_run_row(row)
        if rdir:
            return rdir
        return None

    def get_effective_run_ids_for_accessions(
        self,
        library_id: int,
        accessions: Optional[Sequence[str]] = None,
        *,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ) -> dict[str, int]:
        return self._resolve_busco_runs_for_query(
            int(library_id),
            accessions=accessions,
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            proteome_profile=proteome_profile,
            preferred_proteome_profile=preferred_proteome_profile,
            selection=selection,
            purpose=purpose,
        )

    def get_effective_run_id_for_accession(
        self,
        accession: str,
        library_id: int,
        *,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ) -> Optional[int]:
        run_map = self.get_effective_run_ids_for_accessions(
            int(library_id),
            accessions=[str(accession)],
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            proteome_profile=proteome_profile,
            preferred_proteome_profile=preferred_proteome_profile,
            selection=selection,
            purpose=purpose,
        )
        run_id = run_map.get(str(accession))
        return int(run_id) if run_id is not None else None

    def get_family_results_for_library(
        self,
        library_id,
        accessions=None,
        status=None,
        *,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ):
        requested_accessions = self._normalize_accessions(accessions)
        run_map = self._resolve_busco_runs_for_query(
            int(library_id),
            accessions=requested_accessions,
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            proteome_profile=proteome_profile,
            preferred_proteome_profile=preferred_proteome_profile,
            selection=selection,
            purpose=purpose,
        )
        strict_rows: list[tuple] = []
        if not run_map:
            return []
        return self._fetch_family_rows(
            table="BUSCO_Run_Family_Data",
            library_id=int(library_id),
            accessions=requested_accessions,
            status=status,
            run_ids=list(dict.fromkeys(run_map.values())),
        )

    def get_family_presence_map(
        self,
        library_id: int,
        *,
        accessions: Optional[Sequence[str]] = None,
        family_ids: Optional[Sequence[str]] = None,
        status: Optional[Sequence[int]] = None,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ) -> dict[str, set[str]]:
        rows = self.get_family_results_for_library(
            library_id=library_id,
            accessions=accessions,
            status=status,
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            selection=selection,
            purpose=purpose,
        )
        family_filter = {str(fam) for fam in family_ids or [] if fam is not None}
        presence: dict[str, set[str]] = {}
        for family_id, _lib_id, accession, *_rest in rows or []:
            fam_token = str(family_id)
            if family_filter and fam_token not in family_filter:
                continue
            presence.setdefault(fam_token, set()).add(str(accession))
        return presence

    def get_family_counts_by_accession(
        self,
        library_id: int,
        *,
        accessions: Optional[Sequence[str]] = None,
        status: Optional[Sequence[int]] = None,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ) -> dict[str, int]:
        rows = self.get_family_results_for_library(
            library_id=library_id,
            accessions=accessions,
            status=status,
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            selection=selection,
            purpose=purpose,
        )
        counts: dict[str, set[str]] = {}
        for family_id, _lib_id, accession, *_rest in rows or []:
            counts.setdefault(str(accession), set()).add(str(family_id))
        return {acc: len(families) for acc, families in counts.items()}

    def count_existing_family_locations_by_accession(
        self,
        library_id: int,
        *,
        accessions: Optional[Sequence[str]] = None,
        status: Optional[Sequence[int]] = None,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
        sequence_kind: Optional[str] = None,
    ) -> dict[str, int]:
        rows = self.get_family_results_for_library(
            library_id=library_id,
            accessions=accessions,
            status=status,
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            selection=selection,
            purpose=purpose,
        )
        counts: dict[str, int] = {}
        for family_id, _lib_id, accession, *_rest in rows or []:
            location = self.get_family_location(
                family_id,
                library_id,
                accession,
                sequence_kind=sequence_kind,
                pipeline=pipeline,
                input_mode=input_mode,
                selection=selection,
                purpose=purpose,
            )
            if location and os.path.exists(location):
                counts[str(accession)] = counts.get(str(accession), 0) + 1
        return counts

    @transactional("add legacy BUSCO family data")
    def add_legacy_family_data(self, family_data_list):
        try:
            pairs = {(row[0], row[1]) for row in family_data_list}
            if pairs:
                self.core.executemany(
                    "INSERT OR IGNORE INTO BUSCO_descriptions (family_id, library_id, description, link) VALUES (?, ?, ?, ?)",
                    [(fam, lib, None, "Boo") for fam, lib in pairs],
                )
            self.core.executemany(
                "INSERT OR REPLACE INTO BUSCO_Family_Data (family_id, library_id, accession, status, sequence, score, length) VALUES (?, ?, ?, ?, ?, ?, ?)",
                family_data_list,
            )
            self.conn.commit()
            invalidation: dict[int, set[str]] = {}
            for _fam, lib, acc, *_rest in family_data_list:
                invalidation.setdefault(int(lib), set()).add(str(acc))
            for lib, accs in invalidation.items():
                self.invalidate_adjusted_results_for_busco_scope(lib, accessions=sorted(accs), reason="legacy_busco_family_data_updated")
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding multiple BUSCO family data: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding multiple BUSCO family data: {exc}") from exc

    @transactional("add legacy BUSCO family locations")
    def add_legacy_family_locations(self, family_locations):
        if not family_locations:
            return True
        try:
            pairs = {(fl[0], fl[1]) for fl in family_locations}
            if pairs:
                self.core.executemany(
                    "INSERT OR IGNORE INTO BUSCO_descriptions (family_id, library_id, description, link) VALUES (?, ?, ?, ?)",
                    [(fam, lib, None, "Boo") for fam, lib in pairs],
                )
            self.core.executemany(
                "INSERT OR REPLACE INTO BUSCO_Family_Locations (family_id, library_id, accession, location) VALUES (?, ?, ?, ?)",
                family_locations,
            )
            self.conn.commit()
            invalidation: dict[int, set[str]] = {}
            for _fam, lib, acc, _loc in family_locations:
                invalidation.setdefault(int(lib), set()).add(str(acc))
            for lib, accs in invalidation.items():
                self.invalidate_adjusted_results_for_busco_scope(lib, accessions=sorted(accs), reason="legacy_busco_family_locations_updated")
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding multiple BUSCO family locations: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding multiple BUSCO family locations: {exc}") from exc

    @transactional("add BUSCO run family data")
    def add_run_family_data(self, run_id: int, family_data_list):
        if not family_data_list:
            return True
        rows = []
        for item in family_data_list:
            if len(item) != 7:
                continue
            family_id, library_id, accession, status, sequence, score, length = item
            rows.append((int(run_id), family_id, int(library_id), accession, status, sequence, score, length))
        if not rows:
            return True
        try:
            family_pairs = sorted({(row[1], row[2]) for row in rows})
            self.core.executemany(
                """
                INSERT OR IGNORE INTO BUSCO_descriptions
                    (family_id, library_id, description, link)
                VALUES (?, ?, NULL, NULL)
                """,
                family_pairs,
            )
            self.core.executemany(
                """
                INSERT OR REPLACE INTO BUSCO_Run_Family_Data
                    (run_id, family_id, library_id, accession, status, sequence, score, length)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                rows,
            )
            self.conn.commit()
            self.invalidate_adjusted_results_for_busco_run(int(run_id), reason="busco_run_family_data_updated")
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding BUSCO run family data: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding BUSCO run family data: {exc}") from exc

    @transactional("add BUSCO run family locations")
    def add_run_family_locations(self, run_id: int, family_locations):
        if not family_locations:
            return True
        try:
            pairs = {(fl[0], fl[1]) for fl in family_locations if len(fl) >= 2}
            if pairs:
                self.core.executemany(
                    "INSERT OR IGNORE INTO BUSCO_descriptions (family_id, library_id, description, link) VALUES (?, ?, ?, ?)",
                    [(fam, lib, None, "Boo") for fam, lib in pairs],
                )
            rows = []
            for item in family_locations:
                if len(item) != 4:
                    continue
                family_id, library_id, accession, location = item
                rows.append((int(run_id), family_id, int(library_id), accession, location))
            if not rows:
                return True
            self.core.executemany(
                """
                INSERT OR REPLACE INTO BUSCO_Run_Family_Locations
                    (run_id, family_id, library_id, accession, location)
                VALUES (?, ?, ?, ?, ?)
                """,
                rows,
            )
            self.conn.commit()
            self.invalidate_adjusted_results_for_busco_run(int(run_id), reason="busco_run_family_locations_updated")
            return True
        except sqlite3.Error as exc:
            raise RepositoryWriteError(f"Error adding BUSCO run family locations: {exc}") from exc
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Unexpected error adding BUSCO run family locations: {exc}") from exc

    @transactional("project BUSCO runs to legacy tables")
    def project_legacy(self, accession=None, library_id=None, **kwargs):
        if accession is not None and library_id is not None:
            accessions = [str(accession)]
            lib_id = int(library_id)
        else:
            accessions = [str(a) for a in kwargs.get("accessions", [])]
            lib_id = int(kwargs["library_id"])
        pipeline = kwargs.get("pipeline")
        input_mode = kwargs.get("input_mode")
        run_ids = kwargs.get("run_ids")
        preferred_pipeline = kwargs.get("preferred_pipeline")
        preferred_input_mode = kwargs.get("preferred_input_mode")
        selection = kwargs.get("selection")
        purpose = kwargs.get("purpose", "default")
        if not accessions:
            return True
        run_map = self._resolve_busco_runs_for_query(
            lib_id,
            accessions=accessions,
            run_ids=run_ids,
            pipeline=pipeline,
            input_mode=input_mode,
            preferred_pipeline=preferred_pipeline,
            preferred_input_mode=preferred_input_mode,
            selection=selection,
            purpose=purpose,
        )
        try:
            for acc in accessions:
                run_id = run_map.get(str(acc))
                if run_id is None:
                    continue
                self._sync_legacy_from_run(accession=str(acc), library_id=lib_id, run_id=int(run_id))
            self.conn.commit()
            return True
        except Exception as exc:  # boundary: repository write failures are wrapped as typed errors.
            raise RepositoryWriteError(f"Error projecting BUSCO runs to legacy tables: {exc}") from exc

    def get_family_artifact(self, *, run_id: int, accession: str, family_id: str, sequence_kind: Optional[str] = None):
        clauses = ["run_id = ?", "accession = ?", "family_id = ?"]
        params: list[Any] = [int(run_id), str(accession), str(family_id)]
        if sequence_kind:
            clauses.append("sequence_kind = ?")
            params.append(str(sequence_kind))
        sql = (
            "SELECT run_id, family_id, library_id, accession, artifact_id, sequence_kind, location, metadata_json "
            "FROM BUSCO_Run_Family_Artifacts WHERE "
            + " AND ".join(clauses)
            + " ORDER BY CASE COALESCE(sequence_kind, '') WHEN 'prot' THEN 0 WHEN 'nucl' THEN 1 ELSE 2 END, artifact_id ASC"
        )
        rows = self.core.fetchall(sql, tuple(params))
        return rows[0] if rows else None

    def get_export_family_rows(
        self,
        *,
        run_ids: Sequence[int],
        sequence_kind: str,
        status: Optional[Sequence[int]] = None,
        accessions: Optional[Sequence[str]] = None,
    ) -> list[dict[str, Any]]:
        run_vals = self._normalize_run_ids(run_ids)
        if not run_vals:
            return []
        status_vals = self._normalize_statuses(status)
        accession_vals = self._normalize_accessions(accessions)
        results: list[dict[str, Any]] = []
        seq_kind = str(sequence_kind or "").strip().lower()
        for chunk in self._chunked(run_vals):
            placeholders = ",".join("?" for _ in chunk)
            clauses = [f"d.run_id IN ({placeholders})"]
            params: list[Any] = list(chunk)
            if status_vals:
                status_placeholders = ",".join("?" for _ in status_vals)
                clauses.append(f"d.status IN ({status_placeholders})")
                params.extend(status_vals)
            if accession_vals:
                accession_placeholders = ",".join("?" for _ in accession_vals)
                clauses.append(f"d.accession IN ({accession_placeholders})")
                params.extend(accession_vals)
            clauses.append("COALESCE(a.sequence_kind, '') = ?")
            params.append(seq_kind)
            rows = self.core.fetchall(
                f"""
                SELECT
                    d.run_id,
                    d.family_id,
                    d.library_id,
                    d.accession,
                    d.status,
                    d.sequence,
                    d.score,
                    d.length,
                    a.artifact_id,
                    a.location,
                    a.metadata_json,
                    l.location
                FROM BUSCO_Run_Family_Data d
                LEFT JOIN BUSCO_Run_Family_Artifacts a
                  ON a.run_id = d.run_id
                 AND a.family_id = d.family_id
                 AND a.library_id = d.library_id
                 AND a.accession = d.accession
                LEFT JOIN BUSCO_Run_Family_Locations l
                  ON l.run_id = d.run_id
                 AND l.family_id = d.family_id
                 AND l.library_id = d.library_id
                 AND l.accession = d.accession
                WHERE {" AND ".join(clauses)}
                ORDER BY d.run_id, d.family_id, d.accession, d.sequence
                """,
                tuple(params),
            ) or []
            for row in rows:
                metadata = None
                if row[10]:
                    metadata = self.manager._json_load(row[10])
                resolved_artifact_path = None
                if row[8] is not None:
                    resolved_artifact_path = self.manager.artifacts.resolve_path(int(row[8]))
                results.append(
                    {
                        "run_id": int(row[0]),
                        "family_id": str(row[1]),
                        "library_id": int(row[2]),
                        "accession": str(row[3]),
                        "status": int(row[4]) if row[4] is not None else None,
                        "sequence_id": row[5],
                        "score": row[6],
                        "length": row[7],
                        "artifact_id": int(row[8]) if row[8] is not None else None,
                        "artifact_location": str(row[9]) if row[9] else None,
                        "artifact_path": str(resolved_artifact_path) if resolved_artifact_path else None,
                        "artifact_metadata": metadata if isinstance(metadata, dict) else {},
                        "legacy_location": str(row[11]) if row[11] else None,
                        "sequence_kind": seq_kind,
                    }
                )
        return results

    @transactional("link BUSCO family artifact")
    def link_family_artifact(
        self,
        *,
        run_id: int,
        family_id: str,
        library_id: int,
        accession: str,
        artifact_id: Optional[int],
        sequence_kind: Optional[str],
        location: Optional[str],
        metadata: Optional[dict] = None,
    ) -> bool:
        self.core.execute(
            """
            INSERT INTO BUSCO_Run_Family_Artifacts (
                run_id, family_id, library_id, accession, artifact_id, sequence_kind, location, metadata_json
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id, family_id, accession, sequence_kind) DO UPDATE SET
                artifact_id = excluded.artifact_id,
                location = excluded.location,
                metadata_json = excluded.metadata_json
            """,
            (
                int(run_id),
                str(family_id),
                int(library_id),
                str(accession),
                int(artifact_id) if artifact_id is not None else None,
                sequence_kind or "",
                location,
                self.manager._json_dump(metadata or {}),
            ),
        )
        self.conn.commit()
        return True

    def _resolve_busco_runs_for_query(
        self,
        library_id: int,
        accessions: Optional[Sequence[str]] = None,
        *,
        run_ids: Optional[Sequence[int]] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        preferred_pipeline: Optional[str] = None,
        preferred_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ) -> dict[str, int]:
        ctx = self.manager._get_busco_context()
        pipeline_val = self._normalize_pipeline(pipeline if pipeline is not None else ctx.get("pipeline"))
        mode_val = self._normalize_input_mode(input_mode if input_mode is not None else ctx.get("input_mode"))
        preferred_pipeline_val = self._normalize_pipeline(
            preferred_pipeline if preferred_pipeline is not None else ctx.get("prefer_pipeline")
        )
        preferred_mode_val = self._normalize_input_mode(
            preferred_input_mode if preferred_input_mode is not None else ctx.get("prefer_input_mode")
        )
        profile_val = self._normalize_proteome_profile(
            proteome_profile if proteome_profile is not None else ctx.get("proteome_profile")
        )
        preferred_profile_val = self._normalize_proteome_profile(
            preferred_proteome_profile if preferred_proteome_profile is not None else ctx.get("prefer_proteome_profile")
        )
        policy = (selection if selection is not None else ctx.get("selection") or "primary").strip().lower()
        params: list[Any] = []
        clauses = []
        effective_run_ids = run_ids if run_ids is not None else ctx.get("run_ids")
        run_id_vals = self._normalize_run_ids(effective_run_ids)
        if accessions:
            vals = [str(a) for a in accessions]
            placeholders = ",".join("?" for _ in vals)
            clauses.append(f"r.accession IN ({placeholders})")
            params.extend(vals)
        if pipeline_val:
            if pipeline_val == "proteome":
                clauses.append("LOWER(r.input_mode) = 'protein'")
            elif pipeline_val == "busco":
                clauses.append("LOWER(r.pipeline) IN ('miniprot', 'metaeuk', 'augustus')")
            else:
                clauses.append("LOWER(r.pipeline) = ?")
                params.append(str(pipeline_val).lower())
        if mode_val:
            clauses.append("LOWER(r.input_mode) = ?")
            params.append(str(mode_val).lower())
        if profile_val and profile_val != DEFAULT_CLEAN_PROFILE:
            clauses.append("pp.profile_name = ?")
            params.append(profile_val)
        if run_id_vals:
            placeholders = ",".join("?" for _ in run_id_vals)
            clauses.append(f"r.run_id IN ({placeholders})")
            params.extend(run_id_vals)

        resolved_from_primary: dict[str, int] = {}
        if policy in {"primary", "default"}:
            sql = f"""
                SELECT p.accession, p.run_id, r.pipeline, r.input_mode, pp.profile_name
                FROM BUSCO_Primary p
                JOIN BUSCO_Runs r ON r.run_id = p.run_id
                LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
                WHERE p.library_id = ? AND p.purpose = ? AND r.library_id = ?
                  AND COALESCE(r.status, 'completed') = 'completed'
                {" AND " + " AND ".join(clauses) if clauses else ""}
            """
            rows = self.core.fetchall(sql, tuple([int(library_id), purpose, int(library_id)] + params))
            for acc, run_id, row_pipeline, row_input_mode, row_profile in rows or []:
                if acc is None or run_id is None:
                    continue
                acc_token = str(acc)
                if not (
                    preferred_pipeline_val is None
                    and preferred_mode_val is None
                    and preferred_profile_val is None
                ):
                    if preferred_pipeline_val is not None and not self._matches_pipeline_alias(
                        preferred_pipeline_val,
                        row_pipeline=row_pipeline,
                        row_input_mode=row_input_mode,
                    ):
                        continue
                    if preferred_mode_val is not None and self._normalize_input_mode(row_input_mode) != preferred_mode_val:
                        continue
                    if preferred_profile_val is not None and not self.manager.proteomes.profile_matches_selector(
                        acc_token,
                        self._normalize_proteome_profile(row_profile),
                        preferred_profile_val,
                    ):
                        continue
                resolved_from_primary[acc_token] = int(run_id)
            if resolved_from_primary and (
                accessions is None
                or len(resolved_from_primary) == len({str(acc) for acc in accessions if acc is not None})
            ):
                return resolved_from_primary

        sql = f"""
            SELECT r.accession, r.run_id, r.pipeline, r.input_mode, r.completed_at, pp.profile_name
            FROM BUSCO_Runs r
            LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
            WHERE r.library_id = ?
              AND COALESCE(r.status, 'completed') = 'completed'
            {" AND " + " AND ".join(clauses) if clauses else ""}
            ORDER BY r.accession, r.completed_at DESC, r.run_id DESC
        """
        rows = self.core.fetchall(sql, tuple([int(library_id)] + params))
        grouped: dict[str, list[tuple[Any, ...]]] = {}
        for row in rows or []:
            acc = str(row[0])
            if acc in resolved_from_primary:
                continue
            row_profile = self._normalize_proteome_profile(row[5])
            if profile_val and not self.manager.proteomes.profile_matches_selector(acc, row_profile, profile_val):
                continue
            grouped.setdefault(acc, []).append(row)
        latest: dict[str, int] = dict(resolved_from_primary)
        for acc, acc_rows in grouped.items():
            preferred_rows = [
                row
                for row in acc_rows
                if (
                    (preferred_pipeline_val is None or self._matches_pipeline_alias(preferred_pipeline_val, row_pipeline=row[2], row_input_mode=row[3]))
                    and (preferred_mode_val is None or self._normalize_input_mode(row[3]) == preferred_mode_val)
                    and (preferred_profile_val is None or self.manager.proteomes.profile_matches_selector(acc, self._normalize_proteome_profile(row[5]), preferred_profile_val))
                )
            ]
            pool = preferred_rows or acc_rows
            if pool and pool[0][1] is not None:
                latest[str(acc)] = int(pool[0][1])
        return latest

    def get_family_location(
        self,
        family_id,
        library_id,
        accession,
        *,
        sequence_kind: Optional[str] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        preferred_proteome_profile: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ):
        run_map = self._resolve_busco_runs_for_query(
            int(library_id),
            accessions=[str(accession)],
            pipeline=pipeline,
            input_mode=input_mode,
            proteome_profile=proteome_profile,
            preferred_proteome_profile=preferred_proteome_profile,
            selection=selection,
            purpose=purpose,
        )
        run_id = run_map.get(str(accession)) if run_map else None
        kinds = [str(sequence_kind)] if sequence_kind else ["prot", "nucl", None]
        if run_id is not None:
            for kind in kinds:
                art = self.get_family_artifact(
                    run_id=int(run_id),
                    accession=str(accession),
                    family_id=str(family_id),
                    sequence_kind=kind,
                )
                if art:
                    artifact_id = art[4]
                    fallback_location = art[6]
                    if artifact_id is not None:
                        resolved = self.manager.artifacts.resolve_path(int(artifact_id))
                        if resolved:
                            return resolved
                    if fallback_location:
                        return str(fallback_location)
            row = self.core.fetchone(
                """
                SELECT location
                FROM BUSCO_Run_Family_Locations
                WHERE run_id = ? AND family_id = ? AND library_id = ? AND accession = ?
                """,
                (int(run_id), family_id, int(library_id), accession),
            )
            if row and row[0]:
                return row[0]
        return None

    def get_family_location_legacy(
        self,
        family_id,
        library_id,
        accession,
        *,
        sequence_kind: Optional[str] = None,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        selection: Optional[str] = None,
        purpose: str = "default",
    ):
        run_map = self._resolve_busco_runs_for_query(
            int(library_id),
            accessions=[str(accession)],
            pipeline=pipeline,
            input_mode=input_mode,
            selection=selection,
            purpose=purpose,
        )
        run_id = run_map.get(str(accession)) if run_map else None
        if run_id is not None:
            row = self.core.fetchone(
                """
                SELECT location
                FROM BUSCO_Run_Family_Locations
                WHERE run_id = ? AND family_id = ? AND library_id = ? AND accession = ?
                """,
                (int(run_id), family_id, int(library_id), accession),
            )
            if row and row[0]:
                return row[0]
        result = self.core.fetchone(
            "SELECT location FROM BUSCO_Family_Locations WHERE family_id = ? AND library_id = ? AND accession = ?",
            (family_id, int(library_id), accession),
        )
        return result[0] if result else None

    def register_run_artifact(self, run_id: int, artifact_type: str, path: str, **kwargs):
        return self.manager.artifacts.register(
            owner_type="busco_run",
            owner_id=int(run_id),
            artifact_type=str(artifact_type),
            path=path,
            **kwargs,
        )

    def register_family_artifact(
        self,
        *,
        run_id: int,
        family_id: str,
        library_id: int,
        accession: str,
        path: str,
        sequence_kind: Optional[str],
        role: Optional[str] = None,
        format: Optional[str] = None,
        metadata: Optional[dict] = None,
    ) -> Optional[int]:
        artifact_id = self.manager.artifacts.register(
            owner_type="busco_run",
            owner_id=int(run_id),
            artifact_type="busco_family_sequence",
            path=path,
            role=role or family_id,
            format=format,
            sequence_kind=sequence_kind,
            metadata=metadata or {"family_id": family_id, "accession": accession},
        )
        self.link_family_artifact(
            run_id=int(run_id),
            family_id=str(family_id),
            library_id=int(library_id),
            accession=str(accession),
            artifact_id=artifact_id,
            sequence_kind=sequence_kind,
            location=os.path.abspath(str(path)),
            metadata=metadata,
        )
        return artifact_id
