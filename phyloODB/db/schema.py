from __future__ import annotations

import hashlib
import json

from ..variable_kinds import infer_variable_kind

CURRENT_SCHEMA_VERSION = 5


def _legacy_paralog_run_id(target_library_id, busco_library_id) -> str:
    payload = {
        "target_library_id": int(target_library_id),
        "busco_library_id": int(busco_library_id),
    }
    digest = hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:12]
    return f"paralog_legacy_{target_library_id}_{busco_library_id}_{digest}"


def _ensure_paralog_filtering_run_schema(manager) -> None:
    cursor = manager.cursor
    conn = manager.conn
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS Paralog_Filtering_Runs (
            run_id TEXT PRIMARY KEY,
            target_library_id INT NOT NULL,
            busco_library_id INT NOT NULL,
            targets_json TEXT,
            accessions_json TEXT,
            ref_accessions_json TEXT,
            selection_mode TEXT,
            selection_params_json TEXT,
            config_signature TEXT,
            run_label TEXT,
            report_dir TEXT,
            date DATETIME DEFAULT (datetime('now')),
            FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
            FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id)
        )
        """
    )

    if manager._table_exists("Paralog_Filtering"):
        columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Paralog_Filtering)").fetchall() or [])}
        if "run_id" not in columns:
            cursor.execute("ALTER TABLE Paralog_Filtering RENAME TO Paralog_Filtering_Legacy")
            cursor.execute(
                """
                CREATE TABLE Paralog_Filtering (
                    family_id VARCHAR(20),
                    library_id INT,
                    target_library_id INT,
                    accession VARCHAR(50),
                    run_id TEXT,
                    clean BOOLEAN,
                    selected_ref_count INT,
                    selection_threshold FLOAT,
                    reused INTEGER DEFAULT 0,
                    reason_code TEXT,
                    selection_signature TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (family_id, target_library_id, accession, run_id),
                    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (run_id) REFERENCES Paralog_Filtering_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT family_id, library_id, target_library_id, accession, clean, date
                FROM Paralog_Filtering_Legacy
                """
            ).fetchall() or []
            run_pairs = sorted({(int(row[2]), int(row[1])) for row in rows if row[2] is not None and row[1] is not None})
            for target_library_id, busco_library_id in run_pairs:
                run_id = _legacy_paralog_run_id(target_library_id, busco_library_id)
                cursor.execute(
                    """
                    INSERT OR IGNORE INTO Paralog_Filtering_Runs (
                        run_id, target_library_id, busco_library_id, selection_mode,
                        selection_params_json, config_signature, run_label, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, datetime('now'))
                    """,
                    (
                        run_id,
                        target_library_id,
                        busco_library_id,
                        "legacy",
                        json.dumps({}),
                        run_id,
                        "legacy-migrated",
                    ),
                )
            for family_id, library_id, target_library_id, accession, clean, date in rows:
                run_id = _legacy_paralog_run_id(target_library_id, library_id)
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Paralog_Filtering (
                        family_id, library_id, target_library_id, accession, run_id, clean,
                        selected_ref_count, selection_threshold, reused, reason_code,
                        selection_signature, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        family_id,
                        library_id,
                        target_library_id,
                        accession,
                        run_id,
                        clean,
                        None,
                        None,
                        0,
                        "legacy",
                        run_id,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Paralog_Filtering_Legacy")
        else:
            if "selected_ref_count" not in columns:
                cursor.execute("ALTER TABLE Paralog_Filtering ADD COLUMN selected_ref_count INT")
            if "selection_threshold" not in columns:
                cursor.execute("ALTER TABLE Paralog_Filtering ADD COLUMN selection_threshold FLOAT")
            if "reused" not in columns:
                cursor.execute("ALTER TABLE Paralog_Filtering ADD COLUMN reused INTEGER DEFAULT 0")
            if "reason_code" not in columns:
                cursor.execute("ALTER TABLE Paralog_Filtering ADD COLUMN reason_code TEXT")
            if "selection_signature" not in columns:
                cursor.execute("ALTER TABLE Paralog_Filtering ADD COLUMN selection_signature TEXT")

    if manager._table_exists("Paralog_Filtering_Copy"):
        copy_columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Paralog_Filtering_Copy)").fetchall() or [])}
        if "run_id" not in copy_columns:
            cursor.execute("ALTER TABLE Paralog_Filtering_Copy RENAME TO Paralog_Filtering_Copy_Legacy")
            cursor.execute(
                """
                CREATE TABLE Paralog_Filtering_Copy (
                    family_id VARCHAR(20),
                    library_id INT,
                    target_library_id INT,
                    accession VARCHAR(50),
                    run_id TEXT,
                    query_id TEXT,
                    query_header TEXT,
                    query_status INT,
                    clean BOOLEAN,
                    selected_ref_count INT,
                    reused INTEGER DEFAULT 0,
                    reason_code TEXT,
                    selection_signature TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (family_id, target_library_id, accession, run_id, query_id),
                    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (run_id) REFERENCES Paralog_Filtering_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT family_id, library_id, target_library_id, accession, query_id,
                       query_header, query_status, clean, date
                FROM Paralog_Filtering_Copy_Legacy
                """
            ).fetchall() or []
            run_pairs = sorted({(int(row[2]), int(row[1])) for row in rows if row[2] is not None and row[1] is not None})
            for target_library_id, busco_library_id in run_pairs:
                run_id = _legacy_paralog_run_id(target_library_id, busco_library_id)
                cursor.execute(
                    """
                    INSERT OR IGNORE INTO Paralog_Filtering_Runs (
                        run_id, target_library_id, busco_library_id, selection_mode,
                        selection_params_json, config_signature, run_label, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, datetime('now'))
                    """,
                    (
                        run_id,
                        target_library_id,
                        busco_library_id,
                        "legacy",
                        json.dumps({}),
                        run_id,
                        "legacy-migrated",
                    ),
                )
            for family_id, library_id, target_library_id, accession, query_id, query_header, query_status, clean, date in rows:
                run_id = _legacy_paralog_run_id(target_library_id, library_id)
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Paralog_Filtering_Copy (
                        family_id, library_id, target_library_id, accession, run_id, query_id,
                        query_header, query_status, clean, selected_ref_count, reused,
                        reason_code, selection_signature, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        family_id,
                        library_id,
                        target_library_id,
                        accession,
                        run_id,
                        query_id,
                        query_header,
                        query_status,
                        clean,
                        None,
                        0,
                        "legacy",
                        run_id,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Paralog_Filtering_Copy_Legacy")
        else:
            if "selected_ref_count" not in copy_columns:
                cursor.execute("ALTER TABLE Paralog_Filtering_Copy ADD COLUMN selected_ref_count INT")
            if "reused" not in copy_columns:
                cursor.execute("ALTER TABLE Paralog_Filtering_Copy ADD COLUMN reused INTEGER DEFAULT 0")
            if "reason_code" not in copy_columns:
                cursor.execute("ALTER TABLE Paralog_Filtering_Copy ADD COLUMN reason_code TEXT")
            if "selection_signature" not in copy_columns:
                cursor.execute("ALTER TABLE Paralog_Filtering_Copy ADD COLUMN selection_signature TEXT")


def ensure_selector_preset_schema(manager) -> None:
    manager.cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS Selector_Presets (
            preset_name TEXT PRIMARY KEY,
            selector_json TEXT NOT NULL,
            description TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            updated_at DATETIME DEFAULT (datetime('now'))
        )
        """
    )


def _best_busco_run_for_families(manager, *, accession, busco_library_id, family_ids):
    family_vals = sorted({str(fam) for fam in (family_ids or []) if fam is not None})
    if not family_vals:
        return None
    placeholders = ",".join("?" for _ in family_vals)
    rows = manager.cursor.execute(
        f"""
        SELECT d.run_id, COUNT(DISTINCT d.family_id) AS overlap_cnt
        FROM BUSCO_Run_Family_Data d
        JOIN BUSCO_Runs r ON r.run_id = d.run_id
        WHERE d.accession = ?
          AND d.library_id = ?
          AND d.family_id IN ({placeholders})
          AND COALESCE(r.status, 'completed') = 'completed'
        GROUP BY d.run_id
        ORDER BY overlap_cnt DESC, d.run_id DESC
        """,
        tuple([str(accession), int(busco_library_id), *family_vals]),
    ).fetchall() or []
    if not rows:
        return None
    best_run_id, best_overlap = rows[0]
    if best_overlap is None or int(best_overlap) <= 0:
        return None
    if len(rows) > 1 and int(rows[1][1] or 0) == int(best_overlap):
        return None
    return int(best_run_id)


def _resolve_summary_busco_run_id(manager, *, accession, target_library_id, busco_library_id, run_id):
    vote_rows = manager.cursor.execute(
        """
        SELECT busco_run_id, COUNT(*) AS vote_cnt
        FROM Decontamination_Busco_Votes
        WHERE accession = ?
          AND target_library_id = ?
          AND busco_library_id = ?
          AND run_id = ?
          AND busco_run_id IS NOT NULL
        GROUP BY busco_run_id
        ORDER BY vote_cnt DESC, busco_run_id DESC
        """,
        (str(accession), int(target_library_id), int(busco_library_id), str(run_id)),
    ).fetchall() or []
    if not vote_rows:
        return None
    best_run_id, best_count = vote_rows[0]
    if best_count is None or int(best_count) <= 0:
        return None
    if len(vote_rows) > 1 and int(vote_rows[1][1] or 0) == int(best_count):
        return None
    return int(best_run_id)


def _ensure_analysis_busco_run_link_schema(manager) -> None:
    cursor = manager.cursor
    conn = manager.conn
    cache_invalidated = False

    if manager._table_exists("Paralog_Filtering"):
        columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Paralog_Filtering)").fetchall() or [])}
        if "busco_run_id" not in columns:
            cache_invalidated = True
            cursor.execute("ALTER TABLE Paralog_Filtering RENAME TO Paralog_Filtering_PreBuscoLink")
            cursor.execute(
                """
                CREATE TABLE Paralog_Filtering (
                    family_id VARCHAR(20),
                    library_id INT,
                    target_library_id INT,
                    accession VARCHAR(50),
                    run_id TEXT,
                    busco_run_id INT,
                    clean BOOLEAN,
                    selected_ref_count INT,
                    selection_threshold FLOAT,
                    reused INTEGER DEFAULT 0,
                    reason_code TEXT,
                    selection_signature TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id),
                    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (run_id) REFERENCES Paralog_Filtering_Runs(run_id) ON DELETE CASCADE,
                    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT family_id, library_id, target_library_id, accession, run_id, clean,
                       selected_ref_count, selection_threshold, reused, reason_code,
                       selection_signature, date
                FROM Paralog_Filtering_PreBuscoLink
                """
            ).fetchall() or []
            run_family_cache = {}
            for row in rows:
                family_id, library_id, target_library_id, accession, run_id, clean, selected_ref_count, selection_threshold, reused, reason_code, selection_signature, date = row
                key = (str(accession), int(target_library_id), int(library_id), str(run_id))
                if key not in run_family_cache:
                    fam_rows = cursor.execute(
                        """
                        SELECT family_id
                        FROM Paralog_Filtering_PreBuscoLink
                        WHERE accession = ? AND target_library_id = ? AND library_id = ? AND run_id = ?
                        """,
                        key,
                    ).fetchall() or []
                    run_family_cache[key] = [str(item[0]) for item in fam_rows if item and item[0] is not None]
                busco_run_id = _best_busco_run_for_families(
                    manager,
                    accession=accession,
                    busco_library_id=library_id,
                    family_ids=run_family_cache[key],
                )
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Paralog_Filtering (
                        family_id, library_id, target_library_id, accession, run_id, busco_run_id, clean,
                        selected_ref_count, selection_threshold, reused, reason_code, selection_signature, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        family_id,
                        library_id,
                        target_library_id,
                        accession,
                        run_id,
                        busco_run_id,
                        clean,
                        selected_ref_count,
                        selection_threshold,
                        reused,
                        reason_code,
                        selection_signature,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Paralog_Filtering_PreBuscoLink")

    if manager._table_exists("Paralog_Filtering_Copy"):
        copy_columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Paralog_Filtering_Copy)").fetchall() or [])}
        if "busco_run_id" not in copy_columns:
            cache_invalidated = True
            cursor.execute("ALTER TABLE Paralog_Filtering_Copy RENAME TO Paralog_Filtering_Copy_PreBuscoLink")
            cursor.execute(
                """
                CREATE TABLE Paralog_Filtering_Copy (
                    family_id VARCHAR(20),
                    library_id INT,
                    target_library_id INT,
                    accession VARCHAR(50),
                    run_id TEXT,
                    busco_run_id INT,
                    query_id TEXT,
                    query_header TEXT,
                    query_status INT,
                    clean BOOLEAN,
                    selected_ref_count INT,
                    reused INTEGER DEFAULT 0,
                    reason_code TEXT,
                    selection_signature TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id, query_id),
                    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (run_id) REFERENCES Paralog_Filtering_Runs(run_id) ON DELETE CASCADE,
                    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT family_id, library_id, target_library_id, accession, run_id, query_id,
                       query_header, query_status, clean, selected_ref_count, reused,
                       reason_code, selection_signature, date
                FROM Paralog_Filtering_Copy_PreBuscoLink
                """
            ).fetchall() or []
            for row in rows:
                family_id, library_id, target_library_id, accession, run_id, query_id, query_header, query_status, clean, selected_ref_count, reused, reason_code, selection_signature, date = row
                linked = cursor.execute(
                    """
                    SELECT busco_run_id
                    FROM Paralog_Filtering
                    WHERE family_id = ? AND library_id = ? AND target_library_id = ? AND accession = ? AND run_id = ?
                    ORDER BY rowid DESC
                    LIMIT 1
                    """,
                    (family_id, library_id, target_library_id, accession, run_id),
                ).fetchone()
                busco_run_id = int(linked[0]) if linked and linked[0] is not None else None
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Paralog_Filtering_Copy (
                        family_id, library_id, target_library_id, accession, run_id, busco_run_id,
                        query_id, query_header, query_status, clean, selected_ref_count, reused,
                        reason_code, selection_signature, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        family_id,
                        library_id,
                        target_library_id,
                        accession,
                        run_id,
                        busco_run_id,
                        query_id,
                        query_header,
                        query_status,
                        clean,
                        selected_ref_count,
                        reused,
                        reason_code,
                        selection_signature,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Paralog_Filtering_Copy_PreBuscoLink")

    if manager._table_exists("Decontamination_Busco_Votes"):
        columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Decontamination_Busco_Votes)").fetchall() or [])}
        if "busco_run_id" not in columns:
            cache_invalidated = True
            cursor.execute("ALTER TABLE Decontamination_Busco_Votes RENAME TO Decontamination_Busco_Votes_PreBuscoLink")
            cursor.execute(
                """
                CREATE TABLE Decontamination_Busco_Votes (
                    family_id VARCHAR(20),
                    busco_library_id INT,
                    target_library_id INT,
                    accession VARCHAR(50),
                    run_id TEXT,
                    busco_run_id INT,
                    expected_taxid INT,
                    best_taxid INT,
                    runner_taxid INT,
                    rank VARCHAR(50),
                    best_bitscore FLOAT,
                    delta_bitscore FLOAT,
                    decision VARCHAR(20),
                    top_hits_json TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id),
                    FOREIGN KEY (family_id, busco_library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
                    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (best_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (runner_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT family_id, busco_library_id, target_library_id, accession, run_id,
                       expected_taxid, best_taxid, runner_taxid, rank, best_bitscore,
                       delta_bitscore, decision, top_hits_json, date
                FROM Decontamination_Busco_Votes_PreBuscoLink
                """
            ).fetchall() or []
            run_family_cache = {}
            for row in rows:
                family_id, busco_library_id, target_library_id, accession, run_id, expected_taxid, best_taxid, runner_taxid, rank, best_bitscore, delta_bitscore, decision, top_hits_json, date = row
                key = (str(accession), int(target_library_id), int(busco_library_id), str(run_id))
                if key not in run_family_cache:
                    fam_rows = cursor.execute(
                        """
                        SELECT family_id
                        FROM Decontamination_Busco_Votes_PreBuscoLink
                        WHERE accession = ? AND target_library_id = ? AND busco_library_id = ? AND run_id = ?
                        """,
                        key,
                    ).fetchall() or []
                    run_family_cache[key] = [str(item[0]) for item in fam_rows if item and item[0] is not None]
                busco_run_id = _best_busco_run_for_families(
                    manager,
                    accession=accession,
                    busco_library_id=busco_library_id,
                    family_ids=run_family_cache[key],
                )
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Decontamination_Busco_Votes (
                        family_id, busco_library_id, target_library_id, accession, run_id, busco_run_id,
                        expected_taxid, best_taxid, runner_taxid, rank, best_bitscore, delta_bitscore,
                        decision, top_hits_json, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        family_id,
                        busco_library_id,
                        target_library_id,
                        accession,
                        run_id,
                        busco_run_id,
                        expected_taxid,
                        best_taxid,
                        runner_taxid,
                        rank,
                        best_bitscore,
                        delta_bitscore,
                        decision,
                        top_hits_json,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Decontamination_Busco_Votes_PreBuscoLink")

    if manager._table_exists("Decontamination_Busco_Copy_Votes"):
        copy_columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Decontamination_Busco_Copy_Votes)").fetchall() or [])}
        if "busco_run_id" not in copy_columns:
            cache_invalidated = True
            cursor.execute("ALTER TABLE Decontamination_Busco_Copy_Votes RENAME TO Decontamination_Busco_Copy_Votes_PreBuscoLink")
            cursor.execute(
                """
                CREATE TABLE Decontamination_Busco_Copy_Votes (
                    family_id VARCHAR(20),
                    busco_library_id INT,
                    target_library_id INT,
                    accession VARCHAR(50),
                    run_id TEXT,
                    busco_run_id INT,
                    query_id TEXT,
                    query_header TEXT,
                    query_status INT,
                    expected_taxid INT,
                    best_taxid INT,
                    runner_taxid INT,
                    rank VARCHAR(50),
                    best_bitscore FLOAT,
                    delta_bitscore FLOAT,
                    decision VARCHAR(20),
                    top_hits_json TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id, query_id),
                    FOREIGN KEY (family_id, busco_library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
                    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (best_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (runner_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT family_id, busco_library_id, target_library_id, accession, run_id, query_id,
                       query_header, query_status, expected_taxid, best_taxid, runner_taxid, rank,
                       best_bitscore, delta_bitscore, decision, top_hits_json, date
                FROM Decontamination_Busco_Copy_Votes_PreBuscoLink
                """
            ).fetchall() or []
            for row in rows:
                family_id, busco_library_id, target_library_id, accession, run_id, query_id, query_header, query_status, expected_taxid, best_taxid, runner_taxid, rank, best_bitscore, delta_bitscore, decision, top_hits_json, date = row
                linked = cursor.execute(
                    """
                    SELECT busco_run_id
                    FROM Decontamination_Busco_Votes
                    WHERE family_id = ? AND busco_library_id = ? AND target_library_id = ? AND accession = ? AND run_id = ?
                    ORDER BY rowid DESC
                    LIMIT 1
                    """,
                    (family_id, busco_library_id, target_library_id, accession, run_id),
                ).fetchone()
                busco_run_id = int(linked[0]) if linked and linked[0] is not None else None
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Decontamination_Busco_Copy_Votes (
                        family_id, busco_library_id, target_library_id, accession, run_id, busco_run_id,
                        query_id, query_header, query_status, expected_taxid, best_taxid, runner_taxid,
                        rank, best_bitscore, delta_bitscore, decision, top_hits_json, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        family_id,
                        busco_library_id,
                        target_library_id,
                        accession,
                        run_id,
                        busco_run_id,
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
                        top_hits_json,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Decontamination_Busco_Copy_Votes_PreBuscoLink")

    if manager._table_exists("Decontamination_Summary"):
        columns = {row[1] for row in (cursor.execute("PRAGMA table_info(Decontamination_Summary)").fetchall() or [])}
        if "busco_run_id" not in columns:
            cache_invalidated = True
            cursor.execute("ALTER TABLE Decontamination_Summary RENAME TO Decontamination_Summary_PreBuscoLink")
            cursor.execute(
                """
                CREATE TABLE Decontamination_Summary (
                    accession VARCHAR(50),
                    target_library_id INT,
                    busco_library_id INT,
                    run_id TEXT,
                    busco_run_id INT,
                    expected_taxid INT,
                    majority_taxid INT,
                    rank VARCHAR(50),
                    buscos_tested INT,
                    buscos_supporting INT,
                    buscos_outside INT,
                    off_clade_fraction FLOAT,
                    decision VARCHAR(20),
                    params_json TEXT,
                    date DATETIME DEFAULT (datetime('now')),
                    PRIMARY KEY (accession, target_library_id, busco_library_id, run_id, busco_run_id),
                    FOREIGN KEY (accession) REFERENCES Genome(accession),
                    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
                    FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (majority_taxid) REFERENCES Taxonomy(taxid),
                    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
                )
                """
            )
            rows = cursor.execute(
                """
                SELECT accession, target_library_id, busco_library_id, run_id, expected_taxid, majority_taxid,
                       rank, buscos_tested, buscos_supporting, buscos_outside, off_clade_fraction,
                       decision, params_json, date
                FROM Decontamination_Summary_PreBuscoLink
                """
            ).fetchall() or []
            for row in rows:
                accession, target_library_id, busco_library_id, run_id, expected_taxid, majority_taxid, rank, buscos_tested, buscos_supporting, buscos_outside, off_clade_fraction, decision, params_json, date = row
                busco_run_id = _resolve_summary_busco_run_id(
                    manager,
                    accession=accession,
                    target_library_id=target_library_id,
                    busco_library_id=busco_library_id,
                    run_id=run_id,
                )
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO Decontamination_Summary (
                        accession, target_library_id, busco_library_id, run_id, busco_run_id,
                        expected_taxid, majority_taxid, rank, buscos_tested, buscos_supporting,
                        buscos_outside, off_clade_fraction, decision, params_json, date
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        accession,
                        target_library_id,
                        busco_library_id,
                        run_id,
                        busco_run_id,
                        expected_taxid,
                        majority_taxid,
                        rank,
                        buscos_tested,
                        buscos_supporting,
                        buscos_outside,
                        off_clade_fraction,
                        decision,
                        params_json,
                        date,
                    ),
                )
            cursor.execute("DROP TABLE IF EXISTS Decontamination_Summary_PreBuscoLink")

    if cache_invalidated and manager._table_exists("BUSCO_Adjusted_Results"):
        cursor.execute("DELETE FROM BUSCO_Adjusted_Results")


def ensure_task_queue_schema(manager) -> None:
    if not manager._table_exists("Tasks"):
        return
    cursor = manager.cursor
    conn = manager.conn
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS TaskBlocks (
            block_id INTEGER PRIMARY KEY AUTOINCREMENT,
            task_id INTEGER NOT NULL,
            block_type TEXT NOT NULL,
            condition TEXT NOT NULL,
            message TEXT,
            block_set TEXT,
            block_group TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            satisfied_at DATETIME,
            last_checked_at DATETIME,
            FOREIGN KEY (task_id) REFERENCES Tasks(task_id)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS TaskDependencies (
            dependency_id INTEGER PRIMARY KEY AUTOINCREMENT,
            task_id INTEGER NOT NULL,
            depends_on_task_id INTEGER NOT NULL,
            required_state TEXT NOT NULL,
            block_id INTEGER,
            allow_failed INTEGER DEFAULT 0,
            created_at DATETIME DEFAULT (datetime('now')),
            satisfied_at DATETIME,
            FOREIGN KEY (task_id) REFERENCES Tasks(task_id),
            FOREIGN KEY (depends_on_task_id) REFERENCES Tasks(task_id),
            FOREIGN KEY (block_id) REFERENCES TaskBlocks(block_id)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS TaskTimeConstraints (
            constraint_id INTEGER PRIMARY KEY AUTOINCREMENT,
            task_id INTEGER NOT NULL,
            mode TEXT NOT NULL,
            not_before TEXT NOT NULL,
            block_id INTEGER,
            created_at DATETIME DEFAULT (datetime('now')),
            satisfied_at DATETIME,
            FOREIGN KEY (task_id) REFERENCES Tasks(task_id),
            FOREIGN KEY (block_id) REFERENCES TaskBlocks(block_id)
        )
        """
    )
    if not manager._column_exists("Tasks", "is_barrier"):
        cursor.execute("ALTER TABLE Tasks ADD COLUMN is_barrier INTEGER DEFAULT 0")
    if not manager._column_exists("Tasks", "status_updated_at"):
        cursor.execute("ALTER TABLE Tasks ADD COLUMN status_updated_at DATETIME")
        cursor.execute(
            """
            UPDATE Tasks
            SET status_updated_at = COALESCE(end_time, start_time, queue_time, datetime('now'))
            WHERE status_updated_at IS NULL
            """
        )
    if manager._table_exists("TaskBlocks"):
        if not manager._column_exists("TaskBlocks", "block_set"):
            cursor.execute("ALTER TABLE TaskBlocks ADD COLUMN block_set TEXT")
        if not manager._column_exists("TaskBlocks", "block_group"):
            cursor.execute("ALTER TABLE TaskBlocks ADD COLUMN block_group TEXT")
    if manager._table_exists("Environment_Variables") and not manager._column_exists("Environment_Variables", "updated_at"):
        cursor.execute("ALTER TABLE Environment_Variables ADD COLUMN updated_at DATETIME DEFAULT (datetime('now'))")
        manager._env_updated_at = True
    if manager._table_exists("Assembly") and not manager._column_exists("Assembly", "assembly_status"):
        cursor.execute("ALTER TABLE Assembly ADD COLUMN assembly_status TEXT")
    if manager._table_exists("Assembly") and not manager._column_exists("Assembly", "origin"):
        cursor.execute("ALTER TABLE Assembly ADD COLUMN origin TEXT")
    if manager._table_exists("Genome") and not manager._column_exists("Genome", "assembly_properties"):
        cursor.execute("ALTER TABLE Genome ADD COLUMN assembly_properties TEXT")
    if manager._table_exists("Genome") and not manager._column_exists("Genome", "isoforms_cleaned"):
        cursor.execute("ALTER TABLE Genome ADD COLUMN isoforms_cleaned INTEGER DEFAULT 0")
    if manager._table_exists("Libraries") and not manager._column_exists("Libraries", "status"):
        cursor.execute("ALTER TABLE Libraries ADD COLUMN status TEXT DEFAULT 'ready'")
    if manager._table_exists("OrthoFinder_Results") and not manager._column_exists("OrthoFinder_Results", "status"):
        cursor.execute("ALTER TABLE OrthoFinder_Results ADD COLUMN status TEXT DEFAULT 'ready'")
    if manager._table_exists("OrthoFinder_Results") and not manager._column_exists("OrthoFinder_Results", "mcl_inflation"):
        cursor.execute("ALTER TABLE OrthoFinder_Results ADD COLUMN mcl_inflation REAL")
    if manager._table_exists("OrthoFinder_Results") and not manager._column_exists("OrthoFinder_Results", "command_line"):
        cursor.execute("ALTER TABLE OrthoFinder_Results ADD COLUMN command_line TEXT")
    if manager._table_exists("Libraries"):
        cursor.execute("DROP VIEW IF EXISTS Libraries_View")
        cursor.execute(
            """
            CREATE VIEW Libraries_View AS
            SELECT l.library_id, l.library_name, p.library_name AS parent_name, COALESCE(l.status, 'ready') AS status,
            (
            SELECT GROUP_CONCAT(ra2.accession)
            FROM Reference_Assemblies ra2
            WHERE ra2.library_id = l.library_id
            ) AS accessions
            FROM Libraries l
            LEFT JOIN Libraries p ON p.library_id = l.parent_id
            ORDER BY l.library_id ASC
            """
        )
    if manager._table_exists("OrthoFinder_Results"):
        cursor.execute("DROP VIEW IF EXISTS OrthoFinder_Results_View")
        cursor.execute(
            """
            CREATE VIEW OrthoFinder_Results_View AS
            SELECT ofr.orthofinder_id, ofr.library_id, ofr.location, COALESCE(ofr.status, 'ready') AS status, ofr.date, ofr.mcl_inflation, ofr.command_line, GROUP_CONCAT(oa.accession) AS accessions
            FROM OrthoFinder_Results ofr
            LEFT JOIN OrthoFinder_Accessions oa ON ofr.orthofinder_id = oa.orthofinder_id
            GROUP BY ofr.orthofinder_id
            """
        )
    if not manager._table_exists("Hidden_Genomes"):
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Hidden_Genomes (
                accession TEXT PRIMARY KEY,
                status TEXT,
                reason TEXT,
                created_at DATETIME DEFAULT (datetime('now'))
            )
            """
        )
    _ensure_paralog_filtering_run_schema(manager)
    _ensure_analysis_busco_run_link_schema(manager)
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS Decontamination_Busco_Copy_Votes (
            family_id VARCHAR(20),
            busco_library_id INT,
            target_library_id INT,
            accession VARCHAR(50),
            run_id TEXT,
            busco_run_id INT,
            query_id TEXT,
            query_header TEXT,
            query_status INT,
            expected_taxid INT,
            best_taxid INT,
            runner_taxid INT,
            rank VARCHAR(50),
            best_bitscore FLOAT,
            delta_bitscore FLOAT,
            decision VARCHAR(20),
            top_hits_json TEXT,
            date DATETIME DEFAULT (datetime('now')),
            PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id, query_id),
            FOREIGN KEY (family_id, busco_library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
            FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
            FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
            FOREIGN KEY (best_taxid) REFERENCES Taxonomy(taxid),
            FOREIGN KEY (runner_taxid) REFERENCES Taxonomy(taxid),
            FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
        )
        """
    )
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taskblocks_task ON TaskBlocks(task_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taskdeps_task ON TaskDependencies(task_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_taskdeps_dep ON TaskDependencies(depends_on_task_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_tasktime_task ON TaskTimeConstraints(task_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_runs_target_date ON Paralog_Filtering_Runs(target_library_id, date)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_target_run_acc ON Paralog_Filtering(target_library_id, run_id, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_target_busco_run_acc ON Paralog_Filtering(target_library_id, accession, busco_run_id, run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_busco_run ON Paralog_Filtering(busco_run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_signature ON Paralog_Filtering(target_library_id, selection_signature, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_copy_target_run_acc ON Paralog_Filtering_Copy(target_library_id, run_id, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_paralog_copy_target_busco_run_acc ON Paralog_Filtering_Copy(target_library_id, accession, busco_run_id, run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_decont_copy_target_run_acc ON Decontamination_Busco_Copy_Votes(target_library_id, run_id, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_decont_votes_target_busco_run_acc ON Decontamination_Busco_Votes(target_library_id, accession, busco_run_id, run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_decont_votes_busco_run ON Decontamination_Busco_Votes(busco_run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_decont_summary_target_busco_run_acc ON Decontamination_Summary(target_library_id, accession, busco_run_id, run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_decont_summary_busco_run ON Decontamination_Summary(busco_run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_decont_copy_busco_run ON Decontamination_Busco_Copy_Votes(busco_run_id)")


def ensure_busco_run_schema(manager) -> None:
    if not manager._table_exists("Genome") or not manager._table_exists("Libraries"):
        return
    cursor = manager.cursor
    conn = manager.conn
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS BUSCO_Runs (
            run_id INTEGER PRIMARY KEY AUTOINCREMENT,
            accession VARCHAR(50) NOT NULL,
            library_id INT NOT NULL,
            lineage_name TEXT,
            input_mode TEXT NOT NULL,
            pipeline TEXT NOT NULL,
            pipeline_params_effective_json TEXT,
            pipeline_params_source_json TEXT,
            busco_cli_args_json TEXT,
            busco_version TEXT,
            result_dir TEXT,
            status TEXT DEFAULT 'pending',
            no_sc_complete INT,
            no_duplicated_complete INT,
            no_fragmented INT,
            no_missing INT,
            started_at DATETIME,
            completed_at DATETIME,
            created_at DATETIME DEFAULT (datetime('now')),
            updated_at DATETIME DEFAULT (datetime('now')),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS BUSCO_Run_Family_Data (
            run_id INTEGER,
            family_id VARCHAR(20),
            library_id INT,
            accession VARCHAR(50),
            status INT,
            sequence TEXT,
            score FLOAT,
            length INT,
            PRIMARY KEY (run_id, family_id, accession, sequence),
            FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE,
            FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS BUSCO_Run_Family_Locations (
            run_id INTEGER,
            family_id VARCHAR(20),
            library_id INT,
            accession VARCHAR(50),
            location TEXT,
            PRIMARY KEY (run_id, family_id, accession),
            FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE,
            FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
            FOREIGN KEY (accession) REFERENCES Genome(accession)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS BUSCO_Primary (
            accession VARCHAR(50) NOT NULL,
            library_id INT NOT NULL,
            purpose TEXT NOT NULL,
            run_id INTEGER NOT NULL,
            policy TEXT,
            updated_at DATETIME DEFAULT (datetime('now')),
            updated_by TEXT,
            PRIMARY KEY (accession, library_id, purpose),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (library_id) REFERENCES Libraries(library_id),
            FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS BUSCO_Adjusted_Results (
            cache_key TEXT PRIMARY KEY,
            library_id INT NOT NULL,
            accession VARCHAR(50) NOT NULL,
            species TEXT,
            library_name TEXT,
            effective_busco_run_id INTEGER,
            effective_decont_run_id TEXT,
            effective_decont_library_id INT,
            effective_decont_decision TEXT,
            include_paralog INTEGER NOT NULL DEFAULT 1,
            include_decontam INTEGER NOT NULL DEFAULT 1,
            allow_ambiguous_contaminants INTEGER NOT NULL DEFAULT 0,
            strict_decontamination INTEGER NOT NULL DEFAULT 0,
            rescue_duplicates INTEGER NOT NULL DEFAULT 0,
            complete FLOAT,
            single_copy_complete FLOAT,
            duplicated FLOAT,
            fragmented FLOAT,
            missing FLOAT,
            hidden_paralog FLOAT,
            contaminated FLOAT,
            has_paralog INTEGER NOT NULL DEFAULT 0,
            has_decont INTEGER NOT NULL DEFAULT 0,
            status TEXT NOT NULL DEFAULT 'ready',
            updated_at DATETIME DEFAULT (datetime('now')),
            invalidated_at DATETIME,
            invalidation_reason TEXT,
            FOREIGN KEY (library_id) REFERENCES Libraries(library_id),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (effective_busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE SET NULL
        )
        """
    )
    if manager._table_exists("BUSCO_Adjusted_Results") and not manager._column_exists("BUSCO_Adjusted_Results", "effective_decont_decision"):
        cursor.execute("ALTER TABLE BUSCO_Adjusted_Results ADD COLUMN effective_decont_decision TEXT")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_runs_acc_lib ON BUSCO_Runs(accession, library_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_runs_lib_status ON BUSCO_Runs(library_id, status)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_runs_pipeline_mode ON BUSCO_Runs(pipeline, input_mode)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_run_fam_data_acc ON BUSCO_Run_Family_Data(library_id, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_run_fam_loc_acc ON BUSCO_Run_Family_Locations(library_id, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_primary_run ON BUSCO_Primary(run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_adjusted_results_acc ON BUSCO_Adjusted_Results(library_id, accession, status)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_adjusted_results_run ON BUSCO_Adjusted_Results(effective_busco_run_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_adjusted_results_decont ON BUSCO_Adjusted_Results(effective_decont_run_id)")


def ensure_storage_schema(manager) -> None:
    if not manager._table_exists("Genome") or not manager._table_exists("Libraries"):
        return
    cursor = manager.cursor
    conn = manager.conn
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS StorageRoots (
            storage_root_id INTEGER PRIMARY KEY AUTOINCREMENT,
            logical_kind TEXT NOT NULL,
            label TEXT,
            base_path TEXT NOT NULL,
            writable INTEGER DEFAULT 1,
            is_active INTEGER DEFAULT 1,
            metadata_json TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            updated_at DATETIME DEFAULT (datetime('now')),
            UNIQUE(logical_kind, base_path)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS Artifacts (
            artifact_id INTEGER PRIMARY KEY AUTOINCREMENT,
            owner_type TEXT NOT NULL,
            owner_id TEXT NOT NULL,
            artifact_type TEXT NOT NULL,
            role TEXT,
            status TEXT DEFAULT 'ready',
            storage_root_id INTEGER,
            relative_path TEXT,
            absolute_path TEXT,
            is_dir INTEGER DEFAULT 0,
            format TEXT,
            sequence_kind TEXT,
            checksum TEXT,
            size_bytes INTEGER,
            metadata_json TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            updated_at DATETIME DEFAULT (datetime('now')),
            FOREIGN KEY (storage_root_id) REFERENCES StorageRoots(storage_root_id),
            UNIQUE(owner_type, owner_id, artifact_type, role, sequence_kind, relative_path)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS FilesystemOperations (
            operation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            operation_type TEXT NOT NULL,
            source_path TEXT,
            staging_path TEXT,
            destination_path TEXT,
            status TEXT NOT NULL DEFAULT 'preparing',
            payload_json TEXT,
            error_message TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            updated_at DATETIME DEFAULT (datetime('now'))
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS BUSCO_Run_Family_Artifacts (
            run_id INTEGER,
            family_id VARCHAR(20),
            library_id INT,
            accession VARCHAR(50),
            artifact_id INTEGER,
            sequence_kind TEXT,
            location TEXT,
            metadata_json TEXT,
            PRIMARY KEY (run_id, family_id, accession, sequence_kind),
            FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE,
            FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (artifact_id) REFERENCES Artifacts(artifact_id) ON DELETE SET NULL
        )
        """
    )
    if not manager._column_exists("Genome", "storage_root_id"):
        cursor.execute("ALTER TABLE Genome ADD COLUMN storage_root_id INTEGER")
    if not manager._column_exists("Genome", "relative_path"):
        cursor.execute("ALTER TABLE Genome ADD COLUMN relative_path TEXT")
    if not manager._column_exists("Libraries", "storage_root_id"):
        cursor.execute("ALTER TABLE Libraries ADD COLUMN storage_root_id INTEGER")
    if not manager._column_exists("Libraries", "relative_path"):
        cursor.execute("ALTER TABLE Libraries ADD COLUMN relative_path TEXT")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_storage_roots_kind ON StorageRoots(logical_kind, is_active)")
    cursor.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS idx_storage_roots_label_unique
        ON StorageRoots(label)
        WHERE label IS NOT NULL AND trim(label) <> ''
        """
    )
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_artifacts_owner ON Artifacts(owner_type, owner_id, artifact_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_artifacts_root ON Artifacts(storage_root_id, relative_path)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_run_family_artifacts_run ON BUSCO_Run_Family_Artifacts(run_id, accession)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_filesystem_operations_status ON FilesystemOperations(status, operation_id)")


def ensure_proteome_schema(manager) -> None:
    if not manager._table_exists("Genome"):
        return
    cursor = manager.cursor
    conn = manager.conn
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS Proteome_Profiles (
            proteome_profile_id INTEGER PRIMARY KEY AUTOINCREMENT,
            accession VARCHAR(50) NOT NULL,
            profile_name TEXT NOT NULL,
            kind TEXT NOT NULL,
            parent_profile_id INTEGER,
            artifact_id INTEGER,
            status TEXT DEFAULT 'ready',
            sequence_count INTEGER,
            checksum TEXT,
            is_default INTEGER DEFAULT 0,
            metadata_json TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            updated_at DATETIME DEFAULT (datetime('now')),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (parent_profile_id) REFERENCES Proteome_Profiles(proteome_profile_id) ON DELETE SET NULL,
            FOREIGN KEY (artifact_id) REFERENCES Artifacts(artifact_id) ON DELETE SET NULL,
            UNIQUE(accession, profile_name)
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS Proteome_Preparations (
            preparation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            accession VARCHAR(50) NOT NULL,
            input_profile_id INTEGER NOT NULL,
            output_profile_id INTEGER NOT NULL,
            preparation_type TEXT NOT NULL,
            used_gff INTEGER DEFAULT 0,
            gff_artifact_id INTEGER,
            skip_gff INTEGER DEFAULT 0,
            skip_cdhit INTEGER DEFAULT 0,
            gff_priority INTEGER DEFAULT 0,
            cdhit_identity FLOAT,
            cdhit_threads INTEGER,
            input_count INTEGER,
            output_count INTEGER,
            gff_removed INTEGER,
            cdhit_removed INTEGER,
            total_removed INTEGER,
            status TEXT DEFAULT 'completed',
            params_json TEXT,
            created_at DATETIME DEFAULT (datetime('now')),
            FOREIGN KEY (accession) REFERENCES Genome(accession),
            FOREIGN KEY (input_profile_id) REFERENCES Proteome_Profiles(proteome_profile_id) ON DELETE CASCADE,
            FOREIGN KEY (output_profile_id) REFERENCES Proteome_Profiles(proteome_profile_id) ON DELETE CASCADE,
            FOREIGN KEY (gff_artifact_id) REFERENCES Artifacts(artifact_id) ON DELETE SET NULL
        )
        """
    )
    if manager._table_exists("BUSCO_Runs") and not manager._column_exists("BUSCO_Runs", "proteome_profile_id"):
        cursor.execute("ALTER TABLE BUSCO_Runs ADD COLUMN proteome_profile_id INTEGER")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_proteome_profiles_accession ON Proteome_Profiles(accession, profile_name)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_proteome_profiles_default ON Proteome_Profiles(accession, is_default)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_proteome_preparations_output ON Proteome_Preparations(output_profile_id, preparation_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_busco_runs_profile ON BUSCO_Runs(proteome_profile_id)")


_SETUP_SCHEMA_SQL = """
CREATE TABLE Assembly (
    accession TEXT PRIMARY KEY,
    uid INT,
    assembly_method TEXT,
    assembly_type VARCHAR(50),
    assembly_status TEXT,
    origin TEXT,
    release_date DATE,
    warnings TEXT,
    bioproject_accession VARCHAR(50),
    biosample_accession VARCHAR(50),
    comments TEXT,
    diploid_role VARCHAR(50),
    refseq_category VARCHAR(50),
    sequencing_tech VARCHAR(100),
    submitter TEXT,
    contig_l50 INT,
    contig_n50 INT,
    gc_count BIGINT,
    gc_percent FLOAT,
    genome_coverage VARCHAR(50),
    number_of_component_sequences INT,
    number_of_contigs INT,
    number_of_organelles INT,
    number_of_scaffolds INT,
    scaffold_l50 INT,
    scaffold_n50 INT,
    total_number_of_chromosomes INT,
    total_sequence_length BIGINT,
    total_ungapped_length BIGINT
);

CREATE TABLE Taxonomy (
    taxid INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    rank VARCHAR(50),
    parent_taxid INTEGER
);

CREATE INDEX IF NOT EXISTS idx_taxonomy_parent ON Taxonomy(parent_taxid);
CREATE INDEX IF NOT EXISTS idx_taxonomy_name_rank ON Taxonomy(name, rank);
CREATE INDEX IF NOT EXISTS idx_taxonomy_name_nocase ON Taxonomy(name COLLATE NOCASE);

CREATE TABLE Genome (
    accession VARCHAR(50) PRIMARY KEY,
    taxid INT,
    assembly_level TEXT,
    assembly_properties TEXT,
    assembly_name VARCHAR(100),
    protein BOOLEAN DEFAULT 0,
    isoforms_cleaned BOOLEAN DEFAULT 0,
    comments TEXT,
    dl_date DATETIME,
    location TEXT,
    status INT,
    FOREIGN KEY (accession) REFERENCES Assembly(accession),
    FOREIGN KEY (taxid) REFERENCES Taxonomy(taxid)
);

CREATE TABLE Tasks (
    task_id INTEGER PRIMARY KEY AUTOINCREMENT,
    job_type INT,
    status BOOLEAN,
    priority INT,
    parent_id INT,
    checkpoint INT,
    data TEXT,
    queue_time DATETIME DEFAULT (datetime('now')),
    start_time DATETIME,
    end_time DATETIME,
    error_message TEXT,
    error_stack TEXT,
    is_barrier INTEGER DEFAULT 0,
    status_updated_at DATETIME DEFAULT (datetime('now')),
    FOREIGN KEY (parent_id) REFERENCES Tasks(task_id),
    CHECK (parent_id IS NULL OR parent_id != task_id)
);

CREATE TABLE TaskBlocks (
    block_id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id INTEGER NOT NULL,
    block_type TEXT NOT NULL,
    condition TEXT NOT NULL,
    message TEXT,
    block_set TEXT,
    block_group TEXT,
    created_at DATETIME DEFAULT (datetime('now')),
    satisfied_at DATETIME,
    last_checked_at DATETIME,
    FOREIGN KEY (task_id) REFERENCES Tasks(task_id)
);

CREATE TABLE TaskDependencies (
    dependency_id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id INTEGER NOT NULL,
    depends_on_task_id INTEGER NOT NULL,
    required_state TEXT NOT NULL,
    block_id INTEGER,
    allow_failed INTEGER DEFAULT 0,
    created_at DATETIME DEFAULT (datetime('now')),
    satisfied_at DATETIME,
    FOREIGN KEY (task_id) REFERENCES Tasks(task_id),
    FOREIGN KEY (depends_on_task_id) REFERENCES Tasks(task_id),
    FOREIGN KEY (block_id) REFERENCES TaskBlocks(block_id)
);

CREATE TABLE TaskTimeConstraints (
    constraint_id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id INTEGER NOT NULL,
    mode TEXT NOT NULL,
    not_before TEXT NOT NULL,
    block_id INTEGER,
    created_at DATETIME DEFAULT (datetime('now')),
    satisfied_at DATETIME,
    FOREIGN KEY (task_id) REFERENCES Tasks(task_id),
    FOREIGN KEY (block_id) REFERENCES TaskBlocks(block_id)
);

CREATE INDEX IF NOT EXISTS idx_taskblocks_task ON TaskBlocks(task_id);
CREATE INDEX IF NOT EXISTS idx_taskdeps_task ON TaskDependencies(task_id);
CREATE INDEX IF NOT EXISTS idx_taskdeps_dep ON TaskDependencies(depends_on_task_id);
CREATE INDEX IF NOT EXISTS idx_tasktime_task ON TaskTimeConstraints(task_id);

CREATE TABLE Download_parameters (
    taxid INT,
    member_clade TEXT,
    member_clade_rank VARCHAR(50),
    partition_rank VARCHAR(50),
    frequency INT,
    PRIMARY KEY (taxid)
);

CREATE TABLE Hidden_Genomes (
    accession TEXT PRIMARY KEY,
    status TEXT,
    reason TEXT,
    created_at DATETIME DEFAULT (datetime('now'))
);

CREATE VIEW TaxonomyAssemblySummary AS
SELECT t.taxid, t.name, g.accession AS assembly_accession, g.assembly_level, a.release_date, g.dl_date, g.location, g.status, g.protein
FROM Taxonomy t
JOIN Genome g ON t.taxid = g.taxid
JOIN Assembly a ON g.accession = a.accession;

CREATE VIEW Genome_quick_view AS
SELECT g.accession, t.taxid, t.name, g.assembly_level, a.release_date, g.dl_date, g.location
FROM Taxonomy t
JOIN Genome g ON t.taxid = g.taxid
JOIN Assembly a ON g.accession = a.accession
WHERE g.status > 0
ORDER BY g.dl_date DESC;

CREATE TABLE Search_parameters (
    taxid INT,
    search VARCHAR(50),
    taxonomic_rank VARCHAR(50),
    frequency INT,
    PRIMARY KEY (taxid)
);

CREATE TABLE Reference_Assemblies (
    accession VARCHAR(50),
    library_id INT,
    PRIMARY KEY (accession, library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE Proteome_BlastDBs (
    blastdb_id INTEGER PRIMARY KEY AUTOINCREMENT,
    library_id INT,
    accession VARCHAR(50),
    location TEXT,
    date DATETIME DEFAULT (datetime('now')),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession)
);

CREATE TABLE Libraries (
    library_id INTEGER PRIMARY KEY AUTOINCREMENT,
    library_name VARCHAR(100),
    odb_version VARCHAR(20),
    taxid INT,
    size INT,
    location TEXT,
    parent_id INT,
    status TEXT DEFAULT 'ready',
    UNIQUE(library_name),
    FOREIGN KEY (parent_id) REFERENCES Libraries(library_id)
);

CREATE TRIGGER IF NOT EXISTS set_odb_version_from_parent_after_insert
AFTER INSERT ON Libraries
FOR EACH ROW
WHEN NEW.odb_version IS NULL AND NEW.parent_id IS NOT NULL
BEGIN
    UPDATE Libraries
    SET odb_version = (SELECT odb_version FROM Libraries WHERE library_id = NEW.parent_id)
    WHERE library_id = NEW.library_id;
END;

CREATE TRIGGER IF NOT EXISTS set_odb_version_from_parent_after_update
AFTER UPDATE ON Libraries
FOR EACH ROW
WHEN NEW.odb_version IS NULL AND NEW.parent_id IS NOT NULL
BEGIN
    UPDATE Libraries
    SET odb_version = (SELECT odb_version FROM Libraries WHERE library_id = NEW.parent_id)
    WHERE library_id = NEW.library_id;
END;

CREATE VIEW Libraries_View AS
SELECT l.library_id, l.library_name, p.library_name AS parent_name, COALESCE(l.status, 'ready') AS status,
(
SELECT GROUP_CONCAT(ra2.accession)
FROM Reference_Assemblies ra2
WHERE ra2.library_id = l.library_id
) AS accessions
FROM Libraries l
LEFT JOIN Libraries p ON p.library_id = l.parent_id
ORDER BY l.library_id ASC;

CREATE TABLE BUSCO_Results (
    accession VARCHAR(50),
    library_id INT,
    no_sc_complete INT,
    no_duplicated_complete INT,
    no_fragmented INT,
    no_missing INT,
    date DATETIME DEFAULT (datetime('now')),
    PRIMARY KEY (accession, library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE BUSCO_Family_Data (
    family_id VARCHAR(20),
    library_id INT,
    accession VARCHAR(50),
    status INT,
    sequence TEXT,
    score FLOAT,
    length INT,
    PRIMARY KEY (family_id, library_id, accession, sequence),
    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE BUSCO_Family_Locations (
    family_id VARCHAR(20),
    library_id INT,
    accession VARCHAR(50),
    location TEXT,
    PRIMARY KEY (family_id, library_id, accession),
    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession)
);

CREATE TABLE BUSCO_Runs (
    run_id INTEGER PRIMARY KEY AUTOINCREMENT,
    accession VARCHAR(50) NOT NULL,
    library_id INT NOT NULL,
    lineage_name TEXT,
    input_mode TEXT NOT NULL,
    pipeline TEXT NOT NULL,
    pipeline_params_effective_json TEXT,
    pipeline_params_source_json TEXT,
    busco_cli_args_json TEXT,
    busco_version TEXT,
    result_dir TEXT,
    status TEXT DEFAULT 'pending',
    no_sc_complete INT,
    no_duplicated_complete INT,
    no_fragmented INT,
    no_missing INT,
    started_at DATETIME,
    completed_at DATETIME,
    created_at DATETIME DEFAULT (datetime('now')),
    updated_at DATETIME DEFAULT (datetime('now')),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE BUSCO_Run_Family_Data (
    run_id INTEGER,
    family_id VARCHAR(20),
    library_id INT,
    accession VARCHAR(50),
    status INT,
    sequence TEXT,
    score FLOAT,
    length INT,
    PRIMARY KEY (run_id, family_id, accession, sequence),
    FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE,
    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE BUSCO_Run_Family_Locations (
    run_id INTEGER,
    family_id VARCHAR(20),
    library_id INT,
    accession VARCHAR(50),
    location TEXT,
    PRIMARY KEY (run_id, family_id, accession),
    FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE,
    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession)
);

CREATE TABLE BUSCO_Primary (
    accession VARCHAR(50) NOT NULL,
    library_id INT NOT NULL,
    purpose TEXT NOT NULL,
    run_id INTEGER NOT NULL,
    policy TEXT,
    updated_at DATETIME DEFAULT (datetime('now')),
    updated_by TEXT,
    PRIMARY KEY (accession, library_id, purpose),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
);

CREATE TABLE BUSCO_Adjusted_Results (
    cache_key TEXT PRIMARY KEY,
    library_id INT NOT NULL,
    accession VARCHAR(50) NOT NULL,
    species TEXT,
    library_name TEXT,
    effective_busco_run_id INTEGER,
    effective_decont_run_id TEXT,
    effective_decont_library_id INT,
    effective_decont_decision TEXT,
    include_paralog INTEGER NOT NULL DEFAULT 1,
    include_decontam INTEGER NOT NULL DEFAULT 1,
    allow_ambiguous_contaminants INTEGER NOT NULL DEFAULT 0,
    strict_decontamination INTEGER NOT NULL DEFAULT 0,
    rescue_duplicates INTEGER NOT NULL DEFAULT 0,
    complete FLOAT,
    single_copy_complete FLOAT,
    duplicated FLOAT,
    fragmented FLOAT,
    missing FLOAT,
    hidden_paralog FLOAT,
    contaminated FLOAT,
    has_paralog INTEGER NOT NULL DEFAULT 0,
    has_decont INTEGER NOT NULL DEFAULT 0,
    status TEXT NOT NULL DEFAULT 'ready',
    updated_at DATETIME DEFAULT (datetime('now')),
    invalidated_at DATETIME,
    invalidation_reason TEXT,
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (effective_busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE SET NULL
);

CREATE TABLE Paralog_Filtering (
    family_id VARCHAR(20),
    library_id INT,
    target_library_id INT,
    accession VARCHAR(50),
    run_id TEXT,
    busco_run_id INT,
    clean BOOLEAN,
    selected_ref_count INT,
    selection_threshold FLOAT,
    reused INTEGER DEFAULT 0,
    reason_code TEXT,
    selection_signature TEXT,
    date DATETIME DEFAULT (datetime('now')),
    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id),
    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (run_id) REFERENCES Paralog_Filtering_Runs(run_id) ON DELETE CASCADE,
    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
);

CREATE TABLE Paralog_Filtering_Copy (
    family_id VARCHAR(20),
    library_id INT,
    target_library_id INT,
    accession VARCHAR(50),
    run_id TEXT,
    busco_run_id INT,
    query_id TEXT,
    query_header TEXT,
    query_status INT,
    clean BOOLEAN,
    selected_ref_count INT,
    reused INTEGER DEFAULT 0,
    reason_code TEXT,
    selection_signature TEXT,
    date DATETIME DEFAULT (datetime('now')),
    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id, query_id),
    FOREIGN KEY (family_id, library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (run_id) REFERENCES Paralog_Filtering_Runs(run_id) ON DELETE CASCADE,
    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
);

CREATE TABLE Decontamination_Busco_Votes (
    family_id VARCHAR(20),
    busco_library_id INT,
    target_library_id INT,
    accession VARCHAR(50),
    run_id TEXT,
    busco_run_id INT,
    expected_taxid INT,
    best_taxid INT,
    runner_taxid INT,
    rank VARCHAR(50),
    best_bitscore FLOAT,
    delta_bitscore FLOAT,
    decision VARCHAR(20),
    top_hits_json TEXT,
    date DATETIME DEFAULT (datetime('now')),
    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id),
    FOREIGN KEY (family_id, busco_library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (best_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (runner_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
);

CREATE TABLE Decontamination_Busco_Copy_Votes (
    family_id VARCHAR(20),
    busco_library_id INT,
    target_library_id INT,
    accession VARCHAR(50),
    run_id TEXT,
    busco_run_id INT,
    query_id TEXT,
    query_header TEXT,
    query_status INT,
    expected_taxid INT,
    best_taxid INT,
    runner_taxid INT,
    rank VARCHAR(50),
    best_bitscore FLOAT,
    delta_bitscore FLOAT,
    decision VARCHAR(20),
    top_hits_json TEXT,
    date DATETIME DEFAULT (datetime('now')),
    PRIMARY KEY (family_id, target_library_id, accession, run_id, busco_run_id, query_id),
    FOREIGN KEY (family_id, busco_library_id) REFERENCES BUSCO_descriptions(family_id, library_id),
    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (best_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (runner_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
);

CREATE TABLE Decontamination_Summary (
    accession VARCHAR(50),
    target_library_id INT,
    busco_library_id INT,
    run_id TEXT,
    busco_run_id INT,
    expected_taxid INT,
    majority_taxid INT,
    rank VARCHAR(50),
    buscos_tested INT,
    buscos_supporting INT,
    buscos_outside INT,
    off_clade_fraction FLOAT,
    decision VARCHAR(20),
    params_json TEXT,
    date DATETIME DEFAULT (datetime('now')),
    PRIMARY KEY (accession, target_library_id, busco_library_id, run_id, busco_run_id),
    FOREIGN KEY (accession) REFERENCES Genome(accession),
    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (expected_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (majority_taxid) REFERENCES Taxonomy(taxid),
    FOREIGN KEY (busco_run_id) REFERENCES BUSCO_Runs(run_id) ON DELETE CASCADE
);

CREATE TABLE Decontamination_Runs (
    run_id TEXT PRIMARY KEY,
    target_library_id INT,
    busco_library_id INT,
    targets_json TEXT,
    refs_json TEXT,
    params_json TEXT,
    config_signature TEXT,
    run_label TEXT,
    date DATETIME DEFAULT (datetime('now')),
    FOREIGN KEY (target_library_id) REFERENCES Libraries(library_id),
    FOREIGN KEY (busco_library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE BUSCO_descriptions (
    family_id VARCHAR(20),
    library_id INTEGER,
    description TEXT,
    link TEXT,
    PRIMARY KEY (family_id, library_id),
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE OrthoFinder_Results (
    orthofinder_id INTEGER PRIMARY KEY AUTOINCREMENT,
    library_id INT,
    location TEXT,
    status TEXT DEFAULT 'ready',
    date DATETIME DEFAULT (datetime('now')),
    mcl_inflation REAL,
    command_line TEXT,
    FOREIGN KEY (library_id) REFERENCES Libraries(library_id)
);

CREATE TABLE OrthoFinder_Accessions (
    orthofinder_id INT,
    accession VARCHAR(50),
    PRIMARY KEY (orthofinder_id, accession),
    FOREIGN KEY (orthofinder_id) REFERENCES OrthoFinder_Results(orthofinder_id)
);

CREATE VIEW OrthoFinder_Results_View AS
SELECT ofr.orthofinder_id, ofr.library_id, ofr.location, COALESCE(ofr.status, 'ready') AS status, ofr.date, ofr.mcl_inflation, ofr.command_line, GROUP_CONCAT(oa.accession) AS accessions
FROM OrthoFinder_Results ofr
LEFT JOIN OrthoFinder_Accessions oa ON ofr.orthofinder_id = oa.orthofinder_id
GROUP BY ofr.orthofinder_id;

CREATE TABLE Environment_Variables (
    var_name VARCHAR(100) PRIMARY KEY,
    var_value TEXT,
    var_kind TEXT NOT NULL DEFAULT 'env',
    updated_at DATETIME DEFAULT (datetime('now'))
);

CREATE TABLE Selector_Presets (
    preset_name TEXT PRIMARY KEY,
    selector_json TEXT NOT NULL,
    description TEXT,
    created_at DATETIME DEFAULT (datetime('now')),
    updated_at DATETIME DEFAULT (datetime('now'))
);

CREATE VIEW BUSCO_Results_Percentages AS
SELECT br.accession, t.name AS species, l.library_name,
       ROUND(100.0 * (br.no_sc_complete + br.no_duplicated_complete) / l.size, 2) AS complete,
       ROUND(100.0 * br.no_sc_complete / l.size, 2) AS single_copy_complete,
       ROUND(100.0 * br.no_duplicated_complete / l.size, 2) AS duplicated,
       ROUND(100.0 * br.no_fragmented / l.size, 2) AS fragmented,
       ROUND(100.0 * br.no_missing / l.size, 2) AS missing
FROM BUSCO_Results br
JOIN Genome g ON br.accession = g.accession
JOIN Taxonomy t ON g.taxid = t.taxid
JOIN Libraries l ON br.library_id = l.library_id;
"""


def setup_database(manager) -> None:
    manager.cursor.executescript(_SETUP_SCHEMA_SQL)
    ensure_environment_variable_schema(manager)
    ensure_taxonomy_schema(manager)
    ensure_task_queue_schema(manager)
    ensure_busco_run_schema(manager)
    ensure_storage_schema(manager)
    ensure_proteome_schema(manager)
    manager.storage.ensure_default_roots_from_env()
    manager.conn.commit()


def ensure_environment_variable_schema(manager) -> None:
    cursor = manager.cursor
    if not manager._table_exists("Environment_Variables"):
        return
    if not manager._column_exists("Environment_Variables", "updated_at"):
        cursor.execute("ALTER TABLE Environment_Variables ADD COLUMN updated_at DATETIME DEFAULT (datetime('now'))")
        manager._env_updated_at = True
    if not manager._column_exists("Environment_Variables", "var_kind"):
        cursor.execute("ALTER TABLE Environment_Variables ADD COLUMN var_kind TEXT")

    try:
        known_accessions = {
            str(row[0])
            for row in cursor.execute("SELECT accession FROM Genome").fetchall()
            if row and row[0] is not None
        }
    except Exception:
        known_accessions = set()
    rows = cursor.execute(
        """
        SELECT var_name, var_value
        FROM Environment_Variables
        WHERE var_kind IS NULL OR TRIM(var_kind) = ''
        """
    ).fetchall() or []
    for var_name, var_value in rows:
        try:
            decoded = json.loads(var_value)
        except Exception:
            decoded = var_value
        cursor.execute(
            "UPDATE Environment_Variables SET var_kind = ? WHERE var_name = ?",
            (infer_variable_kind(var_name, decoded, known_accessions), var_name),
        )
    cursor.execute(
        """
        UPDATE Environment_Variables
        SET var_kind = 'env'
        WHERE var_kind NOT IN ('env', 'assemblies', 'busco_runs')
        """
    )


def ensure_taxonomy_schema(manager) -> None:
    manager.cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxonomy_parent ON Taxonomy(parent_taxid)")
    manager.cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxonomy_name_rank ON Taxonomy(name, rank)")
    manager.cursor.execute("CREATE INDEX IF NOT EXISTS idx_taxonomy_name_nocase ON Taxonomy(name COLLATE NOCASE)")
