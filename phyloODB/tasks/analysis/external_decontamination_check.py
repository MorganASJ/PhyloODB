import json
import os
import tempfile
import subprocess
import concurrent.futures

from ..task import Task


class ExternalDecontaminationCheckTask(Task):
    """Run external BLAST checks for BUSCOs flagged outside in an internal decontamination run."""

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=8):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.run_id = self.data.get("run_id")
        self.library_id = self.data.get("library_id")
        self.library_name = self.data.get("library_name")
        self.blast_db_path = self.data.get("blast_db_path") or self.data.get("external_blastdb_path")
        self.blast_db_type = self.data.get("blast_db_type") or self.data.get("external_blastdb_type")
        self.blast_program = self.data.get("blast_program") or self.data.get("external_blast_program")
        self.output_dir = self.data.get("output_dir") or self.data.get("external_blast_output_dir")
        self.reuse_blast_results = self.data.get("reuse_blast_results") or self.data.get("external_reuse_blast_results")
        self.max_target_seqs = self.data.get("max_target_seqs") or self.data.get("external_max_target_seqs")
        self.hit_window = self.data.get("hit_window") or self.data.get("external_hit_window")
        self.force = bool(self.data.get("force", False))
        max_concurrent_raw = self.data.get("max_concurrent")
        if max_concurrent_raw in (None, ""):
            self.max_concurrent = 1
        else:
            try:
                self.max_concurrent = max(1, int(max_concurrent_raw))
            except (TypeError, ValueError):
                self.max_concurrent = 1
        threads_raw = self.data.get("threads")
        if threads_raw in (None, ""):
            threads_raw = self.data.get("required_threads")
        if threads_raw in (None, ""):
            self.threads_budget = max(1, int(required_threads))
        else:
            try:
                self.threads_budget = max(1, int(threads_raw))
            except (TypeError, ValueError):
                self.threads_budget = max(1, int(required_threads))
        blast_threads_raw = self.data.get("blast_threads")
        if blast_threads_raw in (None, ""):
            self.blast_threads = None
        else:
            try:
                self.blast_threads = max(1, int(blast_threads_raw))
            except (TypeError, ValueError):
                self.blast_threads = None
        self.blast_path = None
        self.blast_engine = "blast"

    def _detect_db_type(self):
        if self.blast_db_type in ("prot", "nucl"):
            return self.blast_db_type
        if self.blast_program:
            token = str(self.blast_program).lower()
            if token == "blastn":
                return "nucl"
            if token == "blastp":
                return "prot"
        if self.blast_db_path:
            if any(os.path.exists(f"{self.blast_db_path}.{ext}") for ext in ("psq", "pin", "phr")):
                return "prot"
            if any(os.path.exists(f"{self.blast_db_path}.{ext}") for ext in ("nsq", "nin", "nhr")):
                return "nucl"
        return "prot"

    def _resolve_blast_program(self, db_type: str):
        if self.blast_program:
            token = str(self.blast_program).lower()
            if token == "blastn":
                self.blast_engine = "blast"
                return self.db_manager.env.get("BLASTN_PATH")
            if token == "blastp":
                self.blast_engine = "blast"
                return self.db_manager.env.get("BLASTP_PATH")
            if token == "diamond":
                self.blast_engine = "diamond"
                return self.db_manager.env.get("DIAMOND_PATH") or "diamond"
            return self.blast_program
        if db_type == "nucl":
            self.blast_engine = "blast"
            return self.db_manager.env.get("BLASTN_PATH")
        self.blast_engine = "blast"
        return self.db_manager.env.get("BLASTP_PATH")

    def _blast_cache_path(self, base_dir: str, accession: str, family_id: str) -> str:
        safe_acc = str(accession).replace(os.sep, "_")
        safe_fam = str(family_id).replace(os.sep, "_")
        return os.path.join(base_dir, safe_acc, f"{safe_fam}.blast6")

    def _ensure_output_dir(self):
        if not self.output_dir:
            if self.reuse_blast_results:
                self.output_dir = self.reuse_blast_results
            else:
                self.output_dir = tempfile.mkdtemp(prefix=f"ext_blast_{self.run_id}_")
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except OSError as exc:
            raise RuntimeError(f"Failed to create external BLAST output directory {self.output_dir}: {exc}") from exc
        self.data["output_dir"] = self.output_dir
        try:
            self.db_manager.tasks.update_data(self.task_id, data=self.data)
        except Exception as exc:  # boundary: persist enrichment only; external check can still proceed.
            self.log(f"Failed to persist external decontamination output directory: {exc}", "WARNING")

    def _run_blast_to_file(self, query_faa: str, db_path: str, max_targets: int, threads: int, out_path: str) -> bool:
        if not self.blast_path:
            self.error("BLAST program not set before attempting BLAST.")
            return False
        if self.blast_engine == "diamond":
            # DIAMOND supports blastp for protein DBs; use qcovhsp in place of qcovs.
            command = [
                self.blast_path,
                "blastp",
                "--query",
                query_faa,
                "--db",
                db_path,
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qcovhsp",
                "evalue",
                "bitscore",
                "staxids",
                "--out",
                out_path,
            ]
            if max_targets:
                command.extend(["--max-target-seqs", str(max_targets)])
            if threads and threads > 1:
                command.extend(["--threads", str(threads)])
        else:
            command = [
                self.blast_path,
                "-query",
                query_faa,
                "-db",
                db_path,
                "-outfmt",
                "6 qseqid sseqid pident length qcovs evalue bitscore staxids",
                "-out",
                out_path,
            ]
            if max_targets:
                command.extend(["-max_target_seqs", str(max_targets)])
            if threads and threads > 1:
                command.extend(["-num_threads", str(threads)])
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            stderr_text = (result.stderr or "").strip()
            lines = [ln.strip() for ln in stderr_text.splitlines() if ln.strip()]
            err_line = ""
            for ln in reversed(lines):
                if "error:" in ln.lower():
                    err_line = ln
                    break
            if not err_line:
                err_line = lines[-1] if lines else "no stderr captured"
            self.error(
                f"External BLAST failed (exit={result.returncode}): {err_line}"
            )
            return False
        return True

    def _iter_fasta_records(self, path: str):
        header = None
        seq_chunks = []
        try:
            with open(path, "r") as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header is not None:
                            yield header, "".join(seq_chunks)
                        header = line[1:].strip()
                        seq_chunks = []
                    else:
                        seq_chunks.append(line.strip())
                if header is not None:
                    yield header, "".join(seq_chunks)
        except (OSError, UnicodeError):
            return

    def _write_batched_query(self, accession: str, acc_jobs: list[tuple[str, str, str]]):
        batch_dir = os.path.join(self.output_dir, "_batch")
        os.makedirs(batch_dir, exist_ok=True)
        batch_path = os.path.join(batch_dir, f"{accession}_queries.faa")
        qseqid_to_out = {}
        written = 0
        with open(batch_path, "w") as fout:
            for family_id, query_path, out_path in acc_jobs:
                seq_idx = 0
                for _header, seq in self._iter_fasta_records(query_path):
                    if not seq:
                        continue
                    seq_idx += 1
                    qseqid = f"{family_id}__{seq_idx}"
                    qseqid_to_out[qseqid] = out_path
                    fout.write(f">{qseqid}\n")
                    # write sequence in 60-char lines
                    for i in range(0, len(seq), 60):
                        fout.write(seq[i : i + 60] + "\n")
                    written += 1
                if seq_idx == 0:
                    # No usable sequence: write empty output so downstream can skip
                    try:
                        os.makedirs(os.path.dirname(out_path), exist_ok=True)
                        with open(out_path, "w") as fh:
                            fh.write("")
                    except OSError as exc:
                        self.log(f"Failed to create empty external-check output {out_path}: {exc}", "WARNING")
        return batch_path, qseqid_to_out, written

    def _split_blast_output(self, blast_path: str, qseqid_to_out: dict):
        handles = {}
        try:
            with open(blast_path, "r") as fh:
                for line in fh:
                    if not line:
                        continue
                    qseqid = line.split("\t", 1)[0]
                    out_path = qseqid_to_out.get(qseqid)
                    if not out_path:
                        continue
                    handle = handles.get(out_path)
                    if handle is None:
                        os.makedirs(os.path.dirname(out_path), exist_ok=True)
                        handle = open(out_path, "w")
                        handles[out_path] = handle
                    handle.write(line)
        finally:
            for handle in handles.values():
                try:
                    handle.close()
                except OSError as exc:
                    self.log(f"Failed to close external BLAST split output handle: {exc}", "WARNING")

    def _pending_items(self, run_id: str):
        votes = self.db_manager.filtering.get_decontamination_votes(run_id=run_id)
        pending = []
        for row in votes or []:
            family_id = row[0]
            accession = row[3]
            decision = row[11]
            if str(decision).lower() != "outside":
                continue
            payload = {}
            try:
                payload = json.loads(row[12]) if row[12] else {}
            except (TypeError, json.JSONDecodeError):
                payload = {}
            if payload.get("external_check_run"):
                continue
            if payload.get("external_check_pending") is False:
                continue
            pending.append((accession, family_id))
        return pending

    def run(self):
        if not self.run_id:
            return self.handle_exception("run_id is required for external decontamination check.", {})

        if self.library_name and not self.library_id:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)

        run_row = self.db_manager.filtering.get_decontamination_run(self.run_id)
        if not run_row:
            return self.handle_exception("Decontamination run not found for run_id.", {"run_id": self.run_id})

        _run_id, target_library_id, busco_library_id, _targets_json, _refs_json, params_json, *_rest = run_row
        if self.library_id and target_library_id and int(self.library_id) != int(target_library_id):
            self.log(
                f"External check: target library mismatch (payload={self.library_id}, run={target_library_id}). Using run.",
                "WARNING",
            )
        self.library_id = target_library_id or self.library_id
        self.busco_library_id = busco_library_id

        params = {}
        try:
            params = json.loads(params_json) if params_json else {}
        except (TypeError, json.JSONDecodeError):
            params = {}
        if not self.hit_window and params.get("hit_window"):
            self.hit_window = params.get("hit_window")
        if not self.max_target_seqs:
            self.max_target_seqs = params.get("external_max_target_seqs") or params.get("max_target_seqs")
        if not self.blast_db_path:
            self.blast_db_path = params.get("external_blast_db_path")
        if not self.blast_db_type:
            self.blast_db_type = params.get("external_blast_db_type")
        if not self.blast_program:
            self.blast_program = params.get("external_blast_program")
        if not self.output_dir:
            self.output_dir = params.get("external_blast_output_dir")
        if not self.reuse_blast_results:
            self.reuse_blast_results = params.get("external_reuse_blast_results")
        if not self.data.get("max_concurrent"):
            run_mc = params.get("external_max_concurrent")
            if run_mc not in (None, ""):
                try:
                    self.max_concurrent = max(1, int(run_mc))
                except (TypeError, ValueError):
                    pass
        if not self.data.get("threads"):
            run_thr = params.get("external_threads")
            if run_thr not in (None, ""):
                try:
                    self.threads_budget = max(1, int(run_thr))
                except (TypeError, ValueError):
                    pass

        pending = self._pending_items(self.run_id)
        if not pending:
            self.log("External check: no pending BUSCOs flagged for external BLAST.", "INFO")
            return True

        self._ensure_output_dir()

        db_type = self._detect_db_type()
        if str(self.blast_program or "").lower() == "diamond" and db_type != "prot":
            return self.handle_exception(
                "diamond requires a protein database (db_type=prot).",
                {"blast_db_type": db_type},
            )
        self.blast_path = self._resolve_blast_program(db_type)
        if not self.blast_path and self.blast_db_path:
            if str(self.blast_program or "").lower() == "diamond":
                missing_var = "DIAMOND_PATH"
            else:
                missing_var = "BLASTN_PATH" if db_type == "nucl" else "BLASTP_PATH"
            return self.handle_exception(
                f"{missing_var} is not set in environment variables.",
                {"variable": missing_var},
            )

        max_targets = None
        if self.max_target_seqs not in (None, ""):
            try:
                max_targets = max(1, int(self.max_target_seqs))
            except (TypeError, ValueError):
                max_targets = None
        if max_targets is None:
            hit_window = 1
            try:
                hit_window = max(1, int(self.hit_window or 1))
            except (TypeError, ValueError):
                hit_window = 1
            max_targets = max(1, hit_window + 1)

        jobs = []
        skipped_existing = 0
        for accession, family_id in pending:
            out_path = self._blast_cache_path(self.output_dir, accession, family_id)
            if not self.force and os.path.exists(out_path):
                skipped_existing += 1
                continue
            if self.reuse_blast_results:
                reuse_path = self._blast_cache_path(self.reuse_blast_results, accession, family_id)
                if os.path.exists(reuse_path) and not self.force:
                    skipped_existing += 1
                    continue
            query_path = self.db_manager.busco.get_family_location(
                family_id,
                self.busco_library_id,
                accession,
                sequence_kind=db_type if db_type in ("prot", "nucl") else None,
            )
            if not query_path or not os.path.exists(query_path):
                self.log(
                    f"External check: missing BUSCO query for {accession} family {family_id}; writing empty output.",
                    "WARNING",
                )
                try:
                    os.makedirs(os.path.dirname(out_path), exist_ok=True)
                    with open(out_path, "w") as fh:
                        fh.write("")
                except OSError as exc:
                    self.log(f"External check: failed to write empty output for missing query {out_path}: {exc}", "WARNING")
                continue
            jobs.append((accession, family_id, query_path, out_path))

        if not jobs:
            self.log(
                f"External check: no BLAST jobs queued (skipped existing={skipped_existing}).",
                "INFO",
            )
            return True

        if not self.blast_db_path:
            return self.handle_exception("blast_db_path is required for external BLAST.", {"run_id": self.run_id})

        # Run per-assembly batches, with optional concurrency across assemblies.
        jobs_by_acc = {}
        for acc, family_id, query_path, out_path in jobs:
            jobs_by_acc.setdefault(acc, []).append((family_id, query_path, out_path))

        effective_max_concurrent = max(1, min(int(self.max_concurrent), len(jobs_by_acc)))
        if self.blast_threads is None:
            blast_threads = max(1, int(self.threads_budget // effective_max_concurrent))
        else:
            blast_threads = max(1, int(self.blast_threads))
        self.log(
            f"External check: concurrency={effective_max_concurrent}, threads_budget={self.threads_budget}, "
            f"blast_threads_per_process={blast_threads}, engine={self.blast_engine}",
            "INFO",
        )

        def _run_accession_batch(acc, acc_jobs):
            self.log(
                f"External check: running {len(acc_jobs)} BUSCOs for accession {acc}",
                "INFO",
            )
            batch_path, qseqid_to_out, _written = self._write_batched_query(acc, acc_jobs)
            if not qseqid_to_out:
                return acc, len(acc_jobs), True
            out_path = os.path.join(self.output_dir, "_batch", f"{acc}_blast6.tsv")
            ok = self._run_blast_to_file(batch_path, self.blast_db_path, max_targets, blast_threads, out_path)
            if ok:
                self._split_blast_output(out_path, qseqid_to_out)
                # Ensure every family has a file, even if no hits were returned.
                for fam_out in set(qseqid_to_out.values()):
                    if not os.path.exists(fam_out):
                        try:
                            os.makedirs(os.path.dirname(fam_out), exist_ok=True)
                            with open(fam_out, "w") as fh:
                                fh.write("")
                        except OSError as exc:
                            self.log(f"External check: failed to create empty result file {fam_out}: {exc}", "WARNING")
            else:
                # mark outputs as empty to avoid repeated retries
                for _q, fam_out in qseqid_to_out.items():
                    try:
                        os.makedirs(os.path.dirname(fam_out), exist_ok=True)
                        with open(fam_out, "w") as fh:
                            fh.write("")
                    except OSError as exc:
                        self.log(f"External check: failed to mark failed result as empty {fam_out}: {exc}", "WARNING")
            return acc, len(acc_jobs), ok

        completed = 0
        total = len(jobs)
        with concurrent.futures.ThreadPoolExecutor(max_workers=effective_max_concurrent) as executor:
            futures = [
                executor.submit(_run_accession_batch, acc, acc_jobs)
                for acc, acc_jobs in jobs_by_acc.items()
            ]
            for fut in concurrent.futures.as_completed(futures):
                try:
                    acc, acc_done, ok = fut.result()
                except Exception as exc:  # boundary: one worker failure is aggregated while independent accession batches continue.
                    self.error(f"External check: worker failed unexpectedly: {exc}")
                    continue
                completed += acc_done
                if not ok:
                    self.log(f"External check: accession batch failed ({acc}); wrote empty outputs.", "WARNING")
                self.log(
                    f"External check: {completed}/{total} BUSCOs processed ({acc})",
                    "INFO",
                )

        self.log(
            f"External check: completed BLAST for {len(jobs)} BUSCOs (skipped existing={skipped_existing}).",
            "INFO",
        )
        return True
