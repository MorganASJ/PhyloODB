import os
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List

from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE, RAW_PROFILE, derive_profile_name_from_recipe
from ...proteome_state import summarize_proteome_state
from ..task import Task
from ..utilities import clean_proteome_in_genome_path, prepare_proteome_profile


def _write_gzip_from_plain(src_plain: str, dst_gz: str) -> None:
    with open(src_plain, "rb") as f_in, gzip.open(dst_gz, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def _normalize_gff_files(genome_path: str) -> str | None:
    preferred_gff = None
    gz_by_base: Dict[str, str] = {}
    plain_by_base: Dict[str, str] = {}

    for fname in sorted(os.listdir(genome_path)):
        low = fname.lower()
        path = os.path.join(genome_path, fname)
        if low.endswith((".gff.gz", ".gff3.gz")):
            base = fname[:-3]
            gz_by_base[base] = path
        elif low.endswith((".gff", ".gff3")):
            plain_by_base[fname] = path

    for base, plain_path in plain_by_base.items():
        gz_path = gz_by_base.get(base)
        if gz_path is None:
            gz_path = plain_path + ".gz"
            _write_gzip_from_plain(plain_path, gz_path)
            gz_by_base[base] = gz_path
        if os.path.exists(gz_path) and os.path.exists(plain_path):
            os.remove(plain_path)

    if gz_by_base:
        preferred_gff = gz_by_base[sorted(gz_by_base.keys())[0]]
    elif plain_by_base:
        preferred_gff = plain_by_base[sorted(plain_by_base.keys())[0]]
    return preferred_gff


def _refresh_genome_artifacts(task: Task, accession: str, genome_path: str) -> None:
    proteome_state = summarize_proteome_state(genome_path)
    preferred_faa = proteome_state.active_faa
    archive_faa = proteome_state.archive_faa
    preferred_gff = _normalize_gff_files(genome_path)
    if preferred_faa and os.path.exists(preferred_faa):
        task.db_manager.artifacts.register(
            owner_type="genome",
            owner_id=accession,
            artifact_type="genome_faa",
            path=preferred_faa,
            metadata={"accession": accession},
        )
    else:
        for row in task.db_manager.artifacts.find(owner_type="genome", owner_id=accession, artifact_type="genome_faa"):
            task.db_manager.artifacts.set_status(int(row[0]), "stale")
    if preferred_gff and os.path.exists(preferred_gff):
        task.db_manager.artifacts.register(
            owner_type="genome",
            owner_id=accession,
            artifact_type="genome_gff",
            path=preferred_gff,
            metadata={"accession": accession},
        )
    else:
        for row in task.db_manager.artifacts.find(owner_type="genome", owner_id=accession, artifact_type="genome_gff"):
            task.db_manager.artifacts.set_status(int(row[0]), "stale")
    if archive_faa and os.path.exists(archive_faa):
        task.db_manager.artifacts.register(
            owner_type="genome",
            owner_id=accession,
            artifact_type="genome_faa_archive",
            path=archive_faa,
            metadata={"accession": accession},
        )
    else:
        for row in task.db_manager.artifacts.find(owner_type="genome", owner_id=accession, artifact_type="genome_faa_archive"):
            task.db_manager.artifacts.set_status(int(row[0]), "stale")


class PrepareProteomeTask(Task):
    """Create immutable proteome profile artifacts without mutating the raw proteome."""

    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads=required_threads)
        self.skip_gff = self.payload_bool("skip_gff", not self.env_bool("DEFAULT_PROTEOME_USE_GFF", True))
        self.skip_cdhit = self.payload_bool("skip_cdhit", not self.env_bool("DEFAULT_PROTEOME_USE_CDHIT", True))
        self.gff_priority = self.payload_bool("gff_priority", self.env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False))
        self.downloaded_only = bool(self.data.get("downloaded_only", True))
        self.max_concurrent = int(self.data.get("max_concurrent", self.env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1)) or 1)
        self.threads_per_job = int(self.data.get("threads_per_job", self.env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1)) or 1)
        self.parent_threads_budget = int(self.data.get("parent_threads_budget", 0) or 0)
        self.profile_name = str(self.data.get("profile_name") or DEFAULT_CLEAN_PROFILE).strip() or DEFAULT_CLEAN_PROFILE
        self.profile_is_default_alias = self.profile_name == DEFAULT_CLEAN_PROFILE
        self.input_profile = str(self.data.get("input_profile") or self.env_str("DEFAULT_PROTEOME_INPUT_PROFILE", RAW_PROFILE)).strip() or RAW_PROFILE
        self.replace_existing = bool(self.data.get("replace_existing", False))
        self.set_default = self.payload_bool("set_default", self.env_bool("DEFAULT_PROTEOME_SET_DEFAULT", True))
        self.cdhit_identity = self._resolve_cdhit_identity()

    def _resolve_cdhit_identity(self) -> float:
        explicit = self.data.get("cdhit_identity")
        if explicit is not None:
            return float(explicit)
        env_value = self.db_manager.get_environment_variable("DEFAULT_PROTEOME_CDHIT_IDENTITY")
        if env_value is not None and str(env_value).strip() != "":
            try:
                return float(env_value)
            except (TypeError, ValueError):
                pass
        return 0.96

    def _effective_parallelism(self, total_jobs: int) -> tuple[int, int]:
        budget = self.parent_threads_budget if self.parent_threads_budget > 0 else int(self.REQUIRED_THREADS or 1)
        budget = min(8, max(1, budget))
        requested_threads = max(1, self.threads_per_job)
        effective_threads = min(requested_threads, budget)
        requested_workers = max(1, self.max_concurrent)
        max_workers_by_budget = max(1, budget // effective_threads)
        effective_workers = min(requested_workers, max_workers_by_budget, 8, max(1, total_jobs))
        return effective_workers, effective_threads

    def _resolve_targets(self) -> List[str]:
        return self.prepare_selectors(
            taxid=self.data.get("taxid"),
            additional=self.data.get("accessions"),
            downloaded_only=self.downloaded_only,
            released_after=self.data.get("after"),
            released_before=self.data.get("before"),
            level=self.data.get("level"),
            primary_only=self.data.get("primary_only"),
            require_candidates=True,
        )

    def _resolve_input_profile(self, accession: str, genome_path: str):
        if self.input_profile == RAW_PROFILE:
            raw_id = self.db_manager.proteomes.ensure_raw_profile(accession, is_default=False)
            if raw_id is None:
                raise FileNotFoundError(f"No raw proteome found for accession '{accession}'.")
            row = self.db_manager.proteomes.get(int(raw_id))
        else:
            row = self.db_manager.proteomes.get_profile(accession, self.input_profile)
            if row is None:
                raise FileNotFoundError(f"Proteome profile '{self.input_profile}' does not exist for accession '{accession}'.")
        path = self.db_manager.proteomes.resolve_path(row)
        if not path or not os.path.exists(path):
            raise FileNotFoundError(f"Proteome profile '{self.input_profile}' has no readable artifact for accession '{accession}'.")
        return row, path

    def run(self):
        try:
            selected = self._resolve_targets()
        except ValueError as exc:
            return self.handle_exception(str(exc), {})

        if not selected:
            return self.handle_exception("No assemblies selected for proteome preparation.", {})

        workers, cdhit_threads = self._effective_parallelism(len(selected))
        self.log(
            f"Preparing profile '{self.profile_name}' for {len(selected)} assemblies with max_concurrent={workers}, threads_per_job={cdhit_threads}",
            "INFO",
        )
        results: Dict[str, Dict[str, object]] = {}
        jobs: Dict[str, Dict[str, object]] = {}
        for accession in selected:
            genome_path = self.db_manager.genomes.get_path(accession)
            if not genome_path or not os.path.isdir(genome_path):
                results[accession] = {"ok": False, "status": "missing_genome_path", "accession": accession}
                continue
            existing = self.db_manager.proteomes.get_profile(accession, self.profile_name) if not self.profile_is_default_alias else None
            if existing and not self.replace_existing:
                results[accession] = {
                    "ok": True,
                    "status": "skipped_existing_profile",
                    "accession": accession,
                    "profile_id": int(existing[0]),
                }
                continue
            try:
                input_row, input_path = self._resolve_input_profile(accession, genome_path)
            except FileNotFoundError as exc:
                results[accession] = {"ok": False, "status": "missing_input_profile", "error": str(exc), "accession": accession}
                continue
            profiles_dir = os.path.join(genome_path, "proteome_profiles")
            temp_profile_name = self.profile_name if not self.profile_is_default_alias else f".tmp_{self.task_id}_{accession}"
            output_path = os.path.join(profiles_dir, f"{temp_profile_name}.faa.gz")
            jobs[accession] = {
                "accession": accession,
                "genome_path": genome_path,
                "input_profile_id": int(input_row[0]),
                "input_path": input_path,
                "output_path": output_path,
            }

        def _prepare_one(job: Dict[str, object]) -> Dict[str, object]:
            outcome = prepare_proteome_profile(
                str(job["input_path"]),
                str(job["genome_path"]),
                str(job["output_path"]),
                skip_gff=self.skip_gff,
                skip_cdhit=self.skip_cdhit,
                gff_priority=self.gff_priority,
                cdhit_identity=self.cdhit_identity,
                cdhit_threads=cdhit_threads,
                silent=True,
            )
            outcome["accession"] = str(job["accession"])
            return outcome

        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="PrepProt") as executor:
            future_map = {executor.submit(_prepare_one, job): accession for accession, job in jobs.items()}
            for future in as_completed(future_map):
                acc = future_map[future]
                try:
                    result = future.result()
                except Exception as exc:  # boundary: isolate one proteome-preparation worker failure
                    result = {"ok": False, "status": "exception", "error": str(exc), "accession": acc}
                results[acc] = result

        for accession, res in results.items():
            if str(res.get("status")) != "prepared":
                continue
            job = jobs.get(accession)
            if not job:
                continue
            actual_profile_name = self.profile_name
            if self.profile_is_default_alias:
                actual_profile_name = derive_profile_name_from_recipe(
                    used_gff=bool(res.get("used_gff")),
                    used_cdhit=bool(res.get("used_cdhit")),
                    cdhit_identity=self.cdhit_identity if bool(res.get("used_cdhit")) else None,
                    fallback=DEFAULT_CLEAN_PROFILE,
                )
                final_output_path = os.path.join(str(job["genome_path"]), "proteome_profiles", f"{actual_profile_name}.faa.gz")
                existing_actual = self.db_manager.proteomes.get_profile(accession, actual_profile_name)
                if existing_actual and not self.replace_existing:
                    try:
                        if os.path.exists(str(job["output_path"])):
                            os.remove(str(job["output_path"]))
                    except OSError:
                        pass
                    if self.set_default:
                        self.db_manager.proteomes.set_default_profile(accession, profile_name=actual_profile_name)
                    res["status"] = "skipped_existing_profile"
                    res["profile_id"] = int(existing_actual[0])
                    res["profile_name"] = actual_profile_name
                    continue
                os.makedirs(os.path.dirname(final_output_path), exist_ok=True)
                if os.path.abspath(str(job["output_path"])) != os.path.abspath(final_output_path):
                    if os.path.exists(final_output_path):
                        os.remove(final_output_path)
                    os.replace(str(job["output_path"]), final_output_path)
                job["output_path"] = final_output_path
            res["profile_name"] = actual_profile_name
            profile_id = self.db_manager.proteomes.register_profile(
                accession=accession,
                profile_name=actual_profile_name,
                path=str(job["output_path"]),
                kind="raw" if actual_profile_name == RAW_PROFILE else "derived",
                parent_profile_id=None if actual_profile_name == RAW_PROFILE else int(job["input_profile_id"]),
                is_default=self.set_default,
                metadata={
                    "source": "prepare-proteome",
                    "requested_profile_name": self.profile_name,
                    "actual_profile_name": actual_profile_name,
                    "input_profile": self.input_profile,
                    "skip_gff": self.skip_gff,
                    "skip_cdhit": self.skip_cdhit,
                    "gff_priority": self.gff_priority,
                    "cdhit_identity": self.cdhit_identity,
                },
            )
            gff_artifacts = self.db_manager.artifacts.find(owner_type="genome", owner_id=accession, artifact_type="genome_gff")
            gff_artifact_id = int(gff_artifacts[0][0]) if gff_artifacts else None
            self.db_manager.proteomes.record_preparation(
                accession=accession,
                input_profile_id=int(job["input_profile_id"]),
                output_profile_id=int(profile_id),
                preparation_type="isoform_clean",
                used_gff=bool(res.get("used_gff")),
                gff_artifact_id=gff_artifact_id,
                skip_gff=self.skip_gff,
                skip_cdhit=self.skip_cdhit,
                gff_priority=self.gff_priority,
                cdhit_identity=self.cdhit_identity,
                cdhit_threads=cdhit_threads,
                input_count=res.get("input_count"),
                output_count=res.get("output_count"),
                gff_removed=res.get("gff_removed"),
                cdhit_removed=res.get("cdhit_removed"),
                total_removed=res.get("total_removed"),
                params={
                    "input_profile": self.input_profile,
                    "requested_profile_name": self.profile_name,
                    "profile_name": actual_profile_name,
                    "gff_file": res.get("gff_file"),
                    "cdhit_skipped_due_gff_priority": res.get("cdhit_skipped_due_gff_priority"),
                },
            )
            res["profile_id"] = int(profile_id)
            try:
                _refresh_genome_artifacts(self, accession, str(job["genome_path"]))
            except Exception as exc:  # boundary: collect one failed post-profile artifact refresh
                self.collect_batch_failure(accession, "refresh proteome artifacts", exc)
            try:
                clean_ready = actual_profile_name != RAW_PROFILE
                self.db_manager.genomes.set_isoforms_cleaned(accession, clean_ready)
            except Exception as exc:  # boundary: collect one failed post-profile genome status update
                self.collect_batch_failure(accession, "set isoforms-cleaned status", exc)

        created = sum(1 for r in results.values() if r.get("status") == "prepared")
        skipped = sum(1 for r in results.values() if str(r.get("status", "")).startswith("skipped"))
        failed = sum(1 for r in results.values() if not r.get("ok"))
        for acc, res in results.items():
            status = str(res.get("status") or "")
            if status == "prepared":
                self.log(
                    f"{acc}: prepared profile '{res.get('profile_name', self.profile_name)}' from '{self.input_profile}' "
                    f"({int(res.get('input_count', 0) or 0)}->{int(res.get('output_count', 0) or 0)} proteins).",
                    "DEBUG",
                )
            else:
                self.log(f"{acc}: {status}", "DEBUG" if res.get("ok") else "ERROR")

        if failed:
            failed_items = [f"{acc}:{res.get('status')}" for acc, res in results.items() if not res.get("ok")]
            return self.handle_exception(
                f"Proteome preparation failed for {failed}/{len(selected)} assemblies.",
                {"failures": failed_items[:50], "total_failures": failed},
            )
        if self._batch_failures:
            return self.fail_if_batch_failures("Proteome preparation database update failed")

        self.log(
            f"Proteome preparation summary: created={created}, skipped={skipped}, failed={failed}, profile={self.profile_name}",
            "INFO",
        )
        return True
