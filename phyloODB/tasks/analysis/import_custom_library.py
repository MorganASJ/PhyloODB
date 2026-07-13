import json
import os
import shutil
from datetime import datetime
from json import JSONDecodeError
from typing import Dict, List, Optional, Tuple

from ..task import Task
from ...selector_utils import normalize_accessions, resolve_clade_to_taxid


class ImportCustomLibraryTask(Task):
    """Import a pre-selected BUSCO family list as a custom library."""

    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.library_name = self.data.get("library_name")
        self.parent_library_id = self.data.get("parent_library_id")
        self.parent_library_name = self.data.get("parent_library_name")
        self.coverage = self.data.get("coverage")
        self.coverage_taxid = self.data.get("coverage_taxid")
        self.busco_ids = self.data.get("busco_ids") or []
        self.ref_accessions = self.data.get("ref_accessions") or []
        self.force = bool(self.data.get("force", False))
        self.location = self.data.get("location")
        self.library_id = self.data.get("library_id")
        self.coverage_label = self.data.get("coverage_label", self.coverage)

    def _load_busco_ids_file(self, path: str) -> List[str]:
        ids: List[str] = []
        if path.lower().endswith(".json"):
            try:
                with open(path, "r") as handle:
                    payload = json.load(handle)
                if isinstance(payload, dict):
                    ids = [str(k).strip() for k in payload.keys()]
                elif isinstance(payload, list):
                    ids = [str(item).strip() for item in payload]
                else:
                    raise ValueError("JSON must be a list or object of BUSCO ids.")
            except (OSError, UnicodeError, JSONDecodeError, ValueError) as exc:
                raise ValueError(f"Failed to parse BUSCO ids from JSON file: {exc}") from exc
        else:
            try:
                with open(path, "r") as handle:
                    for line in handle:
                        token = line.strip()
                        if not token or token.startswith("#"):
                            continue
                        token = token.split()[0]
                        if token:
                            ids.append(token)
            except (OSError, UnicodeError) as exc:
                raise ValueError(f"Failed to read BUSCO ids from file: {exc}") from exc

        ids = [str(item) for item in ids if str(item).strip()]
        return list(dict.fromkeys(ids))

    def _resolve_busco_ids(self) -> List[str]:
        raw = self.busco_ids
        if raw is None:
            return []
        if isinstance(raw, str):
            raw_values: List[str] = [raw]
        elif isinstance(raw, (list, tuple, set)):
            raw_values = [str(item) for item in raw if item is not None]
        else:
            raw_values = [str(raw)]

        if len(raw_values) == 1:
            candidate = str(raw_values[0]).strip()
            if candidate:
                if os.path.isfile(candidate):
                    return self._load_busco_ids_file(candidate)
                looks_like_file = any(
                    token in candidate for token in ("/", "\\")
                ) or candidate.lower().endswith((".txt", ".tsv", ".csv", ".json"))
                if looks_like_file:
                    raise ValueError(f"BUSCO ids file not found: {candidate}")

        ids: List[str] = []
        for value in raw_values:
            for token in str(value).split(","):
                token = token.strip()
                if token:
                    ids.append(token)
        ids = [str(item) for item in ids if str(item).strip()]
        return list(dict.fromkeys(ids))

    def _resolve_ref_accessions(self) -> List[str]:
        raw = self.ref_accessions
        if raw is None:
            return []
        if isinstance(raw, str):
            raw_values: List[str] = [raw]
        elif isinstance(raw, (list, tuple, set)):
            raw_values = [str(item) for item in raw if item is not None]
        else:
            raw_values = [str(raw)]

        if len(raw_values) == 1:
            candidate = str(raw_values[0]).strip()
            if candidate:
                if os.path.isfile(candidate):
                    ids = self._load_busco_ids_file(candidate)
                    return list(dict.fromkeys(normalize_accessions(ids)))
                looks_like_file = any(
                    token in candidate for token in ("/", "\\")
                ) or candidate.lower().endswith((".txt", ".tsv", ".csv"))
                if looks_like_file:
                    raise ValueError(f"Reference accession file not found: {candidate}")

        ids: List[str] = []
        for value in raw_values:
            for token in str(value).split(","):
                token = token.strip()
                if token:
                    ids.append(token)
        return list(dict.fromkeys(normalize_accessions(ids)))

    def _fetch_parent_descriptions(self, parent_id: int, families: List[str]) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
        desc_rows: List[Tuple[str, int, Optional[str], Optional[str]]] = []
        for i in range(0, len(families), 900):
            chunk = families[i:i + 900]
            desc_rows.extend(self.db_manager.libraries.get_busco_descriptions(parent_id, chunk) or [])
        return {str(fid): (desc, link) for fid, _lib, desc, link in desc_rows}

    def run(self):
        if not self.library_name:
            return self.handle_exception("Library name is not specified.", {"library_name": self.library_name})
        self.log(
            f"Importing custom library '{self.library_name}' into {self.location or '[default libraries root]'}.",
            "INFO",
        )
        if not self.parent_library_id and not self.parent_library_name:
            return self.handle_exception(
                "Parent library id or name is required.",
                {"parent_library_id": self.parent_library_id, "parent_library_name": self.parent_library_name},
            )

        if self.parent_library_id:
            try:
                self.parent_library_id = int(self.parent_library_id)
            except (TypeError, ValueError):
                return self.handle_exception("Parent library id is invalid.", {"parent_library_id": self.parent_library_id})
            parent_rows = self.db_manager.libraries.get(self.parent_library_id) or []
            if not parent_rows:
                return self.handle_exception(
                    "Parent library id not found in database.",
                    {"parent_library_id": self.parent_library_id},
                )
            self.parent_library_name = parent_rows[0][1]
        else:
            self.parent_library_id = self.db_manager.libraries.get_id(self.parent_library_name)
            if not self.parent_library_id:
                return self.handle_exception(
                    f"Parent library '{self.parent_library_name}' not found in database.",
                    {"parent_library_name": self.parent_library_name},
                )

        if not self.coverage_taxid:
            if self.coverage is None:
                return self.handle_exception("Coverage taxid is not specified.", {"coverage": self.coverage})

            inferred_taxid: Optional[int]
            if isinstance(self.coverage, int):
                inferred_taxid = self.coverage
            else:
                candidate = None
                try:
                    if isinstance(self.coverage, str) and self.coverage.strip():
                        candidate = int(self.coverage.strip())
                except ValueError:
                    candidate = None
                if candidate is not None:
                    inferred_taxid = candidate
                else:
                    inferred_taxid = resolve_clade_to_taxid(self.db_manager, str(self.coverage))

            if not inferred_taxid:
                return self.handle_exception(
                    f"Could not infer coverage taxid from '{self.coverage}'.",
                    {"coverage": self.coverage},
                )
            self.coverage_taxid = inferred_taxid

        try:
            coverage_taxid = int(self.coverage_taxid)
        except (TypeError, ValueError):
            return self.handle_exception("Coverage taxid is not specified or invalid.", {"coverage_taxid": self.coverage_taxid})

        if coverage_taxid <= 0:
            return self.handle_exception("Coverage taxid is not specified or invalid.", {"coverage_taxid": coverage_taxid})

        self.coverage_taxid = coverage_taxid
        self.coverage = coverage_taxid

        try:
            families = self._resolve_busco_ids()
            ref_accessions = self._resolve_ref_accessions()
        except Exception as exc:  # boundary: convert selector/input resolution failures into task error state.
            return self.handle_exception(str(exc), {"busco_ids": self.busco_ids, "ref_accessions": self.ref_accessions})

        if not families:
            return self.handle_exception("BUSCO id list is empty.", {"busco_ids": self.busco_ids})

        libraries_dir = self.db_manager.storage.get_root_base("libraries")
        if not libraries_dir:
            return self.handle_exception("Libraries directory is not configured.", {})
        libraries_dir = str(libraries_dir)

        existing_id = self.db_manager.libraries.get_id(self.library_name)
        existing_record = None
        if existing_id:
            records = self.db_manager.libraries.get(existing_id) or []
            existing_record = records[0] if records else None
            if not self.force:
                return self.handle_exception(
                    f"Library '{self.library_name}' already exists. Re-run with --force to rebuild.",
                    {"library_name": self.library_name},
                )
            self.library_id = existing_id

        if self.location:
            location = self.location
        elif existing_record and existing_record[5]:
            location = existing_record[5]
        else:
            location = os.path.join(libraries_dir, self.library_name)

        if self.force and self.library_id:
            if not self.db_manager.libraries.purge(self.library_id):
                return self.handle_exception(
                    "Failed to purge existing library data before rebuild.",
                    {"library_id": self.library_id},
                )
            if location and os.path.isdir(location):
                try:
                    shutil.rmtree(location)
                except OSError as exc:
                    return self.handle_exception(
                        "Failed to remove existing library directory before rebuild.",
                        {"location": location, "error": str(exc)},
                    )

        try:
            os.makedirs(location, exist_ok=True)
        except OSError as exc:
            return self.handle_exception(
                "Failed to create library directory.",
                {"location": location, "error": str(exc)},
            )

        new_id = self.db_manager.libraries.add(
            library_name=self.library_name,
            taxid=self.coverage_taxid,
            size=len(families),
            location=location,
            parent_id=self.parent_library_id,
            ref_accessions=None,
        )
        if not new_id:
            return self.handle_exception(
                "Failed to create or update library record.",
                {"library_name": self.library_name, "parent_library_id": self.parent_library_id},
            )
        self.library_id = new_id
        self.location = location

        desc_map = self._fetch_parent_descriptions(self.parent_library_id, families)
        missing = [fam for fam in families if str(fam) not in desc_map]
        if missing:
            sample = missing[:25]
            return self.handle_exception(
                f"{len(missing)} BUSCO families not found in parent library '{self.parent_library_name}'.",
                {"missing_sample": sample, "missing_count": len(missing)},
            )

        subset_rows = []
        for fam in families:
            desc, link = desc_map.get(str(fam), (None, None))
            subset_rows.append((str(fam), self.library_id, desc, link))
        if subset_rows:
            if not self.db_manager.libraries.add_busco_descriptions(subset_rows):
                return self.handle_exception(
                    "Failed to record BUSCO subset for custom library.",
                    {"library_id": self.library_id},
                )
            self.db_manager.libraries.update_size(self.library_id, len(families))
        else:
            return self.handle_exception("No BUSCO families provided.", {"library_id": self.library_id})

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            busco_families_file = os.path.join(location, f"cleaned_busco_families_{timestamp}.json")
            with open(busco_families_file, "w") as handle:
                json.dump(families, handle, indent=4)
            canonical_path = os.path.join(location, "cleaned_busco_families.json")
            with open(canonical_path, "w") as handle:
                json.dump(families, handle, indent=4)
        except (OSError, TypeError, ValueError) as exc:
            return self.handle_exception(
                "Failed to save BUSCO families file.",
                {"location": location, "error": str(exc)},
            )
        self.db_manager.artifacts.register(
            owner_type="library",
            owner_id=self.library_id,
            artifact_type="library_root",
            path=location,
            is_dir=True,
            format="directory",
            metadata={"library_id": self.library_id, "library_name": self.library_name},
        )
        self.db_manager.artifacts.register(
            owner_type="library",
            owner_id=self.library_id,
            artifact_type="library_core_set_json",
            path=os.path.join(location, "cleaned_busco_families.json"),
            format="json",
            metadata={"library_id": self.library_id, "library_name": self.library_name},
        )

        if ref_accessions:
            present_rows = self.db_manager.genomes.get_many(ref_accessions) or []
            present = [row[0] for row in present_rows]
            missing = [acc for acc in ref_accessions if acc not in present]
            if missing:
                self.log(
                    f"{len(missing)} reference accessions not found in Genome; they will be skipped.",
                    "WARNING",
                )
            if present:
                if not self.db_manager.libraries.add_reference_assemblies(self.library_id, present):
                    return self.handle_exception(
                        "Failed to attach reference accessions to library.",
                        {"library_id": self.library_id, "ref_accessions": present},
                    )

        self.data.update(
            {
                "library_id": self.library_id,
                "location": self.location,
                "coverage_taxid": self.coverage_taxid,
                "coverage": self.coverage,
                "coverage_label": self.coverage_label,
                "parent_library_id": self.parent_library_id,
                "parent_library_name": self.parent_library_name,
                "ref_accessions": ref_accessions,
            }
        )
        try:
            self.db_manager.tasks.update_data(self.task_id, data=self.data)
        except Exception as exc:  # boundary: persist enrichment only; task result is still scientifically complete.
            self.log(f"Failed to persist imported custom library task metadata: {exc}", "WARNING")

        self.log(
            f"Imported custom library '{self.library_name}' (ID {self.library_id}) with {len(families)} BUSCO families.",
            "INFO",
        )
        return True
