import os
import re
import json
import gzip
import hashlib
from datetime import datetime
from Bio import Entrez
import urllib.request
from xml.etree import ElementTree as ET
import time

MAX_ATTEMPTS = 5
LIBRARIES = "./libraries"
GENOMES = "./genomes"


def infer_ncbi_origin(accession):
    token = str(accession or "").strip().upper()
    if token.startswith("GCF_"):
        return "refseq"
    if token.startswith("GCA_"):
        return "genbank"
    return None

class FnaDownloadError(Exception):
    """Exception raised when failing to download the .fna file."""
    def __init__(self, message="Failed to download .fna file"):
        self.message = message
        super().__init__(self.message)

class FaaDownloadError(Exception):
    """Exception raised when failing to download the .faa file."""
    def __init__(self, message="Failed to download .faa file"):
        self.message = message
        super().__init__(self.message)

class GffDownloadError(Exception):
    """Exception raised when failing to download the .gff file."""
    def __init__(self, message="Failed to download .gff file"):
        self.message = message
        super().__init__(self.message)


def _parse_ncbi_md5checksums(path):
    """Return {filename: md5} from an NCBI md5checksums.txt file."""
    checksums = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            digest = parts[0].strip().lower()
            filename = os.path.basename(parts[-1].lstrip("./"))
            if re.fullmatch(r"[0-9a-f]{32}", digest) and filename:
                checksums[filename] = digest
    return checksums


def _md5_file(path):
    hasher = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _validate_gzip_to_eof(path):
    with gzip.open(path, "rb") as handle:
        for _chunk in iter(lambda: handle.read(1024 * 1024), b""):
            pass


def _download_verified_file(url, destination, *, expected_md5=None, gzip_file=True):
    """Download to a staging path, verify, then atomically move into place."""
    part_path = f"{destination}.part"
    if os.path.exists(part_path):
        os.remove(part_path)
    urllib.request.urlretrieve(url, part_path)
    try:
        if expected_md5:
            actual_md5 = _md5_file(part_path)
            if actual_md5.lower() != str(expected_md5).lower():
                raise ValueError(
                    f"MD5 mismatch for {os.path.basename(destination)}: "
                    f"expected {expected_md5}, observed {actual_md5}"
                )
        if gzip_file:
            _validate_gzip_to_eof(part_path)
        os.replace(part_path, destination)
    except Exception:  # boundary: staged download verification cleans partial files and re-raises original failure.
        try:
            os.remove(part_path)
        except OSError:
            pass
        raise


class NCBIHelper:
    """A class that handles NCBI API calls"""
    def __init__(self, email, db_manager=None, api_key=None):
        self.email = email
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        # Optional DB manager to allow checking for existing taxonomy entries
        self.db_manager = db_manager

    def safe_entrez_call(self, func, *args, max_attempts=MAX_ATTEMPTS, **kwargs):
        """Generic retry wrapper for Entrez requests."""
        attempt = 0
        while attempt < max_attempts:
            try:
                result = func(*args, **kwargs)
                return result
            except Exception as e:  # boundary: retry wrapper handles transient Entrez/network failures.
                attempt += 1
                print(f"Error in {func.__name__} attempt {attempt}: {e}")
                time.sleep(1)  # brief pause before retry
                if attempt == max_attempts:
                    raise RuntimeError(f"Maximum attempts reached for {func.__name__}. Error: {e}") from e

    def get_assembly_summary(self, id):
        """Get esummary for an entrez id"""
        Entrez.email = self.email
        esummary_handle = self.safe_entrez_call(Entrez.esummary, db="assembly", id=id, report="full")
        esummary_record = self.safe_entrez_call(Entrez.read, esummary_handle)
        return esummary_record

    def parse_date(self, date_str):
        """Parse date strings from the assembly summary; accept date or datetime."""
        if not date_str:
            return None
        s = str(date_str).strip()
        # NCBI sometimes uses this sentinel or empty values
        if not s or s == '1/01/01 00:00':
            return None
        # Try common formats, prefer returning ISO date (YYYY-MM-DD)
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%Y/%m/%d %H:%M", "%Y/%m/%d %H:%M:%S"):
            try:
                dt = datetime.strptime(s, fmt)
                return dt.strftime("%Y-%m-%d")
            except ValueError:
                continue
        # Fallback: try ISO parser and return just the date part
        try:
            dt = datetime.fromisoformat(s)
            return dt.date().isoformat()
        except ValueError:
            # Last resort: leave unparsed
            return None

    def parse_datetime(self, dt_str):
        """Parse datetime strings, preserving date and time."""
        if not dt_str:
            return None
        s = str(dt_str).strip()
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"):
            try:
                return datetime.strptime(s, fmt).strftime("%Y-%m-%d %H:%M:%S")
            except ValueError:
                continue
        try:
            return datetime.fromisoformat(s).strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            return None

    def parse_meta(self, meta_str):
        """Extract and parse the first <Stats> block from the metadata."""
        match = re.search(r"(<Stats>.*?</Stats>)", meta_str, re.DOTALL)
        if not match:
            print(f"Could not locate <Stats> block in metadata:\n{meta_str}")
            return {}
        stats_xml = match.group(1)
        try:
            root = ET.fromstring(stats_xml)
        except ET.ParseError as e:
            print(f"Error parsing Stats XML: {e}\n{stats_xml}")
            return {}
        stats = {}
        for stat in root.findall('.//Stat'):
            category = stat.get('category')
            seq_tag = stat.get('sequence_tag', 'all')
            key = f"{category}_{seq_tag}"
            stats[key] = stat.text.strip() if stat.text else '0'
        return stats

    def extract_data_from_esummary(self, summary, uid=None):
        """Parse the assembly summary and insert data into the database."""
        try:
            doc_summary = summary['DocumentSummarySet']['DocumentSummary'][0]
        except (KeyError, IndexError):
            return
        return self.extract_data_from_doc_summary(doc_summary, uid=uid)

    def _suppression_reason(self, doc_summary) -> str | None:
        """Return a suppression reason string if the assembly summary indicates suppression."""
        doc_lower = {str(k).lower(): v for k, v in doc_summary.items()} if isinstance(doc_summary, dict) else {}

        def _get_key(key: str):
            if key in doc_summary:
                return doc_summary.get(key)
            return doc_lower.get(str(key).lower())

        props = doc_summary.get("PropertyList") or []
        if not props:
            props = _get_key("PropertyList") or []
        if isinstance(props, (list, tuple)):
            for prop in props:
                token = str(prop).strip().lower()
                if "suppress" in token or "withdrawn" in token or "removed" in token or "killed" in token:
                    return f"PropertyList={prop}"

        status_keys = (
            "Status",
            "AssemblyStatus",
            "AsmStatus",
            "assemblyStatus",
            "RefSeq_category",
            "Suppressed",
            "SuppressionStatus",
            "GBStatus",
            "GenbankStatus",
            "GB_Status",
        )
        for key in status_keys:
            val = _get_key(key)
            if val is None:
                continue
            if isinstance(val, bool):
                if val:
                    return f"{key}=true"
                continue
            text = str(val).strip().lower()
            if not text:
                continue
            if "suppress" in text or "withdrawn" in text or "removed" in text or "killed" in text:
                return f"{key}={val}"
        return None

    def extract_data_from_doc_summary(self, doc_summary, uid=None):
        """Parse a document summary record into assembly/taxonomy/genome datasets."""
        # Allow doc_summary keys in any case (Entrez JSON often lowercases keys)
        if isinstance(doc_summary, dict):
            doc_lower = {str(k).lower(): v for k, v in doc_summary.items()}
        else:
            doc_lower = {}

        def _get_key(key: str):
            if key in doc_summary:
                return doc_summary.get(key)
            return doc_lower.get(str(key).lower())

        suppression_reason = self._suppression_reason(doc_summary)
        if suppression_reason:
            if self.db_manager:
                acc = doc_summary.get("AssemblyAccession")
                status_val = doc_summary.get("AssemblyStatus") or doc_summary.get("assemblyStatus") or doc_summary.get("Status")
                props = doc_summary.get("PropertyList") or []
                if isinstance(props, (list, tuple)):
                    props_val = ",".join([str(p) for p in props if p is not None])
                else:
                    props_val = str(props) if props else None
                if acc:
                    try:
                        self.db_manager.genomes.hide(
                            str(acc),
                            status=props_val or (str(status_val) if status_val is not None else None),
                            reason=suppression_reason,
                        )
                    except Exception as exc:  # boundary: suppressed-assembly marking is best-effort metadata.
                        print(f"Failed to mark suppressed assembly {acc} hidden: {exc}")
            return None
        # Prefer the UID passed from the caller; fall back to what is present in the summary
        uid_val = str(uid) if uid is not None else str(
            _get_key('Uid') or _get_key('Id') or _get_key('UID') or ''
        )

        meta_stats = self.parse_meta(_get_key('Meta') or '')
        def _first_present(keys):
            for key in keys:
                val = _get_key(key)
                if val is not None and str(val).strip() != "":
                    return val
            return None

        assembly_data = {
            'accession': _get_key('AssemblyAccession'),
            'uid': uid_val,
            'assembly_method': _get_key('AssemblyMethod'),
            'assembly_type': _get_key('AssemblyType'),
            'assembly_status': _first_present(("AssemblyStatus", "assemblyStatus", "Status")),
            'origin': infer_ncbi_origin(_get_key('AssemblyAccession')),
            # Prefer GenBank release date; fall back to RefSeq/other dates if absent.
            'release_date': self.parse_date(
                _get_key('AsmReleaseDate_GenBank')
                or _get_key('AsmReleaseDate_RefSeq')
                or _get_key('AsmReleaseDate')
                or _get_key('SubmissionDate')
                or _get_key('LastUpdateDate')
                or _get_key('SeqReleaseDate')
                or _get_key('SeqReleaseDate_RefSeq')
            ),
            'warnings': _get_key('Warnings'),
            'bioproject_accession': (_get_key('GB_BioProjects')[0].get('BioprojectAccn')
                                    if _get_key('GB_BioProjects') else None),
            'biosample_accession': _get_key('BioSampleAccn'),
            'comments': _get_key('AssemblyDescription'),
            'diploid_role': 'alternate pseudohaplotype' if 'diploid' in str(_get_key('AssemblyType') or '') else None,
            'refseq_category': _get_key('RefSeq_category'),
            'sequencing_tech': _get_key('SequencingTech'),
            'submitter': _get_key('SubmitterOrganization'),
            'contig_l50': int(meta_stats.get('contig_l50_all', 0)),
            'contig_n50': int(_get_key('ContigN50') or 0),
            'gc_count': int(meta_stats.get('gc_count_all', 0)) if meta_stats.get('gc_count_all') else None,
            'gc_percent': float(_get_key('GCPercent') or _get_key('GCpercent') or 0) if (_get_key('GCPercent') or _get_key('GCpercent')) else None,
            'genome_coverage': f"{_get_key('Coverage') or ''}x",
            'number_of_component_sequences': int(meta_stats.get('component_count_all', 0)) if meta_stats.get('component_count_all') else None,
            'number_of_contigs': int(meta_stats.get('contig_count_all', 0)),
            'number_of_organelles': int(meta_stats.get('organelle_count_all', 0)) if meta_stats.get('organelle_count_all') else None,
            'number_of_scaffolds': int(meta_stats.get('scaffold_count_all', 0)),
            'scaffold_l50': int(meta_stats.get('scaffold_l50_all', 0)),
            'scaffold_n50': int(_get_key('ScaffoldN50') or 0),
            'total_number_of_chromosomes': int(meta_stats.get('chromosome_count_all', 0)),
            'total_sequence_length': int(meta_stats.get('total_length_all', 0)),
            'total_ungapped_length': int(meta_stats.get('ungapped_length_all', 0)),
        }

        # cursor = db.cursor
        # cursor.execute("""
        #     INSERT OR IGNORE INTO Assembly (
        #         accession, assembly_method, assembly_type, release_date, warnings,
        #         bioproject_accession, biosample_accession, comments, diploid_role,
        #         refseq_category, sequencing_tech, submitter, contig_l50, contig_n50,
        #         gc_count, gc_percent, genome_coverage, number_of_component_sequences,
        #         number_of_contigs, number_of_organelles, number_of_scaffolds, scaffold_l50,
        #         scaffold_n50, total_number_of_chromosomes, total_sequence_length, total_ungapped_length
        #     ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        # """, tuple(assembly_data.values()))

        taxid = int(_get_key('Taxid') or 0)
        tax_name = _get_key('Organism') or _get_key('SpeciesName')
        # cursor.execute("""
        #     INSERT OR IGNORE INTO Taxonomy (taxid, name, rank, parent_taxid)
        #     VALUES (?, ?, ?, ?)
        # """, (taxid, tax_name, "species", None))

        properties = _get_key("PropertyList")
        if properties and isinstance(properties, (list, tuple)):
            properties_val = ",".join([str(p) for p in properties if p is not None])
        else:
            properties_val = None

        genome_data = {
            'accession': _get_key('AssemblyAccession'),
            'taxid': taxid,
            'assembly_level': _first_present(("AssemblyLevel", "assemblyLevel", "AssemblyStatus", "assemblyStatus")),
            'assembly_properties': properties_val,
            'assembly_name': _get_key('AssemblyName'),
            'comments': _get_key('AssemblyDescription'),
            'dl_date': None,
            'location': None,
            'status': 0,
        }
        # cursor.execute("""
        #     INSERT OR IGNORE INTO Genome (
        #         accession, taxid, assembly_level, assembly_name, comments, dl_date, location, status
        #     ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        # """, tuple(genome_data.values()))

        # db.conn.commit()
        return (assembly_data, {"taxid": taxid, "tax_name": tax_name}, genome_data)

    def get_lineage_info(self, taxid_list, batch_size=100):
        """
        For the given list of taxids, queries NCBI taxonomy in batches and updates the Taxonomy table.
        For each taxid, this function extracts hierarchical lineage information for these ranks:
        superkingdom, kingdom, phylum, class, order, family, genus, species.
        The parent taxid for a taxon is set to the next higher rank (closest available) within this order,
        with superkingdom having no parent.
        Instead of inserting duplicate rows, the unique set of taxonomy rows is accumulated and inserted/updated once.
        """
        if not taxid_list:
            return

        desired_ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        # Normalize ids to integers, then unique
        unique_taxids_int = list({int(str(tid)) for tid in taxid_list if str(tid).strip()})

        # If we have a DB, filter to only taxids missing from the local Taxonomy table
        if self.db_manager and unique_taxids_int:
            try:
                placeholders = ','.join(['?'] * len(unique_taxids_int))
                self.db_manager.cursor.execute(
                    f"SELECT taxid FROM Taxonomy WHERE taxid IN ({placeholders})",
                    unique_taxids_int,
                )
                present = {row[0] for row in self.db_manager.cursor.fetchall()}
                missing_int = [tid for tid in unique_taxids_int if tid not in present]
            except Exception as exc:  # boundary: local taxonomy prefilter failure falls back to querying all ids.
                # On any DB issue, fall back to querying all provided ids
                print(f"Failed to prefilter existing taxonomy ids locally: {exc}")
                missing_int = unique_taxids_int
        else:
            missing_int = unique_taxids_int

        # Nothing to fetch from NCBI
        if not missing_int:
            return {}

        unique_taxids = [str(tid) for tid in missing_int]
        taxa_update_dict = {}

        for i in range(0, len(unique_taxids), batch_size):
            print(f"Processing taxonomy batch {i // batch_size + 1} of {len(unique_taxids) // batch_size + 1}")
            batch = unique_taxids[i:i+batch_size]
            try:
                handle = self.safe_entrez_call(Entrez.efetch, db="taxonomy", id=",".join(batch), retmode="xml")
                tax_records = self.safe_entrez_call(Entrez.read, handle)
            except Exception as e:  # boundary: one taxonomy batch failure does not block remaining batches.
                print(f"Error fetching taxonomy records for batch: {e}")
                continue

            for record in tax_records:
                lineage_mapping = {}
                for taxon in record.get("LineageEx", []):
                    rank = taxon.get("Rank", "").lower()
                    if rank in desired_ranks:
                        lineage_mapping[rank] = taxon
                rec_rank = record.get("Rank", "").lower()
                if rec_rank in desired_ranks:
                    lineage_mapping[rec_rank] = record

                final_lineage = {}
                for rank in desired_ranks:
                    if rank in lineage_mapping:
                        final_lineage[rank] = lineage_mapping[rank]

                for idx, rank in enumerate(desired_ranks):
                    if rank not in final_lineage:
                        continue
                    taxon = final_lineage[rank]
                    parent_taxid = None
                    for higher_rank in reversed(desired_ranks[:idx]):
                        if higher_rank in final_lineage:
                            parent_taxid = int(final_lineage[higher_rank].get("TaxId", 0))
                            break
                    taxid_val = int(taxon.get("TaxId", 0))
                    scientific_name = taxon.get("ScientificName", None)
                    taxa_update_dict[taxid_val] = (scientific_name, rank, parent_taxid)

        return taxa_update_dict
        # cursor = db.cursor
        # for taxid_val, (sci_name, rank, parent_taxid) in taxa_update_dict.items():
        #     cursor.execute("""
        #         INSERT OR REPLACE INTO Taxonomy (taxid, name, rank, parent_taxid)
        #         VALUES (?, ?, ?, ?)
        #     """, (taxid_val, sci_name, rank, parent_taxid))
        # db.conn.commit()

    def fetch_assemblies(self, term, accessions_only=False, debug_raw_path=None, debug_full_path=None):
        """Fetch assembly information in batches of 500 and insert into the database with safe retry."""
        Entrez.email = EMAIL
        retmax = 500
        retstart = 0
        all_ids = []

        if accessions_only:
            total_ids = len(term)
            query = ' OR '.join(term)
        else:
            query = f"txid{term}[Organism:exp]"
            try:
                count_handle = self.safe_entrez_call(Entrez.esearch, db="assembly", term=query, retmax=0)
                count_result = self.safe_entrez_call(Entrez.read, count_handle)
                total_ids = int(count_result.get("Count", 0))
            except Exception as e:  # boundary: top-level Entrez count failure aborts this fetch.
                print(f"Error fetching total count: {e}")
                return [], [], {}, []
            if total_ids == 0:
                return [], [], {}, []

        while retstart < total_ids:
            try:
                search_handle = self.safe_entrez_call(
                    Entrez.esearch, db="assembly", term=query, retmax=retmax, retstart=retstart
                )
                search_record = self.safe_entrez_call(Entrez.read, search_handle)
                ids_batch = search_record.get("IdList", [])
            except Exception as e:  # boundary: one Entrez search page failure is skipped.
                print(f"Error fetching ids at retstart {retstart}: {e}")
                retstart += retmax
                continue

            all_ids.extend(ids_batch)
            retstart += retmax

        taxids = []
        assembly_dataset = []
        tax_info_dataset = []
        genome_dataset = []
        raw_handle = None
        full_handle = None
        if debug_raw_path:
            try:
                os.makedirs(os.path.dirname(debug_raw_path), exist_ok=True)
                raw_handle = open(debug_raw_path, "w")
            except OSError as exc:
                print(f"Failed to open raw debug output {debug_raw_path}: {exc}")
                raw_handle = None
        if debug_full_path:
            try:
                os.makedirs(os.path.dirname(debug_full_path), exist_ok=True)
                full_handle = open(debug_full_path, "w")
            except OSError as exc:
                print(f"Failed to open full debug output {debug_full_path}: {exc}")
                full_handle = None

        for id in all_ids:
            try:
                summary = self.get_assembly_summary(id)
            except Exception as e:  # boundary: one assembly summary failure does not block remaining ids.
                print(f"Failed to fetch summary for id {id}: {e}")
                continue
            if raw_handle:
                try:
                    doc_summary = summary.get("DocumentSummarySet", {}).get("DocumentSummary", [None])[0]
                    if doc_summary:
                        raw_handle.write(json.dumps(doc_summary, default=str, ensure_ascii=True) + "\n")
                except (OSError, TypeError, ValueError) as exc:
                    print(f"Failed to write raw debug summary for id {id}: {exc}")
            if full_handle:
                try:
                    full_handle.write(json.dumps(summary, default=str, ensure_ascii=True) + "\n")
                except (OSError, TypeError, ValueError) as exc:
                    print(f"Failed to write full debug summary for id {id}: {exc}")
            # Pass the UID through so it is stored with assembly and genome records
            parsed = self.extract_data_from_esummary(summary, uid=id)
            if not parsed:
                continue
            assembly_data, tax_info, genome_data = parsed
            # self.db_manager.genomes.insert_assembly(assembly_data)
            # self.db_manager.genomes.insert_taxonomy_information({tax_info.get("taxid"): (tax_info.get("tax_name"), "species", None)})
            # self.db_manager.genomes.insert(genome_data)
            taxids.append(tax_info.get("taxid"))
            assembly_dataset.append(assembly_data)
            tax_info_dataset.append(tax_info)
            genome_dataset.append(genome_data)

        # print(f"Found {len(taxids)} taxids from assemblies.")

        taxonomy_update_dict = self.get_lineage_info(taxids)
        if raw_handle:
            try:
                raw_handle.close()
            except OSError as exc:
                print(f"Failed to close raw debug output {debug_raw_path}: {exc}")
        if full_handle:
            try:
                full_handle.close()
            except OSError as exc:
                print(f"Failed to close full debug output {debug_full_path}: {exc}")
        
        return assembly_dataset, tax_info_dataset, taxonomy_update_dict, genome_dataset

    def fetch_assemblies_v2(
        self,
        term,
        accessions_only=False,
        id_chunk_size=200,
        accession_chunk_size=200,
        progress_cb=None,
        debug_raw_path=None,
        debug_full_path=None,
    ):
        """Fetch assembly metadata with chunked ID summaries and optional progress callback."""
        Entrez.email = self.email
        all_ids = []

        if accessions_only:
            if not term:
                return [], [], {}, []
            for i in range(0, len(term), accession_chunk_size):
                chunk = term[i:i + accession_chunk_size]
                query = ' OR '.join(chunk)
                try:
                    search_handle = self.safe_entrez_call(
                        Entrez.esearch, db="assembly", term=query, retmax=len(chunk)
                    )
                    search_record = self.safe_entrez_call(Entrez.read, search_handle)
                    ids_batch = search_record.get("IdList", [])
                except Exception as e:  # boundary: one accession search batch failure does not block remaining batches.
                    print(f"Error fetching ids for accessions {i}-{i + len(chunk)}: {e}")
                    continue
                all_ids.extend(ids_batch)
        else:
            query = f"txid{term}[Organism:exp]"
            try:
                count_handle = self.safe_entrez_call(Entrez.esearch, db="assembly", term=query, retmax=0)
                count_result = self.safe_entrez_call(Entrez.read, count_handle)
                total_ids = int(count_result.get("Count", 0))
            except Exception as e:  # boundary: top-level Entrez count failure aborts this fetch.
                print(f"Error fetching total count: {e}")
                return [], [], {}, []
            if total_ids == 0:
                return [], [], {}, []

            retmax = 500
            retstart = 0
            while retstart < total_ids:
                try:
                    search_handle = self.safe_entrez_call(
                        Entrez.esearch, db="assembly", term=query, retmax=retmax, retstart=retstart
                    )
                    search_record = self.safe_entrez_call(Entrez.read, search_handle)
                    ids_batch = search_record.get("IdList", [])
                except Exception as e:  # boundary: one Entrez search page failure is skipped.
                    print(f"Error fetching ids at retstart {retstart}: {e}")
                    retstart += retmax
                    continue
                all_ids.extend(ids_batch)
                retstart += retmax

        if not all_ids:
            return [], [], {}, []

        raw_handle = None
        full_handle = None
        if debug_raw_path:
            try:
                os.makedirs(os.path.dirname(debug_raw_path), exist_ok=True)
                raw_handle = open(debug_raw_path, "w")
            except OSError as exc:
                print(f"Failed to open raw debug output {debug_raw_path}: {exc}")
                raw_handle = None
        if debug_full_path:
            try:
                os.makedirs(os.path.dirname(debug_full_path), exist_ok=True)
                full_handle = open(debug_full_path, "w")
            except OSError as exc:
                print(f"Failed to open full debug output {debug_full_path}: {exc}")
                full_handle = None

        taxids = []
        assembly_dataset = []
        tax_info_dataset = []
        genome_dataset = []

        total_ids = len(all_ids)
        processed = 0
        for i in range(0, total_ids, id_chunk_size):
            batch = all_ids[i:i + id_chunk_size]
            try:
                summary_handle = self.safe_entrez_call(
                    Entrez.esummary, db="assembly", id=",".join(batch), report="full"
                )
                summary_record = self.safe_entrez_call(Entrez.read, summary_handle)
            except Exception as e:  # boundary: one summary batch failure does not block remaining batches.
                print(f"Failed to fetch summary batch {i}-{i + len(batch)}: {e}")
                continue
            if full_handle:
                try:
                    full_handle.write(json.dumps(summary_record, default=str, ensure_ascii=True) + "\n")
                except (OSError, TypeError, ValueError) as exc:
                    print(f"Failed to write full debug summary batch {i}-{i + len(batch)}: {exc}")

            doc_summaries = summary_record.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            if isinstance(doc_summaries, dict):
                doc_summaries = [doc_summaries]

            for doc_summary in doc_summaries:
                if raw_handle:
                    try:
                        raw_handle.write(json.dumps(doc_summary, default=str, ensure_ascii=True) + "\n")
                    except (OSError, TypeError, ValueError) as exc:
                        print(f"Failed to write raw debug summary: {exc}")
                parsed = self.extract_data_from_doc_summary(doc_summary)
                if not parsed:
                    continue
                assembly_data, tax_info, genome_data = parsed
                taxids.append(tax_info.get("taxid"))
                assembly_dataset.append(assembly_data)
                tax_info_dataset.append(tax_info)
                genome_dataset.append(genome_data)

            processed += len(doc_summaries)
            if progress_cb:
                progress_cb(processed, total_ids)

        taxonomy_update_dict = self.get_lineage_info(taxids)
        if raw_handle:
            try:
                raw_handle.close()
            except OSError as exc:
                print(f"Failed to close raw debug output {debug_raw_path}: {exc}")
        if full_handle:
            try:
                full_handle.close()
            except OSError as exc:
                print(f"Failed to close full debug output {debug_full_path}: {exc}")

        return assembly_dataset, tax_info_dataset, taxonomy_update_dict, genome_dataset

    def download_assembly(self, accession, location, uid=None, protein=False):
        """Download a single assembly from NCBI."""
        # Create the directory if it doesn't exist
        if not os.path.exists(location):
            os.makedirs(location)

        # Get UID if not provided
        if uid is None:
            esearch_handle = self.safe_entrez_call(Entrez.esearch, db="assembly", term=accession, retmax=1)
            esearch_record = self.safe_entrez_call(Entrez.read, esearch_handle)
            id_list = esearch_record.get("IdList", [])
            if not id_list:
                raise Exception(f"No UID found for accession {accession}")
            uid = id_list[0]

        # Fetch the FTP path for the assembly using UID
        esummary_handle = self.safe_entrez_call(Entrez.esummary, db="assembly", id=uid, report="full")
        esummary_record = self.safe_entrez_call(Entrez.read, esummary_handle)
        doc_summary = esummary_record['DocumentSummarySet']['DocumentSummary'][0]
        ftp_path = doc_summary.get('FtpPath_RefSeq') or doc_summary.get('FtpPath_GenBank')

        # print(f"Downloading files for accession {accession} (UID: {uid}) from {ftp_path}")

        if not ftp_path:
            raise Exception(f"No FTP path found for accession {accession}")

        base_name = os.path.basename(ftp_path)
        fna_name = f"{base_name}_genomic.fna.gz"
        faa_name = f"{base_name}_protein.faa.gz"
        gff_name = f"{base_name}_genomic.gff.gz"
        md5_name = "md5checksums.txt"
        fna_url = f"{ftp_path}/{fna_name}"
        faa_url = f"{ftp_path}/{faa_name}"
        gff_url = f"{ftp_path}/{gff_name}"
        md5_url = f"{ftp_path}/{md5_name}"

        # print(faa_url)

        md5_file_path = os.path.join(location, md5_name)
        try:
            _download_verified_file(md5_url, md5_file_path, gzip_file=False)
            md5s = _parse_ncbi_md5checksums(md5_file_path)
        except Exception as e:  # boundary: required NCBI checksum manifest failure is converted to typed download error.
            raise FnaDownloadError(f"Failed to download checksum manifest from {md5_url}: {e}") from e

        fna_file_path = os.path.join(location, fna_name)
        try:
            _download_verified_file(fna_url, fna_file_path, expected_md5=md5s.get(fna_name), gzip_file=True)
        except Exception as e:  # boundary: required FNA download failure is converted to typed download error.
            raise FnaDownloadError(f"Failed to download file from {fna_url} to {fna_file_path}: {e}") from e

        if protein:
            try:
                faa_file_path = os.path.join(location, faa_name)
                _download_verified_file(faa_url, faa_file_path, expected_md5=md5s.get(faa_name), gzip_file=True)
                
            except Exception as e:  # boundary: required protein download failure is converted to typed download error.
                raise FaaDownloadError(f"Failed to download protein files for {accession}: {e}") from e

            try:
                gff_file_path = os.path.join(location, gff_name)
                _download_verified_file(gff_url, gff_file_path, expected_md5=md5s.get(gff_name), gzip_file=True)
            except Exception as e:  # boundary: required GFF download failure is converted to typed download error.
                raise GffDownloadError(f"Failed to download GFF files for {accession}: {e}") from e

        return {"ftp_path": ftp_path, "md5_path": md5_file_path, "md5": md5s}

        

    def download_assemblies(self, accessions, location, protein=False):
        """Download assemblies from NCBI using the provided accessions."""
        # Max 10 per batch
        if not os.path.exists(location):
            os.makedirs(location)

        # # Use entrez to fetch the .fna files and store in a folder of the accession name
        # batch_size = 10
        # for i in range(0, len(accessions), batch_size):
        #     batch = accessions[i:i + batch_size]
        for accession in accessions:
            # try:
            # Fetch the FTP path for the assembly
            esummary_handle = self.safe_entrez_call(Entrez.esummary, db="assembly", id=accession, report="full")
            # print("Got handle")
            esummary_record = self.safe_entrez_call(Entrez.read, esummary_handle)
            # print("Got record")
            doc_summary = esummary_record['DocumentSummarySet']['DocumentSummary'][0]
            ftp_path = doc_summary.get('FtpPath_GenBank') or doc_summary.get('FtpPath_RefSeq')

            # print(f"Downloading files for accession {accession} from {ftp_path}")

            if not ftp_path:
                raise Exception(f"No FTP path found for accession {accession}")
                # continue

            # Construct file URLs
            base_name = os.path.basename(ftp_path)
            fna_name = f"{base_name}_genomic.fna.gz"
            faa_name = f"{base_name}_protein.faa.gz"
            gff_name = f"{base_name}_genomic.gff.gz"
            fna_url = f"{ftp_path}/{fna_name}"
            faa_url = f"{ftp_path}/{faa_name}"
            gff_url = f"{ftp_path}/{gff_name}"
            md5_url = f"{ftp_path}/md5checksums.txt"

            md5_file_path = os.path.join(location, "md5checksums.txt")
            try:
                _download_verified_file(md5_url, md5_file_path, gzip_file=False)
                md5s = _parse_ncbi_md5checksums(md5_file_path)
            except Exception as e:  # boundary: required NCBI checksum manifest failure is converted to typed download error.
                raise FnaDownloadError(f"Failed to download checksum manifest from {md5_url}: {e}") from e

            # Download .fna file
            fna_file_path = os.path.join(location, fna_name)
            try:
                _download_verified_file(fna_url, fna_file_path, expected_md5=md5s.get(fna_name), gzip_file=True)
            except Exception as e:  # boundary: required FNA download failure is converted to typed download error.
                raise FnaDownloadError(f"Failed to download file from {fna_url} to {fna_file_path}: {e}") from e
            # print(f"Downloaded: {fna_file_path}")

            if protein:
                try:
                # Download .faa file
                    faa_file_path = os.path.join(location, faa_name)
                    _download_verified_file(faa_url, faa_file_path, expected_md5=md5s.get(faa_name), gzip_file=True)
                    # print(f"Downloaded: {faa_file_path}")
                except Exception as e:  # boundary: required protein download failure is converted to typed download error.
                    raise FaaDownloadError(f"Failed to download protein files for {accession}: {e}") from e

                try:
                    # Download .gff file
                    gff_file_path = os.path.join(location, gff_name)
                    _download_verified_file(gff_url, gff_file_path, expected_md5=md5s.get(gff_name), gzip_file=True)
                    # print(f"Downloaded: {gff_file_path}")
                except Exception as e:  # boundary: required GFF download failure is converted to typed download error.
                    raise GffDownloadError(f"Failed to download protein files for {accession}: {e}") from e
            # except Exception as e:
            #     print(f"Failed to download files for accession {accession}: {e}")

        return True
