import argparse
from collections import defaultdict
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, Optional

SILENT = False

def parse_gff_attributes(attribute_string):
    attributes = {}
    for field in attribute_string.strip().split(';'):
        if not field or '=' not in field:
            continue
        key, value = field.split('=', 1)
        attributes[key] = value
    return attributes

def _protein_id_aliases(protein_id):
    aliases = set()
    if not protein_id:
        return aliases

    token = protein_id.strip()
    if not token:
        return aliases

    aliases.add(token)
    parts = [part.strip() for part in token.split("|") if part.strip()]
    aliases.update(parts)
    if len(parts) >= 2:
        aliases.add(parts[-1])
        aliases.add(parts[1])

    # Add versionless aliases (e.g. XP_123.1 -> XP_123) for tolerant matching.
    for alias in list(aliases):
        stripped = re.sub(r"\.[0-9]+$", "", alias)
        if stripped and stripped != alias:
            aliases.add(stripped)

    return aliases

def normalize_feature_id(feature_id):
    if not feature_id:
        return feature_id
    for prefix in ('gene-', 'rna-', 'cds-'):
        if feature_id.startswith(prefix):
            return feature_id[len(prefix):]
    return feature_id

def pick_gene_identifier(attributes, default=None):
    for key in ('gene', 'locus_tag', 'Name'):
        if attributes.get(key):
            return attributes[key]
    if attributes.get('ID'):
        return normalize_feature_id(attributes['ID'])
    return default

def parse_gff(gff_file):
    gene_to_proteins = defaultdict(set)
    gene_id_to_label = {}
    transcript_to_gene = {}

    with open(gff_file, 'r') as file:
        if not SILENT:
            print("Searching",gff_file,"for gene info...")
        for line in file:
            if line.startswith('#'):
                continue
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            feature_type = columns[2]
            attributes = parse_gff_attributes(columns[8])

            if feature_type == 'gene':
                if attributes.get('ID'):
                    gene_id_to_label[attributes['ID']] = pick_gene_identifier(attributes)
                continue

            if feature_type in ('mRNA', 'transcript'):
                transcript_id = attributes.get('ID')
                if not transcript_id:
                    continue
                gene_label = None
                parent_field = attributes.get('Parent')
                if parent_field:
                    for parent in parent_field.split(','):
                        parent = parent.strip()
                        if not parent:
                            continue
                        gene_label = gene_id_to_label.get(parent, normalize_feature_id(parent))
                        if gene_label:
                            break
                if not gene_label:
                    gene_label = pick_gene_identifier(attributes)
                if gene_label:
                    transcript_to_gene[transcript_id] = gene_label
                continue

            if feature_type != 'CDS':
                continue

            protein_id = attributes.get('protein_id')
            if not protein_id:
                continue

            gene_id = attributes.get('gene') or attributes.get('locus_tag')
            if not gene_id:
                parent_field = attributes.get('Parent')
                if parent_field:
                    for parent in parent_field.split(','):
                        parent = parent.strip()
                        if not parent:
                            continue
                        mapped_gene = transcript_to_gene.get(parent)
                        if mapped_gene:
                            gene_id = mapped_gene
                            break
            if not gene_id:
                gene_id = pick_gene_identifier(attributes)
            if gene_id:
                gene_to_proteins[gene_id].add(protein_id)
    return gene_to_proteins

def parse_faa(faa_file):
    proteins = {}
    protein_lookup = {}
    with open(faa_file, 'r') as file:
        protein_id = None
        sequence = []
        for line in file:
            if line.startswith('>'):
                if protein_id:
                    proteins[protein_id] = ''.join(sequence)
                    for alias in _protein_id_aliases(protein_id):
                        protein_lookup.setdefault(alias, protein_id)
                protein_id = line.split()[0][1:]
                # print(protein_id)
                sequence = []
            else:
                sequence.append(line.strip())
        if protein_id:
            proteins[protein_id] = ''.join(sequence)
            for alias in _protein_id_aliases(protein_id):
                protein_lookup.setdefault(alias, protein_id)
    return proteins, protein_lookup

def count_faa_headers(faa_file):
    count = 0
    with open(faa_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                count += 1
    return count

def filter_proteins(gene_to_proteins, proteins, protein_lookup):
    filtered_proteins = {}
    for gene_id, protein_ids in gene_to_proteins.items():
        # Gene id and then a list of the proteins corrosponding to that gene
        present_proteins = [pid for pid in protein_ids if pid in protein_lookup]
        if not present_proteins:
            continue
        longest_protein = max(
            present_proteins,
            key=lambda pid: len(proteins.get(protein_lookup.get(pid, ""), "")),
        )
        longest_faa_id = protein_lookup.get(longest_protein)
        if not longest_faa_id:
            continue
        if not SILENT:
            print("longest protein for", gene_id, "is", longest_faa_id)
        filtered_proteins[longest_faa_id] = proteins[longest_faa_id]
    return filtered_proteins

def write_faa(filtered_proteins, output_file, prefix):
    with open(output_file, 'w') as file:
        for protein_id, sequence in filtered_proteins.items():
            if prefix:
                file.write(f'>{prefix}_{protein_id}\n')
            else:
                file.write(f'>{protein_id}\n')
            for i in range(0, len(sequence), 60):
                file.write(sequence[i:i+60] + '\n')

def filter_isoforms_using_gff(faa_file, gff_file, output_file, prefix=None, silent=False):
    global SILENT
    SILENT = silent

    if not os.path.exists(faa_file):
        print(f"Error: The file {faa_file} does not exist.")
        return

    if not os.path.exists(gff_file):
        print(f"Error: The file {gff_file} does not exist.")
        return

    proteins, protein_lookup = parse_faa(faa_file)

    if not proteins:
        print("Failure could not parse protein file.")
    # print(proteins)

    gene_to_proteins = parse_gff(gff_file)

    if not gene_to_proteins:
        print("Failure could not parse genes from .gff file")
        return
    # print(gene_to_proteins)

    filtered_proteins = filter_proteins(gene_to_proteins, proteins, protein_lookup)
    write_faa(filtered_proteins, output_file, prefix)
    return output_file

def filter_isoforms_using_cdhit(faa_file, output_file, identity=0.95, prefix=None, silent=False, threads=1):
    global SILENT
    SILENT = silent

    if identity <= 0 or identity > 1.0:
        raise ValueError("CD-HIT identity must be between 0 and 1.0.")

    if not os.path.exists(faa_file):
        raise FileNotFoundError(f"CD-HIT input file does not exist: {faa_file}")

    # Run CD-HIT to cluster sequences at the specified identity
    command = [
        "cd-hit",
        "-i", faa_file,
        "-o", output_file,
        "-c", str(identity),
        "-n", "5",
        "-d", "0",
        "-T", str(max(1, int(threads or 1))),
    ]
    try:
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f"Unable to run CD-HIT executable 'cd-hit': {exc}") from exc

    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        suffix = f": {detail}" if detail else ""
        raise RuntimeError(f"CD-HIT failed with exit code {result.returncode}{suffix}")

    if not os.path.isfile(output_file):
        raise RuntimeError(f"CD-HIT completed without creating its output file: {output_file}")

    # Optionally rename headers in the output file
    if prefix:
        with open(output_file, 'r') as infile, open(output_file + '.tmp', 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    protein_id = line[1:].strip()
                    outfile.write(f'>{prefix}_{protein_id}\n')
                else:
                    outfile.write(line)
        os.replace(output_file + '.tmp', output_file)

    return output_file

def _find_primary_faa(genome_path: str) -> Optional[str]:
    faa_candidates = []
    faa_gz_candidates = []
    for fname in os.listdir(genome_path):
        full_path = os.path.join(genome_path, fname)
        if not os.path.isfile(full_path):
            continue
        if fname.endswith(".faa"):
            faa_candidates.append(full_path)
        elif fname.endswith(".faa.gz"):
            faa_gz_candidates.append(full_path)
    if faa_candidates:
        return sorted(faa_candidates)[0]
    if faa_gz_candidates:
        return sorted(faa_gz_candidates)[0]
    return None

def _find_gff_file(genome_path: str, temp_dir: str) -> Optional[str]:
    gff_candidates = []
    gff_gz_candidates = []
    for fname in os.listdir(genome_path):
        full_path = os.path.join(genome_path, fname)
        if not os.path.isfile(full_path):
            continue
        if fname.endswith(".gff") or fname.endswith(".gff3"):
            gff_candidates.append(full_path)
        elif fname.endswith(".gff.gz"):
            gff_gz_candidates.append(full_path)
    if gff_candidates:
        return sorted(gff_candidates)[0]
    if gff_gz_candidates:
        gz_file = sorted(gff_gz_candidates)[0]
        extracted_path = os.path.join(temp_dir, os.path.basename(gz_file[:-3]))
        with gzip.open(gz_file, 'rb') as f_in, open(extracted_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return extracted_path
    return None

def _prepare_working_faa(faa_file: str, temp_dir: str) -> str:
    working_faa = os.path.join(temp_dir, "input.faa")
    if faa_file.endswith(".gz"):
        with gzip.open(faa_file, 'rb') as f_in, open(working_faa, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy2(faa_file, working_faa)
    return working_faa

def _archive_path_for_faa(faa_file: str) -> str:
    if faa_file.endswith(".faa.gz"):
        base_faa = faa_file[:-3]
    elif faa_file.endswith(".faa"):
        base_faa = faa_file
    else:
        base_faa = faa_file + ".faa"
    return base_faa + ".archive.gz"

def _find_archive_faa(genome_path: str) -> Optional[str]:
    archives = sorted(
        os.path.join(genome_path, fname)
        for fname in os.listdir(genome_path)
        if fname.endswith(".faa.archive.gz") or fname.endswith(".faa.archive")
    )
    if not archives:
        return None
    return archives[0]

def _normalize_archive_to_gz(archive_file: str) -> str:
    """Convert legacy .faa.archive files to .faa.archive.gz in place."""
    if archive_file.endswith(".archive.gz"):
        return archive_file
    if not archive_file.endswith(".archive"):
        return archive_file

    gz_archive = archive_file + ".gz"
    if os.path.exists(gz_archive):
        try:
            os.remove(archive_file)
        except OSError as exc:
            if not SILENT:
                print(f"Warning: failed to remove intermediate archive {archive_file}: {exc}")
        return gz_archive

    _write_gzip_from_plain(archive_file, gz_archive)
    try:
        os.remove(archive_file)
    except OSError as exc:
        if not SILENT:
            print(f"Warning: failed to remove intermediate archive {archive_file}: {exc}")
    return gz_archive

def _write_gzip_from_plain(src_plain: str, dst_gz: str) -> None:
    with open(src_plain, 'rb') as f_in, gzip.open(dst_gz, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def _archive_original_faa(faa_file: str) -> str:
    archive_path = _archive_path_for_faa(faa_file)
    if faa_file.endswith(".gz"):
        shutil.copy2(faa_file, archive_path)
    else:
        _write_gzip_from_plain(faa_file, archive_path)
    return archive_path

def _write_cleaned_to_canonical_paths(faa_file: str, cleaned_plain_file: str) -> Dict[str, str]:
    written = {}
    directory = os.path.dirname(faa_file)
    basename = os.path.basename(faa_file)

    if faa_file.endswith(".faa.gz"):
        canonical_gz = faa_file
        _write_gzip_from_plain(cleaned_plain_file, canonical_gz)
        written["canonical_gz"] = canonical_gz

        sibling_plain = os.path.join(directory, basename[:-3])  # strip .gz
        if os.path.exists(sibling_plain):
            shutil.copy2(cleaned_plain_file, sibling_plain)
            written["sibling_plain"] = sibling_plain
    else:
        canonical_plain = faa_file
        shutil.copy2(cleaned_plain_file, canonical_plain)
        written["canonical_plain"] = canonical_plain

        sibling_gz = canonical_plain + ".gz"
        if os.path.exists(sibling_gz):
            _write_gzip_from_plain(cleaned_plain_file, sibling_gz)
            written["sibling_gz"] = sibling_gz
    return written

def revert_proteome_from_archive(genome_path: str) -> Dict[str, object]:
    archive_file = _find_archive_faa(genome_path)
    if not archive_file:
        return {"ok": False, "status": "no_archive"}
    archive_file = _normalize_archive_to_gz(archive_file)
    if archive_file.endswith(".archive.gz"):
        base_faa = archive_file[:-11]  # strip .archive.gz
        target_path = base_faa + ".gz"
        shutil.copy2(archive_file, target_path)
    elif archive_file.endswith(".archive"):
        base_faa = archive_file[:-8]  # strip .archive
        target_path = base_faa
        shutil.copy2(archive_file, target_path)
    else:
        return {"ok": False, "status": "unsupported_archive_suffix", "archive": archive_file}

    plain_target = base_faa
    if target_path != plain_target and os.path.exists(plain_target):
        try:
            os.remove(plain_target)
        except OSError as exc:
            if not SILENT:
                print(f"Warning: failed to remove decompressed target {plain_target}: {exc}")

    return {
        "ok": True,
        "status": "reverted",
        "archive": archive_file,
        "restored": target_path,
    }

def clean_proteome_in_genome_path(
    genome_path: str,
    skip_clean_isoforms: bool = False,
    skip_gff: bool = False,
    skip_cdhit: bool = False,
    gff_priority: bool = False,
    cdhit_identity: float = 0.96,
    cdhit_threads: int = 1,
    silent: bool = True,
    revert_from_archive: bool = False,
    force_reclean: bool = False,
) -> Dict[str, object]:
    if revert_from_archive:
        return revert_proteome_from_archive(genome_path)

    if skip_clean_isoforms:
        return {"ok": True, "status": "skipped_cleaning"}

    existing_archive = _find_archive_faa(genome_path)
    if existing_archive:
        existing_archive = _normalize_archive_to_gz(existing_archive)
        if not force_reclean:
            return {
                "ok": True,
                "status": "skipped_already_cleaned",
                "archive": existing_archive,
            }
        reverted = revert_proteome_from_archive(genome_path)
        if not reverted.get("ok"):
            return {
                "ok": False,
                "status": "failed_revert_before_force_reclean",
                "archive": existing_archive,
            }

    faa_file = _find_primary_faa(genome_path)
    if not faa_file:
        return {"ok": False, "status": "no_faa"}

    with tempfile.TemporaryDirectory(prefix="clean_isoforms_") as temp_dir:
        working_faa = _prepare_working_faa(faa_file, temp_dir)
        input_count = count_faa_headers(working_faa)
        current_file = working_faa

        used_gff = False
        gff_reduced = False
        gff_file = None
        gff_input_count = None
        gff_output_count = None
        gff_applied = False
        if not skip_gff:
            try:
                gff_file = _find_gff_file(genome_path, temp_dir)
            except (OSError, EOFError, gzip.BadGzipFile):
                gff_file = None
            if gff_file:
                gff_input_count = count_faa_headers(current_file)
                gff_out = os.path.join(temp_dir, "filtered_by_gff.faa")
                filtered = filter_isoforms_using_gff(
                    current_file,
                    gff_file,
                    output_file=gff_out,
                    prefix=None,
                    silent=silent,
                )
                if filtered and os.path.exists(filtered):
                    used_gff = True
                    gff_count = count_faa_headers(filtered)
                    if gff_count > 0 and gff_input_count is not None and gff_count <= gff_input_count:
                        current_file = filtered
                        gff_output_count = gff_count
                        gff_applied = True
                        if gff_count < gff_input_count:
                            gff_reduced = True

        used_cdhit = False
        cdhit_input_count = None
        cdhit_output_count = None
        cdhit_skipped_due_gff_priority = bool(gff_priority and gff_reduced)
        if not skip_cdhit and not cdhit_skipped_due_gff_priority:
            cdhit_input_count = count_faa_headers(current_file)
            cdhit_out = os.path.join(temp_dir, "filtered_by_cdhit.faa")
            clustered = filter_isoforms_using_cdhit(
                current_file,
                output_file=cdhit_out,
                identity=cdhit_identity,
                prefix=None,
                silent=silent,
                threads=cdhit_threads,
            )
            if clustered and os.path.exists(clustered):
                used_cdhit = True
                cdhit_count = count_faa_headers(clustered)
                if cdhit_count > 0:
                    cdhit_output_count = cdhit_count
                    current_file = clustered

        output_count = count_faa_headers(current_file)
        if output_count <= 0:
            return {"ok": False, "status": "empty_output", "faa": faa_file}

        archive_path = _archive_original_faa(faa_file)
        written_paths = _write_cleaned_to_canonical_paths(faa_file, current_file)
        gff_removed = (
            max(0, int(gff_input_count) - int(gff_output_count))
            if gff_applied and gff_input_count is not None and gff_output_count is not None
            else 0
        )
        cdhit_removed = (
            max(0, int(cdhit_input_count) - int(cdhit_output_count))
            if used_cdhit and cdhit_input_count is not None and cdhit_output_count is not None
            else 0
        )

        return {
            "ok": True,
            "status": "cleaned",
            "faa": faa_file,
            "archive": archive_path,
            "written": written_paths,
            "input_count": input_count,
            "output_count": output_count,
            "used_gff": used_gff,
            "used_cdhit": used_cdhit,
            "gff_reduced": gff_reduced,
            "gff_file": gff_file,
            "gff_input_count": gff_input_count,
            "gff_output_count": gff_output_count,
            "gff_removed": gff_removed,
            "cdhit_input_count": cdhit_input_count,
            "cdhit_output_count": cdhit_output_count,
            "cdhit_removed": cdhit_removed,
            "cdhit_skipped_due_gff_priority": cdhit_skipped_due_gff_priority,
            "total_removed": max(0, int(input_count) - int(output_count)),
        }


def prepare_proteome_profile(
    source_faa: str,
    genome_path: str,
    output_path: str,
    *,
    skip_gff: bool = False,
    skip_cdhit: bool = False,
    gff_priority: bool = False,
    cdhit_identity: float = 0.96,
    cdhit_threads: int = 1,
    silent: bool = True,
) -> Dict[str, object]:
    if not source_faa or not os.path.exists(source_faa):
        return {"ok": False, "status": "no_source_faa", "faa": source_faa}

    with tempfile.TemporaryDirectory(prefix="prepare_proteome_") as temp_dir:
        working_faa = _prepare_working_faa(source_faa, temp_dir)
        input_count = count_faa_headers(working_faa)
        current_file = working_faa

        used_gff = False
        gff_reduced = False
        gff_file = None
        gff_input_count = None
        gff_output_count = None
        gff_applied = False
        if not skip_gff:
            try:
                gff_file = _find_gff_file(genome_path, temp_dir)
            except (OSError, EOFError, gzip.BadGzipFile):
                gff_file = None
            if gff_file:
                gff_input_count = count_faa_headers(current_file)
                gff_out = os.path.join(temp_dir, "filtered_by_gff.faa")
                filtered = filter_isoforms_using_gff(
                    current_file,
                    gff_file,
                    output_file=gff_out,
                    prefix=None,
                    silent=silent,
                )
                if filtered and os.path.exists(filtered):
                    used_gff = True
                    gff_count = count_faa_headers(filtered)
                    if gff_count > 0 and gff_input_count is not None and gff_count <= gff_input_count:
                        current_file = filtered
                        gff_output_count = gff_count
                        gff_applied = True
                        if gff_count < gff_input_count:
                            gff_reduced = True

        used_cdhit = False
        cdhit_input_count = None
        cdhit_output_count = None
        cdhit_skipped_due_gff_priority = bool(gff_priority and gff_reduced)
        if not skip_cdhit and not cdhit_skipped_due_gff_priority:
            cdhit_input_count = count_faa_headers(current_file)
            cdhit_out = os.path.join(temp_dir, "filtered_by_cdhit.faa")
            clustered = filter_isoforms_using_cdhit(
                current_file,
                output_file=cdhit_out,
                identity=cdhit_identity,
                prefix=None,
                silent=silent,
                threads=cdhit_threads,
            )
            if clustered and os.path.exists(clustered):
                used_cdhit = True
                cdhit_count = count_faa_headers(clustered)
                if cdhit_count > 0:
                    cdhit_output_count = cdhit_count
                    current_file = clustered

        output_count = count_faa_headers(current_file)
        if output_count <= 0:
            return {"ok": False, "status": "empty_output", "faa": source_faa}

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        if str(output_path).lower().endswith(".gz"):
            _write_gzip_from_plain(current_file, output_path)
        else:
            shutil.copy2(current_file, output_path)

        gff_removed = (
            max(0, int(gff_input_count) - int(gff_output_count))
            if gff_applied and gff_input_count is not None and gff_output_count is not None
            else 0
        )
        cdhit_removed = (
            max(0, int(cdhit_input_count) - int(cdhit_output_count))
            if used_cdhit and cdhit_input_count is not None and cdhit_output_count is not None
            else 0
        )

        return {
            "ok": True,
            "status": "prepared",
            "faa": source_faa,
            "written": output_path,
            "input_count": input_count,
            "output_count": output_count,
            "used_gff": used_gff,
            "used_cdhit": used_cdhit,
            "gff_reduced": gff_reduced,
            "gff_file": gff_file,
            "gff_input_count": gff_input_count,
            "gff_output_count": gff_output_count,
            "gff_removed": gff_removed,
            "cdhit_input_count": cdhit_input_count,
            "cdhit_output_count": cdhit_output_count,
            "cdhit_removed": cdhit_removed,
            "cdhit_skipped_due_gff_priority": cdhit_skipped_due_gff_priority,
            "total_removed": max(0, int(input_count) - int(output_count)),
        }

def main():
    parser = argparse.ArgumentParser(description='Filter proteins to retain only the longest one per gene.')
    parser.add_argument('-f', '--faa', required=True, help='Path to the protein .faa file')
    parser.add_argument('-gf', '--gff', required=True, help='Path to the genomic .gff file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output .faa file')
    parser.add_argument('-r', '--rename', required=False, help='Rename fasta headers to Genus_species_accession where argument is Genus_species')
    parser.add_argument('-s', '--silent', action='store_true', help='Suppress output messages')
    args = parser.parse_args()

    global SILENT
    SILENT = args.silent

    if args.rename:
        prefix = args.rename
    else:
        prefix = None

    print("Input Proteome:", args.faa)
    print("Input Genome Annotation File:", args.gff)
    print("Output location:", args.output)

    if not os.path.exists(args.faa):
        print(f"Error: The file {args.faa} does not exist.")
        return

    if not os.path.exists(args.gff):
        print(f"Error: The file {args.gff} does not exist.")
        return

    proteins, protein_lookup = parse_faa(args.faa)

    if not proteins:
        print("Failure could not parse protein file.")
    # print(proteins)

    gene_to_proteins = parse_gff(args.gff)

    if not gene_to_proteins:
        print("Failure could not parse genes from .gff file")
        return
    # print(gene_to_proteins)

    filtered_proteins = filter_proteins(gene_to_proteins, proteins, protein_lookup)
    write_faa(filtered_proteins, args.output, prefix)

if __name__ == '__main__':
    main()
