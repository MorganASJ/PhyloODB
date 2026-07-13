# Exporting a Dataset

When to read this page: read this when you are ready to turn project state into actual dataset files.

There are two export modes users usually care about:

- a quick preview export
- a study-quality export that keeps and documents the active filtering state

## Quick preview export

If you want a deliberately simple exploratory dataset from a parent BUSCO lineage, you can tell export to ignore paralog and decontamination filtering even when those results exist:

```bash
phyloODB mammal_workshop.db run export \
  --library-name mammalia_odb12 \
  --accessions @PRIMATE_PANEL \
  --out-dir wiki_example/quick_export \
  --disable-paralog-filter \
  --disable-decont-filter
```

Use this when:

- you are learning the system
- you want a provisional dataset quickly
- you have BUSCO results but have not yet built a derived library or filtering runs

Do not mistake this for the recommended final study path.

## Study-quality export

Once you have a derived library plus filtering results, use the normal export path:

```bash
phyloODB mammal_workshop.db run export \
  --library-name mammal_core_odb12 \
  --accessions @MAMMAL_TARGETS \
  --out-dir wiki_example/final_export \
  --write-lineage-csv \
  --write-busco-report \
  --write-busco-family-matrix \
  --busco-report-extended
```

You do not need a separate “strict mode” flag for this. Export uses paralog-filtering and decontamination results automatically when they are available for the chosen library and taxa.

## What the filter flags really mean

- `--disable-paralog-filter` means “ignore paralog-filtering results during export, even if they exist”
- `--disable-decont-filter` means “ignore decontamination results during export, even if they exist”
- `--require-paralog-filtering` means “fail rather than export if paralog-filtering results are missing”
- `--require-decontamination` means “fail rather than export if decontamination results are missing”

The default behaviour is:

- if filtering results exist, export uses them
- if they do not exist, export still runs
- if you want missing filtering to be an error, say so explicitly with `--require-paralog-filtering` and/or `--require-decontamination`

That default is intentional. It keeps exploratory work easy without hiding the fact that study-quality work should usually include those filters.

## What export writes

Every export creates a directory containing at least:

- `busco_families/`: one FASTA per retained BUSCO family after export-side filtering and occupancy checks
- `export_filter_report.tsv`: per-family summary of what was kept, removed, and why
- `export_parameters.txt`: resolved export settings, selected BUSCO runs, and report paths
- `export_task.log`: copy of the task log for that export run

Depending on flags, export can also write:

- `lineage.csv`: accession-to-lineage table for the retained taxa
- `busco_report.tsv`: per-accession BUSCO summary in the export context
- `taxa_occupancy.tsv`: which taxa survived the minimum taxa-occupancy rule
- `busco_family_matrix.tsv`: per-accession, per-family status matrix

These files are part of the dataset definition. They explain what the exported FASTAs include and what was filtered out.

## How to interpret the output

You should expect:

- family FASTAs that already reflect the chosen library, available filtering results, required-clause filters, and occupancy thresholds
- report files that explain the retained taxa, retained families, and removal reasons
- sequence sourcing from the selected BUSCO runs, not just whatever FASTA happens to be present on disk

If the export is unexpectedly small, the most common causes are:

- the target panel is too strict
- BUSCO thresholds were too strict
- available paralog filtering or decontamination removed more sequences than expected
- occupancy filtering removed more taxa or families than expected

## Gene trees from an export

If you want gene trees for exactly the exported families, use `build-busco-trees`. That workflow writes:

- `alignments/`: one MAFFT alignment per exported family
- `trees/`: one IQ-TREE result directory per exported family
- `manifest.tsv`: table linking each family to its raw FASTA, alignment, and final tree file

This is different from the gene-tree work inside `add-library`. `build-busco-trees` operates on the final exported family set rather than on the broader orthogroup analysis used to define a derived library.

## Practical pattern

Use this progression:

1. Quick export if you need a preview.
2. Derived library plus filtering if the study matters.
3. Export with the default opportunistic filtering behaviour, or add explicit `--require-*` flags if you want export to refuse unfinished state.

## Next step

If you want the whole workflow in one place, read [End-to-End Example Mammal Dataset](End-to-End-Example-Mammal-Dataset).
