# Command Cheat Sheet

When to read this page: read this when you know the workflow already and just need the most common commands in one place.

This is a task-oriented cheat sheet, not a full command reference.

## Create a project

```bash
phyloODB mammal_workshop.db create
```

## Register metadata

```bash
phyloODB mammal_workshop.db run add --clade Primates
phyloODB mammal_workshop.db run add --clade Rodentia
```

## Inspect candidate assemblies

```bash
phyloODB mammal_workshop.db list assemblies --clade Primates --tidy
phyloODB mammal_workshop.db count assemblies --clade Primates
phyloODB mammal_workshop.db list assemblies --clade Primates --level chromosome --tidy
```

## Build reusable panels

```bash
phyloODB mammal_workshop.db assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_REFS
phyloODB mammal_workshop.db assemblies --clade Primates --ranks family,genus --quantities 2,1 --store PRIMATE_TARGETS
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_REFS --tidy
```

## Download assemblies

```bash
phyloODB mammal_workshop.db run download --accessions @PRIMATE_REFS --protein
phyloODB mammal_workshop.db run verify-downloads --accessions @PRIMATE_REFS --downloaded-only
phyloODB mammal_workshop.db run clean-isoforms --accessions @PRIMATE_REFS --downloaded-only
```

## Download the BUSCO lineage

```bash
phyloODB mammal_workshop.db run download-busco-library --lineage mammalia_odb12 --coverage 1
```

## Run BUSCO

```bash
phyloODB mammal_workshop.db run batch-busco --accessions @PRIMATE_REFS --lineage mammalia_odb12 --format protein
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_REFS --library-name mammalia_odb12 --busco --tidy
```

## Build a derived library

```bash
phyloODB mammal_workshop.db run add-library \
  --name mammal_core_odb12 \
  --coverage Mammalia \
  --parent-library-name mammalia_odb12 \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS
```

## Run hidden-paralog filtering

```bash
phyloODB mammal_workshop.db run paralog-removal \
  --library-name mammal_core_odb12 \
  --ref-accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --targets @MAMMAL_TARGETS \
  --report-dir wiki_example/paralog_removal
```

## Run contamination screening

```bash
phyloODB mammal_workshop.db run internal-decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --rank order \
  --hit-window 8 \
  --off-clade-fraction 0.05 \
  --report-path wiki_example/internal_decontamination/mammal_internal
```

## Quick export

```bash
phyloODB mammal_workshop.db run export \
  --library-name mammalia_odb12 \
  --accessions @PRIMATE_REFS \
  --out-dir wiki_example/quick_export \
  --disable-paralog-filter \
  --disable-decont-filter
```

## Strict final export

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

## Good pages to keep nearby

- [Home](Home)
- [Choosing Taxa and Building Panels](Choosing-Taxa-and-Building-Panels)
- [Running BUSCO](Running-BUSCO)
- [Exporting a Dataset](Exporting-a-Dataset)
