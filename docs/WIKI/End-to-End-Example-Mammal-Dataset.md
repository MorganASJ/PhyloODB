# End-to-End Example Mammal Dataset

When to read this page: read this when you want one coherent copy-pasteable example from database creation through export.

This example builds a primate-rich mammal dataset with a curated reference panel and a broader target panel. It mirrors the real workflow used in the repository tutorial, but it is written for GitHub wiki use and emphasizes the practical decisions a user needs to make.

## Goal

We want:

- a small mammalian reference panel for library construction
- a broader primate-rich target panel for export
- a final export based on a derived library, paralog filtering, and contamination screening

## 1. Create the project

```bash
phyloODB mammal_workshop.db create
```

## 2. Register metadata for the focal clades

```bash
phyloODB mammal_workshop.db run add --clade Primates
phyloODB mammal_workshop.db run add --clade Rodentia
phyloODB mammal_workshop.db run add --clade Carnivora
phyloODB mammal_workshop.db run add --clade Artiodactyla
phyloODB mammal_workshop.db run add --clade Marsupialia
phyloODB mammal_workshop.db run add --clade Monotremata
```

This stage only registers what assemblies exist. It does not download files yet.

## 3. Register the BUSCO lineage

```bash
phyloODB mammal_workshop.db run download-busco-library --lineage mammalia_odb12 --coverage 1
```

## 4. Inspect the primate candidate pool

```bash
phyloODB mammal_workshop.db list assemblies --clade Primates --rank family --quantity 1 --tidy
phyloODB mammal_workshop.db list assemblies --clade Primates --ranks family,genus --quantities 2,1 --tidy
phyloODB mammal_workshop.db list assemblies --clade Primates --level chromosome --tidy
```

## 5. Store the reference and target panels

Reference-oriented primates:

```bash
phyloODB mammal_workshop.db assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_REFS
```

Mammalian outgroups:

```bash
phyloODB mammal_workshop.db assemblies --clade Rodentia --rank order --quantity 1 --store RODENT_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Carnivora --rank family --quantity 1 --store CARNIVORE_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Artiodactyla --rank family --quantity 1 --store UNGULATE_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Marsupialia --quantity 1 --store MARSUPIAL_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Monotremata --quantity 1 --store MONOTREME_OUTGROUPS
```

Broader primate target panel:

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --ranks family,genus \
  --quantities 2,1 \
  --store PRIMATE_TARGETS
```

Combined mammal target panel:

```bash
phyloODB mammal_workshop.db assemblies \
  --accessions @PRIMATE_TARGETS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --store MAMMAL_TARGETS
```

## 6. Download the reference panel, then the target panel

```bash
phyloODB mammal_workshop.db run download \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --protein

phyloODB mammal_workshop.db run download --accessions @MAMMAL_TARGETS --protein
```

## 7. Verify and clean downloads

```bash
phyloODB mammal_workshop.db run verify-downloads --accessions @PRIMATE_REFS,@MAMMAL_TARGETS --downloaded-only
phyloODB mammal_workshop.db run clean-isoforms --accessions @PRIMATE_REFS,@MAMMAL_TARGETS --downloaded-only
```

## 8. Run BUSCO on the references and targets

```bash
phyloODB mammal_workshop.db run batch-busco \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --lineage mammalia_odb12 \
  --format protein

phyloODB mammal_workshop.db run batch-busco \
  --accessions @MAMMAL_TARGETS \
  --lineage mammalia_odb12 \
  --format protein
```

Inspect the target panel afterward:

```bash
phyloODB mammal_workshop.db list assemblies \
  --accessions @MAMMAL_TARGETS \
  --library-name mammalia_odb12 \
  --busco \
  --busco-complete-min 90 \
  --busco-single-min 80 \
  --tidy
```

## 9. Build the derived mammal library

```bash
phyloODB mammal_workshop.db run add-library \
  --name mammal_core_odb12 \
  --coverage Mammalia \
  --parent-library-name mammalia_odb12 \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS
```

## 10. Run hidden-paralog filtering

```bash
phyloODB mammal_workshop.db run paralog-removal \
  --library-name mammal_core_odb12 \
  --ref-accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --targets @MAMMAL_TARGETS \
  --report-dir wiki_example/paralog_removal
```

## 11. Run internal contamination screening

```bash
phyloODB mammal_workshop.db run internal-decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --rank order \
  --hit-window 8 \
  --off-clade-fraction 0.05 \
  --report-path wiki_example/internal_decontamination/mammal_internal
```

Optional reference-based contamination screen:

```bash
phyloODB mammal_workshop.db run decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --refs @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --rank order \
  --off-clade-fraction 0.10 \
  --min-buscos 20 \
  --report-path wiki_example/decontamination/mammal_reference
```

## 12. Export the final filtered dataset

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

## If you only need a quick preview export

Before derived-library and filtering work, you can still do:

```bash
phyloODB mammal_workshop.db run export \
  --library-name mammalia_odb12 \
  --accessions @MAMMAL_TARGETS \
  --out-dir wiki_example/quick_export \
  --disable-paralog-filter \
  --disable-decont-filter
```

Use that as a preview, not as the final study dataset.

## What this example teaches

- Reference and target panels should be separated.
- Metadata registration is not downloading.
- BUSCO is necessary but usually not sufficient.
- The derived library is where the study becomes specific.
- Filtering is what turns a plausible panel into a safer export.
