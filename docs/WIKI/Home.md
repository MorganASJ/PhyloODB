# PhyloODB Wiki

When to read this page: start here if you want the shortest route through the wiki.

PhyloODB is a database-centred system for building phylogenomic datasets from public or local assemblies. Instead of treating download, BUSCO, filtering, and export as disconnected scripts, it stores the state of your project in one SQLite database so you can inspect, revise, and rerun the same study without losing track of what happened.

## What PhyloODB is good for

- Building a repeatable taxon panel from large public clades.
- Downloading assemblies and preparing proteomes for BUSCO-based work.
- Turning a broad BUSCO lineage into a more conservative study-specific library.
- Screening likely hidden paralogs and contamination before export.
- Exporting a dataset that reflects the current filtering state of the project.

## Recommended reading order

1. [Getting Started](Getting-Started)
2. [Core Concepts](Core-Concepts)
3. [Choosing Taxa and Building Panels](Choosing-Taxa-and-Building-Panels)
4. [Downloading Assemblies](Downloading-Assemblies)
5. [Running BUSCO](Running-BUSCO)
6. [Building a Derived Library](Building-a-Derived-Library)
7. [Filtering Paralogs and Contamination](Filtering-Paralogs-and-Contamination)
8. [Exporting a Dataset](Exporting-a-Dataset)
9. [End-to-End Example Mammal Dataset](End-to-End-Example-Mammal-Dataset)
10. [Troubleshooting](Troubleshooting)

## Fast paths

- First-time setup: [Getting Started](Getting-Started)
- “Show me how to choose taxa”: [Choosing Taxa and Building Panels](Choosing-Taxa-and-Building-Panels)
- “I just want a first export”: [Exporting a Dataset](Exporting-a-Dataset)
- Full worked example: [End-to-End Example Mammal Dataset](End-to-End-Example-Mammal-Dataset)
- Common confusion points: [Troubleshooting](Troubleshooting)
- Short command lookup: [Command Cheat Sheet](Command-Cheat-Sheet)

## The shortest successful path

If you want one compact first run, the minimum practical path is:

```bash
phyloODB mammal_workshop.db create
phyloODB mammal_workshop.db run add --clade Primates
phyloODB mammal_workshop.db assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_PANEL
phyloODB mammal_workshop.db run download --accessions @PRIMATE_PANEL --protein
phyloODB mammal_workshop.db run download-busco-library --lineage mammalia_odb12 --coverage 1
phyloODB mammal_workshop.db run batch-busco --accessions @PRIMATE_PANEL --lineage mammalia_odb12 --format protein
phyloODB mammal_workshop.db run export \
  --library-name mammalia_odb12 \
  --accessions @PRIMATE_PANEL \
  --out-dir wiki_example/quick_export \
  --disable-paralog-filter \
  --disable-decont-filter
```

That first export is intentionally simple. For a real study-quality dataset, continue through the library-building and filtering pages before you treat the export as final.

## Important distinctions

- Registering assembly metadata is not the same as downloading sequence files.
- A parent BUSCO lineage is not the same thing as a study-specific derived library.
- A reference panel and a target panel have different jobs.
- `run` is the easiest way to learn the system; `queue` is better once your jobs get larger.

## Existing longer-form docs

The repository also includes:

- `docs/QUICKSTART.md` for a compact local quick start.
- `docs/MANUAL.md` for the full handbook.
- `docs/COMMAND_REFERENCE.md` for dense command lookup.
- `tutorial/README.md` for the source version of the mammal example.
