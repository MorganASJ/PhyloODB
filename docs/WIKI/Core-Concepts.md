# Core Concepts

When to read this page: read this after setup if the CLI feels like a list of commands instead of a coherent workflow.

## The basic model

PhyloODB is built around a project database:

```bash
phyloODB mammal_workshop.db <command> ...
```

That database is not just bookkeeping. It is where PhyloODB stores:

- Taxonomic information
- Assembly metadata
- Downloaded-file bindings
- BUSCO runs
- Derived libraries
- Filtering results
- Reusable variables such as stored accession sets

The database is the project. The files on disk are important, but the database is what makes the workflow inspectable and repeatable.

## Assemblies vs metadata vs downloaded files

These are different states.

- Metadata means PhyloODB knows an assembly exists and has recorded its attributes.
- Downloaded means the sequence files for that assembly have been fetched and bound in the project.
- BUSCO-ready means you have both the downloaded data and the BUSCO lineage context needed for BUSCO analysis.

This is why the usual sequence starts with `run add` and only later uses `run download`.

## Parent BUSCO lineage vs derived library

A parent BUSCO lineage is the broad lineage dataset you download from BUSCO, for example:

```bash
phyloODB mammal_workshop.db run download-busco-library --lineage mammalia_odb12 --coverage 1
```

That lineage is useful, but it is often too broad to be your final study definition.

A derived library is a study-specific library built from reference accessions:

```bash
phyloODB mammal_workshop.db run add-library \
  --name mammal_core_odb12 \
  --coverage Mammalia \
  --parent-library-name mammalia_odb12 \
  --accessions @REFERENCE_PANEL
```

In plain language:

- Parent lineage: broad BUSCO family universe
- Derived library: the subset and interpretation that actually fits your study

## Reference panel vs target panel

These are not interchangeable.

### Reference panel

This is the smaller, curated set used to define the derived library and support filtering decisions. It should be high quality and taxonomically informative.

### Target panel

This is the broader set you eventually want to export. It can be larger, denser, and more exploratory than the reference set.

Good workflow usually looks like this:

- Build a careful reference panel first
- Build or validate the library against that panel
- Screen the broader target panel
- Export the target panel

## Stored variables

PhyloODB lets you save accession selections into named variables:

```bash
phyloODB mammal_workshop.db assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_REFS
```

You can reuse that stored set later:

```bash
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_REFS --tidy
```

This is one of the most useful habits in PhyloODB. It keeps your workflow reproducible and reduces copy-paste mistakes.

## `run` vs `queue`

Both use the same task catalogue.

- `run` executes now in the foreground.
- `queue` records the task for queued execution.

Use `run` when:

- You are learning the workflow
- You want the shortest path through a task
- You are trying a small panel or quick test

Use `queue` when:

- You are building a larger study
- You want to coordinate multiple jobs
- You want to lean into the daemon/scheduling model

## Why the workflow feels staged

The system is intentionally staged because each phase answers a different question:

- Metadata: what is available?
- Selection: which taxa do I actually want?
- Download: do I have the files?
- BUSCO: how complete and single-copy are they?
- Derived library: which BUSCO families are trustworthy for this study?
- Filtering: which sequences or taxa still look suspicious?
- Export: what should I actually keep?

If you keep those questions separate, the CLI becomes much easier to reason about.
