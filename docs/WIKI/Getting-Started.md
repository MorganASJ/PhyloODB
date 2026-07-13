# Getting Started

When to read this page: read this before your first real command, especially on a new machine.

This page is about getting to the point where `phyloODB` runs cleanly and you understand what the very first commands do.

## What you need

At minimum, you need:

- The repository checked out locally.
- The Python package dependencies installed.
- A Python 3.10+ environment. The bundled `environment.yml` creates a
  ready-to-use environment named `phyloodb`.

The codebase also expects several external bioinformatics tools for the later analysis stages. In practice, you should expect to need these before a full end-to-end study:

- BUSCO
- OrthoFinder
- BLAST+
- MAFFT
- IQ-TREE

Some tasks can run before all of those are available. For example, database creation, metadata registration, taxon selection, and downloads do not require the full downstream toolchain.

## Recommended local setup

Create and activate the recommended environment first:

```bash
mamba env create -f environment.yml
conda activate phyloodb
```

Install the package from the repository root:

```bash
pip install .
```

Check that the CLI is visible:

```bash
phyloODB --help
```

If you prefer not to install it globally into the environment, the repository also works via module invocation:

```bash
python -m phyloODB.cli.main --help
```

For development work, install the editable package and test tooling instead:

```bash
pip install -e .[dev]
```

## Create your first database

Every project is anchored to one SQLite database file:

```bash
phyloODB mammal_workshop.db create --email you@example.org --api-key YOUR_NCBI_KEY
```

That database becomes the memory of the project. PhyloODB stores metadata, downloads, BUSCO runs, filtering results, and export-relevant state against it.

The `--email` and `--api-key` flags are optional, but you should set them
before running NCBI-backed add/download tasks.

## The first thing most users do next

After creating the database, most users should register assembly metadata for a clade:

```bash
phyloODB mammal_workshop.db run add --clade Primates
```

This does not download genomes or proteomes. It only tells PhyloODB what assemblies exist and records their metadata so you can choose a sensible panel.

That distinction matters:

- `run add --clade Primates` means “learn what is available”.
- `run download --accessions ...` means “fetch sequence files to disk”.

## A safe first-day workflow

If you are learning the system, this is the right order:

1. Create the database.
2. Register metadata for your clade.
3. Inspect the candidate assemblies.
4. Store a reusable panel in a variable.
5. Download that panel.
6. Download a BUSCO lineage.
7. Run BUSCO.
8. Decide whether you are ready for direct export or should build a derived library.

## Why `run` is recommended first

Use `run` when you are learning:

```bash
phyloODB mammal_workshop.db run add --clade Primates
```

Use `queue` when you are ready to manage a larger study through queued jobs:

```bash
phyloODB mammal_workshop.db queue add --clade Primates
```

The task catalogue is the same. The difference is whether you want immediate foreground execution or queued execution.

## Sanity checks

After your first metadata-registration task, these are useful checks:

```bash
phyloODB mammal_workshop.db count assemblies --clade Primates
phyloODB mammal_workshop.db list assemblies --clade Primates --rank genus --quantity 1 --tidy
```

If those work, you are ready for [Choosing Taxa and Building Panels](Choosing-Taxa-and-Building-Panels).

## Common first-run mistakes

- Creating the database and assuming genomes were downloaded. They were not.
- Going straight to export before any BUSCO work exists.
- Mixing `run` and `queue` examples without noticing which mode you are in.
- Trying advanced filtering before you have a stable reference panel.
