# PhyloODB

Phylogenetic Ortholog Database tooling with a registry-driven task system.

- [Manual](docs/MANUAL.md) – comprehensive user guide covering concepts,
  queueing, storage, recovery, filtering, and export.
- [Quickstart](docs/QUICKSTART.md) – first-run handbook for getting from an
  empty database to a basic export quickly.

RTFM

## Installation

PhyloODB is a Python package, but full end-to-end analyses also require
external bioinformatics tools including BUSCO, OrthoFinder, BLAST+, MAFFT, and
IQ-TREE. The bundled `environment.yml` installs both the Python dependencies
and those command-line tools.

### Install from a release checkout

Clone or download the repository, then create the recommended environment:

```bash
mamba env create -f environment.yml
conda activate phyloodb
```

Install PhyloODB into that environment:

```bash
pip install .
```

Check that the command-line entry points are available:

```bash
phyloODB --help
phyloODB-daemon --help
```

### Install directly from GitHub

If you do not need a local editable checkout, create an environment with the
runtime toolchain and install the Python package directly from GitHub:

```bash
mamba create -n phyloodb -c conda-forge -c bioconda \
  python=3.11 pip blast busco mafft orthofinder iqtree
conda activate phyloodb
pip install "git+https://github.com/MorganASJ/PhyloODB_.git"
```

For development work from a local clone, use an editable install instead:

```bash
pip install -e .[dev]
```

Build release artifacts with:

```bash
python -m build
```

## Documentation

- [Manual](docs/MANUAL.md) – comprehensive user guide covering concepts,
  queueing, storage, recovery, filtering, and export.
- [Quickstart](docs/QUICKSTART.md) – first-run handbook for getting from an
  empty database to a basic export quickly.
- [Command reference](docs/COMMAND_REFERENCE.md) – compact command and task
  lookup reference.
- [Mammal tutorial](tutorial/README.md) – worked end-to-end mammal example with
  serial `run` commands.
