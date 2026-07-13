# PhyloODB

Phylogenetic Ortholog Database tooling with a registry-driven task system.

## Installation

Create the recommended workstation environment:

```bash
mamba env create -f environment.yml
conda activate phyloodb
```

If you already maintain an environment with the same dependencies, activate
that environment instead. For example, this checkout is commonly used from:

```bash
conda activate podb2
```

Install the Python package into that environment:

```bash
pip install .
```

For development work, use an editable install instead:

```bash
pip install -e .[dev]
```

Build release artifacts with:

```bash
python -m build
```

The packaged Python dependencies are only part of the full toolchain. A full
end-to-end PhyloODB analysis setup also expects external tools such as BUSCO,
OrthoFinder, BLAST+, MAFFT, and IQ-TREE. The bundled `environment.yml`
includes those tools.

## Documentation

- `docs/WIKI/Home.md` – GitHub-wiki-ready end-user documentation organised by topic.
- `docs/QUICKSTART.md` – first-run handbook for getting from an empty database
  to a basic export quickly.
- `docs/MANUAL.md` – comprehensive user guide covering concepts, queueing,
  storage, recovery, filtering, and export.
- `docs/COMMAND_REFERENCE.md` – compact command and task lookup reference.
- `tutorial/README.md` – worked end-to-end mammal example with queue-driven
  commands.
