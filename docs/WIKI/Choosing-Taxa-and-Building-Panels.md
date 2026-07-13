# Choosing Taxa and Building Panels

When to read this page: read this once metadata for your clade is in the database and you need to choose what to keep.

This page is about selectors. The main idea is simple: start with a clade, inspect what exists, then store reusable panels instead of rebuilding the same selection from memory.

## Start by inspecting, not downloading

Register metadata first:

```bash
phyloODB mammal_workshop.db run add --clade Primates
```

Now inspect the candidate assemblies:

```bash
phyloODB mammal_workshop.db list assemblies --clade Primates --tidy
phyloODB mammal_workshop.db count assemblies --clade Primates
```

The goal here is to understand the candidate pool before you commit to a download.

## The most useful selector patterns

### One assembly per genus

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --rank genus \
  --quantity 1
```

### One assembly per family

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --rank family \
  --quantity 1
```

### Top five assemblies from the resolved set

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --quantity 5
```

### Multi-stage subsampling

This is useful when you want broader taxonomic coverage without hand-curating every accession:

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --ranks family,genus \
  --quantities 2,1
```

That pattern means: choose two families, then one genus within each chosen family.

## Restrict to better assemblies

If you want to inspect chromosome-level assemblies only:

```bash
phyloODB mammal_workshop.db list assemblies \
  --clade Primates \
  --level chromosome \
  --tidy
```

You can also layer filters later, once BUSCO results exist.

## Store panels as variables

Do this early. It makes the entire rest of the workflow easier.

Store a reference panel:

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --rank genus \
  --quantity 1 \
  --store PRIMATE_REFS
```

Store a denser target panel:

```bash
phyloODB mammal_workshop.db assemblies \
  --clade Primates \
  --ranks family,genus \
  --quantities 2,1 \
  --store PRIMATE_TARGETS
```

Inspect a stored variable:

```bash
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_REFS --tidy
```

## A practical mammal example

For a mammal-wide project, it is common to separate a small reference set from a broader target set.

Reference-oriented outgroups:

```bash
phyloODB mammal_workshop.db assemblies --clade Rodentia --rank order --quantity 1 --store RODENT_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Carnivora --rank family --quantity 1 --store CARNIVORE_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Artiodactyla --rank family --quantity 1 --store UNGULATE_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Marsupialia --quantity 1 --store MARSUPIAL_OUTGROUPS
phyloODB mammal_workshop.db assemblies --clade Monotremata --quantity 1 --store MONOTREME_OUTGROUPS
```

Then combine them with the primate target panel:

```bash
phyloODB mammal_workshop.db assemblies \
  --accessions @PRIMATE_TARGETS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --store MAMMAL_TARGETS
```

## Mental model for good panel design

- Reference panel: smaller, cleaner, taxonomically strategic
- Target panel: broader, export-oriented
- Stored variables: treat them as named project decisions

If you skip this separation, later steps become harder to interpret.

## Copy-paste recipes

### One panel per genus for a clade

```bash
phyloODB mammal_workshop.db assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_PANEL
```

### Count how many candidates a selector would resolve

```bash
phyloODB mammal_workshop.db count assemblies --clade Primates --ranks family,genus --quantities 2,1
```

### Inspect a stored set before download

```bash
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_PANEL --tidy
```

## Next step

Once your panels look sensible, continue to [Downloading Assemblies](Downloading-Assemblies).
