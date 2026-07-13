# Running BUSCO

When to read this page: read this after your assemblies are downloaded and you are ready to assess completeness and single-copy behavior.

BUSCO is the bridge between “I have files” and “I can reason about dataset quality”.

## 1. Download the BUSCO lineage

Start by registering the BUSCO lineage you want to use:

```bash
phyloODB mammal_workshop.db run download-busco-library --lineage mammalia_odb12 --coverage 1
```

For larger cross-clade projects you may keep more than one lineage in the same database, but one well-chosen lineage is enough for a first pass.

## 2. Run BUSCO on your panel

For a reference panel:

```bash
phyloODB mammal_workshop.db run batch-busco \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --lineage mammalia_odb12 \
  --format protein
```

For the broader target panel:

```bash
phyloODB mammal_workshop.db run batch-busco \
  --accessions @MAMMAL_TARGETS \
  --lineage mammalia_odb12 \
  --format protein
```

The command name is `batch-busco` even though some help output still uses the internal task-class label.

## 3. Inspect BUSCO-aware assembly summaries

After the runs complete, inspect the results:

```bash
phyloODB mammal_workshop.db list assemblies \
  --accessions @MAMMAL_TARGETS \
  --library-name mammalia_odb12 \
  --busco \
  --tidy
```

If you want to focus on stronger assemblies:

```bash
phyloODB mammal_workshop.db list assemblies \
  --accessions @MAMMAL_TARGETS \
  --library-name mammalia_odb12 \
  --busco \
  --busco-complete-min 90 \
  --busco-single-min 80 \
  --tidy
```

## What to look for

At this stage, users usually want to answer three questions:

- Are some assemblies obviously weak and worth excluding?
- Does the reference panel still look like a good basis for a derived library?
- Does the broader target set include taxa that are too incomplete for export?

BUSCO is not the final truth about orthology, but it is the main early quality filter.

## A good practical habit

Do not move straight from raw BUSCO results to final export unless this is only a quick exploratory run.

BUSCO tells you:

- how complete an assembly looks
- how much single-copy signal you have

It does not by itself solve:

- hidden paralogs
- study-specific family selection
- contamination

## Common confusion points

### “I ran BUSCO, but export still wants more work.”

That can be normal, but not because export is strict by default. BUSCO is only one layer of evidence. A stronger export may still benefit from:

- a derived library
- hidden-paralog filtering
- decontamination

Export uses those results automatically when they exist, and can be told to require them explicitly.

### “Why use the same lineage on references and targets?”

Because the lineage gives you a comparable basis for ranking and filtering across the whole study.

## Next step

If this is just an exploratory pass, you can skim [Exporting a Dataset](Exporting-a-Dataset).

If you want a stronger study dataset, continue to [Building a Derived Library](Building-a-Derived-Library).
