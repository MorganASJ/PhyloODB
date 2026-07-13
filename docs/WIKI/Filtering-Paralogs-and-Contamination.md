# Filtering Paralogs and Contamination

When to read this page: read this after you have a stable derived library and want a study-quality dataset rather than a quick exploratory export.

This page covers two different problems that are easy to blur together:

- hidden paralogs
- contamination

They are related, but they are not the same step.

## Recommended order

For most studies, a good order is:

1. Build the derived library.
2. Run hidden-paralog filtering.
3. Run contamination screening.
4. Export with the filtering results required.

## Hidden paralog removal

This step asks whether a BUSCO-selected sequence still behaves like the expected ortholog when compared against trusted references.

For the mammal example:

```bash
phyloODB mammal_workshop.db run paralog-removal \
  --library-name mammal_core_odb12 \
  --ref-accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --targets @MAMMAL_TARGETS \
  --report-dir wiki_example/paralog_removal
```

Interpret it like this:

- BUSCO says a family is present.
- Paralog filtering checks whether the chosen copy looks like the expected ortholog.
- The output becomes part of later export decisions.

## Internal contamination screening

Internal decontamination uses the target set itself as the main context:

```bash
phyloODB mammal_workshop.db run internal-decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --rank order \
  --hit-window 8 \
  --off-clade-fraction 0.05 \
  --report-path wiki_example/internal_decontamination/mammal_internal
```

This is often the most natural first contamination screen for a coherent study panel.

## Reference-based contamination screening

If you want an explicit trusted reference panel in the contamination logic:

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

## Which of these are optional?

For a first quick export:

- derived library creation can be skipped if you knowingly export from the parent lineage
- paralog filtering can be skipped
- contamination filtering can be skipped

For a study-quality export:

- derived library is strongly recommended
- paralog filtering is strongly recommended
- contamination screening is strongly recommended

PhyloODB export uses those results automatically when they are available, but it does not require them unless you ask it to with `--require-paralog-filtering` and/or `--require-decontamination`.

## How to think about the results

- Paralog filtering is about orthology confidence.
- Decontamination is about taxonomic plausibility.

A sequence can look complete in BUSCO terms and still fail one of these later screens.

## Common mistakes

- Treating BUSCO completeness alone as enough for final export.
- Using an uncurated reference panel for hidden-paralog screening.
- Forgetting that later export may use these filtering results automatically unless you deliberately ignore them with `--disable-paralog-filter` and/or `--disable-decont-filter`.

## Next step

Continue to [Exporting a Dataset](Exporting-a-Dataset).
