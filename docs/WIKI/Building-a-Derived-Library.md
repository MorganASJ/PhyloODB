# Building a Derived Library

When to read this page: read this once BUSCO runs exist and you want a more conservative, study-specific dataset definition.

This is the page that separates a quick BUSCO-based export from a serious project dataset.

## Why not export directly from the parent BUSCO lineage?

You can export directly from a parent BUSCO lineage for a first-pass result, but many studies should not stop there.

Why:

- The parent lineage is broad by design.
- Not every BUSCO family will behave equally well in your focal clade.
- Your study may need a more conservative family set than the raw lineage provides.

A derived library lets you define that stricter family universe from a curated reference panel.

## Build the right reference panel first

A good reference panel should be:

- smaller than your target panel
- relatively high quality
- taxonomically informative across the study clade

For a mammal-wide project with dense primate sampling, that usually means:

- several primates
- a few strategic non-primate mammals
- enough spread to help expose unstable families and hidden paralogs

## The core command

```bash
phyloODB mammal_workshop.db run add-library \
  --name mammal_core_odb12 \
  --coverage Mammalia \
  --parent-library-name mammalia_odb12 \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS
```

## What `add-library` orchestrates

At a high level, this step can coordinate several things:

- ensuring the reference panel is available
- reusing or running BUSCO as needed
- reusing or running OrthoFinder as needed
- comparing BUSCO families and orthogroups
- storing a derived library for later export/filtering steps

You do not need to treat it as a black box, but the important user-facing point is that this is where your study-specific library definition is created.

## What you should expect when it finishes

After a successful derived-library build, you should have:

- a named library you can refer to with `--library-name mammal_core_odb12`
- a clearer study-specific family universe than the raw BUSCO lineage alone
- a better basis for hidden-paralog screening and final export

## Reuse vs rerun

The CLI exposes flags such as:

- `--rerun-busco`
- `--rerun-orthofinder`
- `--force`

You usually do not need those on the first clean run. Use them when you intentionally want to rebuild state rather than reuse matching existing results.

## Practical advice

- Keep the reference panel curated.
- Avoid building the library from a huge, noisy target set.
- Do not confuse “more taxa” with “better references”.

If your reference panel is weak, the derived library will be harder to trust.

## Next step

Once the derived library exists, continue to [Filtering Paralogs and Contamination](Filtering-Paralogs-and-Contamination).
