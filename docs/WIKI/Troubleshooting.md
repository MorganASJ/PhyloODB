# Troubleshooting

When to read this page: read this when the commands run but the project state still feels confusing.

This page focuses on the most common user-side misunderstandings in the current workflow.

## “Why do I see assemblies but no downloaded files?”

Because metadata registration and download are separate stages.

This command only registers metadata:

```bash
phyloODB mammal_workshop.db run add --clade Primates
```

This command actually downloads files:

```bash
phyloODB mammal_workshop.db run download --accessions @PRIMATE_PANEL --protein
```

Check downloaded state explicitly:

```bash
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_PANEL --downloaded-only --tidy
```

## “Why does export require filtering results?”

It usually does not. By default, `export` uses paralog-filtering and decontamination results when they exist, but it does not require them.

Export only refuses missing filtering state when you explicitly ask for that with `--require-paralog-filtering` and/or `--require-decontamination`.

For a quick preview export where you want to ignore available filtering on purpose:

```bash
phyloODB mammal_workshop.db run export \
  --library-name mammalia_odb12 \
  --accessions @PRIMATE_PANEL \
  --out-dir wiki_example/quick_export \
  --disable-paralog-filter \
  --disable-decont-filter
```

## “When should I use `run` vs `queue`?”

Use `run` when:

- you are learning
- you want immediate feedback
- the panel is small

Use `queue` when:

- the study is larger
- you want scheduling and deferred execution
- you are chaining many jobs

The task names are the same in both modes.

## “What is a parent library vs derived library?”

- Parent library: the downloaded BUSCO lineage, such as `mammalia_odb12`
- Derived library: your study-specific library, such as `mammal_core_odb12`

If you are doing a quick first pass, the parent library may be enough.

If you are building a serious dataset, the derived library is usually the better basis for export.

## “Why did my selector return fewer taxa than expected?”

Common reasons:

- You applied a rank-based rule that intentionally subsampled the clade.
- You filtered to a better assembly level such as `--level chromosome`.
- You used BUSCO filters like `--busco-complete-min` or `--busco-single-min`.
- You stored a panel earlier and forgot what exact selector created it.

Inspect the candidate pool without extra filters first:

```bash
phyloODB mammal_workshop.db list assemblies --clade Primates --tidy
```

Then reapply the selector gradually.

## “Why did BUSCO or export pick a different run than I expected?”

PhyloODB has BUSCO run-selection logic rather than blindly using any matching run.

Things that can affect run choice include:

- `--library-name`
- `--busco-pipeline`
- `--prefer-busco-pipeline`
- `--format`
- `--prefer-format`
- `--proteome-profile`
- `--busco-run-selection`

If the result seems surprising, inspect the available BUSCO runs explicitly:

```bash
phyloODB mammal_workshop.db list busco-runs --accessions @PRIMATE_PANEL --tidy
```

## “I ran BUSCO, but I still cannot export the final dataset I want.”

That usually means one of these is still missing:

- a suitable derived library
- hidden-paralog filtering
- decontamination

BUSCO is the start of the quality-assessment phase, not the end of it.

## “Should I rebuild everything if I change my panel?”

Not always. But if you change the reference panel that defines the study library, you should assume downstream interpretation changes too. In practice:

- changing the target panel is often cheaper
- changing the reference panel is more consequential

## “Which page should I go back to?”

- Setup problems: [Getting Started](Getting-Started)
- Selection confusion: [Choosing Taxa and Building Panels](Choosing-Taxa-and-Building-Panels)
- Download state confusion: [Downloading Assemblies](Downloading-Assemblies)
- BUSCO confusion: [Running BUSCO](Running-BUSCO)
- Export behaviour confusion: [Exporting a Dataset](Exporting-a-Dataset)
