# Downloading Assemblies

When to read this page: read this after you have stored one or more accession panels and are ready to fetch sequence files.

This page covers the transition from “PhyloODB knows these assemblies exist” to “the files are actually on disk and ready for BUSCO”.

## The workflow in one sentence

The normal order is:

metadata registration -> inspection -> stored panel -> download

## Download a stored panel

For the first pass, protein-focused downloads are the simplest path:

```bash
phyloODB mammal_workshop.db run download --accessions @PRIMATE_REFS --protein
```

For the broader target set:

```bash
phyloODB mammal_workshop.db run download --accessions @MAMMAL_TARGETS --protein
```

## What `download` actually does

The download task fetches the assembly data needed for later analysis and binds the downloaded content into the project state.

Important: this is the first step that actually puts sequence files on disk. Earlier metadata steps do not.

## Check what was downloaded

Use selectors to confirm the state:

```bash
phyloODB mammal_workshop.db count assemblies --accessions @PRIMATE_REFS --downloaded-only
phyloODB mammal_workshop.db list assemblies --accessions @PRIMATE_REFS --downloaded-only --tidy
```

If those counts are zero, you have metadata but not downloaded content.

## Verify downloads

Before running expensive downstream work, verify what you downloaded:

```bash
phyloODB mammal_workshop.db run verify-downloads --accessions @PRIMATE_REFS,@MAMMAL_TARGETS --downloaded-only
```

This is a good habit even if the earlier steps appeared to succeed.

## Isoform cleaning

The current CLI defaults are already oriented toward isoform cleaning during download, but it is still useful to know that explicit cleanup exists:

```bash
phyloODB mammal_workshop.db run clean-isoforms --accessions @PRIMATE_REFS,@MAMMAL_TARGETS --downloaded-only
```

Why this matters:

- BUSCO and downstream comparisons are easier to interpret when redundant isoforms are reduced.
- Reference panels especially benefit from cleaner proteomes.

## A simple, safe download sequence

If you want a copy-pasteable pattern:

```bash
phyloODB mammal_workshop.db run download --accessions @PRIMATE_REFS --protein
phyloODB mammal_workshop.db run verify-downloads --accessions @PRIMATE_REFS --downloaded-only
phyloODB mammal_workshop.db run download --accessions @MAMMAL_TARGETS --protein
phyloODB mammal_workshop.db run verify-downloads --accessions @MAMMAL_TARGETS --downloaded-only
```

Many users prefer downloading the references first, because the reference panel drives the later library-building stage.

## Common gotchas

### “I can list assemblies, so why can’t BUSCO run?”

Because listing assemblies only proves metadata exists. BUSCO needs actual downloaded data.

### “Should I download everything before deciding on a reference panel?”

Usually no. Decide on a reference panel first. It keeps the study focused and reduces wasted work.

### “Should I use `run` or `queue` here?”

- Use `run` for a small or first-time workflow.
- Use `queue` when the study is bigger or you want job scheduling.

## Next step

After your files are downloaded, continue to [Running BUSCO](Running-BUSCO).
