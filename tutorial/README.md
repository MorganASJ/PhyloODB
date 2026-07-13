# Tutorial: Building a Mammal Core Dataset with Chimpanzee-Centred Sampling

This tutorial is a worked example of how one might use PhyloODB to construct a conservative mammalian dataset with a strong primate component. It is based conceptually on `tests/chimp_example.sh`, but it is written as a user-facing walkthrough rather than as a smoke test. The aim is not merely to obtain chimpanzee data, but to build a reusable and inspectable mammal project in which primates are represented densely and more distant mammalian outgroups are included to stabilise orthology and downstream export.

Use this tutorial for: a worked end-to-end example. If you want the shortest first-run guide, use [docs/QUICKSTART.md](../docs/QUICKSTART.md). If you want the full handbook, use [docs/MANUAL.md](../docs/MANUAL.md). If you want exact command lookup, use [docs/COMMAND_REFERENCE.md](../docs/COMMAND_REFERENCE.md).

This tutorial uses `queue` on purpose so the project state stays explicit and reusable. That also means queued tasks do not run by themselves. Start a daemon in another terminal before or during the queueing sequence:

```bash
phyloODB-daemon mammal_tutorial.db start --here --log-console --log-level INFO
```

Or run it in the background:

```bash
phyloODB-daemon mammal_tutorial.db start --background
```

When you are done queueing work, stop it cleanly with:

```bash
phyloODB-daemon mammal_tutorial.db stop --drain
```

[Figure placeholder: tutorial overview showing primate core, mammalian outgroups, derived library, paralog removal, decontamination, and export.]

## 1. Aim of the tutorial

By the end of this tutorial, the user should understand how to:

- create a new database for a focused study;
- populate it with assembly metadata for primates and selected mammalian outgroups;
- use selectors to choose a reference panel and a target panel;
- download and prepare the selected assemblies;
- run BUSCO on those assemblies;
- create a mammalian core library from a curated set of references;
- apply hidden paralog removal;
- apply decontamination, including the internal method;
- export a filtered dataset suitable for downstream phylogenomic work.

The tutorial uses chimpanzee and close primates as the centre of the target design, but it does not make the mistake of treating the reference set as a primate-only dataset. For a mammalian core dataset, the reference panel should include strategically chosen non-primate mammals. This improves the chances that the derived library represents Mammalia rather than only a narrow primate subset.

## 2. Conceptual design of the dataset

We distinguish two panels.

### Reference panel

This is the panel used to build the derived library and, later, to support hidden paralog removal. It should contain a small number of high-quality assemblies spanning the intended study clade.

For this tutorial, a sensible reference panel might include:

- several apes or other primates, including chimpanzee;
- one or more rodents;
- one carnivoran;
- one artiodactyl;
- one marsupial;
- one monotreme.

This is not the only valid choice, but it reflects a useful principle: references should be broad enough to catch hidden paralogs and unstable BUSCO families.

### Target panel

This is the set from which sequences will eventually be exported. It can be broader or denser than the reference panel. Here the goal is a primate-rich mammalian dataset, so we sample primates more densely but retain several mammalian outgroups.

## 3. Create the database

```bash
phyloODB mammal_tutorial.db create
```

At this point you have an empty project database. Everything that follows is attached to this file.

## 4. Populate assembly metadata

The first active stage is to tell PhyloODB what assemblies exist for the clades of interest. We begin with Primates and then add a small set of broader mammalian outgroups.

```bash
phyloODB mammal_tutorial.db queue add --clade Primates
phyloODB mammal_tutorial.db queue add --clade Rodentia
phyloODB mammal_tutorial.db queue add --clade Carnivora
phyloODB mammal_tutorial.db queue add --clade Artiodactyla
phyloODB mammal_tutorial.db queue add --clade Marsupialia
phyloODB mammal_tutorial.db queue add --clade Monotremata
```

This does not download sequence files. It populates the database with assembly metadata so that selectors can be used intelligently.

## 5. Register the BUSCO lineages needed for the tutorial

We want a broad metazoan lineage available, but for this tutorial the main working lineage is mammalian.

```bash
phyloODB mammal_tutorial.db queue download-busco-library --lineage metazoa_odb12 --coverage 1
phyloODB mammal_tutorial.db queue --schedule succeeded:LAST download-busco-library --lineage mammalia_odb12 --coverage 1
```

The scheduling here is not strictly required if both libraries are already present, but it illustrates the intended queue model.

At this point it is reasonable to start the daemon if it is not already running. The rest of the tutorial assumes queued tasks will actually be processed.

## 6. Inspect primate assemblies

Before choosing any references, inspect what the database contains.

```bash
phyloODB mammal_tutorial.db list assemblies -c Primates
phyloODB mammal_tutorial.db list assemblies -c Primates --level chromosome
phyloODB mammal_tutorial.db list assemblies -c Primates -r family -q 1
phyloODB mammal_tutorial.db list assemblies -c Primates -r family,genus -q 2,1
```

Use these commands to identify whether the clade is rich enough for dense primate sampling and whether chromosome-level assemblies exist where expected.

## 7. Store candidate panels as variables

A practical PhyloODB workflow stores important selections rather than reconstructing them from memory.

### 7.1 A small primate reference panel

```bash
phyloODB mammal_tutorial.db list assemblies \
  -c Primates \
  -r genus \
  -q 1 \
  -s PRIMATE_REFS
```

This produces one representative per primate genus according to the current ranking rules. In a real project you may refine this by inspecting the resulting list and perhaps replacing a weak assembly with an explicit accession.

### 7.2 Broader mammalian outgroups

Store a few strategic outgroup sets.

```bash
phyloODB mammal_tutorial.db list assemblies --clade Rodentia --rank order --quantity 1 --store RODENT_OUTGROUPS
phyloODB mammal_tutorial.db list assemblies --clade Carnivora --rank family --quantity 1 --store CARNIVORE_OUTGROUPS
phyloODB mammal_tutorial.db list assemblies --clade Artiodactyla --rank family --quantity 1 --store UNGULATE_OUTGROUPS
phyloODB mammal_tutorial.db list assemblies --clade Marsupialia --quantity 1 --store MARSUPIAL_OUTGROUPS
phyloODB mammal_tutorial.db list assemblies --clade Monotremata --quantity 1 --store MONOTREME_OUTGROUPS
```

### 7.3 A denser primate target panel

```bash
phyloODB mammal_tutorial.db list assemblies \
  --clade Primates \
  --ranks family,genus \
  --quantities 2,1 \
  --store PRIMATE_TARGETS
```

This target panel is broader than the reference panel and is intended for downstream export.

### 7.4 Combine the primate-rich target set with mammalian outgroups

To turn the primate target set into a true mammalian dataset, combine it with the stored outgroup panels. This is also a natural place to demonstrate `--append-to`, which unions a newly resolved accession set into an existing stored variable without duplicating accessions.

```bash
phyloODB mammal_tutorial.db list assemblies \
  --accessions @PRIMATE_TARGETS \
  --store MAMMAL_TARGETS

phyloODB mammal_tutorial.db list assemblies \
  --accessions @RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --append-to MAMMAL_TARGETS
```

High-level behavior:

- `--store MAMMAL_TARGETS` creates or replaces the variable;
- `--append-to MAMMAL_TARGETS` unions the new accessions into that variable while preserving existing order and avoiding duplicates.

From this point onward, `@MAMMAL_TARGETS` is the main target panel for BUSCO, filtering, and export.

## 8. Download the selected assemblies

Download the reference panel first. This is the panel on which the derived library will depend.

```bash
phyloODB mammal_tutorial.db queue download --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS --protein
```

Then inspect what has been downloaded.

```bash
phyloODB mammal_tutorial.db count assemblies --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS --downloaded-only
```

After the references, download the denser primate target panel.

```bash
phyloODB mammal_tutorial.db queue download --accessions @MAMMAL_TARGETS --protein
```

If proteome hygiene is a concern, isoform cleaning can either be handled during download or explicitly later using `clean-isoforms`.

## 9. Verify and prepare the downloaded data

It is often worthwhile to verify file integrity before large analyses.

```bash
phyloODB mammal_tutorial.db queue verify-downloads --accessions @MAMMAL_TARGETS,@PRIMATE_REFS --downloaded-only
```

If needed, clean isoforms explicitly.

```bash
phyloODB mammal_tutorial.db queue clean-isoforms --accessions @MAMMAL_TARGETS,@PRIMATE_REFS --downloaded-only
```

These steps are not conceptually glamorous, but they reduce the risk that later BUSCO or BLAST stages are distorted by file corruption or redundant isoform structure.

## 10. Run BUSCO on the reference and target panels

The mammalian BUSCO library is now used to profile the panels.

```bash
phyloODB mammal_tutorial.db queue batch-busco \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --lineage mammalia_odb12 \
  --format protein

phyloODB mammal_tutorial.db queue batch-busco \
  --accessions @MAMMAL_TARGETS \
  --lineage mammalia_odb12 \
  --format protein
```

Once BUSCO results exist, selectors can use BUSCO-based thresholds.

```bash
phyloODB mammal_tutorial.db list assemblies \
  -a @MAMMAL_TARGETS \
  -l mammalia_odb12 \
  --busco \
  --busco-complete-min 90 \
  --busco-single-min 80
```

This stage often leads the user to refine the target set. Poor assemblies can be excluded before the final library is built.

## 11. Create a mammalian core library

This is the defining step of the tutorial. We now construct a derived library using the reference panel. The reference panel is not simply a list of taxa to keep later. It is the empirical basis for deciding which BUSCO families behave well enough across Mammalia to serve as a conservative core set.

```bash
phyloODB mammal_tutorial.db queue add-library \
  --name mammal_core_odb12 \
  --coverage Mammalia \
  --parent-library-name mammalia_odb12 \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS
```

At a high level, this task ensures metadata and downloads exist, runs any missing BUSCO analyses, runs OrthoFinder if necessary, and then produces a filtered library based on the concordance between BUSCO families and OrthoFinder orthogroups.

Current tree behavior matters here:

- default `--gene-tree-source iqtree`: build replacement MAFFT alignments and IQ-TREE trees for the accepted orthogroups;
- those canonical core-set IQ-TREE trees are written under the OrthoFinder result directory in `IQ-TREE_Orthogroup_trees`;
- OrthoFinder's own `Resolved_Gene_Trees` is left untouched;
- `--gene-tree-source fasttree` or `--fast-tree`: reuse OrthoFinder `Resolved_Gene_Trees` directly and skip the replacement MAFFT and IQ-TREE stages;
- `--rerun-gene-trees` forces replacement IQ-TREE trees to be rebuilt even if matching trees already exist in `IQ-TREE_Orthogroup_trees`.

For a first serious build, the default IQ-TREE mode is usually the better conservative choice. If the user wants a quicker or lighter first pass, `--fast-tree` is the explicit shortcut.

After a successful build, useful outputs include:

- `library_build_metadata.json`, which records the effective core-set strategy and gene-tree source;
- `orthogroup_tree_manifest.tsv`, which records the tree paths used downstream;
- the accepted BUSCO family list for the derived library itself.

The reason for building a derived library at this stage is that the parent BUSCO lineage alone is often too permissive for phylogenomic export. The custom library defines a more conservative set of families suitable for this particular study.

## 12. Hidden paralog removal

With the library in place, apply hidden paralog filtering to the target panel using the same broad reference panel.

```bash
phyloODB mammal_tutorial.db queue paralog-removal \
  --library-name mammal_core_odb12 \
  --ref-accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --accessions @MAMMAL_TARGETS \
  --out-dir tutorial/results/paralog_removal
```

High-level interpretation:

- BUSCO says a family is present.
- Paralog removal asks whether the selected BUSCO sequence behaves like the expected ortholog when compared with trusted references.
- Families failing this screen remain present in the database, but they need not be accepted for export.

This is what the older project documents describe as hidden paralog removal: not a family-level orthology redefinition, but an accession-level safeguard against the wrong copy slipping through.

## 13. Decontamination approaches

PhyloODB provides two relevant approaches for this tutorial.

### 13.1 Internal decontamination

This is often the most natural first pass for a coherent target set.

```bash
phyloODB mammal_tutorial.db queue internal-decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --rank order \
  --hit-window 8 \
  --off-clade-fraction 0.05 \
  --report-path tutorial/results/internal_decontamination/mammal_internal
```

Conceptually, the target set is compared against itself using BUSCO-derived sequences. The question is whether each assembly behaves like a consistent member of the intended mammalian grouping.

### 13.2 Reference-based decontamination

A second option is explicit reference-based decontamination using the broader mammalian reference panel.

```bash
phyloODB mammal_tutorial.db queue decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --refs @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --rank order \
  --off-clade-fraction 0.10 \
  --min-buscos 20 \
  --report-path tutorial/results/decontamination/mammal_reference
```

This route is useful when the user wants explicit control over the reference panel used to judge off-clade behaviour.

### 13.3 Optional external confirmation

If the internal route identifies questionable BUSCOs and an external BLAST database is available, the user may follow with:

```bash
phyloODB mammal_tutorial.db queue external-decontamination-check \
  --run-id <internal_run_id> \
  --blast-db-path /path/to/external/db \
  --output-dir tutorial/results/external_decontamination

phyloODB mammal_tutorial.db queue external-decontamination-apply \
  --source-run-id <internal_run_id> \
  --run-label mammal_internal_confirmed \
  --external-blast-output-dir tutorial/results/external_decontamination
```

This two-step process is a refinement rather than a mandatory stage.

## 14. Export the final dataset

Once the custom library, paralog filtering, and decontamination runs exist, export the filtered dataset. In this tutorial the exported panel includes both the dense primate sample and the selected mammalian outgroups.

```bash
phyloODB mammal_tutorial.db queue export \
  --library-name mammal_core_odb12 \
  --accessions @MAMMAL_TARGETS \
  --out-dir tutorial/exports/mammal_core \
  --write-lineage-csv \
  --write-busco-report \
  --write-busco-family-matrix \
  --busco-report-extended
```

The default export behaviour assumes that prior filtering is meaningful and should be respected. This is exactly what makes the export suitable for downstream use.

## 15. Practical checks that are easy to forget

Before trusting the exported dataset, it is worth checking a few concrete things:

- confirm that the reference panel stored in `@PRIMATE_REFS,...` is really the one you intended to use;
- inspect `library_build_metadata.json` for the derived library so the effective `gene_tree_source` and build settings are explicit;
- if default IQ-TREE mode was used, confirm that `IQ-TREE_Orthogroup_trees` was created under the reused or newly built OrthoFinder result directory;
- if `--fast-tree` was used, remember that the core-set analysis is reusing OrthoFinder `Resolved_Gene_Trees` rather than creating replacement IQ-TREE trees;
- review paralog-removal and decontamination outputs before export if the project is intended for downstream phylogenomic inference rather than just a first exploratory pass.

If the final dataset is intended to retain only reasonably well-sampled families, occupancy thresholds can be supplied at this stage.

## 15. What to examine after export

The tutorial is complete when you can inspect:

- the derived library itself;
- the paralog filtering output;
- the decontamination summaries and run ids;
- the exported per-family FASTA files;
- the lineage and BUSCO reports written by export.

The essential point is that the final FASTA set is not simply “whatever BUSCO found”. It is the product of explicit choices about taxon sampling, library definition, paralog filtering, and contamination screening.

## 16. Variants on this tutorial

Several natural extensions are possible.

- Increase primate density by using a multi-stage selector with more genera per family.
- Use `count assemblies` before downloading to estimate the size of the project.
- Replace the generic outgroup sampling with explicit accession curation once a preferred set of assemblies is known.
- Create a second derived library focused on Primates alone and compare its behaviour with the broader mammalian library.

The companion shell script `tutorial/chimp_mammal_core_example.sh` provides a commented queue-oriented version of this tutorial.
