# Tutorial: Building a Phylogenomic dataset for Primates

In this tutorial we will go through the steps of building a dataset for primates. This will involve inspecting the information available online and choosing which assemblies we wish to download and analyse. We will also go through the steps to build a primate core set informed by mammalian outgroups although, to avoid us all competing to run the analysis, I have done the long-running steps for us already.

This tutorial uses `run` commands so each stage executes before the next one is
started. That makes the walkthrough easier to follow in a terminal.
For larger production datasets, you may want to adapt the same task sequence to
your preferred batch or scheduler workflow.

## 1. Aim of the tutorial

By the end of this tutorial, the user should understand how to:

- create a new database for a focused study;
- populate it with assembly metadata for primates and selected mammalian outgroups;
- use selectors to choose a reference panel and a target panel;
- download and prepare the selected assemblies;
- run BUSCO on those assemblies;
- create a primate core library from a curated set of references that includes mammalian outgroups;
- apply hidden paralog removal;
- apply decontamination, including the internal method;
- export a filtered dataset suitable for downstream phylogenomic work.

The tutorial uses chimpanzee and close primates as the centre of the target design, but it does not make the mistake of treating the reference set as a primate-only dataset. For a mammalian core dataset, the reference panel should include strategically chosen non-primate mammals. This improves the chances that the derived library represents Mammalia rather than only a narrow primate subset.

## 2. Conceptual design of the dataset

We distinguish two panels of genome assemblies. 

### Reference panel

We need the reference dataset that we will be using to reciprocally validate BUSCO and OrthoFinder. This will be a small dataset spanning a clade of interest that we through PhyloODB run OrthoFinder on to obtain orthogroups. It is easy to obtain the BUSCO results for these taxa but we can then make a more confident BUSCO subset by comparing the BUSCO families with the OrthoFinder orthogroups. For this we will need annotated genomes with proteomes available. The easiest way to filter to these are to look at refseq assemblies.

For this tutorial, a sensible reference panel might include:

- several apes or other primates, including chimpanzee;
- one or more rodents;
- one carnivoran;
- one artiodactyl;
- one marsupial;
- one monotreme.

This is not the only valid choice, but it reflects a useful principle: references should be broad enough to catch hidden paralogs and unstable BUSCO families. We should also aim to ensure these are high-quality chromosome-level assemblies. 

### Target panel

This is the set from which sequences will eventually be exported. It can be broader or denser than the reference panel. Here the goal is a primate dataset so we will work to create this panel. For the actual demonstration we will subsample this to look at a smaller ape-focused group.

## 3. Install and create database

To install phyloODB we will need conda/mamba installed. If you are working on the brown nugget please run the following to install conda:

```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh
```

This will allow you to install into your home directory.

You can then install phyloODB using the following:

```
conda create -n phyloodb -c conda-forge -c bioconda \
  python=3.11 pip blast busco mafft orthofinder iqtree

conda activate phyloodb

pip install "git+https://github.com/MorganASJ/PhyloODB_.git"
```

I recommend running these line by line. Once this is run you should have a conda environment with phyloODB installed. Make sure to switch to it if you log out again.

### Create the database

I have provided an example database there with information already downloaded and analyses already complete to avoid us from killing trees. In `/shared/podb`, treat commands containing `create`, `run`, or `queue` as examples so do not run them. The safe hands-on commands during the workshop are inspection commands such as `list`, `count`, and `tree`.

However we can run the first section yourself in your home directories:
`cd ~/`

If you were doing this from scratch on your own device or on the cluster, you would first need to build the database.

```bash
phyloODB metazoa.db create --email your.email@domain.ac.uk
```

At this point you have an empty project database. Everything that follows is attached to this file. If you have not already you should create an NCBI account. This gives you access to an API key which will speed up downloading data. You can add this as a flag to the above command or set it later:

```bash
phyloODB metazoa.db set var NCBI_API_KEY [YOUR KEY]
```

There are many variables that can change default behaviour, flags, and mechanics. You can view them all by listing them:

```bash
phyloODB metazoa.db list variables -p
```

## 4. Populate assembly metadata

The first active stage is to tell PhyloODB what assemblies exist for the clades of interest. We begin with Primates and then add a small set of broader mammalian outgroups.

```bash
phyloODB metazoa.db run add --clade Primates
phyloODB metazoa.db run add --clade Rodentia

# We can actually just use:
phyloODB metazoa.db run add --clade Mammalia
```
This does not download sequence files. It instead downloads information on assemblies within those groups. This means the database has ingested knowledge about those taxa.

## 5. Register the BUSCO lineages needed for the tutorial


BUSCO has lineage specific datasets. We can download them into the database to use them in our analyses. Without doing this the system will not know they exist. 

```bash
phyloODB metazoa.db run download-busco-library --lineage metazoa_odb12
phyloODB metazoa.db run download-busco-library --lineage mammalia_odb12
```

## 6. Inspect primate assemblies

Before choosing any references, inspect what the database contains. Here we are saying list assemblies and filter to those in the taxonomic group Primates.

```bash
phyloODB metazoa.db list assemblies --clade Primates
```

This is pretty ugly. Let's view it as a table. In the paged view you can use the arrow keys or page keys to move through the results, and press `q` to quit.

```bash
phyloODB metazoa.db list assemblies -c Primates -p
```

If we add taxonomic ranks it is easier to understand what we may be looking at:
```bash
phyloODB metazoa.db list assemblies -c Primates -p --ranks order,family
```

Adding `--group-by-rank` makes it clearer.

```bash
phyloODB metazoa.db list assemblies -c Primates -p --ranks order,family --group-by-rank
```

But there are a lot of assemblies to see. Perhaps we only care about those that are REFSEQ and chromosome level. We can use filters to do this.

```bash
phyloODB metazoa.db list assemblies -c Primates -p --ranks order,family --group-by-rank --level chromosome --filter origin=refseq
```

If you list the available metadata, you can see more filtering or sorting options.

```bash
phyloODB metazoa.db list metadata
```

Let's select the top 2 from each family:
```bash
phyloODB metazoa.db list assemblies -c Primates -p -r family -q 2
```

Without BUSCO results stored this is just using the metadata to roughly rank assemblies. But we will use these to define our wide panel that we can then filter down.

## 7. Store candidate panels as variables

A practical PhyloODB workflow stores important selections rather than reconstructing them from memory. Here we will store the large panel PRIMATES, which contains up to 10 genera per family and up to 4 assemblies per genus.

```bash
phyloODB metazoa.db list assemblies -c Primates -p --ranks family,genus --group-by-rank --filter origin=refseq --quantities 10,4 --store PRIMATES
```

We can then view these easily:
```bash
phyloODB metazoa.db list assemblies -a @PRIMATES -p
```

And if it is helpful even draw a taxonomic tree in the terminal:
```bash
phyloODB metazoa.db tree -a @PRIMATES --show-accession --colour-by-ranks family
```

### 7.1 A small primate reference panel

I have preselected these to use in our reference panel. Normally we would get BUSCO results first but as I have already run the analysis we can define these here. As you can see we can use direct accession numbers or broad selectors to obtain lists of assemblies.

```bash
phyloODB metazoa.db list assemblies \
  -a GCF_009914755.1,GCF_049354715.1,GCF_041146395.1,GCF_049350105.2 \
  -S PRIMATE_REFS
```

This is our primate reference list. I have pre-chosen these four taxa as they span the tree of primates.

We will use these later when we reciprocally validate our BUSCO dataset for primates.

Store a few strategic outgroup sets.

```bash
phyloODB metazoa.db list assemblies --clade Rodentia --quantity 1 --filter origin=refseq --store RODENT_OUTGROUPS

phyloODB metazoa.db list assemblies --clade Carnivora --quantity 1 --filter origin=refseq --store CARNIVORE_OUTGROUPS

phyloODB metazoa.db list assemblies --clade Artiodactyla --quantity 1 --filter origin=refseq --store UNGULATE_OUTGROUPS

phyloODB metazoa.db list assemblies --clade Metatheria --quantity 1 --filter origin=refseq --store MARSUPIAL_OUTGROUPS

phyloODB metazoa.db list assemblies --clade Monotremata --quantity 1 --filter origin=refseq --store MONOTREME_OUTGROUPS
```

### 7.4 Combine the reference dataset

We will combine the hand-picked PRIMATE_REFS with the OUTGROUPS.

High-level behavior:
- `--store PRIMATE_REFS` creates or replaces the variable;
- `--append-to PRIMATE_REFS` unions the new accessions into that variable while preserving existing order and avoiding duplicates.

```bash
phyloODB metazoa.db list assemblies \
  -a @RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --append-to PRIMATE_REFS
```

This should result in PRIMATE_REFS - our reference dataset:

```bash
(phyloodb) ql22514@it037057:/shared/podb$ phyloODB metazoa.db list assemblies -a @PRIMATE_REFS
accession       species
GCF_009914755.1 Homo sapiens
GCF_015852505.1 Tachyglossus aculeatus
GCF_020171115.1 Neogale vison
GCF_036323735.1 Rattus norvegicus
GCF_041146395.1 Eulemur rufifrons
GCF_049350105.2 Macaca mulatta
GCF_049354715.1 Callithrix jacchus
GCF_054371585.1 Giraffa tippelskirchi
GCF_902635505.1 Sarcophilus harrisiis
```

```text
Metazoa
└─ Mammalia
  ├─────────────────────────Tachyglossus aculeatus
  └─ Theria
    ├─ Boreoeutheria
    │ ├─ Euarchontoglires
    │ │ ├───────────────────Rattus norvegicus
    │ │ └─ Primates
    │ │   ├─ Simiiformes
    │ │   │ ├─ Catarrhini
    │ │   │ │ ├─────────────Macaca mulatta
    │ │   │ │ └─────────────Homo sapiens
    │ │   │ └───────────────Callithrix jacchus
    │ │   └─────────────────Eulemur rufifrons
    │ └─ Laurasiatheria
    │   ├───────────────────Giraffa tippelskirchi
    │   └───────────────────Neogale vison
    └───────────────────────Sarcophilus harrisii
```

## 8. Download the selected assemblies

**From this point I would advise not entering commands that start work with `run` or `queue`. Downloading and running the analysis could take several hours and the machine cannot handle everyone doing this at once. Have a look at the commands and we can view the state of the example database instead.**

```bash
cd /shared/podb
```

**Can you derive a command to view the downloaded content?**

**DO NOT RUN COMMANDS THAT CONTAIN `run` OR `queue` DURING THE WORKSHOP UNLESS INSTRUCTED.**

Download the reference panel first. This is the panel on which the derived library will depend.

```bash
phyloODB metazoa.db run download --accessions @PRIMATE_REFS
```

Then inspect what has been downloaded.

```bash
phyloODB metazoa.db count assemblies --accessions @PRIMATE_REFS --downloaded-only
```

After the references, download the denser primate target panel.

```bash
phyloODB metazoa.db run download -a @PRIMATES
```

If proteome hygiene is a concern, isoform cleaning can either be handled during download or explicitly later using `clean-isoforms`. You will find that this has already occurred automatically during the download process. We can tweak automatic settings using variables.

The following command will show these:
```bash
phyloODB metazoa.db list proteome-profiles -a @PRIMATE_REFS -p
```

## 9. Verify and prepare the downloaded data

It is often worthwhile to verify file integrity before large analyses.

```bash
phyloODB metazoa.db run verify-downloads --accessions @PRIMATES,@PRIMATE_REFS --downloaded-only
```

If needed, clean isoforms explicitly.

```bash
phyloODB metazoa.db run clean-isoforms --accessions @PRIMATES,@PRIMATE_REFS --downloaded-only
```

## 10. Run BUSCO on the reference and target panels

The Metazoa BUSCO library is now used to profile the panels.

```bash
phyloODB metazoa.db run batch-busco \
  --accessions @PRIMATE_REFS \
  --lineage metazoa_odb12 \
  --format protein

phyloODB metazoa.db run batch-busco \
  --accessions @PRIMATES \
  --lineage metazoa_odb12 \
  --format protein
```

Once BUSCO results exist, selectors can use BUSCO-based thresholds.

```bash
phyloODB metazoa.db list assemblies \
  -a @PRIMATES \
  -l metazoa_odb12 \
  --busco \
  --busco-complete-min 90 \
  --busco-single-min 80 \
  --store PRIMATE_TARGETS
```

We can use `list results` to automatically filter to assemblies with BUSCO results:

```bash
phyloODB metazoa.db list results -a @PRIMATES -l metazoa_odb12 -r family -p --sort quality
```

This stage often leads the user to refine the target set. Poor assemblies can be excluded before the final library is built.

## 11. Create a primate core library

This is the defining step of the tutorial. We now construct a derived library using the reference panel. The reference panel is not simply a list of taxa to keep later. It is the empirical basis for deciding which BUSCO families behave well enough across the primate-centred reference design to serve as a conservative core set.

```bash
phyloODB metazoa.db run add-library \
  --name primate_core \
  --coverage Primates \
  --parent-library-name metazoa_odb12 \
  --accessions @PRIMATE_REFS \
  --gene-tree-source iqtree \
  --iqtree-threads 8 \
  --iqtree-flags "-m LG+G+F"
```

At a high level, this task ensures metadata and downloads exist, runs any missing BUSCO analyses, runs OrthoFinder if necessary, and then produces a filtered library based on the concordance between BUSCO families and OrthoFinder orthogroups.

Current tree behavior matters here:

- default `--gene-tree-source iqtree`: build replacement MAFFT alignments and IQ-TREE trees for the accepted orthogroups;
- those canonical core-set IQ-TREE trees are written under the OrthoFinder result directory in `IQ-TREE_Orthogroup_trees`;

For a first serious build, the default IQ-TREE mode is usually the better conservative choice. If the user wants a quicker or lighter first pass, `--fast-tree` can be used to just use OrthoFinder's default fast tree options.

After a successful build, useful outputs include:

- `library_build_metadata.json`, which records the effective core-set strategy and gene-tree source;
- `orthogroup_tree_manifest.tsv`, which records the tree paths used downstream;
- the accepted BUSCO family list for the derived library itself.

These can be found at:

`/shared/podb/libraries/primate_core/`

The reason for building a derived library at this stage is that the parent BUSCO lineage alone is often too permissive. This helps us identify BUSCO families that may have duplicated in the early evolution of our target taxa.

## 12. Hidden paralog removal

With the library in place, apply hidden paralog filtering to the target panel using the same broad reference panel.

```bash
phyloODB metazoa.db queue paralog-removal \
  --library-name primate_core \
  --ref-accessions @PRIMATE_REFS \
  --accessions @PRIMATES \
  --report-dir results/paralog_removal
```

- BUSCO says a family is present.
- Paralog removal asks whether the selected BUSCO sequence is best aligned with the same ortholog in trusted references.
- Families that match with a different sequence first count as failed, they remain in the database but are flagged and will not be accepted for export.

## 13. Decontamination approaches

We can also use the decontamination pipeline to compare BUSCO results against the results of references. This can be done to ensure that the top comparisons within a window match what we expect them to. If they match with an off-clade reference, then we may wish to exclude these from the export. E.g. in this example the refs may be Primates mammals, as well as humans and parasites.

```bash
phyloODB metazoa.db run decontamination \
  --library-name primate_core \
  --targets @PRIMATES \
  --refs @PRIMATE_DECONT_REFS \
  --rank phylum \
  --off-clade-fraction 0.10 \
  --min-buscos 20 \
  --report-path results/decontamination/primate_reference
```

In this command we would be comparing PRIMATES to PRIMATE_DECONT_REFS and ensuring that they match references with the correct phylum.

## 14. Export the final dataset

Once the custom library, paralog filtering, and decontamination runs exist, export the filtered dataset. In this tutorial the exported panel includes the PRIMATE targets.

```bash
phyloODB metazoa.db run export \
  --library-name primate_core \
  --accessions @PRIMATE_TARGETS \
  --out-dir exports/primate_core \
  --busco-report-extended
```

The default export behaviour assumes that prior filtering is meaningful and should be respected. This is exactly what makes the export suitable for downstream use.

To view output:

```cd /shared/podb/exports/primate_core```

## 15. Build BUSCO trees directly

We can ask PhyloODB to pass exported BUSCOs directly to IQ-TREE for processing. In this example we export a subset order HOMINOIDEA to build trees.

```text
Metazoa
└─ Hominoidea
  ├─ Hominidae
  │ ├─ Homininae
  │ │ ├────────────Gorilla gorilla gorilla
  │ │ ├────────────Homo sapiens
  │ │ └─ Pan
  │ │   ├──────────Pan paniscus
  │ │   └──────────Pan troglodytes
  │ └─ Pongo
  │   ├────────────Pongo abelii
  │   └────────────Pongo pygmaeus
  └─ Hylobatidae
    ├──────────────Hylobates moloch
    ├──────────────Nomascus leucogenys
    └──────────────Symphalangus syndactylus
```

```bash
phyloODB metazoa.db queue build-busco-trees \
  --library-name primate_core \
  --accessions @HOMINOIDEA \
  --sequence-type protein \
  --iqtree-threads 4 \
  --iqtree-flags "-m LG+G+F"
```

To view output:

```cd /shared/podb/exports/task_1019_20260714_114154_primate_core```


## 15. What to examine after export

The tutorial is complete when you can inspect:

- the derived library itself;
- the paralog filtering output;
- the decontamination summaries and run ids;
- the exported per-family FASTA files;
- the lineage and BUSCO reports written by export.

The essential point is that the final FASTA set is not simply “whatever BUSCO found”. It is the product of explicit choices about taxon sampling, library definition, paralog filtering, and contamination screening.
