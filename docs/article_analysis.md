# Reproducible commands for the Metazoa phyloODB analysis

This file records the commands and parameters used to construct the
`metazoa_core` and `metazoa_full` libraries and to run hidden-paralog
filtering. Run all commands from the project directory. The commands below use
the phyloODB 0.2.3 command syntax.

Install PhyloODB and requirements:
```bash
mamba create -n phyloodb -c conda-forge -c bioconda \
  python=3.11 pip blast busco cd-hit mafft orthofinder iqtree
conda activate phyloodb
pip install "git+https://github.com/MorganASJ/PhyloODB.git"
```

Create a path for the database
```bash
DB=~/path/to/somewhere/metazoa.db
```

## Software versions

The completed analysis within the paper used:

```text
phyloODB    0.2.3
BUSCO       6.0.0
OrthoFinder 3.1.0
MAFFT       7.526 
IQ-TREE     3.1.2 (Manually compiled Linux ARM64 build, 17 May 2026)
```

All analyses for the paper were conducted using the Isambard 3 supercomputing cluster. Isambard 3 is hosted by the University of Bristol and operated by the GW4 Alliance (https://gw4.ac.uk) and is funded by UK Research and Innovation; and the Engineering and Physical Sciences Research Council [EP/X039137/1].

## 1. Create and initialise the database

Create a new database at the desired location. If you have issues with permissions you can set the working directory to be at the same location rather than /tmp/ with `--working-dir`. As part of this process the taxdump will be downloaded. If you already have a taxdump saved you can point PhyloODB at it using `--taxdump PATH`.

```bash
phyloODB "$DB" create \
  --email morgan.jones@bristol.ac.uk \
  --api-key "$NCBI_API_KEY"
```

It is recommended you set the API key at database initiation. Please use your email for access to ENTREZ services. 

If you need to set the API key later edit the variable.

```bash
phyloODB "$DB" set var NCBI_API_KEY YOUR_KEY
phyloODB "$DB" set var LOG_LEVEL INFO
```

Add information about assemblies available on NCBI. You can point directly to the assemblies listed below or just add information about all Metazoan and outgroup taxa.

```bash
phyloODB "$DB" run add --clade Metazoa
phyloODB "$DB" run add --clade Choanoflagellata
phyloODB "$DB" run add --clade Filasterea
```

## 2. Define the two assembly sets

`METAZOA_CORE` contains 12 assemblies. `METAZOA_FULL` contains those 12 plus
28 additional assemblies (40 total). These exact accession versions, not
unversioned accessions or later replacements, were used.

```bash
phyloODB "$DB" list assemblies -a GCA_000328365.1,GCF_000001405.40,GCF_000002075.1,GCF_000002235.5,GCF_000090795.2,GCF_000150275.1,GCF_000151315.2,GCF_000485595.1,GCF_013753865.1,GCF_026151205.1,GCF_031307605.1,GCF_000188695.1 -S METAZOA_CORE

phyloODB "$DB" list assemblies -a GCA_000328365.1,GCF_000001405.40,GCF_000002075.1,GCF_000002235.5,GCF_000090795.2,GCF_000150275.1,GCF_000151315.2,GCF_000188695.1,GCF_000326865.2,GCF_000485595.1,GCF_000517525.1,GCF_001039355.2,GCF_013753865.1,GCF_019649055.1,GCF_021130785.1,GCF_021730395.1,GCF_026151205.1,GCF_029227915.1,GCF_031307605.1,GCF_032362555.1,GCF_037042905.1,GCF_037392515.1,GCF_037975245.1,GCF_040414725.1,GCF_041260155.1,GCF_049306965.1,GCF_051348905.1,GCF_054371585.1,GCF_910592395.1,GCF_947563725.1,GCF_963422355.1,GCF_963576615.1,GCF_963675165.1,GCF_963678975.1,GCF_964187855.1,GCF_964194025.1,GCF_964199315.1,GCF_964204645.1,GCF_964340395.1,GCF_964355755.1 -S METAZOA_FULL
```

This is equivalent to setting the variable:
```bash
METAZOA_CORE='GCA_000328365.1,GCF_000001405.40,GCF_000002075.1,GCF_000002235.5,GCF_000090795.2,GCF_000150275.1,GCF_000151315.2,GCF_000485595.1,GCF_013753865.1,GCF_026151205.1,GCF_031307605.1,GCF_000188695.1'
METAZOA_FULL='GCA_000328365.1,GCF_000001405.40,GCF_000002075.1,GCF_000002235.5,GCF_000090795.2,GCF_000150275.1,GCF_000151315.2,GCF_000188695.1,GCF_000326865.2,GCF_000485595.1,GCF_000517525.1,GCF_001039355.2,GCF_013753865.1,GCF_019649055.1,GCF_021130785.1,GCF_021730395.1,GCF_026151205.1,GCF_029227915.1,GCF_031307605.1,GCF_032362555.1,GCF_037042905.1,GCF_037392515.1,GCF_037975245.1,GCF_040414725.1,GCF_041260155.1,GCF_049306965.1,GCF_051348905.1,GCF_054371585.1,GCF_910592395.1,GCF_947563725.1,GCF_963422355.1,GCF_963576615.1,GCF_963675165.1,GCF_963678975.1,GCF_964187855.1,GCF_964194025.1,GCF_964199315.1,GCF_964204645.1,GCF_964340395.1,GCF_964355755.1'

phyloODB "$DB" set var --kind assemblies METAZOA_CORE "$METAZOA_CORE"
phyloODB "$DB" set var --kind assemblies METAZOA_FULL "$METAZOA_FULL"
```

The database selectors are subsequently referenced as `@METAZOA_CORE` and
`@METAZOA_FULL`. Check their expansion before starting expensive jobs:

```bash
phyloODB "$DB" list assemblies -a @METAZOA_CORE -r phylum -m -y
phyloODB "$DB" list assemblies -a @METAZOA_FULL -r phylum -m -y
```

## 3. Download the assemblies

Download the protein, GFF and genome records for both sets. We will only issue the download for the full set variable as core is a subset of the full 40 taxa dataset.

Proteome preparation should run automatically as a child task after each
download, provided the preparation defaults are set first. For this analysis,
the automatic step used GFF isoform selection followed by CD-HIT clustering at
96% amino-acid identity, named the result `gff_cdhit96`, and made that profile
the assembly default:

```bash
phyloODB "$DB" set var DEFAULT_PROTEOME_INPUT_PROFILE raw
phyloODB "$DB" set var DEFAULT_PROTEOME_USE_GFF true
phyloODB "$DB" set var DEFAULT_PROTEOME_USE_CDHIT true
phyloODB "$DB" set var DEFAULT_PROTEOME_CDHIT_IDENTITY 0.96
phyloODB "$DB" set var DEFAULT_PROTEOME_SET_DEFAULT true
phyloODB "$DB" set var DEFAULT_PROTEOME_MAX_CONCURRENT 40
phyloODB "$DB" set var DEFAULT_PROTEOME_THREADS_PER_JOB 1
```

```bash
phyloODB "$DB" run download \
  --accessions @METAZOA_FULL \
  --protein \
  --clean-isoforms \
  --clean-cdhit-identity 0.96 \
  --clean-max-concurrent 40 \
  --clean-threads-per-job 1
```

During the original full-set download, a GFF could not be obtained for
`GCF_947563725.1`. That accession was therefore prepared using CD-HIT only,
without the GFF isoform-selection step.

## 4. Create the prepared proteomes

The download commands above normally create the prepared proteomes
automatically. If preparation was disabled, interrupted, or needs to be run
separately, the equivalent manual commands are below. They derive an immutable
profile from the raw NCBI proteins, retain GFF-selected isoforms, cluster at
96% identity, and set the resulting profile as the assembly default.

```bash
phyloODB "$DB" run prepare-proteome \
  --accessions @METAZOA_FULL \
  --exclude-accessions GCF_947563725.1 \
  --input-profile raw \
  --profile-name gff_cdhit96 \
  --gff \
  --cdhit \
  --cdhit-identity 0.96 \
  --set-default \
  --max-concurrent 40 \
  --threads-per-job 1 \
  --threads 40
```

Exception: `GCF_947563725.1` used CD-HIT at 96% identity but no GFF filtering;
its selected profile was therefore named `cdhit96`.

```bash
phyloODB "$DB" run prepare-proteome \
  --accessions GCF_947563725.1 \
  --input-profile raw \
  --profile-name cdhit96 \
  --skip-gff \
  --cdhit \
  --cdhit-identity 0.96 \
  --set-default
```

The `--set-default` in these commands will ensure that the new profile becomes the primray profile for that assembly.

Inspect the selected immutable profiles with:

```bash
phyloODB "$DB" list proteome-profiles > proteome_profiles.tsv
```

The completed library logs confirm `gff_cdhit96` for 39 assemblies and
`cdhit96` for `GCF_947563725.1`.

## 5. Set analysis variables

These were the explicit resource, alignment and tree-inference settings used
for both library builds. They must be set before starting the daemon and
queuing `add-library`.

Threading settings were set as below. If you do not wish to set these manually allow `SET_MAX_THREADS_ON_START` to remain at its default setting of TRUE.
```bash
phyloODB "$DB" set var SET_MAX_THREADS_ON_START false
phyloODB "$DB" set var DAEMON_MAX_THREADS 64
phyloODB "$DB" set var DEFAULT_THREADS_BUSCO_RUN 8
phyloODB "$DB" set var DEFAULT_THREADS_ORTHOFINDER_RUN 64
phyloODB "$DB" set var DEFAULT_THREADS_MAFFT_RUN 2
phyloODB "$DB" set var DEFAULT_THREADS_IQTREE_RUN 8
phyloODB "$DB" set var DEFAULT_THREADS_PARALOG_REMOVAL 64
```

IQTREE and MAFFT were configured using below flags.
```bash
phyloODB "$DB" set var MAFFT_FLAGS '"--auto"'
phyloODB "$DB" set var IQTREE_FLAGS \
  '"-mset LG -madd C10,C20,C30,C40,C50,C60 -mrate +G -mfreq +F,C10,C20,C30,C40,C50,C60"'
```

Thus model selection was restricted to LG models, tested the C10--C60 profile
mixtures, gamma rate heterogeneity, and the listed empirical/profile frequency
options. No bootstrap flag was supplied for these library gene trees.

## 6. Download the BUSCO lineage and run BUSCO

The lineage used was `metazoa_odb12`. BUSCO was run in protein mode against the prepared proteomes.

```bash
phyloODB "$DB" run download-busco-library \
  --lineage metazoa_odb12

phyloODB "$DB" run batch-busco \
  --accessions @METAZOA_FULL \
  --lineage metazoa_odb12 \
  --format protein
```

BUSCO completeness results can be viewed directly after the runs finish. `list results` is an alias of `list assemblies --busco --downloaded`.

```bash
phyloODB "$DB" list results \
  -a @METAZOA_FULL \
  -l metazoa_odb12 \

# A compact, sorted table with --pretty / -p:
phyloODB "$DB" list results \
  -a @METAZOA_FULL \
  -l metazoa_odb12 \
  -s quality \
  -p
```

## 7. Construct the custom libraries

Libraries were produced using the queue system rather than running them individually due to the time it takes to build the trees. 

### Using the Task Daemon
Start the daemon and then end it when the tasks have completed. 


```bash
phyloODB-daemon "$DB" start
```

More control is possible such as setting max threads and using `--here` to start the task in the foreground.

```bash
phyloODB-daemon "$DB" start \
  --here \
  --threads 64 \
  --log-console \
  --log-level INFO
```

You can stop the daemon with the command:
```bash
phyloODB-daemon "$DB" stop
```

### Creating the custom libraries

The library commands used:

```bash
# METAZOA_CORE: 12 reference proteomes
phyloODB "$DB" queue add-library \
  --print-id \
  --name metazoa_core \
  --parent-library-name metazoa_odb12 \
  --coverage Metazoa \
  --coverage-taxid 33208 \
  --accessions @METAZOA_CORE \
  --gene-tree-source iqtree \
  --annotate-og-trees \
  --clean-refs-strict

# METAZOA_FULL: 40 reference proteomes
phyloODB "$DB" queue add-library \
  --print-id \
  --name metazoa_full \
  --parent-library-name metazoa_odb12 \
  --coverage Metazoa \
  --coverage-taxid 33208 \
  --accessions @METAZOA_FULL \
  --gene-tree-source iqtree \
  --annotate-og-trees \
  --clean-refs-strict
```

Each `add-library` workflow ran OrthoFinder, aligned the retained orthogroups
with MAFFT `--auto`, inferred gene trees with the IQ-TREE command above, and
used `metazoa_odb12` as the parent BUSCO library. The completed database records
12 references for `metazoa_core` and 40 for `metazoa_full`.

To automate record the task ID printed by each command (when using `--print-id`) and wait for it to finish. 
If the daemon is stopped, restart it and continue following the existing root task rather than queueing a duplicate library:

```bash
phyloODB "$DB" status TASK_ID --wait
```

### Derived busco-runs

Both commands create derived BUSCO runs under the `orthofinder` pipeline. In
these runs, qualifying single-copy BUSCO families are reclassified as
duplicated when the gene-tree evidence identifies duplication for that
accession. Because the commands explicitly use `--clean-refs-strict`, both
in-paralogs and out-paralogs can cause this accession-specific reclassification.
By comparison, `--clean-refs` reclassifies only accessions identified as having
out-paralogs.

When either cleaning mode is explicitly requested, `--set-cleaned-primary` is
the default: the new OrthoFinder-derived run becomes the primary BUSCO run for
the corresponding accession and parent-library context. Use
`--no-set-cleaned-primary` during `add-library` when the profiles should be
created without changing the existing primary selections. The
`set busco-primary` commands below are therefore needed only to restore or
change the intended profile after the builds have completed.

Strict cleaning provides a conservative BUSCO summary in which all
accession-specific duplications detected in the gene trees, including
in-paralogs, are represented in the duplicated category. This is distinct from
Core-OG family membership: the derived `metazoa_core` and `metazoa_full`
libraries exclude complex mappings and families containing out-paralogs, but
retain families containing only in-paralogs because those copies preserve the
same orthology relationship to sequences in the other taxa.

To set the primary run after the fact use:
```bash

phyloODB "$DB" set busco-primary \
  -a @METAZOA_CORE \
  --orthofinder-target-library metazoa_core

phyloODB "$DB" set busco-primary \
  -a @METAZOA_FULL \
  --orthofinder-target-library metazoa_full
```

*Use `--dry` to preview changes*

Note that for the assemblies in METAZOA_CORE you should have three busco-run records for metazoa_odb12 (and derived libraries):
1. Original 'protein' pipeline results.
2. Derived 'orthofinder' pipeline results from the reciprocal validation of METAZOA_CORE (strict)
3. Derived 'orthofinder' pipeline results from the reciprocal validation of METAZOA_FULL (strict)

These cannot be expected to be the same as the orthogroup definitions may change between runs and the gene trees themselves may differ due to the additional taxa.

### Additional output

Each custom library has its own directory beneath the active libraries storage root:
```text
<libraries-root>/
├── metazoa_core/
└── metazoa_full/
```
A completed IQ-TREE build with --annotate-og-trees has approximately the following structure:
```text
<libraries-root>/<library-name>/
├── cleaned_busco_families.json
├── cleaned_busco_families_<timestamp>.json
├── library_build_metadata.json
├── orthogroup_tree_manifest.tsv
├── annotated-og-trees/
│   ├── OG0000001.nex
│   ├── OG0000002.nex
│   └── ...
├── cleaned_reference_proteomes/
│   ├── GCA_....tsv
│   ├── GCF_....tsv
│   └── ...
└── core_set_analysis/
    ├── busco_fastas/
    │   ├── <BUSCO-family>.fasta
    │   └── ...
    ├── <library-name>_busco_to_orthogroup_map.tsv
    ├── <library-name>_busco_to_orthogroup_exact_map.tsv
    ├── <library-name>_busco_to_orthogroup_1to1.tsv
    ├── <library-name>_busco_to_orthogroup_1to1_occupancy.tsv
    ├── <library-name>_busco_to_orthogroup_og_to_busco_families.tsv
    ├── <library-name>_busco_to_orthogroup_busco_family_to_ogs.tsv
    ├── <library-name>_busco_to_orthogroup_map_with_paralog_class.tsv
    ├── <library-name>_busco_to_orthogroup_species_paralog_status.tsv
    ├── <library-name>_busco_to_orthogroup_unmapped_buscos.tsv
    ├── <library-name>_good_busco_families.txt
    └── paralogs/
        ├── og_paralog_summary.tsv
        ├── OG0000001_inparalogs.txt
        ├── OG0000001_outparalogs.txt
        └── ...
```
Some directories are conditional. annotated-og-trees/ is created only when --annotate-og-trees is used. 
cleaned_reference_proteomes/ is created only when reference cleaning is requested with --clean-refs or --clean-refs-strict.

#### Library-level files

`cleaned_busco_families.json` is the canonical list of BUSCO families retained in the custom library. The timestamped copy records the list written by that particular build.

`library_build_metadata.json` is the best starting point when auditing a library. It records the parent library, reference-cleaning mode, minimum species threshold, effective tree source, accepted-family rule, number of retained families, and paths to the principal analysis outputs.

#### BUSCO-to-orthogroup tables

The files under `core_set_analysis/` describe successive stages of family selection.

`*_map.tsv` is the full sequence-level BUSCO-to-OrthoFinder mapping. Each row associates a BUSCO sequence with an orthogroup sequence. It includes the source BUSCO file and identifier, BUSCO family, orthogroup, OrthoFinder sequence identifier and sequence length.

`*_og_to_busco_families.tsv` is a compact orthogroup-oriented mapping. Use it to identify orthogroups containing more than one BUSCO family.

`*_busco_family_to_ogs.tsv` is the reverse mapping. Use it to identify BUSCO families distributed across multiple orthogroups.

`*_1to1.tsv` summarizes BUSCO families that map to only one orthogroup. `Exact_1_to_1=YES` means that the orthogroup maps back to that BUSCO family and contains the expected multiset of BUSCO-matched sequences. `Relaxed_only=YES` means that the family maps to one orthogroup but fails the exact reciprocal criteria.

`*_exact_map.tsv` is the sequence-level subset of `*_map.tsv` containing only exact BUSCO-family-to-orthogroup pairs.

`*_1to1_occupancy.tsv` applies the minimum-species criterion to the exact pairs. It reports the BUSCO family, orthogroup, number of represented reference species and whether the pair was accepted or rejected.

`*_unmapped_buscos.tsv` summarizes BUSCO records that could not be connected to an assigned OrthoFinder orthogroup. BUSCOs found only among OrthoFinder’s unassigned genes are counted as unmapped because the comparison scans `Orthogroup_Sequences`.

#### Paralog tables

`paralogs/og_paralog_summary.tsv` gives one row per examined orthogroup, with the number of in-paralog and out-paralog sets and one of four classifications:

- `No Paralogs`
- `Only In-Paralogs`
- `Only Out-Paralogs`
- `Both In- and Out-Paralogs`

The accompanying `OG..._inparalogs.txt` and `OG..._outparalogs.txt` files list the actual tree-tip sets supporting those classifications.

`*_species_paralog_status.tsv` converts the orthogroup-level evidence into accession-specific calls. Its columns identify the BUSCO family, orthogroup and accession, followed by `Has_In_Paralogs` and `Has_Out_Paralogs`. This is the most useful table for determining why a BUSCO was reclassified as duplicated in one reference profile but not another.

`*_map_with_paralog_class.tsv` adds the overall orthogroup paralog classification to the sequence-level mapping. `Clean_1to1_NoParalogs=YES` indicates an exact, occupancy-passing pair with either no paralogs or only in-paralogs. Orthogroups containing out-paralogs are excluded from the final custom-library family set.

`*_good_busco_families.txt` is the plain-text final accepted-family list used to produce `cleaned_busco_families.json`.

When reference cleaning is enabled, `cleaned_reference_proteomes/<accession>.tsv` provides an accession-centred audit. For every BUSCO family it reports:

- its orthogroup and mapping category;
- original and rewritten BUSCO status;
- source gene name;
- ortholog, in-paralog or out-paralog classification;
- the supporting paralog file;
- the gene tree used.

This is the most direct file for explaining an individual `single-copy → duplicated` change.

#### Alignments and unannotated gene trees

The large OrthoFinder and IQ-TREE products are stored beneath the active OrthoFinder root rather than copied into every library directory. Their exact paths are recorded in `orthogroup_tree_manifest.tsv` and `library_build_metadata.json`.

For an IQ-TREE build, the relevant OrthoFinder result directory normally contains:

```text
<orthofinder-result>/
├── Orthogroup_Sequences/
├── MAFFT_Alignments/
├── IQ-TREE_Orthogroup_trees/
├── IQTREE_Metadata/
└── Resolved_Gene_Trees/
```

`Orthogroup_Sequences/` contains the raw OrthoFinder FASTAs. `MAFFT_Alignments/` contains the replacement alignments. `IQ-TREE_Orthogroup_trees/` contains the replacement trees used by this build. OrthoFinder’s original `Resolved_Gene_Trees/` directory is retained unchanged.

#### Annotated orthogroup trees

The annotated trees requested by `--annotate-og-trees` are written as NEXUS files under:

```text
<libraries-root>/<library-name>/annotated-og-trees/
```

For example:

```text
<libraries-root>/metazoa_core/annotated-og-trees/OG0000001.nex
```

These can be opened in NEXUS-compatible tree viewers such as FigTree. Tree-tip labels have the general form:

```text
<sequence>_<taxon>_<phylum>_[<BUSCO-status>][*HP]
```

The BUSCO status markers are:

- `BUSCO_SC`: single-copy BUSCO;
- `BUSCO_DUP`: duplicated BUSCO;
- `BUSCO_FRAG`: fragmented BUSCO;
- `BUSCO_MISSING`: missing BUSCO, although missing records normally have no corresponding sequence tip;
- `NON_BUSCO`: an orthogroup sequence not identified as the BUSCO sequence.

Terminal branch colours use the following key:

| Colour | Meaning |
|---|---|
| Green | BUSCO sequence with no stronger paralog annotation |
| Blue | In-paralog |
| Red |  Out-paralog |
| Red |  Hidden paralog, additionally marked `[*HP]` (will not see in this stage) |
| Grey |  Orthogroup sequence not identified as a BUSCO sequence |

Annotations have a priority order. A hidden-paralog annotation takes precedence over out-paralog, which takes precedence over in-paralog, non-BUSCO and ordinary BUSCO colouring.

## 8. Run hidden-paralog filtering

The study treated all 40 full-set proteomes as targets and used the 12 core
proteomes as the reference panel. Duplicated BUSCO families were included in the assessment. The
reference sequences for each family were selected by the lower-quartile mode,
and existing compatible results were reused.

First we set the 12 taxa reference to use the reciprocally validated results:
```bash
phyloODB "$DB" set busco-primary \
  -a @METAZOA_CORE \
  --orthofinder-target-library metazoa_core \
```

The phyloODB command was:

```bash
phyloODB "$DB" queue paralog-removal \
  --print-id \
  --library-name metazoa_odb12 \
  --accessions @METAZOA_FULL \
  --targets @METAZOA_FULL \
  --ref-accessions @METAZOA_CORE \
  --mode lower-quartile \
  --include-duplicated \
  --reuse-existing \
  --run-label metazoa_full_vs_metazoa_core \
  --max-concurrent 32
```

If the daemon stops before filtering completes, restart the daemon and follow
the existing task rather than creating another run:

```bash
phyloODB "$DB" status TASK_ID --wait
```

When all analyses are complete stop the daemon.

```bash
phyloODB-daemon "$DB" stop
```

## 9. Review final results windows

Results of base metazoa_odb12 with no alterations:
```bash
phyloODB "$DB" list results -a @METAZOA_FULL \
  -l metazoa_odb12 \
  -r phylum,class \
  -m proteome-profile \
  --busco-pipeline proteome \
  -s quality \
  -p
```

Results of metazoa_core reciprocally validated:
```bash
phyloODB "$DB" list results -a @METAZOA_CORE \
  -l metazoa_core \
  -r phylum,class \
  -m proteome-profile,orthofinder-target-library \
  --busco-pipeline orthofinder \
  --filter 'orthofinder_target_library=metazoa_core' \
  -s quality \
  -p
```

Results of metazoa_full reciprocally validated:
```bash
phyloODB "$DB" list results -a @METAZOA_FULL \
  -l metazoa_full \
  -r phylum,class \
  -m proteome-profile,orthofinder-target-library \
  --busco-pipeline orthofinder \
  --filter 'orthofinder_target_library=metazoa_full' \
  -s quality \
  -p
```

To show source run of hidden paralog filtering add `--extended-decontamination-headers`. 
To create .tsv files remove the `--pretty` (`-p`) and direct output to a file with `>`. 

## 10. Export a dataset for downstream analysis

This export step was not part of the analyses conducted in the paper, but shows
how the completed database can be used to produce a practical phylogenomic
dataset. In this example, the 40-taxon panel is exported through the
`metazoa_core` library. This applies the family-level reciprocal validation
performed with the 12 reference taxa and the sequence-level hidden-paralog
decisions generated in Section 8.

Before exporting, confirm that the 12 reference taxa use the OrthoFinder
profiles produced by `metazoa_core`. The remaining taxa retain their selected
primary protein runs, normally the `metazoa_full` OrthoFinder profiles selected
by the preceding build. Preview the reference change with `--dry` if required:

```bash
# Set the 28 proteomes to their proteome pipeline versions
phyloODB "$DB" set busco-primary \
  --accessions @METAZOA_FULL \
  --exclude-accessions @METAZOA_CORE \
  --busco-pipeline Proteome \
  --library-name metazoa_odb12

# Set the 12 reference proteomes to their reciprocally validated versions
phyloODB "$DB" set busco-primary \
  --accessions @METAZOA_CORE \
  --busco-pipeline Orthofinder \
  --orthofinder-target-library metazoa_core \
  --library-name metazoa_odb12
```

Choose an explicit output directory and export protein sequences:

```bash
EXPORT_DIR="$PWD/exports/metazoa_core_40taxa_protein"

phyloODB "$DB" run export \
  --library-name metazoa_core \
  --accessions @METAZOA_FULL \
  --require-paralog-filtering \
  --min-occupancy 0.3 \
  --header 'ACCESSION:TAXON:PHYLUM:BUSCO' \
  --out-dir "$EXPORT_DIR"
```

The important choices in this command are:

- `--library-name metazoa_core` restricts the export to BUSCO families retained
  by the 12-taxon reciprocal-validation analysis. Use `metazoa_full` instead to
  export the smaller family set retained by the independent 40-taxon Core-OG
  analysis.
- Export uses the curated primary BUSCO run by default. Use
  `--busco-run-selection latest` to change this behaviour.
- `--require-paralog-filtering` makes the export fail if the selected taxa lack
  the hidden-paralog filtering state created in Section 8. Without this flag,
  filtering is used when available but missing filtering results do not prevent
  export.
- Duplicated BUSCO calls are not rescued or included, so this example produces
  a strictly single-copy export after the active hidden-paralog filtering has
  been applied. This means that accession-specific in-paralog calls represented
  as duplicated in the strict OrthoFinder profiles are omitted from the
  exported FASTAs for those accessions.
- `--min-occupancy 0.3` retains a BUSCO gene family only when it is represented in at
  least 30% of the surviving selected taxa. Taxon-occupancy filtering is currently 
  disabled because `--min-taxa-occupancy` is omitted and defaults to `0`.
- `--header 'ACCESSION:TAXON:PHYLUM:BUSCO'` produces stable FASTA headers that
  record the assembly accession, taxon, phylum, and BUSCO family for downstream
  alignment and tree software.

The export directory contains:

```text
metazoa_core_40taxa_protein/
├── busco_families/
│   ├── <BUSCO-family>.fasta
│   └── ...
├── lineage.csv
├── busco_report.tsv
├── busco_family_matrix.tsv
├── taxa_occupancy.tsv
├── export_filter_report.tsv
├── export_parameters.txt
└── export_task.log
```

`busco_families/` contains one protein FASTA per retained Core-OG family.
`lineage.csv` records the exported taxa, `busco_report.tsv` summarizes their
BUSCO results, and `busco_family_matrix.tsv` gives the family-by-accession
status matrix. `taxa_occupancy.tsv` records taxon-level occupancy decisions,
while `export_filter_report.tsv` records which sequences and families survived
the active filtering and occupancy rules. `export_parameters.txt` and
`export_task.log` provide the settings and execution record needed to audit the
export.

The family FASTAs can then be supplied to a preferred alignment and
phylogenetic workflow. PhyloODB can also export the same family set and queue
MAFFT and IQ-TREE tasks directly:

```bash
phyloODB "$DB" queue build-busco-trees \
  --print-id \
  --library-name metazoa_core \
  --accessions @METAZOA_FULL \
  --require-paralog-filtering \
  --min-occupancy 0.3 \
  --header 'ACCESSION:TAXON:PHYLUM:BUSCO' \
  --mafft-flags "--auto" \
  --iqtree-flags "-mset LG -madd C10,C20,C30,C40,C50,C60 -mrate +G -mfreq +F,C10,C20,C30,C40,C50,C60" \
  --out-dir "$PWD/exports/metazoa_core_40taxa_trees"
```

As with the library builds, it is recommended that a large task such as this
is run via the daemon task scheduler. Record the printed task ID and use
`phyloODB "$DB" status TASK_ID --wait` to follow the existing task to completion. 

A nucleotide export cannot be obtained from these protein-only BUSCO runs merely 
by changing `--sequence-type`: it requires nucleotide-capable genome BUSCO runs, 
such as those created using MetaEuk or Augustus pipelines, together with compatible 
filtering and explicit BUSCO-run selection. The examples in this document and the paper
only used proteomes for simplicity when comparing reciprocally validated datasets. In
practice you would likely complement a core library with unannotated genome data.
