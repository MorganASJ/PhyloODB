## 5. Selectors: defining working sets for analysis

Selectors are how PhyloODB turns a biological idea into a concrete set of assemblies.

In practice, a selector answers questions like:

- which assemblies belong in this dataset?
- should they be chosen directly, by clade, by taxid, or by a sampling rule?
- should only downloaded assemblies be used?
- should BUSCO quality, assembly metadata, root location, or proteome profile affect the choice?
- if several BUSCO runs exist for the same assembly, which run should be used by the next task?

Selectors are used throughout PhyloODB. The same selector language can drive exploratory commands and real analysis tasks:

```bash
phyloODB my_project.db list assemblies --clade Primates
phyloODB my_project.db tree --clade Primates --rank genus --quantity 1
phyloODB my_project.db queue download --clade Primates --rank genus --quantity 1
phyloODB my_project.db queue batch-busco --accessions @PRIMATE_REFS --lineage mammalia_odb12
phyloODB my_project.db queue export-library --accessions @CORE_SET
phyloODB my_project.db queue paralog-removal --accessions @CORE_SET
```

The central idea is simple: first define a working set, inspect it, then reuse it for downstream work.

That working set can be supplied in three main ways:

| Need | Use | Example |
|---|---|---|
| One-off query | Direct selector flags | `--clade Primates --rank genus --quantity 1` |
| Reusable live recipe | Selector preset | `--preset primate_refs` |
| Frozen accession panel | Stored accession set | `--accessions @PRIMATE_REFS` |

The distinction matters. A preset stores the recipe. A stored accession set stores the resolved accessions. If new assemblies are added later, a preset can resolve differently; a stored accession set stays fixed until you replace or append to it.

### 5.1 The selector mental model

Most selector workflows have four steps:

1. Build a candidate pool from accessions, clades, taxids, or stored sets.
2. Apply filters such as downloaded-only, assembly metadata, BUSCO thresholds, roots, or exclusions.
3. Optionally rank and sample candidates, for example one assembly per genus.
4. Optionally choose an analysis context, such as which BUSCO run or proteome profile should represent each selected assembly.

A direct selector can be used immediately:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --rank genus \
  --quantity 1
```

The same idea can be saved as a preset:

```bash
phyloODB my_project.db selector save primate_refs \
  --clade Primates \
  --rank genus \
  --quantity 1

phyloODB my_project.db list assemblies --preset primate_refs
```

Or the resolved accessions can be frozen as a stored set:

```bash
phyloODB my_project.db selector resolve primate_refs -S PRIMATE_REFS
phyloODB my_project.db list assemblies --accessions @PRIMATE_REFS
```

These forms are deliberately compatible. You can explore with direct selectors, save a good recipe as a preset, then freeze a particular result as a stored set when you need a stable dataset for publication or repeated downstream analysis.

### 5.2 Basic assembly selectors

The most common selector inputs are:

- `-a`, `--accessions <list>`: explicit accessions or stored accession sets such as `@CORE_SET`.
- `-c`, `--clade <name>`: assemblies under a taxonomic name.
- `-i`, `--taxid <id>`: assemblies under a specific NCBI taxid.
- `-d`, `--downloaded-only`: only assemblies already present locally.
- `--not-downloaded`: only assemblies that are registered but not downloaded.
- `--local-only`: only local/imported assemblies.
- `--not-local`: exclude local/imported assemblies.
- `--primary-only`: restrict to primary assemblies.
- `-rt`, `--root <id-or-label>`: restrict filesystem-backed selectors to a storage root.
- `-af`, `--after YYYY-MM-DD` and `-bf`, `--before YYYY-MM-DD`: release-date filters.
- `--level`: assembly level, such as `chromosome` or `complete genome`.

Examples:

```bash
phyloODB my_project.db list assemblies --clade Primates
phyloODB my_project.db list assemblies --taxid 9443
phyloODB my_project.db list assemblies --accessions GCF_000001405.40,GCF_037993035.2
phyloODB my_project.db list assemblies --accessions @MAMMALS --downloaded-only
phyloODB my_project.db list assemblies --clade Primates --level chromosome --after 2022-01-01
```

The same selectors can be used on tasks:

```bash
phyloODB my_project.db queue download --clade Primates --rank genus --quantity 1
phyloODB my_project.db queue verify-assembly --accessions @MAMMALS --repair
phyloODB my_project.db run verify-busco --clade Primates --downloaded-only
```

### 5.3 Declarative sampling

Selectors become especially useful when they are used declaratively: instead of typing a final accession list, describe the sampling rule.

For example:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --rank family \
  --quantity 1
```

This means “choose one assembly per family within Primates”. PhyloODB builds the candidate pool, ranks candidate assemblies, and selects representatives.

By default, ranking prefers approximately:

1. better BUSCO score when BUSCO-aware ranking is active;
2. RefSeq over GenBank where applicable;
3. better assembly level;
4. higher contig/scaffold N50;
5. newer release date;
6. accession as a deterministic tie-break.

Multi-stage sampling is also possible:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --ranks family,genus \
  --quantities 2,1
```

Read this as “keep two genera per family, then one assembly per genus”. This lets you control breadth and density separately.

The same sampling recipe can be saved:

```bash
phyloODB my_project.db selector save primate_sampling \
  --clade Primates \
  --ranks family,genus \
  --quantities 2,1

phyloODB my_project.db list assemblies --preset primate_sampling
phyloODB my_project.db queue download --preset primate_sampling
```

#### How positive selectors combine

When more than one positive selector is supplied, PhyloODB first builds a single candidate pool and only then applies filtering and sampling.

That means:

- `--accessions`, `--taxid`, and `--clade` contribute to the same candidate pool;
- those positive selectors are combined as a union, not as an intersection;
- filters such as `--downloaded-only`, `--filter`, `--root`, BUSCO thresholds, and exclusions are applied after the pool is built;
- rank/quantity sampling then chooses from the filtered pool.

For example:

```bash
phyloODB my_project.db list assemblies \
  --accessions @BILATERIA \
  --clade Porifera \
  --ranks phylum,class,order \
  --quantities 1,1,1
```

This means:

1. start with the accessions in `@BILATERIA`;
2. add assemblies resolved from `Porifera`;
3. apply any requested filters;
4. sample from that resulting union.

It does not mean “take only the Porifera members of `@BILATERIA`”.

This union behaviour is useful when widening a hand-curated panel with a new clade, or when combining must-have references with a broader sampling rule.

#### Exclusions and negative selection

Selectors are not limited to inclusion. To subtract assemblies from an otherwise useful rule, use:

- `--exclude-accessions`
- `--exclude-clades`
- `--exclude-taxids`

Example:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --rank genus \
  --quantity 1 \
  --exclude-accessions GCF_000001405.40
```

Stored sets can also be excluded:

```bash
phyloODB my_project.db list assemblies \
  --accessions @PRIMATES \
  --exclude-clades "Hominidae" \
  --exclude-accessions @KNOWN_BAD
```

In set-operation terms, these exclusion flags cover subtraction/difference. There is no separate `--difference` flag.

### 5.4 BUSCO-aware and proteome-aware selection

Some workflows need more than an assembly list. They also need to know which analysis result should represent each assembly.

This matters because one accession can have several BUSCO runs:

- different BUSCO libraries;
- different pipelines, such as `miniprot`, `metaeuk`, or `augustus`;
- different input formats, such as protein or genome;
- different proteome profiles, such as raw or cleaned/CD-HIT reduced proteomes.

So selectors can also describe BUSCO/proteome context.

#### BUSCO quality selectors

BUSCO thresholds can filter or rank assemblies using analysis-derived quality metrics:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --library-name mammalia_odb12 \
  --busco-complete-min 90 \
  --busco-single-min 80
```

Common flags:

- `-l`, `--library-name NAME`: BUSCO lineage/library context.
- `--library-id ID`: numeric library id.
- `--busco-complete-min N`: minimum complete BUSCO score.
- `--busco-single-min N`: minimum single-copy BUSCO score.
- `--has-busco-results`: require BUSCO results.
- `--missing-busco-results`: select accessions missing BUSCO results.

These are often used before queueing more work:

```bash
phyloODB my_project.db queue batch-busco \
  --accessions @MAMMALS \
  --lineage mammalia_odb12 \
  --missing-busco-results
```

#### Choosing which BUSCO run to use

When several BUSCO runs exist for an accession, PhyloODB resolves a run context. The main run-context flags are:

- `--busco-pipeline`
- `--require-busco-pipeline`
- `--prefer-busco-pipeline`
- `--format`
- `--require-format`
- `--prefer-format`
- `--busco-run-selection`
- `--run-id` / `--run-ids` on run-aware commands

The practical difference between require and prefer is:

- require filters shrink the valid run set;
- prefer affects ranking among valid runs;
- if no preferred run exists, PhyloODB can fall back to another valid run.

Examples:

```bash
phyloODB my_project.db list results \
  --accessions @CORE_SET \
  --library-name metazoa_core \
  --busco-pipeline metaeuk \
  --format genome \
  --busco-run-selection primary \
  --tidy
```

This asks for the current primary MetaEuk genome BUSCO run for each accession.

```bash
phyloODB my_project.db list results \
  --accessions @CORE_SET \
  --library-name metazoa_core \
  --prefer-busco-pipeline metaeuk \
  --prefer-format genome \
  --tidy
```

This prefers MetaEuk genome runs but allows fallback.

```bash
phyloODB my_project.db list results \
  --accessions @CORE_SET \
  --library-name metazoa_core \
  --all-runs \
  --tidy
```

This expands to one row per BUSCO run, which is useful when auditing or debugging run selection.

The practical rule is:

- use `primary` when you want the curated/default run;
- use `latest` when you want the newest matching run;
- use `--all-runs` when you want to inspect the full run landscape.

#### Proteome profile selectors

Proteome-aware BUSCO work distinguishes between:

- the BUSCO run's own `proteome_profile`;
- the accession's current `default_proteome_profile`.

Useful flags include:

- `--proteome-profile NAME`: require a specific profile.
- `--prefer-proteome-profile NAME`: prefer that profile while allowing fallback.
- `--isoforms-cleaned`: use the accession's current default cleaned proteome profile.
- `--raw-proteome`: shortcut for `--proteome-profile raw`.

Examples:

```bash
phyloODB my_project.db list proteome-profiles \
  --accessions GCA_000516915.1 \
  --tidy

phyloODB my_project.db list assemblies \
  --busco \
  --all-runs \
  --accessions GCA_000516915.1 \
  --meta proteome_profile,default_proteome_profile \
  --tidy
```

An asterisk in user-facing profile displays marks the accession's current default profile, for example `gff,cdhit96*`.

### 5.5 Inspecting selector results

Before using a selector to run expensive work, inspect it.

The main inspection commands are:

- `list assemblies`: assembly-centric view.
- `list results`: alias-style BUSCO summary view for assemblies.
- `list busco-runs`: one row per BUSCO run.
- `list buscos`: BUSCO family-level rows from selected runs.
- `list proteome-profiles`: registered proteome profiles and preparation provenance.
- `count assemblies`: count without printing all rows.
- `tree`: render the selected taxonomy.

Examples:

```bash
phyloODB my_project.db count assemblies --preset primate_sampling
phyloODB my_project.db list assemblies --preset primate_sampling -p
phyloODB my_project.db tree --preset primate_sampling
```

#### Output formats

By default, list commands produce TSV-like output suitable for redirection:

```bash
phyloODB my_project.db list assemblies --clade Primates > primates.tsv
```

Human-readable terminal modes:

- `-y`, `--tidy`: aligned plain terminal table.
- `-p`, `--pretty`: Rich coloured table with automatic pagination in interactive terminals.
- `--no-pager`, `--no-pagination`: print one complete pretty table and let terminal scrollback handle it.

Pretty pagination uses `[` and `]`, left/right, PageUp/PageDown, Home/End, `q`, and Ctrl-C. TSV output is never hidden, truncated, or paginated.

#### Metadata and BUSCO columns

`--meta` / `-m` appends assembly metadata columns. With no value, the defaults are `release_date`, `level`, `n50`, and `comments`.

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --downloaded-only \
  --meta \
  --tidy
```

Explicit metadata fields can also be requested:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --meta release_date,level,n50,origin \
  --tidy
```

`--busco` / `-b` appends BUSCO summary columns:

```bash
phyloODB my_project.db list assemblies \
  --accessions @CORE_SET \
  --busco \
  --library-name metazoa_core \
  --tidy
```

`list results` is a convenience route for BUSCO-oriented assembly review:

```bash
phyloODB my_project.db list results \
  --accessions @CORE_SET \
  --library-name metazoa_core \
  --tidy
```

#### Filtering and sorting output

`--filter` is a general expression mechanism spanning metadata and BUSCO fields.

- `,` means AND.
- `|` means OR.
- repeated `--filter` clauses add more AND conditions.
- BUSCO fields use the `busco.` prefix.

Examples:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --filter "level contains chrom"

phyloODB my_project.db list assemblies \
  --clade Primates \
  --library-name mammalia_odb12 \
  --filter "busco.complete>=90,busco.single_copy_complete>=80" \
  --busco \
  --tidy
```

Use `-s` / `--sort` to order output rows:

```bash
phyloODB my_project.db list assemblies \
  --clade Primates \
  --meta release_date,level \
  --sort release_date:desc

phyloODB my_project.db list results \
  --accessions @CORE_SET \
  --library-name metazoa_core \
  -s busco.complete:desc,busco.single_copy:desc
```

#### BUSCO run and family inspection

For run-centric inspection:

```bash
phyloODB my_project.db list busco-runs \
  --library-name metazoa_odb12 \
  --clade Primates \
  --prefer-busco-pipeline metaeuk \
  --prefer-format genome \
  --tidy
```

For BUSCO family-level inspection:

```bash
phyloODB my_project.db list buscos \
  --library-name metazoa_odb12 \
  --accessions GCF_000001405.40 \
  --run-id 1234 \
  --tidy
```

Run ids can be stored separately from accession sets:

```bash
phyloODB my_project.db list busco-runs \
  --accessions @CORE_SET \
  --busco-pipeline augustus \
  --ids-only \
  --store-results AUGUSTUS_RUNS

phyloODB my_project.db set busco-primary --run-ids @AUGUSTUS_RUNS
```

Stored run-id sets are distinct from stored accession sets.

### 5.6 Saving and reusing selectors

There are two common reuse patterns:

1. save the resolved accessions as a fixed set;
2. save the selector recipe as a preset.

Most projects use both.

#### Stored accession sets

A stored accession set is a fixed list of resolved accessions.

```bash
phyloODB my_project.db list assemblies \
  --clade "Pan paniscus" \
  -S BONOBO
```

The stored set can be reused with `@BONOBO`:

```bash
phyloODB my_project.db list assemblies --accessions @BONOBO
phyloODB my_project.db queue download --accessions @BONOBO
```

Useful flags:

- `-S`, `--store NAME`, `--save-set NAME`: replace a stored set with the newly resolved accessions.
- `-A`, `--append-to NAME`, `--append-set NAME`: add newly resolved accessions to an existing stored set.
- `--intersection SET`: keep only accessions shared with an explicit accession list and/or `@SET`.

Examples:

```bash
phyloODB my_project.db list assemblies --clade Primates -S PRIMATES
phyloODB my_project.db list assemblies --clade Rodentia -A MAMMALS
phyloODB my_project.db list assemblies --accessions @MAMMALS --intersection @REFSEQ_ONLY -S CURATED_MAMMALS
```

Use stored accession sets when you want the panel to remain stable.

#### Selector presets

A selector preset stores the recipe, not the resolved accession list.

Save a preset:

```bash
phyloODB my_project.db selector save primate_refs \
  --clade Primates \
  --rank genus \
  --quantity 1 \
  --description "One representative assembly per primate genus"
```

Inspect presets:

```bash
phyloODB my_project.db selector list
phyloODB my_project.db selector show primate_refs
phyloODB my_project.db selector show primate_refs --json
```

Preview or resolve a preset:

```bash
phyloODB my_project.db selector preview primate_refs
phyloODB my_project.db selector resolve primate_refs -S PRIMATE_REFS
```

Use a preset directly:

```bash
phyloODB my_project.db list assemblies --preset primate_refs
phyloODB my_project.db queue download --preset primate_refs
phyloODB my_project.db tree --preset primate_refs
```

Delete a preset:

```bash
phyloODB my_project.db selector delete primate_refs
```

Explicit selector flags supplied alongside `--preset` refine or override the saved recipe:

```bash
phyloODB my_project.db list assemblies \
  --preset primate_refs \
  --downloaded-only \
  --busco-complete-min 90
```

Use presets when you want a reusable, living selector definition.

### 5.7 Using selectors in tasks

Most expensive work in PhyloODB should be driven by inspected selectors, presets, or stored sets.

#### Download assemblies

Direct:

```bash
phyloODB my_project.db queue download \
  --clade Primates \
  --rank genus \
  --quantity 1
```

Preset:

```bash
phyloODB my_project.db queue download --preset primate_refs
```

Stored set:

```bash
phyloODB my_project.db queue download --accessions @PRIMATE_REFS
```

#### Run BUSCO

```bash
phyloODB my_project.db queue batch-busco \
  --accessions @PRIMATE_REFS \
  --lineage mammalia_odb12
```

To queue BUSCO only where results are missing:

```bash
phyloODB my_project.db queue batch-busco \
  --accessions @PRIMATE_REFS \
  --lineage mammalia_odb12 \
  --missing-busco-results
```

#### Verify and repair

```bash
phyloODB my_project.db queue verify-assembly \
  --accessions @PRIMATE_REFS \
  --repair

phyloODB my_project.db queue verify-busco \
  --preset primate_refs \
  --repair
```

#### Build or export a library

```bash
phyloODB my_project.db queue export-library \
  --accessions @CORE_SET \
  --library-name metazoa_core \
  --busco-pipeline metaeuk \
  --format genome
```

The accession selector decides which assemblies are in scope. The BUSCO selectors decide which BUSCO run context supplies the genes/results used by the export.

#### Tree building

```bash
phyloODB my_project.db tree --preset primate_refs
phyloODB my_project.db queue build-busco-trees --accessions @CORE_SET
```

#### Paralog filtering and decontamination

```bash
phyloODB my_project.db queue paralog-removal \
  --accessions @CORE_SET \
  --library-name metazoa_core

phyloODB my_project.db queue decontamination \
  --preset metazoa_targets \
  --library-name metazoa_core
```

The exact task names and options depend on the workflow, but the pattern is the same: define the assembly set, choose relevant analysis context, inspect if necessary, then queue the work.

### 5.8 Variables, defaults, and manual overrides

Stored accession sets and stored BUSCO run-id sets are kept in the project variable table. You can inspect them with:

```bash
phyloODB my_project.db list variables
```

Project variables include several categories:

- path and tool configuration, such as `GENOME_DIR`, `LIBRARIES_DIR`, `BUSCO_BINARIES_PATH`, and `ORTHOFINDER_OUTPUT_DIR`;
- daemon configuration, such as `DAEMON_MAX_THREADS` and polling intervals;
- selector defaults, such as `SELECTOR_DEFAULT_DOWNLOADED_ONLY`, `SELECTOR_DEFAULT_PRIMARY_ONLY`, and `SELECTOR_SCORE_ORDER`;
- stored accession sets, such as `METAZOA_CORE` or `MAMMALS`;
- stored BUSCO run-id sets, such as `AUGUSTUS_RUNS`;
- convenience pointers, such as `LAST` and task-specific `LAST_*` values.

Manual variable editing is available:

```bash
phyloODB my_project.db set var SELECTOR_DEFAULT_DOWNLOADED_ONLY true
phyloODB my_project.db set var MY_PANEL '["GCF_000001405.40","GCF_037993035.2"]'
```

Most users should prefer selector commands and stored-set flags over hand-editing variables, because they validate the selector and normalize values.

#### Default proteome profiles

Proteome profiles are separate from BUSCO primary runs.

- `set proteome-profile` changes the accession's default proteome profile.
- `set busco-primary` changes which BUSCO run is primary for an accession/library context.

Examples:

```bash
phyloODB my_project.db list proteome-profiles \
  --accessions GCA_000516915.1 \
  --tidy

phyloODB my_project.db set proteome-profile \
  --accessions GCA_000516915.1 \
  --profile-name cdhit96 \
  --dry

phyloODB my_project.db set proteome-profile \
  --accessions GCA_000516915.1 \
  --profile-name cdhit96
```

The default proteome profile is what proteome-aware commands use when no explicit `--proteome-profile` is supplied.

#### Manual BUSCO primary overrides

Primary BUSCO selection is automatic by default. PhyloODB maintains a `primary` run per accession, library, and purpose.

Automatic ranking prefers:

1. higher single-copy complete count;
2. then higher total complete count;
3. then lower duplicated count;
4. then later completion time;
5. then higher run id.

When you need to pin a run manually, use `set busco-primary`.

Examples:

```bash
phyloODB my_project.db set busco-primary \
  --accessions @CORE_SET \
  --busco-pipeline augustus \
  --format genome \
  --dry

phyloODB my_project.db set busco-primary \
  --accessions @CORE_SET \
  --busco-pipeline augustus \
  --format genome

phyloODB my_project.db set busco-primary \
  --accession GCA_000001405.1 \
  --run-id 535

phyloODB my_project.db set busco-primary \
  --accessions @CORE_SET \
  --run-ids @AUGUSTUS_RUNS
```

Important rules:

- if you do not supply `--refresh`, then without `--run-id` or `--run-ids` you must provide at least one run-disambiguating selector such as `--format` or `--busco-pipeline`;
- the selected run updates every primary purpose that the run can support;
- manual overrides are persistent and are not replaced by automatic refreshes.

To recompute automatic primaries rather than pin a manual override, use:

```bash
phyloODB my_project.db set busco-primary --refresh
phyloODB my_project.db set busco-primary --refresh --accessions @CORE_SET --library-name metazoa_core
phyloODB my_project.db set busco-primary --refresh --dry --accessions ACC1,ACC2
```

### 5.9 Practical workflow pattern

A reliable selector-driven workflow usually looks like this:

1. Start with a direct selector and inspect it.
2. Save the recipe as a preset if it is reusable.
3. Resolve the preset into a stored accession set when the dataset should become fixed.
4. Use the stored set or preset to queue downstream work.
5. Inspect BUSCO/proteome context before exports, core-set construction, paralog filtering, decontamination, or tree building.

For example:

```bash
# 1. Explore
phyloODB my_project.db list assemblies \
  --clade Metazoa \
  --rank phylum \
  --quantity 2 \
  --downloaded-only \
  -p

# 2. Save the recipe
phyloODB my_project.db selector save metazoa_broad_sample \
  --clade Metazoa \
  --rank phylum \
  --quantity 2 \
  --downloaded-only

# 3. Freeze today's resolved accession set
phyloODB my_project.db selector resolve metazoa_broad_sample -S METAZOA_BROAD_SAMPLE

# 4. Inspect BUSCO context
phyloODB my_project.db list results \
  --accessions @METAZOA_BROAD_SAMPLE \
  --library-name metazoa_core \
  --prefer-busco-pipeline metaeuk \
  --prefer-format genome \
  --tidy

# 5. Queue downstream work
phyloODB my_project.db queue export-library \
  --accessions @METAZOA_BROAD_SAMPLE \
  --library-name metazoa_core \
  --prefer-busco-pipeline metaeuk \
  --prefer-format genome
```

This keeps the project reproducible without forcing every command to carry a huge accession list or a long selector expression.
