# PhyloODB Command Reference

This file is a compact reference to the current command tree and the main built-in tasks. It is not a substitute for `--help`, but it is intended to give a stable overview of how the program is organised.

Use this reference for: fast command lookup. If you are new to the project, start with the [Quick Start](./QUICKSTART.md). If you want fuller explanations and operational guidance, use the [User Manual](./MANUAL.md). If you want a worked example, use the [Tutorial](../tutorial/README.md).

## 1. Top-level commands

```text
phyloODB <database> {list,watch,tree,storage,discover,count,assemblies,selector,queue,status,run,set,info,create,migrate,clear,purge,reset,kill,cancel}
```

### `create`
Create and initialise a database.

```bash
phyloODB project.db create [--force] [--taxdump PATH] [--retain-taxdump] [--working-dir DIR]
```

### `list`
List tasks, queue/errors, assemblies, libraries, variables, ranks, or metadata.

Common forms:

```bash
phyloODB project.db list tasks
phyloODB project.db list queue
phyloODB project.db list queue --pretty --watch
phyloODB project.db list errors --stack
phyloODB project.db list variables
phyloODB project.db list assemblies --clade Primates --rank family --quantity 1
phyloODB project.db list libraries --parent-name metazoa_odb12
phyloODB project.db list busco-runs --accessions @PANEL --prefer-busco-pipeline metaeuk --prefer-format genome
phyloODB project.db list busco-runs --accessions @PANEL --busco-pipeline augustus --ids-only --store-results AUGUSTUS_RUNS
phyloODB project.db list roots
```

### `watch`
Exact convenience alias for the corresponding list command with `--watch`.

Common forms:

```bash
phyloODB project.db list queue --watch
phyloODB project.db list queue --watch --all --refresh 1
phyloODB project.db list errors --watch
phyloODB project.db list errors --watch --stack --limit 50
```

Pretty terminal lists paginate automatically when they exceed the available
height. Use `[`/`]`, Left/Right, or Page Up/Page Down to move, Home/End to jump,
and `q` to close the pager. Add `--no-pager` to print the complete pretty table
without opening the interactive pager.

### `migrate`

Upgrade an older database explicitly. Normal list, daemon, and task commands
never migrate a database as a side effect.

```bash
phyloODB project.db migrate
```

The current schema is version 4. It includes transactional filesystem-operation
recovery, storage/proteome state, and selector presets.

### `storage`
Manage storage roots and move bound data.

Common forms:

```bash
phyloODB project.db storage roots
phyloODB project.db storage add-root --kind genomes --base-path /mnt/hdd/genomes --label hdd
phyloODB project.db storage rename-root hdd --label archive
phyloODB project.db storage activate-root 7
phyloODB project.db storage deactivate-root 1
phyloODB project.db storage move-genomes --accessions @OLD_SET --to-root 7 --apply
phyloODB project.db storage move-libraries --library-name metazoa_core --to-root 8 --apply
phyloODB project.db storage recover
phyloODB project.db storage recover --operation-id 12 --apply
```

Notes:

- `storage add-root` creates missing parent directories and verifies the resolved directory is readable and writable before changing the database
- `storage rename-root ROOT --label LABEL` accepts a root id or exact current label and changes only its unique label
- `storage move-genomes --apply` queues a genome-move finalization task
- with default `--verify`, that task copies first, updates the binding, runs `verify-assembly --repair --tidy` first, then `verify-busco --repair --reingest`, suspends until they complete, and only then deletes the original source
- if verification fails, the task rolls the genome binding back to the original source and removes the copied destination

### `discover`
Discover assemblies and BUSCO runs from registered genomes roots.

Common forms:

```bash
phyloODB project.db discover
phyloODB project.db discover --root HDD_GENOMES
phyloODB project.db discover --path /mnt/hdd/genomes/subset_a
phyloODB project.db discover --root HDD_GENOMES --dry-run
phyloODB project.db discover --root HDD_GENOMES --overwrite
```

Rules:

- plain `discover` scans all registered genomes roots, including inactive ones
- `--root` accepts a root id or exact unique label and scans one registered genomes root
- `--path` must be inside a registered genomes root and scans only that subtree
- out-of-root paths are rejected; register the root first with `storage add-root`
- discovery already ingests BUSCO runs found under accession folders
- `verify-busco --repair` is still recommended after discovery to reconcile metadata, artifacts, primary assignments, and duplicates

### `assemblies`
Resolve a selector to accessions and optionally store it.

```bash
phyloODB project.db assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_SET
phyloODB project.db assemblies --clade Rodentia --quantity 5 --append-to MAMMAL_SET
phyloODB project.db assemblies --accessions @MAMMAL_SET --intersection @REFSEQ_ONLY --store CURATED_MAMMALS
```

Notes:

- `-S`, `--store NAME`, `--save-set NAME` replaces `NAME` with the final resolved accession set.
- `-A`, `--append-to NAME`, `--append-set NAME` unions the final resolved accession set into `NAME`.
- `--intersection SET` intersects the resolved selector with an explicit set of accessions and/or `@VARIABLE` references.
- subtraction/difference is already covered by `--exclude-accessions`, `--exclude-clades`, and `--exclude-taxids`.
- variable-targeting flags tolerate `@NAME` input, but stored variable names themselves cannot contain `@`.

### `selector`
Save, inspect, and resolve named selector presets.

```bash
phyloODB project.db selector save primate_refs --clade Primates --rank genus --quantity 1
phyloODB project.db selector list
phyloODB project.db selector show primate_refs
phyloODB project.db selector preview primate_refs
phyloODB project.db selector resolve primate_refs --store PRIMATE_REFS
phyloODB project.db selector delete primate_refs
```

Notes:

- selector presets store the selector recipe, not the resolved accession list
- `--preset NAME` reruns the recipe against the current database state
- `--accessions @NAME` reuses a frozen accession panel stored as a variable
- explicit selector flags supplied with `--preset` override the stored preset fields

### `count`
Count assemblies matching a selector.

```bash
phyloODB project.db count assemblies --clade Primates --downloaded-only
```

### `tree`
Render a selector-defined taxonomic tree and optionally write Newick output.

```bash
phyloODB project.db tree --accessions @PANEL --colour-by-ranks phylum,class -o panel_tree.nwk
```

### `queue`
Queue a task for daemon execution.

```bash
phyloODB project.db queue <task> [task options]
```

Useful queue options:

- `--schedule EXPR`: release condition. Forms: `started:<task>`, `finished:<task>`, `succeeded:<task>`, `failed:<task>`, `delay:<Ns|Nm|Nh>`, or `at:HH:MM`; combine alternatives with `|` and requirements with `&`.
- `--parent ID`
- `--as-subtask-of SELECTOR`
- `--print-id`
- `--output-json`
- `--payload-file PATH`
- `--json '{...}'`

For shell workflows, `queue` exits `0` when the task was successfully submitted. Use `--print-id` to capture the queued task id:

```bash
task_id=$(phyloODB project.db queue download --print-id --accessions @PANEL)
```

### `status`
Check a task by id or selector and return a shell-friendly exit code.

```bash
phyloODB project.db status "$task_id"
phyloODB project.db status LAST --wait --quiet
phyloODB project.db status LAST_DOWNLOAD --json
```

Exit codes:

- `0`: completed successfully
- `1`: failed or errored
- `2`: not complete yet, including timeout while waiting
- `3`: task id or selector not found
- `4`: status check failed

Useful status options:

- `--wait`
- `--interval SECONDS`
- `--timeout SECONDS`
- `--quiet`
- `--json`

### `run`
Queue one root task, then follow that task chain in a temporary foreground daemon.

```bash
phyloODB project.db run <task> [task options]
```

Notes:

- `run` now follows suspending parent tasks through any subtasks they queue, then exits when that task chain finishes
- unrelated queued tasks are ignored; `run` does not drain the whole database queue
- `--threads N` is a hard total budget for the whole chain, not just the root task
- `--threads N` must be less than or equal to the detected available thread count
- logs stream to the terminal by default during `run`, but scheduler lifecycle messages are hidden unless you add `--show-scheduler`
- add `--quiet` to suppress console log streaming entirely while preserving normal log-file output

### `set`
Set database variables, proteome profile defaults, or BUSCO primary overrides.

```bash
phyloODB project.db set var VAR value
phyloODB project.db set var --json variables.json
phyloODB project.db set env VAR value
phyloODB project.db set proteome-profile ...
phyloODB project.db set busco-primary ...
```

Notes:

- `set env` is an alias of `set var`
- `set proteome-profile` sets the default proteome profile for matched accessions
- `set proteome-profile --dry` previews the accessions and target profile without writing
- the default proteome profile affects future proteome-aware operations when no explicit `--proteome-profile` selector is supplied
- `set busco-primary` creates persistent manual BUSCO primary overrides
- `set busco-primary --refresh` recomputes automatic BUSCO primaries for all matched accession/library pairs with completed BUSCO runs
- `set busco-primary --dry` previews current versus proposed primary run changes
- without `--refresh`, `set busco-primary` requires either `--run-id`/`--run-ids` or at least one run-disambiguating selector such as `--format` or `--busco-pipeline`
- `--refresh` preserves existing manual overrides (`manual_override`) and rewrites only non-manual primary assignments
- `--refresh` recomputes all three primary purposes: `default`, `export_protein`, and `export_nucleotide`
- `--refresh` cannot be combined with manual pinning flags such as `--run-id`, `--run-ids`, `--format`, or `--busco-pipeline`
- the chosen run updates every primary purpose that the run actually supports; unsupported purposes are left unchanged

### `info`
Show database stats or describe a task.

```bash
phyloODB project.db info
phyloODB project.db info add-library
```

### `clear`, `purge`, `reset`, `kill`, `cancel`
Administrative and destructive or semi-destructive operations. Use deliberately.

Common purge forms:

```bash
phyloODB project.db purge variables --custom-only
phyloODB project.db purge roots --inactive-only
phyloODB project.db purge busco-primary --accessions @PANEL
phyloODB project.db purge busco-primary --accessions @PANEL --apply
```

## 2. Selector grammar in brief

Selectors appear in assembly inspection commands and in many tasks.

Core selector fields:

- `--preset NAME`
- `-a`, `--accessions` (`--accession` remains available)
- `-c`, `--clade`
- `-i`, `--taxid`
- `--exclude-accessions`
- `--exclude-clades`
- `--exclude-taxids`
- `-d`, `--downloaded-only`
- `--not-downloaded`
- `--local-only`
- `--not-local`
- `--primary-only`
- `-af`, `--after YYYY-MM-DD`
- `-bf`, `--before YYYY-MM-DD`
- `-rt`, `--root`
- `--level {complete genome,chromosome,scaffold,contig}`

Rule-based selection:

- `-r`, `--rank` with `-q`, `--quantity`
- `-r`, `--ranks` with `-q`, `--quantities`
- `--sample-strategy {rank,random}`
- `--sample-seed`

BUSCO-aware selection:

- `-li`, `--library-id` or `-l`, `--library-name`
- `--busco-pipeline` or `--require-busco-pipeline`
- `--prefer-busco-pipeline`
- `--format` or `--require-format`
- `--prefer-format`
- `--proteome-profile`
- `--prefer-proteome-profile`
- `--isoforms-cleaned`
- `--raw-proteome`
- `--busco-run-selection`
- `--busco-complete-min`
- `--busco-single-min`
- `--has-busco-results`
- `--missing-busco-results`

Semantics:

- require flags shrink the candidate BUSCO run set
- prefer flags only affect ranking inside the already valid candidate set
- `--busco-pipeline` and `--require-busco-pipeline` are the same hard filter
- `--proteome-profile` is a hard BUSCO/proteome filter; `--prefer-proteome-profile` is a soft preference
- `--isoforms-cleaned` is a shortcut for "use the accession's current default cleaned proteome profile"
- `--raw-proteome` is a shortcut for `--proteome-profile raw`
- hard BUSCO filters do not override `--busco-run-selection`; `--busco-pipeline augustus --busco-run-selection primary` means "use Augustus only where Augustus is already the stored primary", while `--busco-pipeline augustus --busco-run-selection latest` means "use the latest Augustus run if one exists"
- `--format` here means BUSCO input format (`protein` or `genome`)
- `--export-format` means required output capability of the selected BUSCO run (`protein` or `nucleotide`)

Filtering by earlier screening:

- `--paralog-filtered`
- `--not-paralog-filtered`
- `--min-hidden-paralogs`
- `--max-hidden-paralogs`
- `--decontaminated`
- `--not-decontaminated`
- `--contaminated`
- `--decontamination-run`
- `--ignore-contaminated-assemblies`
- `--include-contaminated-assemblies`

Stored sets:

- `-S`, `--store NAME`, `--save-set NAME` stores a resolved accession set.
- `-A`, `--append-to NAME`, `--append-set NAME` appends a resolved accession set to an existing stored set.
- `@NAME` reuses that stored set later.

Notes:

- one-dash multi-letter aliases such as `-mc`, `-li`, `-rt`, `-af`, and `-bf` are intentional in PhyloODB
- `-q` is reused for both `--quantity` and `--quantities`, and `-r` is reused for both `--rank` and `--ranks`

### `list assemblies` output controls

- default output: TSV, suitable for `>` or `>>`
- `--no-header`: omit the column header row, useful when appending TSV rows with `>>`
- `-y`, `--tidy`: aligned terminal table
- `-p`, `--pretty`: colour in interactive terminals
- `-m`, `--meta [FIELDS]`: append metadata columns
- `-b`, `--busco`: append BUSCO summary columns
- `--filter EXPR`: metadata/BUSCO filter expressions. Operators: `=`, `!=`, `<`, `<=`, `>`, `>=`, `~`/`contains`, `!~`/`not contains`, `in`, `not in`, `exists`, `missing`. Use `,` for AND and `|` for OR.
- `-s`, `--sort FIELD[:asc|desc][,...]`: sort output rows before rendering. Examples: `--sort accession`, `--sort latest`, `--sort quality`, `--sort busco.complete:desc,busco.single_copy_complete:desc`.
- `--output-path PATH`: write directly to a file

Examples:

```bash
phyloODB project.db list assemblies -c Primates -d > primates.tsv
phyloODB project.db list assemblies -c Primates -d -y -p
phyloODB project.db list assemblies -c Primates -d -m
phyloODB project.db list assemblies --clade Primates --filter "level contains chrom"
phyloODB project.db list assemblies -c Primates -l metazoa_odb12 --filter "busco.complete>=90,busco.single_copy_complete>=80" -b -y
phyloODB project.db list assemblies -c Primates -l metazoa_odb12 --busco --sort quality
```

The default `--meta` fields are `release_date`, `level`, `n50`, and `comments`. `origin` is also available explicitly and records values such as `local`, `refseq`, or `genbank`.

Use `list metadata` to see a compact accepted-terms list for `--filter` and `--sort`. BUSCO fields use the `busco.` prefix consistently, for example `busco.complete` and `busco.single_copy_complete`; `quality` is a sort/filter alias for BUSCO single-copy completeness, and `latest` sorts by `release_date:desc`.

Proteome-profile-related metadata fields:

- `proteome_profile`: the proteome profile used by the BUSCO run shown on that row when `--busco` is active
- `default_proteome_profile`: the accession's current default proteome profile

Useful examples:

```bash
phyloODB project.db list busco-runs --accessions GCA_000516915.1
phyloODB project.db list proteome-profiles --accessions GCA_000516915.1
phyloODB project.db list assemblies --busco --all-runs --meta proteome_profile,default_proteome_profile
```

### `prepare-proteome` defaults

`prepare-proteome` creates a derived proteome profile without mutating the raw proteome.

When recipe options are omitted, `prepare-proteome` reads the database
`DEFAULT_PROTEOME_*` variables. New databases default to:

- `DEFAULT_PROTEOME_USE_GFF=true`
- `DEFAULT_PROTEOME_USE_CDHIT=false`
- `DEFAULT_PROTEOME_CDHIT_IDENTITY=0.96`

Explicit flags such as `--skip-gff`, `--skip-cdhit`, `--no-skip-cdhit`, and
`--cdhit-identity` override these defaults for one run.

Examples:

```bash
phyloODB project.db run prepare-proteome --accessions GCA_000516915.1 --profile-name cdhit96
phyloODB project.db run prepare-proteome --accessions GCA_000516915.1 --profile-name cdhit98 --cdhit-identity 0.98
phyloODB project.db run prepare-proteome --accessions GCA_000516915.1 --profile-name gff --skip-cdhit
phyloODB project.db set var DEFAULT_PROTEOME_USE_CDHIT true
phyloODB project.db set var DEFAULT_PROTEOME_CDHIT_IDENTITY 0.98
```

Profile naming notes:

- stored profile names should normally reflect the preparation recipe, for example `cdhit96`, `gff`, or `gff_cdhit96`
- in BUSCO-consuming commands, `--isoforms-cleaned` means "use the accession's current default cleaned profile"

Protein-mode BUSCO examples:

```bash
phyloODB project.db queue batch-busco --accessions @METAZOA_CANDIDATES --lineage metazoa_odb12 --format protein
phyloODB project.db queue batch-busco --accessions @METAZOA_CANDIDATES --lineage metazoa_odb12 --format protein --isoforms-cleaned
phyloODB project.db queue batch-busco --accessions @METAZOA_CANDIDATES --lineage metazoa_odb12 --format protein --raw-proteome
phyloODB project.db queue batch-busco --accessions @METAZOA_CANDIDATES --lineage metazoa_odb12 --format protein --proteome-profile gff_cdhit99
```

Selector behavior for protein-mode BUSCO:

- no profile option: use the accession's default proteome profile
- `--isoforms-cleaned`: use the accession's default cleaned profile
- `--raw-proteome`: use `raw`
- `--proteome-profile NAME`: use that exact stored profile name

### Variables and manual configuration

Variables are stored in the database and can be inspected or changed through:

```bash
phyloODB project.db list variables
phyloODB project.db list variables --json > variables.json
phyloODB project.db list variables --kind assemblies --json > panels.json
phyloODB project.db set var VAR value
phyloODB project.db set var --json variables.json
phyloODB project.db set env VAR value
```

Notes:

- keys are uppercased by `set var`
- values are parsed as JSON when possible
- variables are stored with an explicit kind: `env`, `assemblies`, or `busco_runs`
- `list variables --json` exports a kinded JSON document with `environment`, `assemblies`, and `busco_runs` objects
- `set var --json PATH` imports those objects into the matching stored kinds
- `@NAME` is reference syntax for reading a stored variable; it is not part of the variable name itself
- variable-writing flags such as `--store`/`--save-set`, `--append-to`/`--append-set`, and `--store-results` accept either `NAME` or `@NAME`, but stored names may not contain `@`
- automatically maintained variables include `LAST`, `LAST_<TASK>`, and active decontamination run pointers
- configuration variables include `GENOME_DIR`, `LIBRARIES_DIR`, `ORTHOFINDER_OUTPUT_DIR`, `EXPORTS_DIR`, `REPORTS_DIR`, `LOG_DIR`, daemon settings, selector defaults, and binary paths
- thread defaults are task-name variables such as `DEFAULT_THREADS_BUSCO_RUN`; `SET_MAX_THREADS_ON_START=true` refreshes `DAEMON_MAX_THREADS` and task thread defaults at daemon/run startup
- `list variables --kind busco-runs` shows stored BUSCO run-id variables

## 2.1 Storage roots and active write targets

PhyloODB now uses a storage-root registry for filesystem base paths.

Important rules:

- `genomes`, `libraries`, `orthofinder`, `exports`, and `logs` use a single-active-root model
- only the active root for one of those kinds is used for new writes
- inactive roots remain valid for existing bound data
- `reports`, `cache`, and `misc` are non-strict roots and are not part of the single-active switching model
- task and daemon logs default to `phyloodb.log` under the active `logs` root; daemon `--logfile` remains an explicit per-run override

This means a common SSD/HDD workflow is:

1. keep the SSD root active for new downloads and new work
2. add a second HDD root
3. move older genomes or libraries there
4. leave that HDD root inactive until you want new writes to land there

Key commands:

```bash
phyloODB project.db list roots
phyloODB project.db storage add-root --kind genomes --base-path /mnt/hdd/genomes --label HDD
phyloODB project.db storage rename-root HDD --label HDD_ARCHIVE
phyloODB project.db storage activate-root 7
phyloODB project.db storage deactivate-root 7
phyloODB project.db purge roots --inactive-only
phyloODB project.db discover --root HDD
```

Notes:

- `storage add-root` creates and verifies the resolved directory before registration, and creates a non-first strict root inactive by default
- `storage rename-root` preserves the root path, bindings, writable flag, and active state
- if no active root exists for `genomes`, `libraries`, `orthofinder`, or `exports`, tasks that need to create new data for that kind are blocked until a root is activated

## 3. Queue scheduling in brief

Queue dependencies are expressed with `--schedule`.

Accepted patterns include:

- `started:<selector>`
- `finished:<selector>`
- `succeeded:<selector>`
- `failed:<selector>`
- `delay:30s`
- `delay:5m`
- `delay:2h`
- `at:02:00`
- `at:2026-01-15T02:00:00`
- `queued-drained`

Task selectors commonly include:

- numeric task ids
- `LAST`
- task-specific forms such as `LAST_DOWNLOAD_BUSCO_LIBRARY`
- workflow step ids in workflow files

### Daemon control

Queued tasks are executed by the separate daemon process:

```bash
phyloODB-daemon project.db start --background
phyloODB-daemon project.db start --here --log-console --log-level INFO
phyloODB-daemon project.db stop
phyloODB-daemon project.db stop --drain
```

Useful daemon controls:

- `--max-threads` / `--threads`
- `--polling`
- `--blocked-polling`
- `--logfile`
- `--log-level`

Useful queue inspection commands:

```bash
phyloODB project.db list queue --watch
phyloODB project.db list queue --watch --all --refresh 1
phyloODB project.db list errors --watch
phyloODB project.db list errors --watch --stack --limit 50
phyloODB project.db watch queue
phyloODB project.db watch errors
```

`list queue` and `watch queue` support `-s`, `--sort latest|new|old|errors|running|status`.
Aliases are `changed`, `newest`, `oldest`, and `active`. The default `latest`
profile sorts top-level task blocks by the newest status change anywhere in the
task/subtask tree, while preserving hierarchy and sorting children recursively
inside each block. `errors` puts blocks with errored descendants first;
`running` puts blocks with active descendants first; `status` groups by the
root task status.

## 4. Task catalogue

### Bootstrap and data registration

#### `create-taxonomy`
Populate taxonomy tables from an NCBI taxdump.

Key options:

- `--path-to-taxdump`
- `--retain-taxdump`
- `--working-dir`

#### `update-assembly` (`add`)
Fetch assembly metadata for a taxid or explicit accession list.

Key options:

- `--taxid` or `--accessions`
- `--force-update`
- `--after`, `--before`
- `--level`
- `--exclude-accessions`
- `--exclude-taxids`
- `--exclude-clades`
- `--primary-only`

#### `download-assemblies` (`download`)
Download assemblies, optionally using selectors and built-in isoform cleaning.

NCBI downloads are written to staging files, checked against the assembly FTP `md5checksums.txt` manifest where available, fully validated as gzip streams, and only then promoted into place. The checksum manifest is registered as a genome artifact, and expected MD5 values are stored on downloaded file artifacts. Automatic prepare-proteome work is a required subtask; if that child task fails, the parent download/import task reports the failure too.

Key options:

- `--accessions` or `--taxid`
- `--protein`
- `--max-concurrent`
- `--force-redownload`
- `--download-retries`
- `--rank`, `--quantity`
- `--use-busco`
- `--min-completeness`
- `--min-single-copy-complete`
- cleaning flags such as `--skip-clean-isoforms`, `--clean-skip-gff`, and `--no-clean-skip-cdhit`
- project defaults such as `DEFAULT_PROTEOME_USE_GFF`, `DEFAULT_PROTEOME_USE_CDHIT`, and `DEFAULT_PROTEOME_CDHIT_IDENTITY`

#### `import-local-assembly`
Register a local assembly or proteome already present on disk.

At least one of `--fna` or `--faa` is required, together with either `--taxid`, `--taxon-name`, or `--genus` plus `--species`.

Examples:

```bash
phyloODB project.db queue import-local-assembly \
  --faa /data/sample/proteins.faa.gz \
  --gff /data/sample/annotations.gff3.gz \
  --accession SAMPLE_001 \
  --taxid 9606

phyloODB project.db queue import-local-assembly \
  --fna /data/sample/genome.fna.gz \
  --faa /data/sample/proteins.faa.gz \
  --genus Pan \
  --species troglodytes \
  --metadata '{"assembly_level":"chromosome"}'
```

By default local import also enables isoform cleaning with CD-HIT skipped; use `--no-clean-isoforms`, `--no-clean-skip-cdhit`, or related `clean-*` flags to change that behaviour. Local imports set `Assembly.origin=local` automatically, while NCBI metadata fetches populate `origin=refseq` or `origin=genbank`.

Local imports validate gzip inputs and write `phyloodb_md5checksums.txt` beside the imported assembly files. The manifest uses the same format as NCBI `md5checksums.txt`, is registered as a genome artifact, and supplies expected MD5 values for local FNA/FAA/GFF artifacts.

#### `batch-import-local-assembly`
Import many local assemblies from a directory.

Example:

```bash
phyloODB project.db queue batch-import-local-assembly \
  --assembly-dir /data/local_batch \
  --accessions-for-import SAMPLE_001 SAMPLE_002
```

This task also defaults to isoform cleaning with CD-HIT skipped.

#### `verify-downloads`
Check downloaded assemblies for file presence, gzip integrity, and stored NCBI or PhyloODB-local checksum matches where available, optionally reorganising or reacquiring data.

Important options:

- `--discover`
- `--discover-protein`
- `--reaquire`
- `--tidy`
- `--organise`
- `--organise-check-only`
- `--split-isolated-proteomes`
- `--report`

#### `verify-busco`

Supports `--root` for genome-root scoping. The root may be given as a numeric id or exact unique label.

Examples:

```bash
phyloODB project.db queue verify-assembly --root HDD_GENOMES --repair --tidy
phyloODB project.db queue verify-busco --root HDD_GENOMES --repair
phyloODB project.db queue verify-libraries --root HDD_LIBRARIES --repair
phyloODB project.db queue verify-orthofinder --root OF_HDD --repair
```
Check BUSCO outputs on disk against BUSCO records in the database.

Notes:

- verify-task reports default under the shared `reports` root in `verify-reports/task_<task_id>_<timestamp>...`
- `verify-orthofinder --repair` requires an active `orthofinder` storage root so result directories can be registered as rooted artifacts

Important options:

- `--library-id` or `--library-name`
- `--discover`
- `--reingest`
- `--queue-missing`
- `--report`

#### `split-records`
Split multiple FASTA files in a genome folder into new accessions.

#### `prepare-proteome`
Create an immutable derived proteome profile from a raw proteome.

Key options:

- `--accessions`
- `--profile-name`
- `--skip-gff`
- `--skip-cdhit` / `--no-skip-cdhit`
- `--cdhit-identity`
- `--set-default` / `--no-set-default`

### BUSCO and library construction

#### `download-busco-library`
Download and register a BUSCO lineage dataset.

Key options:

- `--lineage`
- `--libraries-dir`
- `--busco-path`
- `--parent-library-name`
- `--coverage`
- `--size`

#### `busco-run` (`busco`)
Run BUSCO for one accession.

Key options:

- `--accession`
- `--lineage`
- `--format {auto,protein,genome,nucleotide}`
- `--output-path`
- `--force`

Notes:

- here `--format` is the input mode for the new BUSCO run being created
- existing-run selector flags such as `--run-ids`, `--export-format`, `--require-format`, and `--prefer-format` do not apply to `busco-run`

#### `batch-busco`
Run BUSCO for multiple accessions.

Key options:

- `--accessions`
- `--lineage`
- `--format`
- `--output-dir`
- `--max-concurrent`
- `--busco-lib-wait-seconds`
- `--busco-lib-retries`

Notes:

- here `--format` is the input mode for the new BUSCO runs being created
- existing-run selector flags such as `--run-ids`, `--export-format`, `--require-format`, and `--prefer-format` do not apply to `batch-busco`

#### `orthofinder-run`
Run OrthoFinder across a panel of proteomes.

Key options:

- `--accessions` or `--input-dir`
- `--library-id` or `--library-name`
- `--proteome-profile`
- `--mcl-inflation`
- `--force`

Reuse semantics:

- OrthoFinder run reuse now depends on both the accession set and the effective MCL inflation value.
- if `--mcl-inflation` is omitted, PhyloODB treats that as the default OrthoFinder clustering setting and reuses only runs recorded with no explicit inflation override;
- if `--mcl-inflation` is supplied, PhyloODB requires the stored run to have the same value before reusing it.
- when an old run is discovered by scanning folders, PhyloODB reads `Log.txt` and parses the `Command Line:` entry to recover `-I` when present.

#### `add-library`
Build a derived study-specific library from a parent BUSCO lineage and reference accessions.

Key options:

- `--name`
- `--coverage` or `--coverage-taxid`
- `--accessions`
- `--parent-library-name` or `--parent-library-id`
- `--clean-refs` / `--clean-refs-strict`
- `--set-cleaned-primary` / `--no-set-cleaned-primary`
- `--rerun-busco`
- `--rerun-orthofinder`
- `--orthofinder-mcl-inflation`
- `--rerun-gene-trees`
- `--skip-paralog-analysis`
- `--gene-tree-source {iqtree,fasttree}`
- `--fast-tree`
- `--force`

Rebuild semantics:

- without `--force`, `add-library` refuses to overwrite an existing derived library of the same name;
- `--force` purges the stored derived-library state for that library name and removes the derived library output directory before rebuilding;
- `--force` does not by itself rerun BUSCO or OrthoFinder;
- `--rerun-busco` forces fresh BUSCO child tasks for the reference accessions;
- `--rerun-orthofinder` forces a fresh OrthoFinder child task for the reference accession set;
- `--orthofinder-mcl-inflation` is forwarded to the child OrthoFinder run and becomes part of the reuse identity for matching existing OrthoFinder results;
- `--rerun-gene-trees` forces fresh replacement IQ-TREE orthogroup trees even when matching trees already exist in `IQ-TREE_Orthogroup_trees`;
- `--skip-paralog-analysis` accepts exact 1:1 BUSCO/orthogroup families directly after occupancy filtering and skips tree building plus paralog classification;
- default `--gene-tree-source iqtree` builds canonical core-set trees in `IQ-TREE_Orthogroup_trees` while leaving OrthoFinder `Resolved_Gene_Trees` untouched;
- `--gene-tree-source fasttree` or `--fast-tree` reuses OrthoFinder `Resolved_Gene_Trees` directly and skips both replacement MAFFT and IQ-TREE subtasks;
- if `--force` is used without `--rerun-orthofinder`, any existing OrthoFinder run already linked to that derived library is preserved and can be reused;
- if `--force --rerun-orthofinder` is used, that existing library-linked OrthoFinder run is dropped and a new one is expected.

Library outputs now include `library_build_metadata.json`, which records the chosen core-set strategy and effective tree source.

#### `import-custom-library`
Register a library from a user-supplied list of BUSCO family ids.

Key options:

- `--library-name`
- `--coverage` or `--coverage-taxid`
- `--parent-library-name` or `--parent-library-id`
- `--busco-ids`
- `--ref-accessions`

### BLAST and paralog screening

#### `create-proteome-blast-db`
Create a BLAST database from one accession’s protein set.

#### `construct-busco-blast-db`
Create a BLAST database from BUSCO sequences for a set of accessions.

#### `paralog-removal` (`paralog-filtering`, `filter-paralogs`)
Screen BUSCO-derived sequences against reference proteomes to identify likely hidden paralogs.

Key options:

- `--ref-accessions`
- `--accessions`
- `--targets`
- `--library-id` or `--library-name`
- `--mode {median,percent,bitscore,lower-quartile,upper-quartile}`
- `--percentile`
- `--bitscore-threshold`
- `--report-dir`
- `--run-label`
- `--max-concurrent`
- `--reuse-existing`
- `--avoid-unclean-buscos`

`--accessions` defines the BUSCO pool used to compute family thresholds for `median`, `lower-quartile`, and `upper-quartile`. `--targets` optionally restricts which accessions receive new paralog-filtering results; if omitted, targets default to the resolved `--accessions` set.

BUSCO selector flags on this command choose which BUSCO run supplies the family rows for both targets and references. In practice:

- `--busco-pipeline augustus` or `--require-busco-pipeline augustus` requires Augustus runs
- the task still obeys `--busco-run-selection`
- so `--busco-pipeline augustus` by itself means "use the Augustus primary run if one exists", not "use any Augustus run"

If the intent is "use Augustus-derived BUSCO families wherever Augustus results exist", add `--busco-run-selection latest`.

Mode summary:

- `median`: default; compare against reference proteomes whose BUSCO bitscore is at or above the family median.
- `lower-quartile`: more permissive than median; include references at or above the family 25th percentile.
- `upper-quartile`: stricter than median; include references at or above the family 75th percentile.
- `percent`: rank reference proteomes by BUSCO bitscore per family and compare against the top `N%` via `--percentile`.
- `bitscore`: compare against references whose BUSCO bitscore is at or above `--bitscore-threshold`.

Reports default to the shared reports root under `paralog-filtering-reports/task_<task_id>_<timestamp>...`. If `--report-dir` is supplied, all paralog-filtering report TSVs are written directly into that directory.

### Decontamination

#### `decontamination` (`decontam`)
Reference-based BUSCO decontamination.

Key options:

- `--targets` or `--accessions`
- `--refs` or `--ref-accessions`
- `--library-id` or `--library-name`
- reference selectors: `--ref-clade`, `--ref-rule-rank`, `--ref-rule-quantity`
- thresholds: `--rank`, `--off-clade-fraction`, `--min-buscos`, `--min-identity`, `--min-coverage`, `--min-delta-bitscore`, `--min-hits`, `--hit-window`
- `--config-path`
- `--run-label`
- `--report-path`

JSON config shape:

```json
{
  "params": {
    "rank": "order",
    "off_clade_fraction": 0.1,
    "min_buscos": 20,
    "ref_clade": "Mammalia",
    "ref_rule_rank": "order",
    "ref_rule_quantity": 1
  },
  "targets": ["GCF_000001405.40"],
  "references": ["Primates", "Rodentia"],
  "groups": [
    {
      "members": ["Homo sapiens"],
      "clades": ["Primates"],
      "blacklist": [],
      "min_hits": 1
    }
  ]
}
```

#### `internal-decontamination` (`internal-decontam`, `idc`)
Internal BUSCO-consistency screen with optional external BLAST confirmation.

Key options:

- `--targets` or `--accessions`
- `--library-id` or `--library-name`
- `--rank`
- `--hit-window`
- `--p-value-threshold`
- `--off-clade-fraction`
- `--save-blast-output`
- `--reuse-blast-results`
- external stage options such as `--external-blast-db-path`, `--external-blast-output-dir`

JSON config shape:

```json
{
  "params": {
    "rank": "order",
    "hit_window": 8,
    "p_value_threshold": 0.05,
    "off_clade_fraction": 0.05,
    "external_blast_db_path": "/db/nr/nr"
  },
  "targets": ["GCF_000001405.40", "GCF_037993035.2"],
  "acceptance_rules": [
    {
      "members": ["Primates"],
      "clades": ["Primates"],
      "blacklist": [],
      "min_hits": 1
    }
  ]
}
```

#### `external-decontamination-check`
Run external BLAST checks for BUSCOs flagged by an internal decontamination run.

#### `external-decontamination-apply`
Apply external BLAST results to a previous internal decontamination run under a new run id.

### Export and reporting

#### `export` (`export-library`)
Export a filtered dataset for a library and target panel.

Key options:

- target definition: `--accessions`, `--accession`, `--taxid`, `--clade`, `--rank`, `--quantity`
- library selection: `--library-id` or `--library-name`
- `--out-dir`
- `--disable-paralog-filter`
- `--disable-decont-filter`
- `--require-paralog-filtering`
- `--require-decontamination`
- BUSCO run selection: `--busco-pipeline` / `--require-busco-pipeline`, `--prefer-busco-pipeline`, `--format` / `--require-format`, `--prefer-format`, `--busco-run-selection`
- report outputs: `--write-lineage-csv`, `--write-busco-report`, `--write-busco-family-matrix`
- `--busco-report-extended`
- `--retain-headers`
- occupancy thresholds: `--min-occupancy`, `--min-taxa-occupancy`

Notes:

- on BUSCO-aware export selection, `--format` means BUSCO input format (`protein` or `genome`)
- export output sequence type is still controlled separately by `--sequence-type`
- if `--out-dir` is omitted, export writes to a default task-named directory under the active `exports` root
- export reports are written alongside the exported data in that export directory rather than under the generic reports root
- default filtering behaviour is opportunistic: if paralog-filtering or decontamination results exist, export uses them; if they do not exist, export still runs
- `--require-paralog-filtering` and `--require-decontamination` make missing filtering state an error
- `--disable-paralog-filter` and `--disable-decont-filter` ignore those filtering results even when they are available

`--require` semantics:

- clauses are ANDed together
- `|` provides OR within a clause
- clade names and taxids are accepted
- exact accessions are accepted via `acc:ACCESSION` or `accession:ACCESSION`

Examples:

```bash
phyloODB project.db queue export --library-name metazoa_core --accessions @MY_PANEL --require Primates Rodentia
phyloODB project.db queue export --library-name metazoa_core --accessions @MY_PANEL --require '(Primates|Glires)' Carnivora
phyloODB project.db queue export --library-name metazoa_core --accessions @MY_PANEL --require '(acc:GCA_000001405.29|Primates)' Carnivora
```

`--header` syntax:

- active only when `--retain-headers` is not used
- valid tokens: `ACCESSION`, `TAXON`, `KINGDOM`, `PHYLUM`, `CLASS`, `ORDER`, `FAMILY`, `GENUS`, `SPECIES`, `RANK`, `BUSCO`, `LENGTH`, `GENE`, `TAXID`, `BITSCORE`
- allowed separators: `.`, `|`, `_`, `-`, `:`, `[`, `]`
- use `--header-rank` when the template includes `RANK`
- `FAMILY` means taxonomic family; `BUSCO` means the BUSCO family id

Export-side auxiliary outputs can include:

- `busco_families/` FASTA directory
- lineage CSV
- BUSCO report
- taxa occupancy TSV
- BUSCO family matrix
- filter report
- export parameters summary
- export task log copy

Example:

```bash
phyloODB project.db queue export \
  --library-name metazoa_core \
  --accessions @MY_PANEL \
  --header 'BUSCO:TAXON:RANK' \
  --header-rank phylum
```

#### `build-busco-trees`
Build MAFFT alignments and IQ-TREE gene trees for an exported BUSCO family set.

Key options:

- all `export` target/library/run-selection options
- `--out-dir`
- `--mafft-flags`
- `--iqtree-flags`
- `--mafft-threads`
- `--iqtree-threads`

Outputs:

- `busco_families/`: export-stage family FASTAs used as tree inputs
- `alignments/`: one MAFFT alignment per family
- `trees/`: one IQ-TREE result directory per family
- `manifest.tsv`: per-family table with `family_id`, `raw_fasta`, `alignment_path`, `tree_dir`, and final `tree_path`

This command wraps `export` first, then builds trees from the exported families that survived export-side filtering and occupancy thresholds.

#### `generate-lineage-csv`
Export lineage information for a selector-defined set of accessions.

## 5. Recommended habits

- Use `phyloODB <db> info <task>` for a task summary.
- Use `phyloODB <db> queue <task> --help` for exact flags.
- Use `list assemblies` or `assemblies` before committing to large downloads.
- Store important selector results with `--store` so later commands remain reproducible.
