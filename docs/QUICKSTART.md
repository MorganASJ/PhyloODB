# PhyloODB Quick Start

Use this quick start for a first run from an empty database to a basic BUSCO-family export. It is intentionally brief. For the full explanation of any step, follow the links to the [User Manual](./MANUAL.md). For exact flag lookup, use the [Command Reference](./COMMAND_REFERENCE.md). For a richer worked example, use the [Tutorial](../tutorial/README.md).

This guide uses `run` so each task is executed immediately in the foreground. For larger studies, use `queue` and the daemon model described in [queueing, running, and the daemon model](./MANUAL.md#7-queueing-running-and-the-daemon-model).

## 1. Create a project database

```bash
phyloODB quickstart.db create --email your_email@domain.com --api-key your_NCBI_api_key
```

The database stores assembly metadata, downloaded file locations, BUSCO runs, variables, task history, and exported results. See [creating a database](./MANUAL.md#4-creating-a-database).

## 2. Add assembly metadata for a clade

```bash
phyloODB quickstart.db run add --clade Primates
```

This records which assemblies are available for the clade. It does not download sequence files yet. See [taxonomy and assembly metadata](./MANUAL.md#811-taxonomy-and-assembly-metadata).

## 3. Inspect and sample assemblies

```bash
phyloODB quickstart.db list assemblies \
  --clade Primates \
  --rank genus \
  --quantity 1 \
  --tidy

phyloODB quickstart.db count assemblies --clade Primates
```

The first command shows a simple one-assembly-per-genus sample. The second reports the size of the available clade. Selectors are reusable rules for choosing assemblies; see [selectors and accession resolution](./MANUAL.md#5-selectors-and-accession-resolution) and [`list assemblies`](./MANUAL.md#55-list-assemblies-metadata-busco-columns-filtering-and-output-modes).

## 4. Store a reusable accession panel

```bash
phyloODB quickstart.db assemblies \
  --clade Primates \
  --rank genus \
  --quantity 1 \
  --store PRIMATE_PANEL

phyloODB quickstart.db list variables --kind assemblies
```

Stored accession panels are referenced with `@NAME`, so this panel can be reused later as `@PRIMATE_PANEL`. See [stored accession sets as variables](./MANUAL.md#52-stored-accession-sets-as-variables). If you want to save the selector recipe rather than the resolved accession list, use [selector presets](./MANUAL.md#53-selector-presets).

## 5. Download the selected assemblies

```bash
phyloODB quickstart.db run download \
  --accessions @PRIMATE_PANEL \
  --protein
```

This downloads genome data and requires protein FASTA files, which are useful for a first protein-mode BUSCO run. New downloads can also trigger automatic proteome preparation using the database defaults. See [downloading assemblies](./MANUAL.md#812-downloading-assemblies) and [automatic proteome preparation after download](./MANUAL.md#813-automatic-proteome-preparation-after-download).

## 6. Download a BUSCO lineage

```bash
phyloODB quickstart.db run download-busco-library \
  --lineage mammalia_odb12
```

BUSCO lineages provide the marker families used for completeness assessment and this simple export. See [BUSCO lineage libraries](./MANUAL.md#821-busco-lineage-libraries).

## 7. Run BUSCO on the panel

```bash
phyloODB quickstart.db run batch-busco \
  --accessions @PRIMATE_PANEL \
  --lineage mammalia_odb12 \
  --format protein

phyloODB quickstart.db list assemblies \
  --accessions @PRIMATE_PANEL \
  --library-name mammalia_odb12 \
  --busco \
  --tidy
```

`batch-busco` creates one BUSCO run per accession. The `--format protein` flag means the new BUSCO runs use protein input. See [running BUSCO on assemblies](./MANUAL.md#822-running-busco-on-assemblies) and [BUSCO input modes and proteome profiles](./MANUAL.md#823-busco-input-modes-and-proteome-profiles).

## 8. Export a first-pass dataset

```bash
phyloODB quickstart.db run export \
  --library-name mammalia_odb12 \
  --accessions @PRIMATE_PANEL \
  --out-dir quickstart_export \
  --disable-paralog-filter \
  --disable-decont-filter
```

This writes a simple BUSCO-family export using the parent BUSCO lineage directly. It deliberately disables paralog and decontamination filters so the first run is easy to complete. Real study-quality exports should normally use an explicit filtering strategy. See [export and reporting](./MANUAL.md#85-export-and-reporting), [hidden paralog filtering](./MANUAL.md#83-hidden-paralog-filtering), and [decontamination](./MANUAL.md#84-decontamination).

## 9. Next steps

After this basic pass, the normal production path is to build a study-specific derived library, run filtering, and export again.

- Build a derived library with `add-library`: [derived library construction](./MANUAL.md#825-derived-library-construction).
- Understand BUSCO run selection before exporting: [choosing the BUSCO run for export](./MANUAL.md#choosing-the-busco-run-for-export).
- Learn how files are placed and moved: [storage roots and recovery](./MANUAL.md#12-database-management-and-recovery).
- Use `queue` and the daemon for larger task sets: [the daemon and queue inspection](./MANUAL.md#73-the-daemon-and-queue-inspection).

The quick start proves that the database, metadata, downloads, BUSCO, and export path are working. The [User Manual](./MANUAL.md) explains how to turn that into a reproducible analysis.
