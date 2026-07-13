#!/usr/bin/env bash
set -euo pipefail

# Example tutorial workflow for a primate-rich mammal project in PhyloODB.
# This script is illustrative. Read and adapt it before running on a real study.
#
# This example uses `run` so each stage completes before the next one starts.
# For larger production datasets, adapt the task sequence to your preferred
# batch or scheduler workflow.

DB="mammal_tutorial.db"
RESULTS_DIR="tutorial/results"
EXPORT_DIR="tutorial/exports"

mkdir -p "$RESULTS_DIR" "$EXPORT_DIR"

# 1. Create the project database.
phyloODB "$DB" create

# 2. Populate metadata for the central clade and a small set of mammalian outgroups.
phyloODB "$DB" run add --clade Primates
phyloODB "$DB" run add --clade Rodentia
phyloODB "$DB" run add --clade Carnivora
phyloODB "$DB" run add --clade Artiodactyla
phyloODB "$DB" run add --clade Marsupialia
phyloODB "$DB" run add --clade Monotremata

# 3. Register BUSCO libraries. Mammalia is the main working lineage here.
phyloODB "$DB" run download-busco-library --lineage metazoa_odb12 --coverage 1
phyloODB "$DB" run download-busco-library --lineage mammalia_odb12 --coverage 1

# 4. Inspect the primate clade before committing to a panel.
phyloODB "$DB" list assemblies --clade Primates --rank family --quantity 1
phyloODB "$DB" list assemblies --clade Primates --ranks family,genus --quantities 2,1
phyloODB "$DB" list assemblies --clade Primates --level chromosome

# 5. Store reference and target panels.
phyloODB "$DB" list assemblies --clade Primates --rank genus --quantity 1 --store PRIMATE_REFS
phyloODB "$DB" list assemblies --clade Rodentia --rank order --quantity 1 --store RODENT_OUTGROUPS
phyloODB "$DB" list assemblies --clade Carnivora --rank family --quantity 1 --store CARNIVORE_OUTGROUPS
phyloODB "$DB" list assemblies --clade Artiodactyla --rank family --quantity 1 --store UNGULATE_OUTGROUPS
phyloODB "$DB" list assemblies --clade Marsupialia --quantity 1 --store MARSUPIAL_OUTGROUPS
phyloODB "$DB" list assemblies --clade Monotremata --quantity 1 --store MONOTREME_OUTGROUPS
phyloODB "$DB" list assemblies --clade Primates --ranks family,genus --quantities 2,1 --store PRIMATE_TARGETS

# Demonstrate variable union with --append-to while building the final target panel.
phyloODB "$DB" list assemblies --accessions @PRIMATE_TARGETS --store MAMMAL_TARGETS
phyloODB "$DB" list assemblies --accessions @RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS --append-to MAMMAL_TARGETS

# 6. Check the stored variables.
phyloODB "$DB" list variables
phyloODB "$DB" list assemblies --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS
phyloODB "$DB" list assemblies --accessions @PRIMATE_TARGETS
phyloODB "$DB" list assemblies --accessions @MAMMAL_TARGETS

# 7. Download the references first, then the denser primate target panel.
phyloODB "$DB" run download --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS --protein
phyloODB "$DB" run download --accessions @MAMMAL_TARGETS --protein

# 8. Verify downloaded content and optionally normalise proteomes.
phyloODB "$DB" run verify-downloads --accessions @PRIMATE_REFS,@MAMMAL_TARGETS --downloaded-only
phyloODB "$DB" run clean-isoforms --accessions @PRIMATE_REFS,@MAMMAL_TARGETS --downloaded-only

# 9. Run BUSCO on the reference panel and the primate-rich target panel.
phyloODB "$DB" run batch-busco \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --lineage mammalia_odb12 \
  --format protein

phyloODB "$DB" run batch-busco \
  --accessions @MAMMAL_TARGETS \
  --lineage mammalia_odb12 \
  --format protein

# 10. Review BUSCO-aware rankings after BUSCO is complete.
phyloODB "$DB" list assemblies \
  --accessions @MAMMAL_TARGETS \
  --library-name mammalia_odb12 \
  --busco \
  --busco-complete-min 90 \
  --busco-single-min 80

# 11. Build a mammal-wide derived core library from the curated references.
# Default behavior builds replacement MAFFT alignments plus IQ-TREE orthogroup
# trees and writes the canonical core-set trees under:
#   <orthofinder results>/IQ-TREE_Orthogroup_trees
# To reuse OrthoFinder's own resolved trees instead, add:
#   --fast-tree
phyloODB "$DB" run add-library \
  --name mammal_core_odb12 \
  --coverage Mammalia \
  --parent-library-name mammalia_odb12 \
  --accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS

# 12. Remove hidden paralogs from the target panel using the same references.
phyloODB "$DB" run paralog-removal \
  --library-name mammal_core_odb12 \
  --ref-accessions @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --accessions @MAMMAL_TARGETS \
  --out-dir "$RESULTS_DIR/paralog_removal"

# 13. Internal decontamination is often the best first-pass contamination screen.
phyloODB "$DB" run internal-decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --rank order \
  --hit-window 8 \
  --off-clade-fraction 0.05 \
  --report-path "$RESULTS_DIR/internal_decontamination/mammal_internal"

# 14. Reference-based decontamination is available when an explicit mammalian
# reference panel is preferred.
phyloODB "$DB" run decontamination \
  --library-name mammal_core_odb12 \
  --targets @MAMMAL_TARGETS \
  --refs @PRIMATE_REFS,@RODENT_OUTGROUPS,@CARNIVORE_OUTGROUPS,@UNGULATE_OUTGROUPS,@MARSUPIAL_OUTGROUPS,@MONOTREME_OUTGROUPS \
  --rank order \
  --off-clade-fraction 0.10 \
  --min-buscos 20 \
  --report-path "$RESULTS_DIR/decontamination/mammal_reference"

# 15. Optional external confirmation after an internal run.
# Replace <internal_run_id> with the run_id reported by the internal decontamination task.
# phyloODB "$DB" run external-decontamination-check \
#   --run-id <internal_run_id> \
#   --blast-db-path /path/to/external/db \
#   --output-dir "$RESULTS_DIR/external_decontamination"
#
# phyloODB "$DB" run external-decontamination-apply \
#   --source-run-id <internal_run_id> \
#   --run-label mammal_internal_confirmed \
#   --external-blast-output-dir "$RESULTS_DIR/external_decontamination"

# 16. Export the final filtered dataset.
phyloODB "$DB" run export \
  --library-name mammal_core_odb12 \
  --accessions @MAMMAL_TARGETS \
  --out-dir "$EXPORT_DIR/mammal_core" \
  --write-lineage-csv \
  --write-busco-report \
  --write-busco-family-matrix \
  --busco-report-extended
