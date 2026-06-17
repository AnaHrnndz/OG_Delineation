# 🧬 OGD: Orthologous Group Delineation

**OGD** is a Python bioinformatics pipeline for **Orthologous Group Delineation** from gene phylogenetic trees. It identifies orthologous groups (OGs) based on the detection and scoring of gene duplication events, using ETE4 for tree manipulation and supporting NCBI and GTDB taxonomy.

---

## 📋 Table of Contents

1. [Main Workflow](#main-workflow)
2. [Requirements and Dependencies](#requirements-and-dependencies)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Command-Line Options](#command-line-options)
   - [Required Arguments](#required-arguments)
   - [Tree Options](#tree-options)
   - [Taxonomy Options](#taxonomy-options)
   - [Species Overlap Thresholds](#species-overlap-thresholds)
   - [Algorithm Parameters](#algorithm-parameters)
   - [Visualization](#visualization)
   - [Optional Modules](#optional-modules)
6. [Output Files](#output-files)
7. [Usage Examples](#usage-examples)

---

## 💡 Main Workflow

The core script (`og_delineation.py`) follows these steps to delineate **Orthologous Groups (OGs)**:

1. **Data Loading** — Imports the gene tree, taxonomy database, reference species tree, and taxonomic counter.
2. **Pre-analysis Setup** — Performs tree adjustments:
   - Polytomy resolution.
   - Tree rooting (*Midpoint* or *MinVar* via FastRoot).
   - Taxonomic annotation (NCBI or GTDB).
3. **Outlier Detection and Score Calculation** — Computes key metrics:
   - Detection of long-branch outliers and taxonomically misplaced sequences.
   - Calculation of **Species Overlap (SO)**, *Duplication Score*, *Inparalogs Rate*, and *Lineage Loss*.
4. **HQ-Duplication Detection** — Selects high-quality duplication events that define OG roots.
5. **OG Generation** — Delineates three OG types: monophyletic (mOGs), paraphyletic (pOGs), and root-level OGs.
6. **Ortholog Pairs** — Generates a table of ortholog pairs.
7. **Output** — Writes the annotated tree, OG info table, sequence-to-OG mapping, and ortholog pairs.

### Optional Modules

#### 🔄 Recovery Pipeline

Attempts to recover sequences initially excluded from OGs (long-branch or taxonomically misplaced) using an HMM-based approach:

1. Generates FASTA files for sequences excluded from OGs.
2. Builds an HMM profile for each OG using HMMER.
3. Runs **HMMscan** to assign excluded sequences to their best-matching OG.
4. Updates the gene tree and OG information with recovered sequences.

#### 📝 eggNOG-mapper Annotation

Enriches sequences with functional annotations (DIAMOND + HMM) using **eggNOG-mapper**, adding GO terms, KEGG pathways, and PFAM domains.

---

## ⚙️ Requirements and Dependencies

### Python

- Python 3.10
- [ETE4](https://github.com/etetoolkit/ete) — tree manipulation and taxonomy

### External Binaries

| Tool | Usage | Notes |
|------|-------|-------|
| **FAMSA** | Multiple sequence alignment | Binary expected at `og_delineation/bin/famsa` |
| **HMMER** (`hmmbuild`, `hmmscan`) | HMM-based sequence recovery | Must be in `$PATH` |
| **FastRoot.py** | MinVar tree rooting | Required only with `--rooting MinVar` |
| **eggNOG-mapper** | Functional annotation | Required only with `--run_emapper` |

### Taxonomy Databases

- **NCBI taxonomy** (default): ETE4 SQLite database (e.g., `e6.taxa.sqlite`)
- **GTDB taxonomy**: via ETE4's `GTDBTaxa`

---

## 🛠️ Installation

```bash
# Clone the repository
git clone https://github.com/AnaHrnndz/OG_Delineation.git
cd OG_Delineation

# Create and activate conda environment
conda env create -f og_delineation/ogd_env2.yml
conda activate ogd_env2
```

---

## 🚀 Quick Start

```bash
# Minimal run with default settings
./og_delineation.py --tree tree.nw --output_path ./ogd_results

# Run and open visualization in browser when done
./og_delineation.py --tree tree.nw --output_path ./ogd_results --open_visualization

# All available options
./og_delineation.py -h
```

---

## 🔧 Command-Line Options

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--tree PATH` | Input gene tree file (Newick format) |
| `--output_path PATH` | Output directory |

---

### Tree Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--rooting` | `Midpoint` | Rooting method: `Midpoint` or `MinVar` (FastRoot) |
| `--sp_delimitator STR` | `.` | Delimiter used to extract species ID from leaf names |

**Leaf name format:** `SpeciesID<delimiter>GeneName`

```bash
# Default: "9606.Gene001" → species ID = "9606"
./og_delineation.py --tree tree.nw --output_path ./output

# Custom delimiter: "9606_Gene001" → species ID = "9606"
./og_delineation.py --tree tree.nw --output_path ./output --sp_delimitator "_"
```

---

### Taxonomy Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--taxonomy_type` | `NCBI` | Taxonomy backend: `NCBI` or `GTDB` |
| `--user_taxonomy PATH` | None | Custom taxonomy SQLite database |
| `--user_taxonomy_counter PATH` | None | JSON file with precomputed species counts per taxonomic node |
| `--reftree PATH` | None | Custom reference species tree |

```bash
# NCBI taxonomy (default)
./og_delineation.py --tree tree.nw --output_path ./output --taxonomy_type NCBI

# GTDB taxonomy (for prokaryotic genomes)
./og_delineation.py --tree tree.nw --output_path ./output --taxonomy_type GTDB

# Custom taxonomy
./og_delineation.py --tree tree.nw --output_path ./output \
  --user_taxonomy /path/to/custom.db \
  --user_taxonomy_counter /path/to/counter.json \
  --reftree /path/to/species_ref.nw
```

---

### Species Overlap Thresholds

**Species Overlap (SO)** measures the fraction of shared species between the two subtrees at a duplication node. A lower SO indicates a cleaner duplication event (less co-speciation). OGD accepts global and per-domain SO thresholds.

| Argument | Default | Description |
|----------|---------|-------------|
| `--sp_ovlap_all FLOAT` | `0.1` | Global SO threshold (all nodes) |
| `--sp_ovlap_euk FLOAT` | None | SO threshold for Eukaryota-only nodes |
| `--sp_ovlap_bact FLOAT` | None | SO threshold for Bacteria-only nodes |
| `--sp_ovlap_arq FLOAT` | None | SO threshold for Archaea-only nodes |

When a domain-specific threshold is set, it overrides `--sp_ovlap_all` for nodes belonging to that domain.

```bash
# Default: 10% overlap allowed globally
./og_delineation.py --tree tree.nw --output_path ./output

# More stringent globally: fewer but cleaner OGs
./og_delineation.py --tree tree.nw --output_path ./output \
  --sp_ovlap_all 0.05

# More permissive globally: more OGs, accepts noisier duplications
./og_delineation.py --tree tree.nw --output_path ./output \
  --sp_ovlap_all 0.20

# Domain-specific thresholds
# (Bacteria more permissive due to HGT; Eukaryota more strict)
./og_delineation.py --tree tree.nw --output_path ./output \
  --sp_ovlap_all 0.1 \
  --sp_ovlap_euk 0.05 \
  --sp_ovlap_bact 0.25 \
  --sp_ovlap_arq 0.10
```

---

### Algorithm Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--lineage_threshold FLOAT` | `0.05` | Min fraction of sequences to consider a taxon represented at a node (below → outlier) |
| `--best_taxa_threshold FLOAT` | `0.9` | Confidence threshold for LCA taxonomic annotation |
| `--species_losses_perct FLOAT` | `0.7` | Max fraction of species loss tolerated at a duplication node |
| `--no_inherit_outliers` | False | Prevent outlier status from propagating to parent nodes |
| `--skip_get_pairs` | False | Skip ortholog pairs table generation |

---

### Visualization

OGD uses the **ETE4 SmartView** web interface for interactive tree visualization. There are two visualization modes:

#### Mode A — Run pipeline then visualize (`--open_visualization`)

Runs the full pipeline and opens SmartView automatically when it finishes.

```bash
# Local access
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --open_visualization

# Remote access (specify server IP)
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --open_visualization --user_IP 192.168.1.100
```

#### Mode B — Visualize only (`--only_visualization`)

Opens SmartView directly on an already-annotated tree, **without re-running the pipeline**. Useful for re-inspecting results from a previous run.

> **Input requirement:** The tree file should be an annotated Newick file (e.g., `tree.tree_annot.nw`) produced by a previous OGD run.

```bash
# Visualize a previously generated annotated tree
./og_delineation.py \
  --tree ./ogd_results/tree.tree_annot.nw \
  --output_path ./dummy \
  --only_visualization

# Include alignment for sequence display in SmartView
./og_delineation.py \
  --tree ./ogd_results/tree.tree_annot.nw \
  --output_path ./dummy \
  --only_visualization \
  --raw_alg seqs.aln

# Remote access
./og_delineation.py \
  --tree ./ogd_results/tree.tree_annot.nw \
  --output_path ./dummy \
  --only_visualization \
  --user_IP 192.168.1.100
```

| Flag | Runs pipeline | Generates output files | Opens SmartView |
|------|:---:|:---:|:---:|
| *(none)* | ✅ | ✅ | ❌ |
| `--open_visualization` | ✅ | ✅ | ✅ at the end |
| `--only_visualization` | ❌ | ❌ | ✅ immediately |

---

### Optional Modules

#### Sequence Recovery (HMMER-based)

Requires a sequence alignment file (`--raw_alg`).

| Argument | Values | Description |
|----------|--------|-------------|
| `--raw_alg PATH` | — | Input alignment file (FASTA) |
| `--run_recovery` | `run-align` / `skip-align` | Enable recovery; `run-align` re-aligns, `skip-align` uses `--raw_alg` directly |

```bash
# Recovery with re-alignment
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --raw_alg seqs.faa --run_recovery run-align

# Recovery skipping alignment step (use provided alignment directly)
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --raw_alg seqs.aln --run_recovery skip-align
```

#### eggNOG-mapper Annotation

Requires a sequence alignment file and eggNOG databases.

| Argument | Description |
|----------|-------------|
| `--run_emapper` | Run eggNOG-mapper (DIAMOND + HMM) |
| `--emapper_dmnd_db PATH` | DIAMOND database for eggNOG-mapper |
| `--emapper_pfam_db PATH` | PFAM HMM database |
| `--path2emapper_main PATH` | Use pre-computed eggNOG-mapper results (main table) |
| `--path2emapper_pfams PATH` | Use pre-computed eggNOG-mapper results (PFAM table) |

```bash
# Run eggNOG-mapper annotation
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --raw_alg seqs.aln \
  --run_emapper \
  --emapper_dmnd_db /path/to/eggnog.dmnd \
  --emapper_pfam_db /path/to/pfam.hmm

# Use pre-computed emapper results (skip re-running emapper)
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --path2emapper_main /path/to/emapper.annotations \
  --path2emapper_pfams /path/to/emapper.pfam_hits
```

---

## 📂 Output Files

| File | Description | Generated when |
|------|-------------|----------------|
| `<tree>.tree_annot.nw` | Annotated gene tree in Newick format | Always |
| `<tree>.ogs_info.tsv` | OG table: taxonomy, scores, members, hierarchy | Always |
| `<tree>.seq2ogs.tsv` | Sequence → OG mapping (TSV) | Always |
| `<tree>.seq2ogs.jsonl` | Sequence → OG mapping (JSONL, one line per sequence) | Always |
| `<tree>_pairs_clean.tsv` | All ortholog pairs | Without `--skip_get_pairs` |
| `<tree>_pairs_strict.tsv` | High-confidence ortholog pairs | Without `--skip_get_pairs` |
| HMM profiles + recovery files | Recovered sequences and HMM assignments | With `--run_recovery` |
| eggNOG-mapper annotations | Functional annotations, PFAM domains | With `--run_emapper` |

### `ogs_info.tsv` columns

| Column | Description |
|--------|-------------|
| `#OG_name` | OG identifier |
| `Lca_Dup` | Taxonomic ID of the LCA at the duplication node |
| `TaxoLevel` | Taxonomic level of the OG |
| `SciName_TaxoLevel` | Scientific name of the taxonomic level |
| `NumSP` | Number of species in the OG |
| `NumSeqs` | Number of sequences in the OG |
| `NumRecoverySeqs` | Number of sequences added by recovery |
| `OG_up` / `OG_down` | Parent / child OGs (hierarchy) |
| `Inparalogs_Rate` | In-paralogs rate at the duplication node |
| `SP_overlap_dup` | Species overlap score at the duplication node |
| `Seqs` | List of sequence IDs |

---

## 📚 Usage Examples

### Minimal run

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results
```

### With MinVar rooting

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results --rooting MinVar
```

### Run and visualize immediately

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results --open_visualization
```

### Visualize a previous result without re-running

```bash
./og_delineation.py \
  --tree ./ogd_results/tree.tree_annot.nw \
  --output_path ./dummy \
  --only_visualization
```

### Custom Species Overlap for mixed-domain tree

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --sp_ovlap_all 0.1 \
  --sp_ovlap_euk 0.05 \
  --sp_ovlap_bact 0.25 \
  --sp_ovlap_arq 0.10
```

### With sequence recovery

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --raw_alg seqs.aln --run_recovery run-align
```

### Full run with recovery and functional annotation

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --rooting MinVar \
  --raw_alg seqs.aln \
  --run_recovery run-align \
  --run_emapper \
  --emapper_dmnd_db /path/to/eggnog.dmnd \
  --emapper_pfam_db /path/to/pfam.hmm \
  --open_visualization
```

### GTDB taxonomy (prokaryotic genomes)

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --taxonomy_type GTDB
```

### Remote visualization

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --open_visualization --user_IP 192.168.1.100
```
