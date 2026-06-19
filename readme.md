## 🧬 OGD: Orthologous Group Delineation

**OGD** is a Python bioinformatics pipeline for **Orthologous Group Delineation** from gene phylogenetic trees. It identifies orthologous groups (OGs) based on the detection and scoring of gene duplication events, using ETE4 for tree manipulation. NCBI taxonomy is the supported backend; GTDB support is experimental (see [Current Limitations](#current-limitations)).

---

## 📋 Table of Contents

1. [Main Workflow](#main-workflow)
2. [Requirements and Dependencies](#requirements-and-dependencies)
3. [Installation](#installation)
4. [Taxonomy Database (required)](#taxonomy-database-required)
5. [Quick Start](#quick-start)
6. [Command-Line Options](#command-line-options)
   - [Required Arguments](#required-arguments)
   - [Tree Options](#tree-options)
   - [Taxonomy Options](#taxonomy-options)
   - [Species Overlap Thresholds](#species-overlap-thresholds)
   - [Algorithm Parameters](#algorithm-parameters)
   - [Visualization](#visualization)
   - [Optional Modules](#optional-modules)
7. [Output Files](#output-files)
8. [Usage Examples](#usage-examples)
9. [Current Limitations](#current-limitations)

---

## 💡 Main Workflow

The core script (`og_delineation.py`) follows these steps to delineate **Orthologous Groups (OGs)**:

1. **Data Loading** — Imports the gene tree, taxonomy database, reference species tree, and taxonomic counter.
2. **Pre-analysis Setup** — Performs tree adjustments:
   - Polytomy resolution.
   - Tree rooting (*Midpoint* or *MinVar* via FastRoot).
   - Taxonomic annotation (NCBI).
3. **Outlier Detection and Score Calculation** — Computes key metrics:
   - Detection of long-branch outliers and taxonomically misplaced sequences.
   - Calculation of **Species Overlap (SO)**, *Duplication Score*, *Inparalogs Rate*, and *Lineage Loss*.
4. **HQ-Duplication Detection** — Selects high-quality duplication events that define OG roots.
5. **OG Generation** — Delineates three OG types: monophyletic (mOGs), paraphyletic (pOGs), and root-level OGs.
6. **Ortholog Pairs** — Generates a table of ortholog pairs.
7. **Output** — Writes the annotated tree, OG info table, sequence-to-OG mapping, and ortholog pairs.

### Optional Module

#### 📝 eggNOG-mapper Annotation

Enriches sequences with functional annotations (DIAMOND + HMM) using **eggNOG-mapper**, adding GO terms, KEGG pathways, and PFAM domains. See [Optional Modules](#optional-modules) below.

---

## ⚙️ Requirements and Dependencies

### Python

- Python 3.10
- [ETE4](https://github.com/etetoolkit/ete) — tree manipulation and taxonomy

### External Binaries

| Tool | Usage | Notes |
|------|-------|-------|
| **FastRoot.py** | MinVar tree rooting | Required only with `--rooting MinVar` |
| **eggNOG-mapper** | Functional annotation (**Linux only**) | Required only with `--run_emapper` (see [Optional Modules](#optional-modules)) |

### Taxonomy Database

- **NCBI taxonomy** via ETE4's `NCBITaxa` (required — see [Taxonomy Database](#taxonomy-database-required)).

---

## 🛠️ Installation

### Option 1 — pip (recommended)

​```bash
git clone https://github.com/AnaHrnndz/OG_Delineation.git
cd OG_Delineation

pip install .              # core OG delineation
pip install ".[emapper]"   # optional: eggNOG-mapper annotation (Linux only)
​```

This pulls the dependencies (`ete4`, `FastRoot`, `numpy`) and installs the **`og-delineation`** command. To install without cloning (e.g. as a dependency of another tool):

​```bash
pip install "git+https://github.com/AnaHrnndz/OG_Delineation.git"
​```

### Option 2 — conda environment

​```bash
git clone https://github.com/AnaHrnndz/OG_Delineation.git
cd OG_Delineation

conda env create -f ogd_env.yml
conda activate ogd_env
pip install .              # install OGD into the environment
​```

The conda environment pins the dependency versions; `pip install .` then adds the `og-delineation` command.

---

## 🗃️ Taxonomy Database (required)

OGD reads the species from each leaf name and resolves it against a taxonomy database through **ETE4's `NCBITaxa`**. You must set this database up **before** running OGD.

ETE4 handles the download, build, and storage of the NCBI taxonomy database, and lets you choose how and where to store it. Follow ETE4's own documentation for this step:

> **ETE4 — NCBITaxa documentation:** `https://etetoolkit.github.io/ete/tutorial/tutorial_taxonomy.html`

Once it is available, you can either rely on ETE4's default location or point OGD to a specific SQLite file with `--user_taxonomy`:

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --user_taxonomy /path/to/taxa.sqlite
```

> **Note — missing taxids stop the run.** Every species identifier in your tree must exist in the taxonomy database. If one does not, OGD aborts with a clear message:
>
> ```
> Failed to build reference species tree. Taxid not found in taxonomy DB: <taxid>.
> Check that all species identifiers in your tree exist in the taxonomy database.
> ```

---

## 🚀 Quick Start

Leaf names must encode the species as `SpeciesID<delimiter>GeneName`, with the delimiter defaulting to `.` (for example `9606.ENSP00000269305`). Change it with `--sp_delimitator` if your names use a different separator.

```bash
# Minimal run on the example tree shipped in data/
./og_delineation.py --tree data/P53.fa.nw --output_path ./ogd_results

# Run and open the interactive ETE SmartView visualization when done
./og_delineation.py --tree data/P53.fa.nw --output_path ./ogd_results --open_visualization

# All available options
./og_delineation.py -h
```

> **Running OGD.** If you installed with pip, use the `og-delineation` command. The examples below use `./og_delineation.py` (running from a clone) — the two are equivalent, so just swap `./og_delineation.py` for `og-delineation`. You can also use `python -m ogd`.

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
| `--taxonomy_type` | `NCBI` | Taxonomy backend: `NCBI` or `GTDB` (GTDB is experimental) |
| `--user_taxonomy PATH` | None | Custom taxonomy SQLite database |
| `--user_taxonomy_counter PATH` | None | JSON file with precomputed species counts per taxonomic node |
| `--reftree PATH` | None | Custom reference species tree |

```bash
# NCBI taxonomy (default)
./og_delineation.py --tree tree.nw --output_path ./output --taxonomy_type NCBI

# Custom taxonomy database / counter / reference tree
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
```

| Flag | Runs pipeline | Generates output files | Opens SmartView |
|------|:---:|:---:|:---:|
| *(none)* | ✅ | ✅ | ❌ |
| `--open_visualization` | ✅ | ✅ | ✅ at the end |
| `--only_visualization` | ❌ | ❌ | ✅ immediately |

---

### Optional Modules

#### eggNOG-mapper Annotation

OGD can enrich the sequences in the tree with functional annotations using **eggNOG-mapper** (emapper). This step is **optional** and requires a sequence alignment (`--raw_alg`).

> **Note — Linux only.** The eggNOG-mapper module currently runs on **Linux only**: the pip-installed eggNOG-mapper bundles a Linux DIAMOND binary that fails on macOS (`cannot execute binary file`). Everything else (the core OG delineation, rooting, visualization) works on both Linux and macOS.

eggNOG-mapper and its databases (DIAMOND + Pfam) are installed and managed separately from OGD. Install emapper and download its databases following the eggNOG-mapper documentation:

> **eggNOG-mapper:** https://github.com/eggnogdb/eggnog-mapper

| Argument | Description |
|----------|-------------|
| `--raw_alg PATH` | Input alignment file (FASTA), required for emapper |
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
| `OG_up` / `OG_down` | Parent / child OGs (hierarchy) |
| `Inparalogs_Rate` | In-paralogs rate at the duplication node |
| `SP_overlap_dup` | Species overlap score at the duplication node |
| `Seqs` | List of sequence IDs |

---

## 📚 Usage Examples

### Minimal run

```bash
./og_delineation.py --tree data/P53.fa.nw --output_path ./ogd_results
```

### With MinVar rooting

```bash
./og_delineation.py --tree data/P53.fa.nw --output_path ./ogd_results --rooting MinVar
```

### Run and visualize immediately

```bash
./og_delineation.py --tree data/P53.fa.nw --output_path ./ogd_results --open_visualization
```

### Visualize a previous result without re-running

```bash
./og_delineation.py \
  --tree ./ogd_results/<tree>.tree_annot.nw \
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

### Run with functional annotation (eggNOG-mapper)

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --rooting MinVar \
  --raw_alg seqs.aln \
  --run_emapper \
  --emapper_dmnd_db /path/to/eggnog.dmnd \
  --emapper_pfam_db /path/to/pfam.hmm \
  --open_visualization
```

### Remote visualization

```bash
./og_delineation.py --tree tree.nw --output_path ./ogd_results \
  --open_visualization --user_IP 192.168.1.100
```

---

## Current Limitations

- The **eggNOG-mapper annotation** module (`--run_emapper`) runs on **Linux only** (the pip-installed eggNOG-mapper bundles a Linux-only DIAMOND binary). The core OG delineation works on Linux and macOS.
- **GTDB taxonomy** (`--taxonomy_type GTDB`) is experimental and not fully supported yet.
- The **sequence recovery** module is not ready for general use yet and is therefore not documented here.