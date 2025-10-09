# 🧬 OG_Delineation: Orthologous Group Delineation (OGD)


The **OG_Delineation** algorithm is designed for **Orthologous Group Delineation (OGD)** from gene phylogenetic trees. It identifies orthologous groups based on the detection and scoring of gene duplication events.

---

## 💡 Main Workflow (The OGD Algorithm)

The core script (`og_delineation.py`) follows these key steps to figure out the **Basal Orthologous Groups (OGs)**:

1.  **Data Loading:** Imports essential files (Gene Tree, Taxonomy Database, Reference Species Tree, etc.).
2.  **Pre-analysis Setup:** Performs fundamental tree adjustments:
    * Polytomy resolution.
    * Tree rooting (e.g., *Midpoint* or *MinVar*).
    * NCBI annotation.
3.  ** Detect outlierts and Score Calculation:** Determines key metrics:
    * Detection of long branches and **taxonomical outliers**.
    * Calculation of: **Species Overlap (SO)**, *Duplication Score*, *Inparalogs Rate*, and *Lineage Lost*.
4.  **HQ-Duplication Detection:** Selects high-quality duplication events that define the root of the **Orthologous Groups (OGs)**.
5.  **OG Generation:** Structures and annotates the final orthologous groups.
6.  **Get orthologs pairs**
6.  **Annotation & Output:** Annotates the tree root, flags sequences outside OGs, and writes output files.

---

## 🧩 Optional Pipelines

The program includes modules to refine and functionally annotate the results:

### 🔄 1. Recovery Pipeline (Sequence Recovery)

This module attempts to recover sequences that were initially excluded from the OGs (termed "left out" sequences) using an HMM-based process:

1.  Generates FASTA files for the excluded sequences and for each OG.
2.  Creates an HMM (*Hidden Markov Model*) profile for each OG.
3.  Uses **HMMscan** to assign excluded sequences to their best-matching OG.
4.  Updates the gene tree and OG information with the recovered sequences.

### 📝 2. eggnog-mapper Annotation

This allows you to enrich the sequences within the tree with functional annotations using the **eggnog-mapper** (*emapper*) tool.

---

## ⚙️ Requirements and Dependencies

This project requires Python 3.x and several specialized bioinformatics and data manipulation libraries.

* **Python 3.x**
* **ete3** (or **ete4**)
* **Biopython**
* Other standard dependencies (e.g., `NumPy`, `Pandas`).

> **Note:** The use of **Apptainer/Singularity** is highly recommended for a reproducible environment, as suggested by the `apptainer_ogd` folder.

---

## 🚀 Usage (Quick Start)

1.  **Clone the repository:**

    ```bash
    git clone [https://github.com/AnaHrnndz/OG_Delineation.git](https://github.com/AnaHrnndz/OG_Delineation.git)
    cd OG_Delineation
    ```

2.  **Install dependencies:**

    *(Add your specific installation steps here, or refer to your Apptainer setup.)*

3.  **Execution:**

    ```bash
    python og_delineation.py <tree_file> <taxonomy_file> [OPTIONS]
    ```

    *To view all available command-line options, run:*
    ```bash
    python og_delineation.py -h
    ```

---
