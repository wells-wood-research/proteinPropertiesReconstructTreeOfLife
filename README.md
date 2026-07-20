# Protein Properties & the Reconstruction of the Tree of Life

Can the biophysical properties of a proteome — computed entirely from AlphaFold2 structural models — carry enough phylogenetic signal to reconstruct an approximation of the Tree of Life? This project tests that hypothesis across 41 organisms spanning six kingdoms.

**Conda environment:** `treeoflife` (Python 3.10.15, R 4.3.1)

---

## Contents

1. [Scientific overview](#scientific-overview)
2. [Environment setup](#environment-setup)
3. [Repository structure](#repository-structure)
4. [Data pipeline](#data-pipeline)
5. [Source module reference](#source-module-reference)
6. [Analyses](#analyses)
7. [Key configuration variables](#key-configuration-variables)
8. [Organism coverage](#organism-coverage)
9. [Antibody subproject](#antibody-subproject)
10. [Development notes](#development-notes)

---

## Scientific overview

Phylogenetic inference has traditionally relied on DNA or protein *sequence* similarity. This project explores an orthogonal source of signal: the **aggregate biophysical properties** of an organism's AlphaFold2-predicted proteome, as computed by the [DE-STRESS](https://pragmaticproteindesign.bio.ed.ac.uk/de-stress/) structural analysis tool.

DE-STRESS produces ~60 metrics per structure — packing density, secondary structure composition, isoelectric point, Rosetta and EvoEF2 energy terms, aggregation propensity (Aggrescan3D), and more. These are averaged per organism across tens of thousands of AF2 structural models, then projected with PCA and fed into hierarchical clustering to produce a dendrogram. The resulting Newick tree is compared quantitatively against the NCBI reference phylogeny using the **Clustering Information Distance (CID)** metric from the R `TreeDist` package.

Secondary analyses explore subcellular location, mitochondrially-encoded proteins, individual FoldSeek structural clusters, and single protein families (e.g. NADH-ubiquinone oxidoreductase chain 4) as alternative phylogenetic lenses.

There is also a self-contained sub-project (`antibodyproduction/`) that applies the same DE-STRESS feature pipeline to a separate problem: predicting the yeast-display expression level of computationally designed antibody scFv fragments.

---

## Environment setup

```bash
conda env create -f environment.yml
conda activate treeoflife
```

The key dependencies are:

| Package | Role |
|---|---|
| `pandas`, `numpy` | Data manipulation |
| `scikit-learn`, `scipy` | PCA, clustering, scaling |
| `matplotlib`, `seaborn`, `plotly` | Visualisation |
| `ete3`, `dendropy` | Newick tree I/O |
| `requests`, `beautifulsoup4` | AFDB / UniProt API queries |
| R `TreeDist`, `ape` | Phylogenetic tree distance (CID) |

> **Hardcoded paths:** Most scripts contain absolute paths set during development. Before running any script, search for path strings near the top and update them for your environment.

---

## Repository structure

Source code is version-controlled. `data/` and `analysis/` are git-ignored (large, reproducible from the pipeline).

```
proteinPropertiesReconstructTreeOfLife/
│
├── src/
│   ├── data_download/          # Download scripts (AFDB, UniProt, rsync)
│   ├── data_prep/              # Feature engineering pipeline
│   ├── dim_red/                # PCA and visualisation
│   ├── clustering/             # Hierarchical + k-means + tree distance (R)
│   ├── pymol_calign_script.pml
│   ├── pymol_calign_trim_structure_script.py
│   └── pymol_save_all_objects.py
│
├── antibodyproduction/         # Independent scFv expression-prediction subproject
│   └── src/
│       ├── data_prep_tools.py / data_prep.py
│       ├── dim_red_tools.py / dim_red.py
│       ├── feature_selection_tools.py / feature_selection.py
│       └── model_building_tools.py / model_building.py
│
├── data/                       # git-ignored — raw and processed data
│   ├── raw_data/
│   └── processed_data/
│       ├── af2/{standard,robust,minmax}/
│       └── pdb/{standard,robust,minmax}/
│
├── analysis/                   # git-ignored — all analysis outputs (PNGs, CSVs, .nwk)
│
├── environment.yml
└── .gitignore
```

---

## Data pipeline

The pipeline runs in four sequential stages.

### Stage 1 — Data acquisition

Download AF2 proteome tar files from the AFDB FTP, extract pLDDT scores from PDB B-factor columns, query the AFDB REST API for organism metadata, query the UniProt REST API for subcellular location and gene encoding type, and filter the 500M-row FoldSeek cluster table to only proteins present in the DE-STRESS dataset.

| Script | Output |
|---|---|
| `download_af2_data.py` | PDB tar files on HPC (`NUM_WORKERS=20`, `wget`) |
| `extract_plddt_score.py` | `raw_data/af2_plddt_scores.csv` |
| `download_af2_org_sci_name.py` | `raw_data/af2db_org_uniprot_desc_results.csv` |
| `download_uniprot_data.py` | `raw_data/uniprot_results_org_subcellloc_uniprotkb.csv` |
| `filtering_af2db_clustering_data.py` | `raw_data/filtered_af2db_clusters_data.csv` |

### Stage 2 — Data preparation & feature engineering

Join DE-STRESS output with organism metadata and FoldSeek cluster assignments. Remove features with >5% missing values, drop constant features, normalize energy terms by sequence length, remove highly correlated features (Spearman |r| > 0.6), and scale using three methods in parallel: Standard (z-score), Robust (IQR), and MinMax.

Produces one **non-redundant** dataset (one structure per organism × FoldSeek cluster, random seed 42) and one full dataset.

| Script | Role |
|---|---|
| `uniprot_af2db_data_prep.py` | Joins UniProt + AFDB metadata → `processed_af2db_uniprot_data.csv` |
| `af2db_clusters_data_prep.py` | Cleans filtered FoldSeek cluster table |
| `random_select_structures_by_org_and_cluster.py` | Generates `af2_structures_non_redundant.csv` |
| `data_prep.py` | Main execution script — calls `data_prep_tools.py` |
| `data_prep_tools.py` | Shared utility library (see module reference below) |

Scaled outputs are written to `data/processed_data/af2/{standard,robust,minmax}/`:
- `processed_destress_data_scaled.csv` — full dataset
- `processed_destress_data_scaled_nonredundant.csv` — non-redundant dataset
- `{method}_scaler.pkl` — fitted scaler for applying to new data
- `labels.csv` / `labels_nonredundant.csv` — metadata columns saved separately

### Stage 3 — Dimensionality reduction & clustering

Run PCA across all dataset × scaling method combinations. Average scaled DE-STRESS metrics per organism, then run hierarchical clustering across all combinations of linkage method and distance metric. Export each dendrogram as a Newick `.nwk` file.

| Script | Role |
|---|---|
| `pca_analysis.py` | Full-proteome PCA (all models) |
| `pca_avg_destress_metrics_by_org.py` | Organism-averaged PCA |
| `pca_all_mitochondrial_proteins.py` | Mito-encoded proteins only |
| `pca_subcellular_location.py` | Membrane / Nucleus / Cytoplasm |
| `single_af2db_cluster_pca_analysis.py` | One FoldSeek cluster at a time |
| `single_protein_pca_analysis.py` | One protein family at a time |
| `hierarchical_clustering_avg_destress_metrics.py` | Tree reconstruction → `.nwk` files |
| `kmeans_avg_destress_metrics.py` | K-means evaluation (k = 2–20) |

### Stage 4 — Tree distance comparison (R)

Read all reconstructed Newick trees and compute the Clustering Information Distance (CID) against the NCBI reference phylogeny. Compare against a pre-computed null distribution from `randomTreeDistances.rda`.

```bash
Rscript src/clustering/tree_dist.R
```

Output: `tree_distances.csv` in the relevant `analysis/` subdirectory.

Two reference trees are available:
- `data/processed_data/ncbi_phylo_tree.phy` — full ~40-organism tree
- `data/processed_data/ncbi_phylo_tree_euk.phy` — eukaryotes only

---

## Source module reference

### `src/data_prep/data_prep_tools.py` — Core utility library

Consumed by all data prep scripts.

| Function | Description |
|---|---|
| `remove_missing_val_features(data, output_path, threshold)` | Drops columns where proportion of missing values exceeds `threshold`. Saves diagnostic CSV. |
| `remove_constant_features(data, constant_features_threshold, output_path)` | Drops columns where the modal value (rounded to 2 d.p.) accounts for more than `threshold` fraction of rows. |
| `remove_highest_correlators(data, corr_coeff_threshold, output_path)` | Greedy iterative removal of the column with the most pairwise Spearman correlations above `threshold`. Saves before/after correlation matrices. |
| `adding_af2db_uniprot_columns(...)` | Joins DE-STRESS data with AFDB organism metadata and FoldSeek cluster assignments. Creates `organism_group` (7 kingdoms) and `organism_group2` (Eukaryotes/Prokaryotes/Other) label columns. |
| `adding_destress_summary_cols(destress_data)` | Adds binned categorical columns: `dssp_bin` (dominant secondary structure), `isoelectric_point_bin`, `packing_density_bin`, `aggrescan3d_avg_bin`. |
| `scale_destress_data_remove_high_corr(...)` | Scales with chosen method, saves fitted scaler as `.pkl`, then removes highly correlated features. |
| `process_af2_data(...)` | Full end-to-end pipeline: missing value removal → dropna → summary cols → join labels → normalize energies by `num_residues` → deduplicate by organism/cluster → filter → scale. |
| `process_pdb_data(...)` | Same pipeline without the organism/cluster joining steps. Used for PDB reference structures. |

### `src/dim_red/dim_red_tools.py` — PCA utilities

| Function | Description |
|---|---|
| `pca_var_explained(data, n_components, ...)` | Fits PCA, saves cumulative variance explained as CSV and line plot, saves fitted model as `.pkl`. |
| `perform_pca(data, labels_df, n_components, ...)` | Fits and transforms; saves per-component top-10 feature loadings and the full transformed DataFrame as CSV. |
| `plot_latent_space_2d(...)` | Seaborn scatterplot of two PC dimensions. Configurable hue, style, alpha, palette, and legend position. |
| `plot_pca_plotly(...)` | Interactive Plotly 2D scatter saved as HTML with hover data. |
| `spectral_plot(...)` | Mean PC value across PC1–PC7 per organism group with 99% CI bands. Compares kingdom-level PC signatures. |
| `distance_to_reference(data, dim_columns, ...)` | All-vs-all pairwise distances in PCA space. |
| `plot_pca_boxplots(...)` | Boxplots of PC values grouped by a categorical variable (e.g. secondary structure type). |

### `src/clustering/clustering_tools.py` — Clustering utilities

| Function | Description |
|---|---|
| `get_newick(node, parent_dist, leaf_names, newick)` | Recursive conversion of a scipy `to_tree()` output to Newick format with branch lengths. |
| `plot_dendrogram(model, **kwargs)` | Converts a sklearn `AgglomerativeClustering` model to a scipy linkage matrix and plots the dendrogram. |
| `adj_rand_ind_wssd_plot(...)` | Dual-axis plot: inertia (left) and adjusted Rand index (right) vs. number of k-means clusters. |

### `src/clustering/tree_dist.R`

Reads all reconstructed `.nwk` trees and the NCBI reference phylogeny. Computes CID for every tree using the `TreeDist` R package. Loads a pre-computed null distribution from `randomTreeDistances.rda` for statistical significance assessment.

### `src/data_download/` — Download scripts

| Script | What it downloads |
|---|---|
| `download_af2_data.py` | AFDB proteome tar files (parses AFDB index HTML with BeautifulSoup, 20 parallel `wget` processes) |
| `download_af2_org_sci_name.py` | Organism name + UniProt description via AFDB REST API (4 worker processes) |
| `download_uniprot_data.py` | Subcellular location, GO codes, lineage class, and gene encoding type via UniProt REST API (6 worker processes) |
| `download_structures_af2_pdb_files.py` | Individual AF2 PDB files from AFDB for a specific set of UniProt IDs |
| `download_af2db_files_from_server.py` | rsync of specific PDB files from HPC to local disk |
| `extract_plddt_score.py` | Parses PDB B-factor columns to extract per-residue pLDDT; computes per-model mean from Cα atoms |
| `filtering_af2db_clustering_data.py` | Filters the 500M-row FoldSeek cluster TSV in chunks of 500,000 |

### PyMOL scripts

- `pymol_calign_script.pml` — Loads PDB files for a cluster and aligns all to the first loaded structure using `cmd.align`. Currently configured for cluster `A0A0G9LKG2`.
- `pymol_calign_trim_structure_script.py` — Runs inside PyMOL. Aligns all structures, finds the intersection of well-aligned residue pairs (distance ≤ 2.0 Å), and creates trimmed selections.
- `pymol_save_all_objects.py` — Runs inside PyMOL. Saves all loaded objects as PDB files to `data/raw_data/superoxide_dismutase_trunc/`.

> PyMOL must be available separately — it is not in `environment.yml`.

---

## Analyses

| Script | Analysis | Input subset | Key outputs |
|---|---|---|---|
| `pca_analysis.py` | Full-proteome PCA | All AF2 models, all 3 scalers | 2D scatter plots by kingdom, spectral plots, boxplots by secondary structure and isoelectric point |
| `pca_avg_destress_metrics_by_org.py` | Organism-averaged PCA | Non-redundant; averaged per organism | PCA of 41 organisms coloured by Domain, shaped by Kingdom |
| `pca_all_mitochondrial_proteins.py` | Mitochondrial proteome PCA | Mitochondrion gene-encoding-type only | PC1×PC2, PC1×PC3, PC2×PC3 scatter plots |
| `pca_subcellular_location.py` | Subcellular location PCA | Non-redundant; Membrane/Nucleus/Cytoplasm only | PCA of organism × location combinations |
| `single_protein_pca_analysis.py` | Single protein family PCA | One UniProt description | PCA of orthologs across organisms, FASTA export |
| `single_af2db_cluster_pca_analysis.py` | Single FoldSeek cluster PCA | One FoldSeek cluster representative | PCA coloured by kingdom |
| `hierarchical_clustering_avg_destress_metrics.py` | Hierarchical tree reconstruction | Organism-averaged features, all scalers | Dendrograms (PNG) and Newick trees (`.nwk`) for all linkage × distance combinations |
| `kmeans_avg_destress_metrics.py` | K-means evaluation | Organism-averaged features | Adjusted Rand index and inertia vs k (k = 2–20) |
| `tree_dist.R` | Phylogenetic accuracy | Reconstructed `.nwk` + NCBI reference tree | `tree_distances.csv` with CID scores |

---

## Key configuration variables

All parameters are set near the top of each execution script. There is no central config file.

### `src/data_prep/data_prep.py`

| Variable | Default | Description |
|---|---|---|
| `scaling_method_list` | `["standard", "robust", "minmax"]` | All three are always run in parallel. |
| `af2_corr_coeff_threshold` | `0.6` | Spearman \|r\| threshold for removing correlated features. |
| `af2_constant_features_threshold` | `0.25` | Drops a feature if its modal value accounts for >25% of rows. |
| `missing_val_threshold_af2` | `0.05` | Drops a feature if >5% of values are missing. |
| `remove_low_quality_af2_models` | `True` | Filters out AF2 models with mean pLDDT < 70. |
| `remove_redundant_af2_models` | `True` | One structure per (organism, FoldSeek cluster), seed 42. |
| `af2_energy_field_list` | 25 Rosetta/EvoEF2 + 4 BuDEFF fields | Normalized by `num_residues` before scaling to remove length bias. |

### `src/dim_red/pca_analysis.py` and related

| Variable | Default | Description |
|---|---|---|
| `n_components` | `7` | Number of principal components retained. |

### `src/clustering/hierarchical_clustering_avg_destress_metrics.py`

| Variable | Default | Description |
|---|---|---|
| Linkage methods | `["single", "average", "complete", "ward"]` | All four are run and compared. Ward linkage requires Euclidean distance. |
| Distance metrics | `["euclidean", "cityblock", "cosine", "correlation"]` | All combinations produce a separate Newick tree and dendrogram. |

### `src/data_download/download_af2db_files_from_server.py`

| Variable | Value | Description |
|---|---|---|
| `server` | `"glutamate.bio.ed.ac.uk"` | HPC server for rsync of AF2 PDB files. Update to your own credentials. |
| `remote_folder` | `"/mnt/scratch/alphafold_model_organisms/"` | Remote path to proteome PDB files. |

---

## Organism coverage

41 organisms across 6 kingdoms. Kingdom assignments are hardcoded in `data_prep_tools.py` and used to create the `organism_group` label column throughout the pipeline.

**Animals (13):** *Homo sapiens*, *Mus musculus*, *Rattus norvegicus*, *Danio rerio*, *Caenorhabditis elegans*, *Drosophila melanogaster*, *Brugia malayi*, *Dracunculus medinensis*, *Onchocerca volvulus*, *Schistosoma mansoni*, *Strongyloides stercoralis*, *Trichuris trichiura*, *Wuchereria bancrofti*

**Bacteria (16):** *Escherichia coli*, *Mycobacterium tuberculosis*, *Mycobacterium leprae*, *Mycobacterium ulcerans*, *Staphylococcus aureus*, *Streptococcus pneumoniae*, *Pseudomonas aeruginosa*, *Klebsiella pneumoniae*, *Helicobacter pylori*, *Campylobacter jejuni*, *Enterococcus faecium*, *Salmonella typhimurium*, *Shigella dysenteriae*, *Haemophilus influenzae*, *Neisseria gonorrhoeae*, *Nocardia brasiliensis*

**Fungi (9):** *Saccharomyces cerevisiae*, *Schizosaccharomyces pombe*, *Candida albicans*, *Ajellomyces capsulatus*, *Paracoccidioides lutzii*, *Cladophialophora carrionii*, *Fonsecaea pedrosoi*, *Madurella mycetomatis*, *Sporothrix schenckii*

**Plants (4):** *Arabidopsis thaliana*, *Glycine max*, *Oryza sativa*, *Zea mays*

**Protozoa (5):** *Plasmodium falciparum*, *Trypanosoma brucei*, *Trypanosoma cruzi*, *Leishmania infantum*, *Dictyostelium discoideum*

**Archaea (1):** *Methanocaldococcus jannaschii*

---

## Antibody subproject

`antibodyproduction/` is an independent pipeline that mirrors the main project's structure (data prep → dim red → feature selection → model building). It applies DE-STRESS features to a classification problem: predicting whether computationally designed scFv antibody fragments (Fleishman lab designs targeting insulin and *M. tuberculosis* ACP) will express in a yeast-display system.

**Task:** Three-class classification — Low / Medium / High expression level.  
**Model:** Gaussian Naive Bayes with 10-repeat 5-fold cross-validation.  
**Feature selection:** Mutual information (`feature_select_mi`) and Random Forest importance (`feature_select_rf`, 1000 trees, balanced class weights).

The pipeline loops over all combinations of scaling method (`standard`, `robust`, `minmax`) × amino acid composition included/excluded × feature selection method (`mi`/`rf`).

### Key parameters (`antibodyproduction/src/data_prep.py`)

| Variable | Value | Description |
|---|---|---|
| `test_size` | `0.25` | Stratified train/test split by expression bin. |
| `random_state` | `42` | |
| `corr_coeff_threshold` | `0.7` | Slightly more permissive than main project (0.6) due to smaller dataset. |
| `constant_features_threshold` | `0.8` | More permissive constant-feature cutoff. |
| `num_cvs` / `num_folds` | `10` / `5` | 50 total validation folds; mean metrics and CI reported. |

### Outputs

| Directory | Contents |
|---|---|
| `antibodyproduction/data/processed_data/{scaler}/{comp_flag}/` | `X_train_scaled.csv`, `X_test_scaled.csv`, `y_train.csv`, `y_test.csv`, `pdb_scaled.csv` |
| `antibodyproduction/feature_selection/{scaler}/{comp_flag}/` | `selected_features_mi.csv`, `selected_features_rf.csv` |
| `antibodyproduction/analysis/dim_red/pca/` | PCA scatter plots and spectral plots coloured by expression level |
| `antibodyproduction/models/` | Confusion matrices, ROC curves, `model_val_master.csv`, `model_test_master.csv` |

---

## Development notes

**Three scaling methods are always run in parallel.** Every analysis loops over `["standard", "robust", "minmax"]`. Results are compared across scaling methods to assess robustness — there is no single canonical output.

**Non-redundant vs. full dataset.** The pipeline produces both a full set (all structures) and a non-redundant set (one structure per organism × FoldSeek cluster). Most clustering and organism-averaged PCA analyses use the non-redundant set.

**Single-cluster analysis is partially templated.** `single_af2db_cluster_pca_analysis.py` has infrastructure to loop over many FoldSeek clusters but is currently overridden to run only against cluster `A0A0G9LKG2`. The PyMOL alignment scripts are also configured for this cluster.

**BeautifulSoup import bug.** `download_af2_data.py` contains the incorrect import `from beautifulsoup4 import BeautifulSoup`. The correct form is `from bs4 import BeautifulSoup`. The package installs as `beautifulsoup4` but imports as `bs4`.

**PyMOL scripts run inside a PyMOL session.** `pymol_calign_trim_structure_script.py` and `pymol_save_all_objects.py` are not standalone scripts — they must be executed from within PyMOL.
