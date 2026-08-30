# Protein Properties & the Reconstruction of the Tree of Life

Code and analysis for the paper: *Large-scale analysis of AlphaFold structures reveals organism-specific physicochemical signatures*

**Conda environment:** `treeoflife` (Python 3.10.15, R 4.3.1)

---

## Environment setup

```bash
conda env create -f environment.yml
conda activate treeoflife
```

| Package | Role |
|---|---|
| `pandas`, `numpy` | Data manipulation |
| `scikit-learn`, `scipy` | PCA, clustering, scaling |
| `matplotlib`, `seaborn`, `plotly` | Visualisation |
| `ete3`, `dendropy` | Newick tree I/O |
| `requests`, `beautifulsoup4` | AFDB / UniProt API queries |
| R `TreeDist`, `ape` | Phylogenetic tree distance (CID) |

> **Note:** PyMOL must be installed separately - it is not in `environment.yml`. PyMOL scripts must be run from within a PyMOL session, not as standalone scripts.

---

## Repository structure

Source code is version-controlled. `data/` and `analysis/` are git-ignored (large, reproducible from the pipeline).

```
proteinPropertiesReconstructTreeOfLife/
│
├── src/
│   ├── data_download/
│   ├── data_prep/
│   ├── dim_red/
│   ├── clustering/
│   ├── pymol_calign_script.pml
│   ├── pymol_calign_hide_nterm_script.pml
│   ├── pymol_calign_trim_structure_script.py
│   └── pymol_save_all_objects.py
│
├── antibodyproduction/
│   └── src/
│
├── data/                       # git-ignored
│   ├── raw_data/
│   └── processed_data/
│       ├── af2/{standard,robust,minmax}/
│       └── pdb/{standard,robust,minmax}/
│
├── analysis/                   # git-ignored
├── environment.yml
└── .gitignore
```

---

## Data pipeline

The pipeline runs in four sequential stages.

### Stage 1 - Data acquisition (`src/data_download/`)

| Script | What it does |
|---|---|
| `download_af2_data.py` | Downloads AFDB proteome tar files (parses AFDB index HTML, 20 parallel `wget` processes). **Note:** contains incorrect import `from beautifulsoup4 import BeautifulSoup` - correct form is `from bs4 import BeautifulSoup`. |
| `extract_plddt_score.py` | Parses PDB B-factor columns to extract per-residue pLDDT; computes per-model mean from Ca atoms |
| `download_af2_org_sci_name.py` | Downloads organism name and UniProt description via AFDB REST API (4 worker processes) |
| `download_uniprot_data.py` | Downloads subcellular location, GO codes, lineage class, and gene encoding type via UniProt REST API (6 worker processes) |
| `download_mitochondrial_encoded_proteins_uniprot.py` | Downloads mitochondrially-encoded protein accessions from UniProt |
| `filtering_af2db_clustering_data.py` | Filters the 500M-row FoldSeek cluster TSV in chunks of 500,000 to proteins present in the DE-STRESS dataset |
| `download_structures_af2_pdb_files.py` | Downloads individual AF2 PDB files from AFDB for a specific set of UniProt IDs |
| `download_af2db_files_from_server.py` | rsync of specific PDB files from HPC to local disk. Update `server` and `remote_folder` at the top of the script for your environment. |

### Stage 2 - Data preparation & feature engineering (`src/data_prep/`)

Joins DE-STRESS output with organism metadata and FoldSeek cluster assignments. Removes features with >5% missing values, drops constant features, normalizes energy terms by sequence length, removes highly correlated features (Spearman |r| > 0.6), and scales using three methods in parallel: Standard (z-score), Robust (IQR), and MinMax.

Produces one **non-redundant** dataset (one structure per organism x FoldSeek cluster, random seed 42) and one full dataset.

Scaled outputs are written to `data/processed_data/af2/{standard,robust,minmax}/`:
- `processed_destress_data_scaled.csv` - full dataset
- `processed_destress_data_scaled_nonredundant.csv` - non-redundant dataset
- `{method}_scaler.pkl` - fitted scaler
- `labels.csv` / `labels_nonredundant.csv` - metadata columns

| Script | Role |
|---|---|
| `data_prep.py` | Main execution script - calls `data_prep_tools.py` |
| `data_prep_tools.py` | Shared utility library (see module reference below) |
| `uniprot_af2db_data_prep.py` | Joins UniProt + AFDB metadata |
| `af2db_clusters_data_prep.py` | Cleans filtered FoldSeek cluster table |
| `random_select_structures_by_org_and_cluster.py` | Generates the non-redundant structure list |
| `data_prep_single_proteins.py` | Data preparation for single protein family analyses |
| `csv_to_fasta_file.py` | Converts a CSV of sequences to FASTA format |
| `destress_features_by_dataset.py` | Summarises which DE-STRESS features are used in each analysis across the paper |
| `dssp_and_aa_proportions.py` | Computes DSSP bin, secondary structure residue, and amino acid proportions across datasets |
| `sequence_length_histogram.py` | Plots protein sequence length distributions across datasets |

### Stage 3 - Dimensionality reduction (`src/dim_red/`)

Runs PCA across dataset x scaling method combinations. Averaged scaled DE-STRESS metrics per organism feed into hierarchical clustering.

| Script | Role |
|---|---|
| `dim_red_tools.py` | Shared utility library (see module reference below) |
| `pca_analysis.py` | Full-proteome PCA (all models) |
| `pca_avg_destress_metrics_by_org.py` | Organism-averaged PCA |
| `pca_all_mitochondrial_proteins.py` | Mitochondrially-encoded proteins only |
| `pca_subcellular_location.py` | Membrane / Nucleus / Cytoplasm subsets |
| `pca_analysis_spectral_all_organisms.py` | Spectral PCA plots coloured by kingdom across all organisms |
| `subcellular_location_analysis.py` | Subcellular location breakdown analysis |
| `single_af2db_cluster_pca_analysis.py` | PCA for a single FoldSeek cluster (currently set to cluster `A0A0G9LKG2`) |
| `single_af2db_cluster_pca_analysis_filtered.py` | PCA for all 11 FoldSeek clusters with sufficient organism coverage (>=40 organisms) |
| `single_protein_pca_analysis.py` | PCA for a single protein family |
| `plddt_pca_correlation_sod.py` | Correlation between pLDDT and PC scores for the superoxide dismutase cluster |
| `plddt_pca_correlation_sod_by_kingdom.py` | Same analysis broken down by kingdom |
| `plddt_vs_feature_pc1_correlations_sod.py` | Correlation between pLDDT and individual DE-STRESS features for the superoxide dismutase cluster |

### Stage 4 - Clustering & tree distance (`src/clustering/`)

Averages scaled DE-STRESS metrics per organism, runs hierarchical clustering across all combinations of linkage method and distance metric, exports dendrograms as Newick `.nwk` files, and computes Clustering Information Distance (CID) against the NCBI reference phylogeny.

Two NCBI reference trees are available:
- `data/processed_data/ncbi_phylo_tree.phy` - full 48-organism tree
- `data/processed_data/ncbi_phylo_tree_euk.phy` - eukaryotes only

| Script | Role |
|---|---|
| `clustering_tools.py` | Shared utility library (see module reference below) |
| `hierarchical_clustering_avg_destress_metrics.py` | Hierarchical clustering across all organisms - produces dendrograms and `.nwk` trees for all linkage x distance combinations |
| `kmeans_avg_destress_metrics.py` | K-means evaluation (k = 2-20) |
| `tree_dist.R` | Computes CID between all reconstructed trees and the NCBI reference phylogeny |
| `hierarchical_clustering_avg_destress_metrics_single_af2db_cluster_filtered.py` | Hierarchical clustering for each of the 11 filtered FoldSeek clusters |
| `plot_cid_summary_single_af2db_cluster_filtered.py` | Plots CID summary across the 11 filtered clusters |
| `tree_dist_single_af2db_cluster_filtered.R` | CID computation for the filtered single-cluster trees |
| `test_tree_dist_single_af2db_cluster_filtered.R` | Validates tip label normalisation and CID computation for filtered cluster trees |
| `summarise_cluster_descriptions.py` | Summarises protein descriptions within a single FoldSeek cluster |
| `summarise_cluster_descriptions_all.py` | Summarises protein descriptions across all 22 FoldSeek clusters |

---

## Source module reference

### `src/data_prep/data_prep_tools.py`

| Function | Description |
|---|---|
| `remove_missing_val_features(data, output_path, threshold)` | Drops columns where proportion of missing values exceeds `threshold`. Saves diagnostic CSV. |
| `remove_constant_features(data, constant_features_threshold, output_path)` | Drops columns where the modal value (rounded to 2 d.p.) accounts for more than `threshold` fraction of rows. |
| `remove_highest_correlators(data, corr_coeff_threshold, output_path)` | Greedy iterative removal of the column with the most pairwise Spearman correlations above `threshold`. Saves before/after correlation matrices. |
| `adding_af2db_uniprot_columns(...)` | Joins DE-STRESS data with AFDB organism metadata and FoldSeek cluster assignments. Creates `organism_group` (6 kingdoms) and `organism_group2` (Eukaryotes/Prokaryotes/Other) label columns. |
| `adding_destress_summary_cols(destress_data)` | Adds binned categorical columns: `dssp_bin` (dominant secondary structure), `isoelectric_point_bin`, `packing_density_bin`, `aggrescan3d_avg_bin`. |
| `scale_destress_data_remove_high_corr(...)` | Scales with chosen method, saves fitted scaler as `.pkl`, then removes highly correlated features. |
| `process_af2_data(...)` | Full end-to-end pipeline: missing value removal -> dropna -> summary cols -> join labels -> normalize energies by `num_residues` -> deduplicate by organism/cluster -> filter -> scale. |
| `process_pdb_data(...)` | Same pipeline without the organism/cluster joining steps. Used for PDB reference structures. |

### `src/dim_red/dim_red_tools.py`

| Function | Description |
|---|---|
| `pca_var_explained(data, n_components, ...)` | Fits PCA, saves cumulative variance explained as CSV and line plot, saves fitted model as `.pkl`. |
| `perform_pca(data, labels_df, n_components, ...)` | Fits and transforms; saves per-component top-10 feature loadings and the full transformed DataFrame as CSV. |
| `plot_latent_space_2d(...)` | Seaborn scatterplot of two PC dimensions. Configurable hue, style, alpha, palette, and legend position. |
| `plot_pca_plotly(...)` | Interactive Plotly 2D scatter saved as HTML with hover data. |
| `spectral_plot(...)` | Mean PC value across PC1-PC7 per organism group with 99% CI bands. |
| `distance_to_reference(data, dim_columns, ...)` | All-vs-all pairwise distances in PCA space. |
| `plot_pca_boxplots(...)` | Boxplots of PC values grouped by a categorical variable. |

### `src/clustering/clustering_tools.py`

| Function | Description |
|---|---|
| `get_newick(node, parent_dist, leaf_names, newick)` | Recursive conversion of a scipy `to_tree()` output to Newick format with branch lengths. |
| `plot_dendrogram(model, **kwargs)` | Converts a sklearn `AgglomerativeClustering` model to a scipy linkage matrix and plots the dendrogram. |
| `adj_rand_ind_wssd_plot(...)` | Dual-axis plot: inertia (left) and adjusted Rand index (right) vs. number of k-means clusters. |

### PyMOL scripts

- `pymol_calign_script.pml` - Loads PDB files for a cluster and aligns all to the first loaded structure using `cmd.align`. Currently configured for cluster `A0A0G9LKG2`.
- `pymol_calign_hide_nterm_script.pml` - Same alignment, but hides low-pLDDT N-terminal tails before aligning.
- `pymol_calign_trim_structure_script.py` - Runs inside PyMOL. Aligns all structures, finds the intersection of well-aligned residue pairs (distance <= 2.0 A), and creates trimmed selections.
- `pymol_save_all_objects.py` - Runs inside PyMOL. Saves all loaded objects as PDB files.

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

---

## Organism coverage

48 organisms across 6 kingdoms. Kingdom assignments are hardcoded in `data_prep_tools.py`.

**Animals (13):** *Homo sapiens*, *Mus musculus*, *Rattus norvegicus*, *Danio rerio*, *Caenorhabditis elegans*, *Drosophila melanogaster*, *Brugia malayi*, *Dracunculus medinensis*, *Onchocerca volvulus*, *Schistosoma mansoni*, *Strongyloides stercoralis*, *Trichuris trichiura*, *Wuchereria bancrofti*

**Bacteria (16):** *Escherichia coli*, *Mycobacterium tuberculosis*, *Mycobacterium leprae*, *Mycobacterium ulcerans*, *Staphylococcus aureus*, *Streptococcus pneumoniae*, *Pseudomonas aeruginosa*, *Klebsiella pneumoniae*, *Helicobacter pylori*, *Campylobacter jejuni*, *Enterococcus faecium*, *Salmonella typhimurium*, *Shigella dysenteriae*, *Haemophilus influenzae*, *Neisseria gonorrhoeae*, *Nocardia brasiliensis*

**Fungi (9):** *Saccharomyces cerevisiae*, *Schizosaccharomyces pombe*, *Candida albicans*, *Ajellomyces capsulatus*, *Paracoccidioides lutzii*, *Cladophialophora carrionii*, *Fonsecaea pedrosoi*, *Madurella mycetomatis*, *Sporothrix schenckii*

**Plants (4):** *Arabidopsis thaliana*, *Glycine max*, *Oryza sativa*, *Zea mays*

**Protozoa (5):** *Plasmodium falciparum*, *Trypanosoma brucei*, *Trypanosoma cruzi*, *Leishmania infantum*, *Dictyostelium discoideum*

**Archaea (1):** *Methanocaldococcus jannaschii*

---

## Antibody subproject

`antibodyproduction/` is an independent pipeline that applies DE-STRESS features to a classification problem: predicting whether computationally designed scFv antibody fragments (Fleishman lab designs targeting insulin and *M. tuberculosis* ACP) will express in a yeast-display system.

**Task:** Three-class classification - Low / Medium / High expression level.
**Model:** Gaussian Naive Bayes with 10-repeat 5-fold cross-validation.
**Feature selection:** Mutual information (`feature_select_mi`) and Random Forest importance (`feature_select_rf`, 1000 trees, balanced class weights).

The pipeline loops over all combinations of scaling method (`standard`, `robust`, `minmax`) x amino acid composition included/excluded x feature selection method (`mi`/`rf`).

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
