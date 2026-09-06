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
│   │   ├── download_af2_data.py
│   │   ├── extract_plddt_score.py
│   │   ├── download_af2_org_sci_name.py
│   │   ├── download_uniprot_data.py
│   │   ├── download_mitochondrial_encoded_proteins_uniprot.py
│   │   ├── filtering_af2db_clustering_data.py
│   │   ├── download_structures_af2_pdb_files.py
│   │   └── download_af2db_files_from_server.py
│   │
│   ├── data_prep/
│   │   ├── data_prep.py
│   │   ├── data_prep_tools.py
│   │   ├── data_prep_aa_composition.py
│   │   ├── data_prep_single_proteins.py
│   │   ├── uniprot_af2db_data_prep.py
│   │   ├── af2db_clusters_data_prep.py
│   │   ├── random_select_structures_by_org_and_cluster.py
│   │   ├── csv_to_fasta_file.py
│   │   ├── destress_features_by_dataset.py
│   │   ├── dssp_and_aa_proportions.py
│   │   ├── sequence_length_histogram.py
│   │   └── correlation_dssp_hb_isoelectric_aggrescan.py
│   │
│   ├── dim_red/
│   │   ├── dim_red_tools.py
│   │   ├── pca_analysis.py
│   │   ├── pca_avg_destress_metrics_by_org.py
│   │   ├── pca_all_mitochondrial_proteins.py
│   │   ├── pca_subcellular_location.py
│   │   ├── pca_analysis_spectral_all_organisms.py
│   │   ├── subcellular_location_analysis.py
│   │   ├── single_af2db_cluster_pca_analysis.py
│   │   ├── single_af2db_cluster_pca_analysis_filtered.py
│   │   ├── single_protein_pca_analysis.py
│   │   ├── plddt_pca_correlation_sod.py
│   │   ├── plddt_pca_correlation_sod_by_kingdom.py
│   │   └── plddt_vs_feature_pc1_correlations_sod.py
│   │
│   ├── clustering/
│   │   ├── clustering_tools.py
│   │   ├── hierarchical_clustering_avg_destress_metrics.py
│   │   ├── hierarchical_clustering_avg_destress_metrics_single_af2db_cluster_filtered.py
│   │   ├── hierarchical_clustering_aa_composition.py
│   │   ├── kmeans_avg_destress_metrics.py
│   │   ├── kmeans_aa_composition.py
│   │   ├── tree_dist.R
│   │   ├── tree_dist_single_af2db_cluster_filtered.R
│   │   ├── tree_dist_aa_composition.R
│   │   ├── test_tree_dist_single_af2db_cluster_filtered.R
│   │   ├── plot_cid_summary_single_af2db_cluster_filtered.py
│   │   ├── compare_cid_destress_vs_aa_composition.py
│   │   ├── plot_cid_destress_vs_aa_composition.py
│   │   ├── format_tree_distances_table.py
│   │   ├── summarise_cluster_descriptions.py
│   │   └── summarise_cluster_descriptions_all.py
│   │
│   ├── pymol_calign_script.pml
│   ├── pymol_calign_hide_nterm_script.pml
│   ├── pymol_calign_trim_structure_script.py
│   └── pymol_save_all_objects.py
│
├── antibodyproduction/
│   └── src/
│       ├── data_prep.py
│       ├── data_prep_tools.py
│       ├── feature_selection.py
│       ├── feature_selection_tools.py
│       ├── dim_red.py
│       ├── dim_red_tools.py
│       ├── model_building.py
│       ├── model_building_tools.py
│       └── pymol_calign_script.pml
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

Joins DE-STRESS output with organism metadata and FoldSeek cluster assignments. Removes features with >5% missing values, drops constant features, normalizes energy terms by sequence length, removes highly correlated features (Spearman |r| > 0.6), and scales using three methods in parallel: Standard (z-score), Robust (IQR), and MinMax. Produces a full dataset and a non-redundant dataset (one structure per organism x FoldSeek cluster, random seed 42).

| Script | Role |
|---|---|
| `data_prep.py` | Main execution script - calls `data_prep_tools.py` |
| `data_prep_tools.py` | Shared utility library |
| `uniprot_af2db_data_prep.py` | Joins UniProt + AFDB metadata |
| `af2db_clusters_data_prep.py` | Cleans filtered FoldSeek cluster table |
| `random_select_structures_by_org_and_cluster.py` | Generates the non-redundant structure list |
| `data_prep_single_proteins.py` | Data preparation for single protein family analyses |
| `data_prep_aa_composition.py` | Prepares amino acid composition baseline data: extracts the 20 standard `composition_*` columns and scales with all three methods |
| `csv_to_fasta_file.py` | Converts a CSV of sequences to FASTA format |
| `destress_features_by_dataset.py` | Summarises which DE-STRESS features are used in each analysis across the paper |
| `dssp_and_aa_proportions.py` | Computes DSSP bin, secondary structure residue, and amino acid proportions across datasets |
| `sequence_length_histogram.py` | Plots protein sequence length distributions across datasets |
| `correlation_dssp_hb_isoelectric_aggrescan.py` | Spearman correlation analysis between selected feature pairs (DSSP/HB energies, isoelectric point/Aggrescan3D, VdW/packing density) |

### Stage 3 - Dimensionality reduction (`src/dim_red/`)

Runs PCA across dataset x scaling method combinations. Averaged scaled DE-STRESS metrics per organism feed into hierarchical clustering.

| Script | Role |
|---|---|
| `dim_red_tools.py` | Shared utility library |
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

**DE-STRESS clustering**

| Script | Role |
|---|---|
| `clustering_tools.py` | Shared utility library |
| `hierarchical_clustering_avg_destress_metrics.py` | Hierarchical clustering across all organisms - dendrograms and `.nwk` trees for all linkage x distance combinations |
| `kmeans_avg_destress_metrics.py` | K-means evaluation (k = 2-20) |
| `tree_dist.R` | Computes CID between all reconstructed trees and the NCBI reference phylogeny |
| `hierarchical_clustering_avg_destress_metrics_single_af2db_cluster_filtered.py` | Hierarchical clustering for each of the 11 filtered FoldSeek clusters |
| `plot_cid_summary_single_af2db_cluster_filtered.py` | Plots CID summary across the 11 filtered clusters |
| `tree_dist_single_af2db_cluster_filtered.R` | CID computation for the filtered single-cluster trees |
| `test_tree_dist_single_af2db_cluster_filtered.R` | Validates tip label normalisation and CID computation for filtered cluster trees |
| `summarise_cluster_descriptions.py` | Summarises protein descriptions within a single FoldSeek cluster |
| `summarise_cluster_descriptions_all.py` | Summarises protein descriptions across all 22 FoldSeek clusters |

**Amino acid composition baseline**

| Script | Role |
|---|---|
| `hierarchical_clustering_aa_composition.py` | Hierarchical clustering using organism-averaged AA composition across all linkage x distance combinations |
| `kmeans_aa_composition.py` | K-means evaluation using AA composition (k = 2-20) |
| `tree_dist_aa_composition.R` | Computes CID between AA composition trees and the NCBI reference phylogeny |
| `compare_cid_destress_vs_aa_composition.py` | Joins DE-STRESS and AA composition CID results and computes per-combination differences |
| `plot_cid_destress_vs_aa_composition.py` | Plots CID difference between DE-STRESS and AA composition clustering |
| `format_tree_distances_table.py` | Reformats a `tree_distances.csv` into a publication-ready table sorted ascending by CID |

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

| Script | Role |
|---|---|
| `data_prep.py` | Main execution script - processes expression and DE-STRESS data into train/test splits; calls `data_prep_tools.py` |
| `data_prep_tools.py` | Shared utility library for data preparation |
| `feature_selection.py` | Runs mutual information and random forest feature selection across all scaling x composition combinations |
| `feature_selection_tools.py` | Shared utility library for feature selection |
| `dim_red.py` | PCA analysis of the antibody expression data |
| `dim_red_tools.py` | Shared utility library for dimensionality reduction |
| `model_building.py` | Trains Gaussian Naive Bayes classifiers across all feature sets and evaluates with 10×5-fold cross-validation |
| `model_building_tools.py` | Shared utility library for model training and evaluation |
| `pymol_calign_script.pml` | Structural alignment of antibody PDB files in PyMOL |

