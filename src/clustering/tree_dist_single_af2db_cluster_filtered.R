if (!requireNamespace("TreeDist", quietly = TRUE)) install.packages("TreeDist", repos = "https://cloud.r-project.org")

library(TreeDist)
library(ape)

# Load pre-computed CID mean and SD for random trees of various sizes.
# Used to contextualise observed CID scores (how much better than random?).
load("~/GitRepos/proteinPropertiesReconstructTreeOfLife/data/raw_data/randomTreeDistances.rda")

# Paths
base_path <- "~/GitRepos/proteinPropertiesReconstructTreeOfLife/"
labels_path <- paste0(base_path, "data/processed_data/af2/labels.csv")

# Full reference tree (all organisms, not eukaryotes-only) since clusters span all domains.
reference_nwk_file <- paste0(base_path, "data/processed_data/ncbi_phylo_tree.phy")

# Output from hierarchical_clustering_avg_destress_metrics_single_af2db_cluster_filtered.py
output_base_path <- paste0(base_path, "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/")

# 1. Build cluster list from pLDDT-filtered processed data --------------------------------

labels_df <- read.csv(labels_path)

# Count unique organisms per cluster from the actual data used in analysis.
org_counts <- aggregate(organism_scientific_name ~ cluster_representative,
                        data = labels_df,
                        FUN = function(x) length(unique(x)))
colnames(org_counts)[2] <- "n_orgs"
af2db_cluster_list_df <- org_counts[org_counts$n_orgs >= 40, ]

# Attach most common UniProt description per cluster for filtering.
get_top_desc <- function(cluster_id) {
  descs <- labels_df$uniprot_description[labels_df$cluster_representative == cluster_id]
  names(sort(table(descs), decreasing = TRUE))[1]
}
af2db_cluster_list_df$top_description <- sapply(
  af2db_cluster_list_df$cluster_representative, get_top_desc
)

# Exclude uncharacterised and domain-containing proteins — same criteria as Python scripts.
exclude_pattern <- "Uncharacterized protein|domain-containing protein"
af2db_cluster_list_df <- af2db_cluster_list_df[
  !grepl(exclude_pattern, af2db_cluster_list_df$top_description), ]
rownames(af2db_cluster_list_df) <- NULL

cat(sprintf("Running tree distance for %d clusters:\n", nrow(af2db_cluster_list_df)))
for (i in seq_len(nrow(af2db_cluster_list_df))) {
  cat(sprintf("  %s: %s (%d orgs)\n",
    af2db_cluster_list_df$cluster_representative[i],
    af2db_cluster_list_df$top_description[i],
    af2db_cluster_list_df$n_orgs[i]))
}

# 2. Load full NCBI reference tree --------------------------------------------------------

full_reference_tree <- ape::read.tree(reference_nwk_file)
# Strip single quotes and spaces from tip labels.
# The .phy file wraps names in quotes which ape preserves, but the clustering .nwk files
# have unquoted labels which ape reads by stripping spaces (e.g. "Homo sapiens" -> "Homosapiens").
# Both are normalised the same way so tip labels match.
full_reference_tree$tip.label <- gsub("'| ", "", full_reference_tree$tip.label)
cat(sprintf("\nLoaded NCBI reference tree with %d tips\n", length(full_reference_tree$tip.label)))

# 3. Loop over clusters -------------------------------------------------------------------

all_distances <- list()

for (i in seq_len(nrow(af2db_cluster_list_df))) {

  af2db_cluster <- af2db_cluster_list_df$cluster_representative[i]
  n_orgs        <- af2db_cluster_list_df$n_orgs[i]
  top_desc      <- af2db_cluster_list_df$top_description[i]

  cat(sprintf("\n--- %s: %s (%d orgs) ---\n", af2db_cluster, top_desc, n_orgs))

  # Organisms in this cluster — spaces stripped to match ape's unquoted label parsing.
  cluster_orgs <- gsub(" ", "", unique(
    labels_df$organism_scientific_name[labels_df$cluster_representative == af2db_cluster]
  ))

  # Prune the reference tree to only organisms present in this cluster.
  # Trees must share the same tip set for TreeDistance to be valid.
  tips_to_drop <- full_reference_tree$tip.label[
    !full_reference_tree$tip.label %in% cluster_orgs
  ]
  pruned_reference_tree <- ape::drop.tip(full_reference_tree, tips_to_drop)
  cat(sprintf("  Reference tree pruned to %d shared tips\n",
              length(pruned_reference_tree$tip.label)))

  # n for expected CID lookup — based on shared tips after pruning
  n_shared <- length(pruned_reference_tree$tip.label)
  expectedCID <- randomTreeDistances["cid", "mean", as.character(n_shared)]
  sdCID       <- randomTreeDistances["cid", "sd",   as.character(n_shared)]

  # Find all .nwk files across scaling method subdirectories for this cluster
  cluster_dir <- paste0(output_base_path, af2db_cluster, "/")
  nwk_files   <- list.files(cluster_dir, pattern = "\\.nwk$",
                             full.names = TRUE, recursive = TRUE)

  if (length(nwk_files) == 0) {
    cat(sprintf("  No .nwk files found in %s — skipping\n", cluster_dir))
    next
  }
  cat(sprintf("  Found %d .nwk files\n", length(nwk_files)))

  distances_list <- list()

  for (nwk_file in nwk_files) {

    current_tree <- ape::read.tree(nwk_file)
    distance     <- TreeDistance(pruned_reference_tree, current_tree)
    file_name <- basename(nwk_file)
    # Include the scaling method subdirectory in the label for clarity
    rel_path  <- sub(cluster_dir, "", nwk_file, fixed = TRUE)
    distances_list[[rel_path]] <- distance
  }

  distances_df <- data.frame(
    file         = names(distances_list),
    tree_dist    = unlist(distances_list),
    expected_cid = expectedCID,
    sd_cid       = sdCID,
    n_shared     = n_shared,
    cluster      = af2db_cluster,
    description  = top_desc,
    row.names    = NULL
  )

  # Save per-cluster CSV alongside the .nwk files
  csv_file <- paste0(cluster_dir, "tree_distances.csv")
  write.csv(distances_df, file = csv_file, row.names = FALSE)
  cat(sprintf("  Saved to %s\n", csv_file))

  all_distances[[af2db_cluster]] <- distances_df
}

# 4. Save combined summary across all clusters --------------------------------------------

all_distances_df <- do.call(rbind, all_distances)
rownames(all_distances_df) <- NULL

summary_csv <- paste0(output_base_path, "tree_distances_all_clusters.csv")
write.csv(all_distances_df, file = summary_csv, row.names = FALSE)
cat(sprintf("\nSaved combined summary to %s\n", summary_csv))
