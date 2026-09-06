library(TreeDist)

load("~/GitRepos/proteinPropertiesReconstructTreeOfLife/data/raw_data/randomTreeDistances.rda")

directory_path <- "~/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/hier_clustering_avg_by_org_aa_composition/af2/"

nwk_files <- list.files(directory_path, pattern = "\\.nwk$", full.names = TRUE)

reference_nwk_file <- "~/GitRepos/proteinPropertiesReconstructTreeOfLife/data/processed_data/ncbi_phylo_tree.phy"

reference_tree <- ape::read.tree(reference_nwk_file)
print(reference_tree)

distances_list <- list()

for (nwk_file in nwk_files) {
  current_tree <- ape::read.tree(nwk_file)
  print(current_tree)

  distance <- TreeDistance(reference_tree, current_tree)

  file_name <- basename(nwk_file)
  distances_list[[file_name]] <- distance
}

print(distances_list)

distances_df <- as.data.frame(distances_list)

csv_file <- "~/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/hier_clustering_avg_by_org_aa_composition/af2/tree_distances.csv"

write.csv(distances_df, file = csv_file, row.names = FALSE)

expectedCID <- randomTreeDistances["cid", "mean", "48"]
sdCID <- randomTreeDistances["cid", "sd", "48"]
print(expectedCID)
print(sdCID)
