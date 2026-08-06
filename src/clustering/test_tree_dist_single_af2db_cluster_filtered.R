library(ape)
library(TreeDist)

BASE        <- "~/GitRepos/proteinPropertiesReconstructTreeOfLife/"
LABELS      <- paste0(BASE, "data/processed_data/af2/labels.csv")
REF_NWK     <- paste0(BASE, "data/processed_data/ncbi_phylo_tree.phy")
OUTPUT_BASE <- paste0(BASE, "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/")

pass <- 0
fail <- 0

check <- function(label, condition, detail = "") {
  if (condition) {
    cat(sprintf("  PASS  %s\n", label))
    pass <<- pass + 1
  } else {
    msg <- sprintf("  FAIL  %s", label)
    if (nchar(detail) > 0) msg <- paste0(msg, " — ", detail)
    cat(msg, "\n")
    fail <<- fail + 1
  }
}

# ---------------------------------------------------------------------------
cat("=== 1. Reference tree loading and tip label normalisation ===\n")

full_ref <- ape::read.tree(REF_NWK)
raw_tips  <- full_ref$tip.label
full_ref$tip.label <- gsub("'| ", "", full_ref$tip.label)
norm_tips <- full_ref$tip.label

check("Reference tree loads",        !is.null(full_ref))
check("Reference tree has tips",     length(norm_tips) > 0)
check("Tip count unchanged",         length(raw_tips) == length(norm_tips))
check("No quotes in normalised tips", !any(grepl("'", norm_tips)))
check("No spaces in normalised tips", !any(grepl(" ", norm_tips)))

cat(sprintf("  INFO  %d tips in reference tree\n", length(norm_tips)))
cat("\n  Tip label transformation (reference tree, first 5):\n")
cat(sprintf("    %-40s  ->  %s\n", "RAW", "NORMALISED"))
for (j in seq_len(min(5, length(raw_tips)))) {
  cat(sprintf("    %-40s  ->  %s\n", raw_tips[j], norm_tips[j]))
}

# ---------------------------------------------------------------------------
cat("\n=== 2. Cluster list construction from labels.csv ===\n")

labels_df <- read.csv(LABELS)

check("labels.csv loads",                        nrow(labels_df) > 0)
check("cluster_representative column present",   "cluster_representative"   %in% colnames(labels_df))
check("organism_scientific_name column present", "organism_scientific_name" %in% colnames(labels_df))
check("uniprot_description column present",      "uniprot_description"      %in% colnames(labels_df))

org_counts <- aggregate(organism_scientific_name ~ cluster_representative,
                        data = labels_df, FUN = function(x) length(unique(x)))
colnames(org_counts)[2] <- "n_orgs"
af2db <- org_counts[org_counts$n_orgs >= 40, ]

get_top_desc <- function(cluster_id) {
  descs <- labels_df$uniprot_description[labels_df$cluster_representative == cluster_id]
  names(sort(table(descs), decreasing = TRUE))[1]
}
af2db$top_description <- sapply(af2db$cluster_representative, get_top_desc)

exclude_pattern <- "Uncharacterized protein|domain-containing protein"
af2db <- af2db[!grepl(exclude_pattern, af2db$top_description), ]
rownames(af2db) <- NULL

check("All clusters have n_orgs >= 40",               all(af2db$n_orgs >= 40))
check("Exactly 11 clusters after exclusion filter",   nrow(af2db) == 11,
      sprintf("found %d", nrow(af2db)))
check("No 'Uncharacterized protein' in included",
      !any(grepl("Uncharacterized protein",  af2db$top_description)))
check("No 'domain-containing protein' in included",
      !any(grepl("domain-containing protein", af2db$top_description)))

cat(sprintf("  INFO  %d clusters included:\n", nrow(af2db)))
for (i in seq_len(nrow(af2db))) {
  cat(sprintf("          %s: %s (%d orgs)\n",
    af2db$cluster_representative[i], af2db$top_description[i], af2db$n_orgs[i]))
}

# ---------------------------------------------------------------------------
cat("\n=== 3. Per-cluster tip labels, pruning, and .nwk files ===\n")

for (i in seq_len(nrow(af2db))) {

  cluster_id <- af2db$cluster_representative[i]
  top_desc   <- af2db$top_description[i]
  n_orgs     <- af2db$n_orgs[i]

  cat(sprintf("\n  --- %s: %s ---\n", cluster_id, top_desc))

  # Organism names with spaces stripped — must match ape's unquoted label parsing
  cluster_orgs <- gsub(" ", "", unique(
    labels_df$organism_scientific_name[labels_df$cluster_representative == cluster_id]
  ))

  check(sprintf("%s: no spaces in cluster_orgs",          cluster_id), !any(grepl(" ", cluster_orgs)))
  check(sprintf("%s: cluster_orgs count matches n_orgs",  cluster_id),
        length(cluster_orgs) == n_orgs,
        sprintf("got %d expected %d", length(cluster_orgs), n_orgs))

  missing_from_ref <- cluster_orgs[!cluster_orgs %in% full_ref$tip.label]
  check(sprintf("%s: all cluster orgs in reference tree", cluster_id),
        length(missing_from_ref) == 0,
        if (length(missing_from_ref) > 0) paste("missing:", paste(missing_from_ref, collapse = ", ")) else "")

  # Prune reference tree to cluster organisms
  tips_to_drop     <- full_ref$tip.label[!full_ref$tip.label %in% cluster_orgs]
  pruned_ref       <- ape::drop.tip(full_ref, tips_to_drop)
  n_shared         <- length(pruned_ref$tip.label)

  check(sprintf("%s: pruned reference has correct tip count", cluster_id),
        n_shared == n_orgs, sprintf("got %d expected %d", n_shared, n_orgs))
  check(sprintf("%s: pruned reference tip set equals cluster_orgs", cluster_id),
        setequal(pruned_ref$tip.label, cluster_orgs))

  # Print tip labels at each stage for first 5 organisms (alphabetical)
  example_orgs <- sort(cluster_orgs)[seq_len(min(5, length(cluster_orgs)))]
  nwk_example  <- list.files(paste0(OUTPUT_BASE, cluster_id, "/standard"),
                              pattern = "\\.nwk$", full.names = TRUE)[1]
  ctree_example <- ape::read.tree(nwk_example)

  # Raw organism name from labels.csv (before space stripping)
  raw_orgs <- unique(
    labels_df$organism_scientific_name[labels_df$cluster_representative == cluster_id]
  )
  raw_orgs_sorted <- sort(raw_orgs)[seq_len(min(5, length(raw_orgs)))]

  cat(sprintf("    Tip label stages (first 5, %s):\n", cluster_id))
  cat(sprintf("    %-35s  %-35s  %-35s  %s\n",
              "labels.csv (raw)", "cluster_orgs (spaces stripped)",
              "pruned ref tip", "clustering tree tip"))
  for (j in seq_along(raw_orgs_sorted)) {
    raw_o    <- raw_orgs_sorted[j]
    norm_o   <- gsub(" ", "", raw_o)
    ref_tip  <- if (norm_o %in% pruned_ref$tip.label)   norm_o else "NOT IN REF"
    clust_tip <- if (norm_o %in% ctree_example$tip.label) norm_o else "NOT IN TREE"
    cat(sprintf("    %-35s  %-35s  %-35s  %s\n", raw_o, norm_o, ref_tip, clust_tip))
  }
  cat("\n")

  # Locate .nwk files for this cluster
  cluster_dir <- paste0(OUTPUT_BASE, cluster_id, "/")
  nwk_files   <- list.files(cluster_dir, pattern = "\\.nwk$", recursive = TRUE, full.names = TRUE)

  check(sprintf("%s: .nwk files exist",                      cluster_id), length(nwk_files) > 0)
  check(sprintf("%s: 39 .nwk files (13 combos x 3 scaling)", cluster_id),
        length(nwk_files) == 39, sprintf("found %d", length(nwk_files)))

  tip_mismatches <- c()
  cid_values     <- c()

  for (nwk in nwk_files) {
    ctree <- ape::read.tree(nwk)

    if (!setequal(ctree$tip.label, pruned_ref$tip.label)) {
      tip_mismatches <- c(tip_mismatches, basename(nwk))
    }

    cid_values <- c(cid_values, TreeDistance(pruned_ref, ctree))
  }

  check(sprintf("%s: all clustering trees match pruned reference tip set", cluster_id),
        length(tip_mismatches) == 0,
        if (length(tip_mismatches) > 0)
          sprintf("%d mismatched: %s", length(tip_mismatches), paste(head(tip_mismatches, 3), collapse = ", "))
        else "")
  check(sprintf("%s: no NA CID values",       cluster_id), !any(is.na(cid_values)))
  check(sprintf("%s: all CID values in [0,1]", cluster_id),
        all(cid_values >= 0 & cid_values <= 1, na.rm = TRUE))
  check(sprintf("%s: all observed CIDs below expected random", cluster_id), {
    load(paste0(BASE, "data/raw_data/randomTreeDistances.rda"))
    exp_cid <- randomTreeDistances["cid", "mean", as.character(n_shared)]
    all(cid_values < exp_cid, na.rm = TRUE)
  })

  cat(sprintf("    INFO  CID: min=%.4f  max=%.4f  mean=%.4f\n",
              min(cid_values, na.rm = TRUE),
              max(cid_values, na.rm = TRUE),
              mean(cid_values, na.rm = TRUE)))
}

# ---------------------------------------------------------------------------
cat("\n=== 4. Summary CSV integrity ===\n")

summary_csv <- paste0(OUTPUT_BASE, "tree_distances_all_clusters.csv")
check("tree_distances_all_clusters.csv exists", file.exists(summary_csv))

if (file.exists(summary_csv)) {
  df <- read.csv(summary_csv)
  expected_cols <- c("file", "tree_dist", "expected_cid", "sd_cid", "n_shared", "cluster", "description")
  check("Expected columns present",
        all(expected_cols %in% colnames(df)),
        paste("missing:", paste(setdiff(expected_cols, colnames(df)), collapse = ", ")))
  check("11 distinct clusters in summary",
        length(unique(df$cluster)) == 11, sprintf("found %d", length(unique(df$cluster))))
  check("39 rows per cluster",
        all(table(df$cluster) == 39),
        paste(names(which(table(df$cluster) != 39)), collapse = ", "))
  check("No NA tree_dist values",      !any(is.na(df$tree_dist)))
  check("All tree_dist in [0, 1]",    all(df$tree_dist >= 0 & df$tree_dist <= 1, na.rm = TRUE))
  check("All expected_cid in [0, 1]", all(df$expected_cid >= 0 & df$expected_cid <= 1, na.rm = TRUE))
  check("All observed CIDs below expected random", all(df$tree_dist < df$expected_cid, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
cat(sprintf("\n=== Results: %d passed, %d failed ===\n", pass, fail))
if (fail > 0) quit(status = 1)
