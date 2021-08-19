## code to prepare pbmc_small_moduli

library(Seurat)
set.seed(123)

data("pbmc_small")

# normalize data
pbmc_small <- SCTransform(pbmc_small, verbose = F)

# cluster genes
pam <- gene_medioid_clustering(pbmc_small, 6)

# compute moduli
pbmc_small_moduli <- get_moduli(pbmc_small, gene.membership = pam$clustering, npcs = 3, verbose = F, filebacked = F)

usethis::use_data(pbmc_small_moduli, overwrite = TRUE)
