library(tidyverse)
library(ggfortify)
library(PRPS)
library(glue)

source("src/libs.R")

# Load metadata and subset to included classes
class_include <- c("BL", "FL", "HGBCL-DH-BCL2", "DLBCL", "HGBCL-DH-BCL6")

md <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    filter(seq_type == "mrna") %>%
    mutate(dlbcl_call = ifelse(is.na(dlbcl_call), "ND", dlbcl_call)) %>% 
    group_by(patient_id, biopsy_id) %>% 
    slice_max(protocol == "Ribodepletion", n=1, with_ties = FALSE) %>% 
    ungroup() %>% 
    filter(ICC_class %in% class_include) %>% 
    mutate(group = case_when(
        ICC_class == "DLBCL" & dlbcl_call != "ND" ~ glue("{dlbcl_call}-DLBCL"), 
        TRUE ~ ICC_class
    ))
    
# Set colour palettes and factor levels
group_levels <- c(
    "BL",
    "FL",
    "HGBCL-DH-BCL2",
    "DLBCL",
    "GCB-DLBCL",
    "UNCLASS-DLBCL",
    "ABC-DLBCL",
    "HGBCL-DH-BCL6"
)

colours_groups <- c(
    colours$ICC_class,
    colours$group,
    "UNCLASS-DLBCL" = unname(get_gambl_colours("coo")["UNCLASS"])
)[group_levels]

# Get a list of COO and DHITsig genes 
load(system.file("extdata", "WrightCOO.rda", package = "PRPS"))
coo_genes <- rownames(WrightCOO)
load(system.file("extdata", "DHITsigENSG.rda", package = "PRPS"))
dhitsig_genes <- names(DHITsigENSG)

# Load expression data
expr <- read_tsv("results/expression/uncorrected_vst_matrix.tsv")

# Convert expression to matrix and subset to COO and DHITsig genes
expr_mat <- expr %>%
    select(-gene_id, -hgnc_symbol) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
    column_to_rownames("ensembl_gene_id")

expr_mat <- as.matrix(expr_mat[
    rownames(expr_mat) %in% unique(c(coo_genes, dhitsig_genes)), 
    colnames(expr_mat) %in% md$sample_id
])

# Match the metadata to the expression matrix
md_mat <- column_to_rownames(md, "sample_id")
md_mat <- md_mat[colnames(expr_mat), ]
md_mat$group <- factor(md_mat$group, levels = group_levels)

identical(colnames(expr_mat), rownames(md_mat))

# Perform PCA
pca_res <- prcomp(t(expr_mat))

# Plot PCA by biological variables and batches
p1 <- autoplot(
    pca_res, 
    data = md_mat, 
    colour = "group"
) + 
scale_colour_manual(name = "", values = colours_groups) + 
xlim(-0.10, 0.10) + 
theme(legend.position = "bottom") + 
ggtitle("Tumor Classification")
p1
ggsave("results/expression/uncorrected_PCA_by_ICC_class.pdf", height = 6, width = 8)

p2 <- autoplot(
    pca_res,
    data = md_mat,
    colour = "protocol"
) +
    xlim(-0.10, 0.10) +
    theme(
        legend.position = "bottom", 
    legend.title = element_blank(), 
    legend.margin = margin()) + 
    ggtitle("Protocol")
p2
ggsave("results/expression/uncorrected_PCA_by_protocol.pdf", height = 6, width = 8)

p3 <- autoplot(
    pca_res,
    data = md_mat,
    colour = "preservation"
) +
xlim(-0.10, 0.10) +
theme(legend.position = "bottom", legend.title = element_blank()) + 
ggtitle("Preservation")
p3
ggsave("results/expression/uncorrected_PCA_by_preservation.pdf", height = 6, width = 8)


cowplot::plot_grid(
    p1 + guides(colour = guide_legend(nrow = 3)), p2, p3, 
    nrow = 1, 
    align = "hv", 
    axis = "bt"
)
ggsave("results/expression/uncorrected_combined_PCA.pdf", height = 5, width = 15)
