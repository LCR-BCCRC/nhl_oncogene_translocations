source("src/libs.R")
library(tidyverse)
library(GAMBLR)
library(glue)
library(ggbeeswarm)
library(rstatix)
library(circlize)

md <- read_tsv("data/metadata/breakpoint_capture_best.tsv")

fct_levels <- c(
    "BL", 
    "FL",
    "HGBCL-DH-BCL2", 
    "GCB-DLBCL", 
    "ABC-DLBCL"
)

md <- md %>% 
    mutate(group = ifelse(ICC_class == "DLBCL", glue("{dlbcl_call}-DLBCL"), ICC_class)) %>% 
    filter(group %in% fct_levels) %>% 
    mutate(group = factor(group, levels = fct_levels))
    


sigprofiler <- read_tsv("data/maf/sigprofiler_proportions.tsv") %>%
    filter(sample_id %in% md$sample_id) %>%
    filter(variable %in% c("SBS1", "SBS5", "SBS5", "SBS84", "SBS9", "SBS96A", "SBS96C", "SBS96F")) %>% 
    complete(sample_id, nesting(region, variable))
    
sigprofiler %>%
    count(region, variable) %>%
    print(n = 80)
    
md_heatmap <- md %>%
    select(sample_id, group) %>%
    filter(sample_id %in% sigprofiler$sample_id) %>%
    column_to_rownames("sample_id")
    
col_fun <- colorRamp2(
    breaks = c(0, 0.3, 0.7, 1),
    col = c("white", "orange", "purple", "darkblue")
)

lgd <- Legend(
    col_fun = col_fun, title = "Signature Exposure",
    legend_width = unit(4, "cm"),
    direction = "horizontal"
)

pdf("results/sigprofiler/heatmap_legend.pdf", height = 1, width = 4)
draw(lgd)
dev.off()

heatmap_by_region <- function(this_region){
    # this_region = "IGHC"

    sp_proportion_matrix <- sigprofiler %>%
        filter(region == this_region) %>%
        select(sample_id, proportion, signature = variable) %>%
        distinct(sample_id, signature, .keep_all = TRUE) %>%
        pivot_wider(
            names_from = "sample_id",
            values_from = "proportion",
            id_cols = signature
        ) %>%
        column_to_rownames("signature") %>%
        as.matrix()

    sp_proportion_matrix <- sp_proportion_matrix[, rownames(md_heatmap)]

    identical(colnames(sp_proportion_matrix), rownames(md_heatmap))


    # Make heatmap

    ha <- HeatmapAnnotation(
        " " = md_heatmap$group,
        col = list(
            " " = colours$group[c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")]
        ),
        show_legend = FALSE
    )

    col_fun <- colorRamp2(
        breaks = c(0, 0.3, 0.7, 1),
        col = c("white", "orange", "purple", "darkblue")
    )

    heatmap_ashm <- Heatmap(sp_proportion_matrix,
        col = col_fun,
        name = " ",
        top_annotation = ha,
        column_split = md_heatmap$group,
        row_labels = rownames(sp_proportion_matrix),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 10), 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        column_names_side = "top",
        show_column_names = FALSE,
        row_title = this_region,
        show_heatmap_legend = FALSE
    )

    ht <- draw(heatmap_ashm)
    
    pdf_height <- 1 + 0.5 * length(rownames(sp_proportion_matrix))
    
    pdf(paste0("results/sigprofiler/", this_region, "_sigprofiler_heatmap.pdf"), height = pdf_height, width = 8)
    print(ht)
    dev.off()
}    

all_regions <- unique(sigprofiler$region)

lapply(all_regions, heatmap_by_region)

# Plot WRCY enrichment

enrich <- list()
for(tsv in dir("data/maf/region_mafs/", pattern = ".tsv", full.names = TRUE)){
    class = str_remove(str_remove(tsv, ".*/"), "_genome.*")
    this_tsv <- read_tsv(tsv) %>% mutate(ICC_class = class)
    enrich[[class]] <- this_tsv
}
enrich <- bind_rows(enrich) %>% 
    filter(ICC_class %in% c("BL", "FL", "HGBCL-DH-BCL2", "DLBCL")) %>% 
    mutate(ICC_class = factor(ICC_class, levels = c("BL", "FL", "HGBCL-DH-BCL2", "DLBCL"))) %>% 
    filter(gene %in% c("IGH", "IGK", "IGL", "MYC", "BCL2", "BCL6", "PAX5"))
    
motif_or_mat <- enrich %>%
    select(ICC_class, gene, motif_odds_ratio) %>% 
    pivot_wider(
        names_from = gene, 
        values_from = motif_odds_ratio
    ) %>% 
    column_to_rownames("ICC_class") %>%
    select(c("BCL2", "MYC", "PAX5", "BCL6", "IGH", "IGK", "IGL"))%>%
    as.matrix()

motif_pval_mat <- enrich %>%
    select(ICC_class, gene, motif_bias) %>%
    mutate(motif_bias = case_when(
        motif_bias < 0.0001 ~ "****", 
        motif_bias < 0.001 ~ "***", 
        motif_bias < 0.01 ~ "**", 
        motif_bias < 0.05 ~ "*", 
        TRUE ~ "ns"
    )) %>% 
    pivot_wider(
        names_from = gene,
        values_from = motif_bias
    ) %>%
    column_to_rownames("ICC_class") %>%
    select(c("BCL2", "MYC", "PAX5", "BCL6", "IGH", "IGK", "IGL")) %>%
    as.matrix()
    
site_or_mat <- enrich %>%
    select(ICC_class, gene, motif_site_odds_ratio) %>%
    pivot_wider(
        names_from = gene,
        values_from = motif_site_odds_ratio
    ) %>%
    column_to_rownames("ICC_class") %>%
    select(c("BCL2", "MYC", "PAX5", "BCL6", "IGH", "IGK", "IGL")) %>%
    as.matrix()

site_pval_mat <- enrich %>%
    select(ICC_class, gene, motif_site_bias) %>%
    mutate(motif_site_bias = case_when(
        motif_site_bias < 0.0001 ~ "****",
        motif_site_bias < 0.001 ~ "***",
        motif_site_bias < 0.01 ~ "**",
        motif_site_bias < 0.05 ~ "*",
        TRUE ~ "ns"
    )) %>%
    pivot_wider(
        names_from = gene,
        values_from = motif_site_bias
    ) %>%
    column_to_rownames("ICC_class") %>%
    select(c("BCL2", "MYC", "PAX5", "BCL6", "IGH", "IGK", "IGL"))%>% 
    as.matrix()
    
col_fun <- colorRamp2(
    breaks = c(0, 1, 4),
    col = c("blue", "white", "red")
)

heatmap_motif <- Heatmap(motif_or_mat,
    col = col_fun,
    name = " ",
    row_labels = rownames(motif_or_mat),
    column_labels = colnames(motif_or_mat),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%s", motif_pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
    },
    rect_gp = gpar(col = "black", lwd = 0.5),
    heatmap_legend_param = list(
        title = "Motif Odds Ratio\n(vs Expected)",
        legend_width = unit(5, "cm"),
        direction = "horizontal"
    ),
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_side = "top",
    show_column_names = TRUE,
    row_title = NA
)

ht <- draw(heatmap_motif, heatmap_legend_side = "bottom")

pdf("results/sigprofiler/wrcy_motif_heatmap.pdf", height = 2.5, width = 4)
ht
dev.off()

heatmap_motif <- Heatmap(site_or_mat,
    col = col_fun,
    name = " ",
    row_labels = rownames(site_or_mat),
    column_labels = colnames(site_or_mat),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%s", site_pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
    },
    rect_gp = gpar(col = "black", lwd = 0.5),
    heatmap_legend_param = list(
        title = "Motif Odds Ratio\n(vs Expected)",
        legend_width = unit(5, "cm"),
        direction = "horizontal"
    ),
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_side = "top",
    show_column_names = TRUE,
    row_title = NA
)

ht <- draw(heatmap_motif, heatmap_legend_side = "bottom")

pdf("results/sigprofiler/wrCy_site_heatmap.pdf", height = 2.5, width = 4)
ht
dev.off()
