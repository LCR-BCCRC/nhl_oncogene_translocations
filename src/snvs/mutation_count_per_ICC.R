source("src/libs.R")
library(tidyverse)
library(GAMBLR)
library(glue)
library(ggbeeswarm)
library(rstatix)
library(circlize)

source("src/snvs/functions/calc_mut_freq_window.R")
source("src/snvs/functions/get_mutation_frequency_bin_matrix.r")
source("src/snvs/functions/process_regions.R")
source("src/snvs/functions/heatmap_mutation_frequency_bin.R")

md <- read_tsv("data/metadata/breakpoint_capture_best.tsv") %>%
    left_join(select(read_tsv("data/metadata/breakpoint_capture_biopsies.tsv"), patient_id, biopsy_id, matches("partner")))

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

maf <- read_tsv("data/maf/genome_capture.hg38.maf")

maf_md <- maf %>%
    left_join(md, by = c("Tumor_Sample_Barcode" = "sample_id")) %>%
    filter(!is.na(group))

regions <- read_tsv(
    "data/region_data/regions_for_mutsig.tsv"
) %>%
    rename(Chromosome = chrom) %>%
    mutate(region = ifelse(name == "ZCCHC7", "PAX5", name))

maf_count <- maf_md %>%
    left_join(regions, by = "Chromosome", relationship = "many-to-many") %>%
    filter(Start_Position >= start & End_Position <= end) %>%
    group_by(group, region) %>%
    count(Tumor_Sample_Barcode) %>%
    ungroup() %>%
    complete(
        region,
        nesting(group, Tumor_Sample_Barcode),
        fill = list(n = 0)
    ) %>%
    left_join(select(
        md,
        Tumor_Sample_Barcode = sample_id,
        patient_id,
        biopsy_id
    ))

pvals <- maf_count %>%
    filter(region %in% c("IGHC", "IGHVDJ", "IGK", "IGL", "BCL6", "PAX5")) %>%
    mutate(region = factor(region, levels = c(
        "BCL6",
        "PAX5",
        "IGHC",
        "IGHVDJ",
        "IGK",
        "IGL"
    ))) %>%
    filter(!(region == "IGK" & n > 200)) %>%
    group_by(region) %>%
    wilcox_test(n ~ group, ref.group = "HGBCL-DH-BCL2", p.adjust.method = "BH") %>%
    add_y_position(step.increase = 0, scales = "free_y")


maf_count %>%
    filter(region %in% c("IGHC", "IGHVDJ", "IGK", "IGL", "BCL6", "PAX5")) %>%
    mutate(region = factor(region, levels = c(
        "BCL6",
        "PAX5",
        "IGHC",
        "IGHVDJ",
        "IGK",
        "IGL"
    ))) %>%
    filter(!(region == "IGK" & n > 200)) %>%
    ggplot(aes(x = group, y = n)) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = group)) +
    facet_wrap(~region, scales = "free_y", nrow = 1) +
    stat_pvalue_manual(pvals, x = "group2", hide.ns = FALSE) +
    scale_colour_manual(values = colours$group[fct_levels], name = "") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none"
    ) +
    xlab("") +
    ylab("Mutation Count")

ggsave("results/shm_heatmaps/mutation_counts_per_region.pdf", height = 4, width = 12)

# Make mini heatmap
maf_count_mat <- maf_count %>%
    filter(region %in% c("IGHC", "IGHVDJ", "IGK", "IGL", "BCL6", "PAX5")) %>%
    mutate(region = factor(region, levels = c(
        "BCL6",
        "PAX5",
        "IGHC",
        "IGHVDJ",
        "IGK",
        "IGL"
    ))) %>%
    group_by(region, group) %>%
    summarize(
        mut_rate = median(n),
    ) %>%
    ungroup() %>%
    pivot_wider(
        names_from = region,
        values_from = mut_rate
    ) %>%
    column_to_rownames("group")

pvals_mat <- pvals %>%
    select(region, group2, p.adj.signif) %>%
    pivot_wider(
        names_from = region,
        values_from = p.adj.signif
    ) %>%
    add_row(
        group2 = "HGBCL-DH-BCL2"
    ) %>%
    mutate(across(-group2, ~ ifelse(!str_detect(.x, "[*]") | is.na(.x), "", .x))) %>%
    column_to_rownames("group2")

pvals_mat <- pvals_mat[rownames(maf_count_mat), ]

col_fun <- colorRamp2(
    breaks = c(0, 5, 40),
    col = c("white", "orange", "purple")
)

heatmap_ashm <- Heatmap(maf_count_mat,
    col = col_fun,
    name = " ",
    row_labels = rownames(maf_count_mat),
    column_labels = colnames(maf_count_mat),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%s", pvals_mat[i, j]), x, y, gp = gpar(fontsize = 10))
    },
    rect_gp = gpar(col = "black", lwd = 0.5),
    heatmap_legend_param = list(
        title = "Median Mutation Count",
        legend_width = unit(4, "cm"),
        direction = "horizontal"
    ),
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_side = "top",
    show_column_names = TRUE,
    row_title = NA
)

ht <- draw(heatmap_ashm, heatmap_legend_side = "bottom")

pdf("results/shm_heatmaps/mutation_counts_per_region_collapsed.pdf", height = 2.2, width = 3.3)
print(ht)
dev.off()

md_heatmap <- md %>%
    mutate(across(matches("partner"), ~ replace_na(.x, "no break"))) %>%
    mutate(myc_partner_group = factor(
        myc_partner_group,
        levels = c(
            "IGH",
            "IGK/L",
            "recurrent non-IG",
            "non-IG",
            "no break"
        )
    ))

regions2 <- regions %>%
    rename(chrom = Chromosome) %>%
    mutate(name = factor(region, levels = c("BCL6", "MYC", "PAX5", "BCL2", "IGK", "IGL", "IGHC", "IGHM", "Emu", "IGHVDJ")))

# Make full heatmap
pdf("results/shm_heatmaps/mutation_per_region_heatmap.pdf", height = 10, width = 15)
heatmap_mutation_frequency_bin(
    regions_df = regions2,
    these_samples_metadata = rename(
        md_heatmap,
        "ICC Class" = "group",
        "MYC Partner" = "myc_partner_group"
    ),
    maf = maf_md,
    metadataColumns = c("ICC Class", "MYC Partner"),
    sortByColumns = c("ICC Class", "MYC Partner"),
    customColours = list(
        `ICC Class` = colours$group,
        `MYC Partner` = colours$myc_partner_group
    ),
    return_heatmap_obj = FALSE,
    slide_by = 500,
    window_size = 500,
    region_padding = 0,
    min_bin_recurrence = 10,
    label_regions_by = "name",
    projection = "hg38"
)
dev.off()

maf_count_dzsig <- maf_md %>%
    filter(group %in% c("HGBCL-DH-BCL2", "GCB-DLBCL")) %>%
    mutate(dz_group = case_when(
        group == "HGBCL-DH-BCL2" ~ factor(group),
        dhitsig_call == "POS" ~ factor("DHITsig+"),
        TRUE ~ factor("GCB-DLBCL")
    )) %>%
    mutate(
        dz_group = factor(dz_group,
            levels = c("GCB-DLBCL", "DHITsig+", "HGBCL-DH-BCL2")
        )
    ) %>%
    left_join(regions, by = "Chromosome", relationship = "many-to-many") %>%
    filter(Start_Position >= start & End_Position <= end) %>%
    group_by(dz_group, region) %>%
    count(Tumor_Sample_Barcode) %>%
    ungroup() %>%
    complete(
        region,
        nesting(dz_group, Tumor_Sample_Barcode),
        fill = list(n = 0)
    )

pvals_dzsig <- maf_count_dzsig %>%
    filter(region %in% c("IGHC", "IGHVDJ", "IGK", "IGL", "BCL6", "PAX5")) %>%
    mutate(region = factor(region, levels = c(
        "BCL6",
        "PAX5",
        "IGHC",
        "IGHVDJ",
        "IGK",
        "IGL"
    ))) %>%
    filter(!(region == "IGK" & n > 200)) %>%
    group_by(region) %>%
    wilcox_test(n ~ dz_group, ref.group = "HGBCL-DH-BCL2", p.adjust.method = "BH") %>%
    add_y_position(step.increase = 0, scales = "free_y")


maf_count_dzsig %>%
    filter(region %in% c("IGHC", "IGHVDJ", "IGK", "IGL", "BCL6", "PAX5")) %>%
    mutate(region = factor(region, levels = c(
        "BCL6",
        "PAX5",
        "IGHC",
        "IGHVDJ",
        "IGK",
        "IGL"
    ))) %>%
    filter(!(region == "IGK" & n > 200)) %>%
    ggplot(aes(x = dz_group, y = n)) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = dz_group)) +
    facet_wrap(~region, scales = "free_y", nrow = 1) +
    stat_pvalue_manual(pvals_dzsig, x = "group2", hide.ns = FALSE) +
    scale_colour_manual(values = colours$group[unique(as.character(maf_count_dzsig$dz_group))], name = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
    xlab("") +
    ylab("Mutation Count")

ggsave("results/shm_heatmaps/mutation_counts_per_region_dzsig.pdf", height = 4, width = 12)

# Does aSHM correlate with MYC-R in DH?

sv_data <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv") %>%
    select(patient_id, biopsy_id, myc_partner, myc_partner_igh_discrete)

pvals_myc_in_region <- maf_count %>%
    filter(group == "HGBCL-DH-BCL2") %>%
    filter(region %in% c("BCL6", "IGHC", "IGK", "IGL", "PAX5")) %>%
    left_join(sv_data) %>%
    mutate(mycr_here = factor(case_when(
        region == myc_partner_igh_discrete ~ region,
        region == "IGHC" & str_detect(myc_partner_igh_discrete, "IGH[AGE]") ~ region,
        region == myc_partner ~ region,
        TRUE ~ "Other"
    ), levels = c("BCL6", "IGHC", "IGK", "IGL", "PAX5", "Other"))) %>%
    filter(!(region == "IGK" & n > 200)) %>%
    group_by(region) %>%
    wilcox_test(n ~ mycr_here, ref.group = "Other", p.adjust.method = "BH") %>%
    add_y_position(step.increase = 0, scales = "free_y")


maf_count %>%
    filter(group == "HGBCL-DH-BCL2") %>%
    filter(region %in% c("BCL6", "IGHC", "IGK", "IGL", "PAX5")) %>%
    left_join(sv_data) %>%
    mutate(mycr_here = factor(case_when(
        region == myc_partner_igh_discrete ~ region,
        region == "IGHC" & str_detect(myc_partner_igh_discrete, "IGH[AGE]") ~ region,
        region == "IGHVDH" & str_detect(myc_partner_igh_discrete, "IGH[VDJ]") ~ region,
        region == myc_partner ~ region,
        TRUE ~ "Other"
    ), levels = c("BCL6", "IGHC", "IGK", "IGL", "PAX5", "Other"))) %>%
    ggplot(aes(
        x = mycr_here,
        y = n
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom() +
    stat_pvalue_manual(pvals_myc_in_region) +
    facet_wrap(~region, scales = "free", nrow = 1) +
    xlab("MYC Partner") +
    ylab("Mutation Count") +
    theme_bw()

ggsave("results/shm_heatmaps/mutation_counts_per_region_by_myc_partner.pdf", height = 4, width = 10)
