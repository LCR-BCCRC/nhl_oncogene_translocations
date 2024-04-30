source("src/libs.R")
library(tidyverse)
library(ComplexHeatmap)
library(GAMBLR)
library(data.table)
library(cowplot)
library(circlize)

source("src/snvs/functions/calc_mut_freq_window.R")
source("src/snvs/functions/get_mutation_frequency_bin_matrix.r")
source("src/snvs/functions/process_regions.R")
source("src/snvs/functions/heatmap_mutation_frequency_bin.R")

gambl_meta <- get_gambl_metadata(seq_type_filter = c("capture", "genome")) %>%
    filter(pathology == "DLBCL")


# Load metadata. Select one representative sample.
md <- read_tsv("data/metadata/breakpoint_capture_best.tsv")

# Load SV data
sv_data_cat <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv") %>%
    left_join(select(md, patient_id, biopsy_id, sample_id)) %>%
    mutate(ICC_class = factor(ICC_class, levels = c("BL", "FL", "DLBCL", "HGBCL-DH-BCL2", "HGBCL-DH-BCL6"))) %>%
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
    )) %>%
    filter(!is.na(ICC_class))

# Load MAF data
maf_data <- read_tsv("data/maf/genome_capture.WRCY.hg38.maf")
maf <- maf_data


panel_regions <- fread(
    "data/region_data/capture_TE99028370--hg38.bed",
    col.names = c("chrom", "start", "end", "name")
)

pdf("results/shm_heatmaps/IGH_shm_heatmap.pdf", height = 10, width = 12)
heatmap_mutation_frequency_bin(
    regions_df = panel_regions[panel_regions$name == "IGH", ],
    these_samples_metadata = sv_data_cat,
    maf = maf_data,
    metadataColumns = c("ICC_class", "myc_partner_group"),
    sortByColumns = c("ICC_class", "myc_partner_group"),
    slide_by = 1000,
    window_size = 1000,
    return_heatmap_obj = FALSE,
    min_bin_recurrence = 25,
    projection = "hg38"
)
dev.off()

pdf("results/shm_heatmaps/HGBCL-DH-BCL2_IGH_shm_heatmap.pdf", height = 8, width = 12)
heatmap_mutation_frequency_bin(
    regions_df = panel_regions[panel_regions$name == "IGH", ],
    these_samples_metadata = sv_data_cat[sv_data_cat$ICC_class == "HGBCL-DH-BCL2", ],
    maf = maf_data,
    metadataColumns = c("ICC_class", "myc_partner_group"),
    sortByColumns = c("ICC_class", "myc_partner_group"),
    slide_by = 1000,
    window_size = 1000,
    return_heatmap_obj = FALSE,
    min_bin_recurrence = 25,
    projection = "hg38"
)
dev.off()

bcl2_mat <- get_mutation_frequency_bins(
    maf = maf_data,
    regions_df = panel_regions[panel_regions$name == "BCL2", ],
    these_samples_metadata = sv_data_cat[sv_data_cat$ICC_class %in% c("FL", "DLBCL", "HGBCL-DH-BCL2"), ],
    region_padding = 1000,
    window_size = 100,
    slide_by = 100,
    projection = "hg38"
) %>%
    column_to_rownames("sample_id")


bcl2_regions <- panel_regions %>%
    filter(name == "BCL2") %>%
    mutate(name = ifelse(start > 63300000, "BCL2_TSS", "BCL2_exon2"))

pdf("results/shm_heatmaps/BCL2_shm_heatmap.pdf", height = 8, width = 8)
heatmap_mutation_frequency_bin(
    regions_df = bcl2_regions,
    these_samples_metadata = sv_data_cat[sv_data_cat$ICC_class %in% c("FL", "DLBCL", "HGBCL-DH-BCL2"), ],
    maf = maf_data,
    metadataColumns = c("ICC_class", "bcl2_partner_group", "bcl2_ba"),
    sortByColumns = c("ICC_class", "bcl2_partner_group", "bcl2_ba"),
    return_heatmap_obj = FALSE,
    slide_by = 250,
    window_size = 250,
    min_bin_recurrence = 100,
    projection = "hg38"
)
dev.off()

igh_discrete <- read_tsv("data/region_data/genes_ig_myc_zoomed.hg38.tsv") %>%
    filter(chrom == "chr14") %>%
    mutate(end = case_when(
        name == "M" ~ 105860400,
        name == "Emu" ~ end,
        TRUE ~ end + 5000
    )) %>%
    mutate(start = ifelse(start == "Emu", 105860401, start)) %>%
    mutate(name = factor(name, levels = name))

sv_data_dlbcl_coo <- sv_data_cat %>%
    filter(!(ICC_class == "DLBCL" & !dlbcl_call %in% c("ABC", "GCB"))) %>%
    mutate(ICC_class = ifelse(ICC_class == "DLBCL", glue::glue("{dlbcl_call}-DLBCL"), as.character(ICC_class))) %>%
    mutate(ICC_class = factor(ICC_class, levels = c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL"))) %>%
    filter(!is.na(ICC_class))

pdf("results/shm_heatmaps/IGHC_shm_heatmap.pdf", height = 5, width = 7)
heatmap_mutation_frequency_bin(
    regions_df = igh_discrete,
    these_samples_metadata = rename(
        sv_data_dlbcl_coo,
        "ICC Class" = "ICC_class",
        "MYC Partner" = "myc_partner_group"
    ),
    maf = maf_data,
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

bl_vs_abc <- sv_data_cat %>%
    filter(
        ICC_class == "BL" |
            (ICC_class == "DLBCL" & dlbcl_call %in% c("ABC", "GCB"))
    ) %>%
    mutate(Group = ifelse(ICC_class == "DLBCL", dlbcl_call, "BL")) %>%
    mutate(Group = factor(Group, levels = c("BL", "ABC", "GCB")))

pdf("results/shm_heatmaps/IGHC_shm_heatmap_BL_ABC.pdf", height = 5, width = 7)
heatmap_mutation_frequency_bin(
    regions_df = igh_discrete,
    these_samples_metadata = bl_vs_abc,
    maf = maf_data,
    metadataColumns = c("Group"),
    sortByColumns = c("Group"),
    return_heatmap_obj = FALSE,
    slide_by = 500,
    window_size = 500,
    region_padding = 0,
    min_bin_recurrence = 10,
    label_regions_by = "name",
    projection = "hg38"
)
dev.off()

maf_igh <- maf_data %>%
    select(-name) %>%
    inner_join(igh_discrete, by = c("Chromosome" = "chrom")) %>%
    filter(Start_Position >= start, End_Position <= end) %>%
    left_join(sv_data_cat, by = c("Tumor_Sample_Barcode" = "sample_id")) %>%
    mutate(VAF = t_alt_count / t_depth) %>%
    left_join(read_tsv("data/ig_rearrangements/svar_master_evidence_of_csr.tsv"),
        by = c("Tumor_Sample_Barcode" = "sample_id")
    )

maf_igh %>%
    filter(!str_detect(name, "Emu|D|RR")) %>%
    filter(!ICC_class %in% c("HGBCL-DH-BCL6", NA)) %>%
    ggplot(aes(x = name, y = VAF, colour = ICC_class)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.8, size = 0.5) +
    geom_hline(yintercept = 0.65) +
    scale_colour_manual(values = colours$ICC_class[1:4]) +
    theme_cowplot() +
    facet_grid(CSR_event ~ ICC_class) +
    theme(legend.position = "none")

ggsave("results/shm_heatmaps/IGH_VAF_vs_C_gene.pdf", height = 10, width = 15)


# TODO: Add zero counts for all name groups
levels <- names(colours$ighc)
igh_mut_counts <- maf_igh %>%
    mutate(ICC_class = ifelse(
        ICC_class == "DLBCL",
        paste0(dlbcl_call, "-DLBCL"),
        as.character(ICC_class)
    )) %>%
    filter(!str_detect(name, "D|RR")) %>%
    filter(!ICC_class %in% c("HGBCL-DH-BCL6", NA)) %>%
    count(Tumor_Sample_Barcode, ICC_class, name) %>%
    full_join(select(
        sv_data_cat,
        Tumor_Sample_Barcode = sample_id,
        ICC_class
    )) %>%
    pivot_wider(
        names_from = name,
        values_from = n,
        values_fill = 0
    ) %>%
    select(-`NA`) %>%
    pivot_longer(
        -c("Tumor_Sample_Barcode", "ICC_class"),
        names_to = "name",
        values_to = "n"
    ) %>%
    mutate(name = factor(name, levels = rev(levels(maf_igh$name)))) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")
        )
    ) %>%
    filter(ICC_class %in% c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL"))

igh_mut_test <- igh_mut_counts %>%
    filter(!ICC_class %in% c("HGBCL-DH-BCL6", NA)) %>%
    group_by(name) %>%
    rstatix::wilcox_test(n ~ ICC_class, ref.group = "HGBCL-DH-BCL2", p.adjust.method = "BH") %>%
    rstatix::add_y_position(scales = "free", step.increase = 0)

igh_mut_counts %>%
    filter(!ICC_class %in% c("HGBCL-DH-BCL6", NA)) %>%
    ggplot(aes(x = ICC_class, y = n)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.5, aes(colour = ICC_class)) +
    ggpubr::stat_pvalue_manual(igh_mut_test, hide.ns = TRUE, x = "group2") +
    scale_colour_manual(
        values = colours$group[c(1, 2, 3, 9, 8)],
        name = ""
    ) +
    theme_bw() +
    facet_wrap(~name, scales = "free_y") +
    theme(legend.position = "bottom") +
    xlab("") +
    ylab("Mutation Count per Tumour") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none"
    )

ggsave("results/shm_heatmaps/count_shm_ighc.pdf", height = 7, width = 7)

# Make mini heatmap
igh_mut_by_group <- igh_mut_counts %>%
    group_by(name, ICC_class) %>%
    summarize(
        mut_rate = median(n),
    ) %>%
    ungroup() %>%
    pivot_wider(
        names_from = name,
        values_from = mut_rate
    ) %>%
    column_to_rownames("ICC_class")


igh_mut_test_mat <- igh_mut_test %>%
    select(name, group2, p.adj.signif) %>%
    pivot_wider(
        names_from = name,
        values_from = p.adj.signif
    ) %>%
    add_row(
        group2 = "HGBCL-DH-BCL2"
    ) %>%
    mutate(across(-group2, ~ ifelse(!str_detect(.x, "[*]") | is.na(.x), "", .x))) %>%
    column_to_rownames("group2")

igh_mut_test_mat <- igh_mut_test_mat[rownames(igh_mut_by_group), ]

col_fun <- colorRamp2(
    breaks = c(0, 5, 40),
    col = c("white", "orange", "purple")
)

heatmap_ashm <- Heatmap(igh_mut_by_group,
    col = col_fun,
    name = " ",
    row_labels = rownames(igh_mut_by_group),
    column_labels = colnames(igh_mut_by_group),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%s", igh_mut_test_mat[i, j]), x, y, gp = gpar(fontsize = 10))
    },
    rect_gp = gpar(col = "black", lwd = 0.5),
    heatmap_legend_param = list(
        title = "Median Mutation Count",
        legend_height = unit(4, "cm"),
        direction = "horizontal"
    ),
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_side = "top",
    show_column_names = TRUE,
    row_title = NA
)

ht <- draw(heatmap_ashm, heatmap_legend_side = "bottom")

pdf("results/shm_heatmaps/ighc_shm_heatmap_collapsed.pdf", height = 2, width = 4)
ht
dev.off()

igh_counts_dh <- maf_igh %>%
    filter(!str_detect(name, "D|RR")) %>%
    filter(ICC_class == "HGBCL-DH-BCL2", myc_partner_group != "no break") %>%
    count(Tumor_Sample_Barcode, myc_partner_group, name) %>%
    left_join(select(
        sv_data_cat,
        Tumor_Sample_Barcode = sample_id,
        myc_partner_group
    )) %>%
    pivot_wider(
        names_from = name,
        values_from = n,
        values_fill = 0
    ) %>%
    # select(-`NA`) %>%
    pivot_longer(
        -c("Tumor_Sample_Barcode", "myc_partner_group"),
        names_to = "name",
        values_to = "n"
    ) %>%
    mutate(name = factor(name, levels = rev(levels(maf_igh$name)))) %>%
    mutate(myc_partner_group = factor(myc_partner_group, levels = names(colours$myc_partner_group[1:4]))) %>%
    filter(!is.na(myc_partner_group))

igh_counts_dh %>%
    ggplot(aes(x = myc_partner_group, y = n)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.5, aes(colour = myc_partner_group)) +
    # ggpubr::stat_pvalue_manual(igh_mut_test, hide.ns = TRUE, x = "group2") +
    scale_colour_manual(
        values = colours$myc_partner_group[1:4],
        name = ""
    ) +
    theme_bw() +
    facet_wrap(~name, scales = "free_y") +
    theme(legend.position = "bottom") +
    xlab("") +
    ylab("Mutation Count per Tumour") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("results/shm_heatmaps/count_shm_ighc_dh_bcl2.pdf", height = 8, width = 8)


igh_counts_dh_igh <- maf_igh %>%
    filter(!str_detect(name, "D|RR")) %>%
    filter(ICC_class == "HGBCL-DH-BCL2", myc_partner_group == "IGH") %>%
    count(Tumor_Sample_Barcode, myc_annotation_partner, name) %>%
    left_join(select(
        sv_data_cat,
        Tumor_Sample_Barcode = sample_id,
        myc_annotation_partner
    )) %>%
    pivot_wider(
        names_from = name,
        values_from = n,
        values_fill = 0
    ) %>%
    # select(-`NA`) %>%
    pivot_longer(
        -c("Tumor_Sample_Barcode", "myc_annotation_partner"),
        names_to = "name",
        values_to = "n"
    ) %>%
    mutate(name = factor(name, levels = rev(levels(maf_igh$name)))) %>%
    # mutate(myc_partner_group = factor(myc_partner_group, levels = names(colours$myc_partner_group[1:4]))) %>%
    mutate(myc_annotation_partner = factor(
        str_remove(myc_annotation_partner, "IGH"),
        levels = levels(igh_counts_dh$name)
    )) %>%
    filter(!is.na(myc_annotation_partner)) %>%
    group_by(Tumor_Sample_Barcode) %>%
    complete(
        nesting(myc_annotation_partner, name),
        fill = list(n = 0)
    ) %>%
    ungroup()

igh_counts_dh_igh %>%
    filter(!myc_annotation_partner %in% c("Emu", "D")) %>%
    ggplot(aes(x = myc_annotation_partner, y = n)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.5, aes(colour = name, alpha = myc_annotation_partner == name)) +
    # ggpubr::stat_pvalue_manual(igh_mut_test, hide.ns = TRUE, x = "group2") +
    # scale_colour_manual(
    #     values = colours$myc_partner_group[1:4],
    #     name = ""
    # ) +
    scale_alpha_manual(
        values = c("TRUE" = 1, "FALSE" = 0.4),
        guide = "none"
    ) +
    theme_bw() +
    facet_wrap(~name, scales = "fixed", nrow = 3) +
    theme(legend.position = "bottom") +
    xlab("MYC Partner") +
    ylab("Mutation Count per Tumour")

ggsave("results/shm_heatmaps/count_shm_ighc_dh_bcl2_igh.pdf", height = 8, width = 8)
