source("src/libs.R")
library(tidyverse)
library(ComplexHeatmap)
library(GAMBLR)
library(data.table)
library(cowplot)
library(circlize)

source("src/circos/prepare_circos.R")
source("src/circos/plot_circos.R")

# Set a colour scheme for plotting breakpoints
colours_links <- c(
    "MYC" = "grey",
    "BCL2" = "grey",
    "IGK" = "#2E86AB",
    "IGL" = "#FC5130",
    "IGH" = "#87BC14",
    "PAX5" = "#C28CAE",
    "RFTN1" = "#52154E",
    "BCL6" = "#E40066",
    "IRAG2" = "#03CEA4",
    "SOCS1" = "#19535F"
)


##### Define regions for plotting #####
# Load a bed file of genes of interest
gene_bed_raw <- read_tsv("data/region_data/genes_hg38.bed")
gene_bed_raw$colour <- colours_links[gene_bed_raw$region]
gene_bed <- gene_bed_raw[!gene_bed_raw$region %in% c("BCL2", "SOCS1"), ]

# Collapse IG genes for zoomed out plots
gene_bed_large <- gene_bed %>%
    mutate(name = case_when(
        name == "IGHD" ~ "IGHC",
        str_detect(name, "IGH[AGDEM]") ~ "IGHC",
        str_detect(name, "IGH[VDJ]") ~ "IGHV",
        str_detect(name, "IG[KL]") ~ str_sub(name, 1, 3),
        TRUE ~ name
    )) %>%
    group_by(name) %>%
    mutate(
        start = min(start),
        end = max(end)
    ) %>%
    ungroup() %>%
    distinct() %>%
    filter(!str_detect(name, "Emu|RR"))

gene_bed_ig <- gene_bed_large %>%
    filter(region %in% c("MYC", "IGH", "IGK", "IGL"))

# Load a file of genes for zoomed IGH plots
colours_links_zoomed <- c(
    "MYC" = "grey",
    "CASC8" = "grey",
    "PVT1" = "grey",
    "Variable" = "#595959",
    "Emu" = "#c41230",
    "M" = "#F9BD1F",
    "D" = "#B76D29",
    "G3" = "#E90C8B",
    "G1" = "#2CACE3",
    "A1" = "#046852",
    "3'RR1" = "#c41230",
    "G2" = "#5C266C",
    "G4" = "#39B54B",
    "E" = "#FE9003",
    "A2" = "#115284",
    "3'RR2" = "#c41230"
)

gene_bed_igh_zoomed <- read_tsv("data/region_data/genes_ig_myc_zoomed.hg38.tsv")
gene_bed_igh_zoomed$colour <- colours_links_zoomed[
    gene_bed_igh_zoomed$name
]
gene_bed_igh_zoomed$name <- ifelse(
    str_detect(gene_bed_igh_zoomed$name, "RR|Emu"),
    "",
    gene_bed_igh_zoomed$name
)

regions_list <- list()
regions_list$all <- gene_bed_large
regions_list$ig <- gene_bed_ig
regions_list$igh <- gene_bed_igh_zoomed
regions_list$BCL2 <- gene_bed_large %>%
    filter(region %in% c("IGH", "IGK", "IGL")) %>%
    bind_rows(filter(gene_bed_raw, region == "BCL2")) %>%
    mutate(end = ifelse(name == "BCL2", 63128759, end)) %>%
    add_row(chrom = "chr18", region = "BCL2", start = 63318082, end = 63318666, name = "BCL2", colour = "grey")
    
##### Load metadata and SV data #####

# Load metadata. Select one representative sample.
md <- read_tsv("data/metadata/breakpoint_capture_best.tsv") 

# Load SV data
myc_data_raw <- read_tsv("data/sv_data/myc_annotated_breaks_unique.tsv")
myc_data <- myc_data_raw %>%
    select(-target) %>%
    rename(target = annotation_target) %>%
    filter(!is.na(start_target)) %>% 
    left_join(select(md, patient_id, biopsy_id, sample_id))
    
sv_data_cat <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv")%>% 
    left_join(select(md, patient_id, biopsy_id, sample_id))

# Load MAF data
maf_data <- read_tsv("data/maf/genome_capture.WRCY.hg38.maf")

panel_regions <- fread(
    "data/region_data/capture_TE99028370--hg38.bed",
    col.names = c("Chromosome", "Start_Position", "End_Position", "Region")
)

setkey(panel_regions, Chromosome, Start_Position, End_Position)

maf_in_region <- foverlaps(data.table(maf_data), panel_regions) %>%
    filter(!is.na(Region)) %>%
    mutate(
        Start_Position = i.Start_Position,
        End_Position = i.End_Position
    ) %>%
    select(-matches("^i[.]"))


# Define sample subsets

ICC_subsets <- split(md$sample_id, md$ICC_class)
ICC_subsets$`HGBCL-DH-BCL2_High_Grade` <- md %>%
    filter(ICC_class == "HGBCL-DH-BCL2" & morphology %in% c("HGBL", "High grade") & cpr_complete) %>%
    pull(sample_id)
ICC_subsets$`HGBCL-DH-BCL2_DLBCL` <- md %>%
    filter(ICC_class == "HGBCL-DH-BCL2" & morphology == "DLBCL" & cpr_complete) %>%
    pull(sample_id)
BL_EBV_subsets <- split(md$sample_id, md$BL_EBV)
BL_age_subsets <- split(md$sample_id, md$BL_agegrp)
BL_subsets <- append(BL_EBV_subsets, BL_age_subsets)
ICC_subsets <- append(ICC_subsets, BL_subsets)

circos_by_ICC <- function(
    subset,
    zoom_level,
    dim = 6,
    angle = 0,
    target_downstream_buffer = 5e5,
    target_upstream_buffer = 0,
    partner_downstream_buffer = 1e5,
    partner_upstream_buffer = 1e5,
    prop_chroms,
    prop_regions,
    region_order = NULL,
    colour_links_by = "partner",
    colours = colours_links) {
    circos <- prepare_circos(
        breakpoints = sv_data,
        genes = regions_list[[zoom_level]],
        sample_ids = ICC_subsets[[subset]],
        maf = maf_in_region,
        target_downstream_buffer = target_downstream_buffer,
        target_upstream_buffer = target_upstream_buffer,
        partner_downstream_buffer = partner_downstream_buffer,
        partner_upstream_buffer = partner_upstream_buffer,
        region_order = region_order
    )
    
    saveRDS(circos, paste0("results/sv_circos/", subset, ".", zoom_level, ".RDS"))

    pdf(paste0("results/sv_circos/", subset, ".", zoom_level, ".pdf"),
        height = dim, width = dim
    )
    plot_circos(
        circos,
        angle = angle,
        prop_chroms = prop_chroms,
        prop_regions = prop_regions,
        nest_colours = colours,
        link_colours = colours,
        colour_links_by = colour_links_by
    )
    dev.off()

}

##### Plot all recurrent MYC partners #####

sv_data <- myc_data

prop_chroms_all <- select(
    regions_list$all,
    chrom,
    region
) %>%
    distinct() %>%
    count(chrom) %>%
    arrange(as.numeric(str_remove(chrom, "chr"))) %>%
    mutate(n = ifelse(chrom == "chr8", 2, n)) %>%
    pull(n)
prop_regions_all <- c(1, 1, 1, 2, 1, 1, 1, 1)

lapply(
    names(ICC_subsets),
    circos_by_ICC,
    zoom_level = "all",
    dim = 6,
    angle = -117,
    target_downstream_buffer = 2e5,
    prop_chroms = prop_chroms_all,
    prop_regions = prop_regions_all
)

# Barplots of partner counts IG/non-IG
make_ig_barplot <- function(class) {
    plot <- sv_data_cat %>%
        filter(str_detect(sample_id, paste0(ICC_subsets[[class]], collapse = "|"))) %>%
        filter(!is.na(myc_partner)) %>% 
        mutate(partner_upper = ifelse(!str_detect(myc_partner, "IG[HKL]"), "Non-IG", str_remove(myc_partner, ".*--"))) %>%
        mutate(partner_upper = factor(partner_upper,
            levels = c("IGH", "IGK", "IGL", "Non-IG")
        )) %>%
        count(partner_upper, .drop = FALSE) %>%
        mutate(Percent = n / sum(n) * 100) %>%
        ggplot(aes(x = partner_upper, y = Percent, fill = partner_upper)) +
        geom_col() +
        scale_fill_manual(values = c(colours_links[c("IGH", "IGK", "IGL")], "Non-IG" = "darkgrey")) +
        xlab("") +
        ggtitle(class) +
        ylim(0, 85) +
        theme_cowplot() +
        theme(legend.position = "none")
    ggsave(paste0("results/sv_circos/", class, "_IG_barplot.pdf"), plot, height = 3, width = 3)
    saveRDS(plot, paste0("results/sv_circos/", class, "_IG_barplot.RDS"))
}

lapply(names(ICC_subsets), make_ig_barplot)

make_ig_barplot_mod <- function(class) {
    plot <- sv_data_cat %>%
        filter(!is.na(myc_partner)) %>% 
        filter(str_detect(sample_id, paste0(ICC_subsets[[class]], collapse = "|"))) %>%
        mutate(partner_upper = ifelse(!str_detect(myc_partner, "IG[HKL]|BCL6"), "Other\nNon-IG", str_remove(myc_partner, ".*--"))) %>%
        mutate(partner_upper = factor(partner_upper,
            levels = c("IGH", "IGK", "IGL", "BCL6", "Other\nNon-IG")
        )) %>%
        group_by(ICC_class) %>%
        count(partner_upper, .drop = FALSE) %>%
        mutate(Percent = n / sum(n) * 100) %>%
        ggplot(aes(x = partner_upper, y = Percent, fill = partner_upper)) +
        geom_col() +
        scale_fill_manual(values = c(colours_links[c("IGH", "IGK", "IGL", "BCL6")], "Other\nNon-IG" = "darkgrey")) +
        xlab("") +
        ggtitle(class) +
        ylim(0, 85) +
        theme_cowplot() +
        theme(legend.position = "none")
    ggsave(paste0("results/sv_circos/", class, "_IG_barplot_BCL6.pdf"), plot, height = 3, width = 3)
    saveRDS(plot, paste0("results/sv_circos/", class, "_IG_barplot_BCL6.RDS"))
}

lapply("HGBCL-DH-BCL6", make_ig_barplot_mod)

##### Focus on IG partners only #####

prop_chroms_ig <- c(2, 1, 1, 2)
prop_regions_ig <- prop_chroms_ig
region_order_ig <- c("chr8", "chr2", "chr22", "chr14")

lapply(
    names(ICC_subsets),
    # "HGBCL-DH-BCL2",
    circos_by_ICC,
    zoom_level = "ig",
    dim = 4,
    angle = 145,
    target_downstream_buffer = 2e5,
    prop_chroms = prop_chroms_ig,
    prop_regions = prop_regions_ig,
    region_order = region_order_ig
)

##### Focus on IGH partners only #####

prop_chroms <- c(1, 1)
prop_regions <- prop_chroms
sv_data <- mutate(
    sv_data,
    annotation_partner = str_remove(annotation_partner, "IGH")
)

lapply(
    names(ICC_subsets),
    # "HGBCL-DH-BCL2",
    circos_by_ICC,
    zoom_level = "igh",
    dim = 4,
    angle = -180,
    target_downstream_buffer = 1e5,
    partner_upstream_buffer = 1e4,
    partner_downstream_buffer = 1e4,
    prop_chroms = c(1, 1),
    prop_regions = c(1, 1),
    colour_links_by = "annotation_partner",
    colours = colours_links_zoomed
)

# Barplots of partner counts IGH regions

igh_levels <- names(colours_links_zoomed)[
    !names(colours_links_zoomed) %in% c("MYC", "CASC8", "PVT1")
]
igh_labels = str_replace(igh_levels, "Variable", "V")
igh_labels <- str_replace(igh_labels, "3'RR[12]", "RR")
make_igh_barplot <- function(class) {
    plot <- sv_data_cat %>%
        filter(!is.na(myc_partner)) %>% 
        filter(str_detect(sample_id, paste0(ICC_subsets[[class]], collapse = "|"))) %>%
        filter(myc_partner == "IGH", !str_detect(myc_annotation_partner, "RR")) %>%
        mutate(partner_upper = case_when(
            str_detect(myc_annotation_partner, "[VDJ]") ~ "Variable",
            TRUE ~ str_remove(myc_annotation_partner, "IGH")
        )) %>% 
        mutate(partner_upper = factor(
            partner_upper,
            levels = igh_levels,
            labels = igh_labels
        )) %>% 
        filter(!is.na(partner_upper)) %>%
        group_by(ICC_class) %>%
        count(partner_upper, .drop = FALSE) %>%
        mutate(Percent = n / sum(n) * 100) %>%
        ggplot(aes(x = partner_upper, y = Percent, fill = partner_upper)) +
        geom_col() +
        scale_fill_manual(values = c(colours_links_zoomed)) +
        xlab("") +
        ggtitle(class) +
        ylim(0, 55) +
        theme_cowplot() +
        theme(legend.position = "none")
    ggsave(paste0("results/sv_circos/", class, "_IGH_barplot.pdf"), plot, height = 3.5, width = 4)
    saveRDS(plot, paste0("results/sv_circos/", class, "_IGH_barplot.RDS"))
}

lapply(names(ICC_subsets), make_igh_barplot)

# Binary constant region vs Eu/Su/IGHVDJ

make_igh_barplot_binary <- function(class) {
    plot <- sv_data_cat %>%
        filter(!is.na(myc_partner)) %>%
        filter(str_detect(sample_id, paste0(ICC_subsets[[class]], collapse = "|"))) %>%
        filter(myc_partner == "IGH") %>%
        mutate(partner_upper = case_when(
            str_detect(myc_annotation_partner, "[VDJMm]") ~ "VDJ/Emu/IGHM",
            TRUE ~ "IGHC"
        )) %>%
        mutate(partner_upper = factor(
            partner_upper,
            # labels = c("VDJ/\nEmu/\nIGHM", "IGHC"),
            levels = c("VDJ/Emu/IGHM", "IGHC")
        )) %>%
        filter(!is.na(partner_upper)) %>%
        group_by(ICC_class) %>%
        count(partner_upper, .drop = FALSE) %>%
        mutate(Percent = n / sum(n) * 100) %>% 
        ggplot(aes(x = partner_upper, y = Percent, fill = partner_upper)) +
        geom_col() +
        scale_fill_manual(
            values = c("VDJ/Emu/IGHM" = "#c41230", "IGHC" = "#115284")
        )+
        scale_x_discrete(
            labels = c("VDJ\nEmu\nIGHM", "IGHC")
        ) +
        xlab("") +
        ylim(0,100) +
        theme_cowplot() +
        theme(legend.position = "none")
    ggsave(paste0("results/sv_circos/", class, "_IGH_barplot_binary.pdf"), plot, height = 5, width = 3)
    saveRDS(plot, paste0("results/sv_circos/", class, "_IGH_barplot_binary.RDS"))
}

lapply(names(ICC_subsets), make_igh_barplot_binary)


# Quantify HGBCL-DH-BCL2 morphology differences
sv_data %>%
    left_join(select(md, patient_id, biopsy_id, morphology, cpr_complete)) %>%
    mutate(morphology = ifelse(morphology == "High grade", "HGBL", morphology)) %>%
    filter(ICC_class == "HGBCL-DH-BCL2", morphology %in% c("DLBCL", "HGBL"), cpr_complete) %>%
    mutate(partner_ig = ifelse(str_detect(partner, "IG[HKL]"), TRUE, FALSE)) %>%
    summarize(table = list(table(morphology, partner_ig))) %>%
    mutate(fish = map(table, fisher.test), test = map(fish, broom::tidy)) %>%
    unnest(test) %>%
    select(-table, -fish) %>%
    write_tsv("results/sv_circos/fish_test_dh_hg_vs_dlbcl_ig.tsv")

for_table <- sv_data %>%
    left_join(select(md, patient_id, biopsy_id, morphology, cpr_complete)) %>%
    mutate(morphology = ifelse(morphology == "High grade", "HGBL", morphology)) %>%
    filter(ICC_class == "HGBCL-DH-BCL2", morphology %in% c("DLBCL", "HGBL"), cpr_complete) %>%
    mutate(partner_ig = ifelse(str_detect(partner, "IG[HKL]"), TRUE, FALSE))

table(morphology = for_table$morphology, IG_partner = for_table$partner_ig)

##### BCL2 IG partners only #####

bcl2_sv_data <- read_tsv("data/sv_data/bcl2_annotated_breaks_unique.tsv") %>%
    filter(chrom_target == "chr18") %>%
    mutate(target = "BCL2") %>% 
    left_join(select(md, patient_id, biopsy_id, sample_id))


sv_data <- bcl2_sv_data

lapply(
    names(ICC_subsets)[3:5],
    # "HGBCL-DH-BCL2",
    circos_by_ICC,
    zoom_level = "BCL2",
    dim = 4,
    angle = 145,
    target_downstream_buffer = 5e4,
    target_upstream_buffer = 5e4,
    prop_chroms = c(2, 1, 1, 2),
    prop_regions = c(2, 1, 1, 2),
    region_order = c("chr18", "chr2", "chr22", "chr14")
)

##### BCL6 rearrangements #####

bcl6_sv_data <- read_tsv("data/sv_data/bcl6_annotated_breaks_unique.tsv") %>%
    filter(chrom_target == "chr3") %>%
    mutate(target = "BCL6") %>%
    left_join(select(md, patient_id, biopsy_id, sample_id))

sv_data <- bcl6_sv_data

# ST6GAL1 chr3:186,930,526-187,078,553
# EIF4A2 chr3:186,783,577-186,789,897
# CIITA chr16:10,877,202-10,936,394
# MBNL1 chr3:152,268,929-152,465,780

gene_bed_bcl6 <- gene_bed_large %>%
    mutate(colour = case_when(
        region == "MYC" ~ "#E40066",
        region == "BCL6" ~ "grey",
        TRUE ~ colour
    )) %>% 
    filter(!region %in% c("IRAG2", "RFTN1")) %>% 
    add_row(
        chrom = c("chr3", "chr3", "chr16", "chr3"), 
        start = c(186930523, 186783577, 10877202, 152268929), 
        end = c(187078553, 186783577, 10936394, 152465780), 
        name = c("", "EIF4A2", "CIITA", "MBNL1"), 
        region = c("BCL6", "BCL6", "CIITA", "MBNL1"), 
        colour = c("grey", "grey", "#03CEA4", "#52154E")
    )

regions_list$BCL6 <- gene_bed_bcl6 %>% 
    mutate(name = ifelse(name %in% c("PVT1", "IGHC", "ZCCHC7"), "", name)) %>% 
    mutate(name = ifelse(name == "IGHV", "IGH", name))
    
colours_bcl6 <- gene_bed_bcl6$colour
names(colours_bcl6) <- gene_bed_bcl6$region
colours_bcl6 <- colours_bcl6[unique(names(colours_bcl6))]

prop_chroms_bcl6 <- select(
    regions_list$BCL6,
    chrom,
    region
) %>%
    distinct() %>%
    count(chrom) %>%
    arrange(as.numeric(str_remove(chrom, "chr"))) %>%
    mutate(n = ifelse(chrom == "chr8", 2, n)) %>%
    mutate(n = ifelse(chrom == "chr3", 3, n)) %>%
    pull(n)
prop_regions_bcl6 <- c(1, 1, 2, 2, 1, 1, 1, 1)

sv_data <- bcl6_sv_data %>% 
    mutate(partner = case_when(
        partner %in% c("LPP", "ST6GAL1", "EIF4A2") ~ "BCL6", 
        TRUE ~ partner
    ))

lapply(
    names(ICC_subsets),
    circos_by_ICC,
    zoom_level = "BCL6",
    dim = 5,
    angle = 190,
    target_downstream_buffer = 5e5,
    target_upstream_buffer = 5e5,
    prop_chroms = prop_chroms_bcl6,
    prop_regions = prop_regions_bcl6,
    colours = colours_bcl6
)


make_bcl6_barplot <- function(class) {
    plot <- sv_data_cat %>%
        filter(!is.na(bcl6_partner)) %>%
        filter(str_detect(sample_id, paste0(ICC_subsets[[class]], collapse = "|"))) %>%
        mutate(partner_upper = case_when(
            str_detect(bcl6_partner, "IG[HKL]") ~ str_remove(bcl6_partner, ".*--"),
            bcl6_partner == "MYC" ~ "MYC",
            !is.na(bcl6_partner) ~ "Other\nnon-IG"
        )) %>%
        mutate(partner_upper = factor(partner_upper,
            levels = c("IGH", "IGK", "IGL", "MYC", "Other\nnon-IG")
        )) %>%
        count(partner_upper, .drop = FALSE) %>%
        mutate(Percent = n / sum(n) * 100) %>%
        ggplot(aes(x = partner_upper, y = Percent, fill = partner_upper)) +
        geom_col() +
        scale_fill_manual(values = c(colours_bcl6[c("IGH", "IGK", "IGL", "MYC")], "Other\nnon-IG" = "darkgrey")) +
        xlab("BCL6 Partner") +
        ggtitle(class) +
        ylim(0, 85) +
        theme_cowplot() +
        theme(legend.position = "none")
    ggsave(paste0("results/sv_circos/", class, "_BCL6_barplot.pdf"), plot, height = 3, width = 3)
    saveRDS(plot, paste0("results/sv_circos/", class, "_BCL6_barplot.RDS"))
}

lapply(names(ICC_subsets), make_bcl6_barplot)

# Zoomed BCL6-IGH plots

regions_list$bcl6_zoomed <- regions_list$igh %>%
    filter(region == "IGH") %>%
    bind_rows(filter(regions_list$BCL6, region == "BCL6"))
    
prop_chroms_bcl6_zoomed <- c(1, 1)
prop_regions_bcl6_zoomed <- c(1, 1)

colours_bcl6_zoomed <- c(
    colours_links_zoomed[4:16]
)

sv_data <- mutate(
    sv_data,
    annotation_partner = str_remove(annotation_partner, "IGH")
)

lapply(
    names(ICC_subsets),
    circos_by_ICC,
    zoom_level = "bcl6_zoomed",
    dim = 5,
    angle = 190,
    target_downstream_buffer = 5e5,
    target_upstream_buffer = 5e5,
    prop_chroms = prop_chroms_bcl6_zoomed,
    prop_regions = prop_regions_bcl6_zoomed,
    colour_links_by = "annotation_partner",
    colours = colours_bcl6_zoomed
)

make_igh_barplot_bcl6 <- function(class) {
    plot <- sv_data_cat %>%
        filter(str_detect(sample_id, paste0(ICC_subsets[[class]], collapse = "|"))) %>%
        filter(bcl6_partner == "IGH", !str_detect(bcl6_annotation_partner, "RR")) %>%
        mutate(partner_upper = case_when(
            str_detect(bcl6_annotation_partner, "[VDJ]") ~ "Variable",
            TRUE ~ str_remove(bcl6_annotation_partner, "IGH")
        )) %>% 
        mutate(partner_upper = factor(
            partner_upper,
            levels = igh_levels,
            labels = igh_labels
        )) %>%
        filter(!is.na(partner_upper)) %>%
        group_by(ICC_class) %>%
        count(partner_upper, .drop = FALSE) %>%
        mutate(Percent = n / sum(n) * 100) %>%
        ggplot(aes(x = partner_upper, y = Percent, fill = partner_upper)) +
        geom_col() +
        scale_fill_manual(values = c(colours_links_zoomed)) +
        xlab("") +
        ggtitle(class) +
        ylim(0, 80) +
        theme_cowplot() +
        theme(legend.position = "none")
    ggsave(paste0("results/sv_circos/", class, "_BCL6_IGH_barplot.pdf"), plot, height = 3.5, width = 4)
    saveRDS(plot, paste0("results/sv_circos/", class, "_BCL6_IGH_barplot.RDS"))
}

lapply(names(ICC_subsets), make_igh_barplot_bcl6)
