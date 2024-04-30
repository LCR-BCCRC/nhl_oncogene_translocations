source("src/libs.R")

library(tidyverse)
library(ComplexHeatmap)
library(GAMBLR)
library(data.table)
library(cowplot)
library(lemon)
library(ggrepel)

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
    left_join(select(md, patient_id, biopsy_id, sample_id)) %>% 
    mutate(myc_partner_group = ifelse(
        ICC_class == "BL" & is.na(myc_partner_group), 
        "not found", myc_partner_group
    ))

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

myc_labels <- data.frame(
    x1 = 127735434, 
    x2 = 127741428, 
    y1 = c(50, 2.5, 15, 1, 20, 6), 
    y2 = c(50, 2.5, 15, 1, 20, 6), 
    labels = c("MYC", "MYC", rep.int(NA, 4)), 
    ICC_class = c("BL", "BL", "DLBCL", "DLBCL", "HGBCL-DH-BCL2", "HGBCL-DH-BCL2"), 
    myc_partner_igh_binary = rep(c("IGH", "non-IGH"), 3)
)

myc_data %>% 
    select(-matches("ba|status")) %>% 
    filter(ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2")) %>% 
    left_join(sv_data_cat) %>% 
    filter(!is.na(myc_partner_group)) %>% 
    ggplot(aes(
        x = start_target, 
        fill = myc_partner_group
    )) + 
    geom_histogram(binwidth = 5000) +
    geom_segment(
        data = myc_labels,
        aes(
            x = x1-1000,
            xend = x2+1000, 
            y = y1,
            yend = y2
        ),
        inherit.aes = FALSE,
        size = 1.5
    ) +
    ggrepel::geom_text_repel(
        data = myc_labels, 
        aes(
            x = x1, 
            y = y1, 
            label = labels
        ), 
        inherit.aes = FALSE
    ) +
    facet_rep_wrap(myc_partner_igh_binary ~ ICC_class, scales = "free_y", ncol = 3) + 
    scale_fill_manual(
        values = colours$myc_partner_group[1:4], 
        name = "MYC Partner"
    ) +
    scale_x_continuous(
        breaks = seq(127600000, 128400000, by = 200000), 
        labels = c("127.6MB", "127.8MB", "128.0MB", "128.2MB", "128.4MB"), 
        limits = c(127600000, 128400000)
    ) +
    xlab("Position (chr8)") + ylab("")+ 
    theme_bw() + 
    theme(
        legend.position = "bottom"
    )
    
ggsave("results/myc_svs/myc_break_by_partner_hist.pdf", height = 6, width = 12)


myc_data %>% 
    select(-matches("ba|status")) %>% 
    filter(ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2")) %>% 
    left_join(sv_data_cat) %>% 
    filter(partner == "IGH") %>%
    mutate(IGH_partner_mechanism = ifelse(start_partner > 105861000, "SHM", "CSR")) %>%  
    ggplot(aes(
        x = start_target, 
        fill = IGH_partner_mechanism
    )) + 
    geom_histogram(binwidth = 100) +
    # geom_density() +
    geom_segment(
        data = filter(myc_labels, myc_partner_igh_binary == "IGH"),
        aes(
            x = x1,
            xend = x2, 
            y = 0,
            yend = 0
        ),
        inherit.aes = FALSE,
        size = 1.5
    ) +
    facet_rep_wrap( ~ ICC_class, scales = "free_y", ncol = 1) +
    scale_fill_manual(
        values = c("SHM" = "#c41230", "CSR" = "#115284"), 
        name = "IGH Mechanism"
    ) +
    scale_x_continuous(
        breaks = seq(127730000, 127745000, by = 5000), 
        labels = c("127.3MB", "127.35MB", "127.4MB", "127.45MB"), 
        limits = c(127730000, 127745000)
    ) +
    xlab("Position (chr8)") + ylab("")+ 
    theme_bw() + 
    theme(
        legend.position = "bottom"
    )
    
ggsave("results/myc_svs/myc_break_by_IGH_mechanism.pdf", height = 6, width = 5)

myc_data %>% 
    left_join(select(
        sv_data_cat, 
        patient_id, 
        biopsy_id, 
        myc_partner_group, 
        BL_EBV
    )) %>% 
    mutate(group = case_when(
        ICC_class == "BL" & BL_EBV == "EBV-negative" ~ "BL EBV-", 
        ICC_class == "BL" & BL_EBV == "EBV-positive" ~ "BL EBV+",
        TRUE ~ ICC_class
    )) %>% 
    mutate(group = factor(group, levels = c(
        "BL EBV+", "BL EBV-", "DLBCL", "HGBCL-DH-BCL2"
    ))) %>% 
    mutate(myc_focal = ifelse(
        start_target > 127730000 & start_target < 127745000, 
        TRUE, 
        FALSE
    )) %>% 
    filter(!is.na(group), myc_partner_group == "IGH") %>% 
    count(group, myc_focal) %>% 
    group_by(group) %>% 
    mutate(percent = round(n/sum(n) * 100)) %>% 
    ggplot(aes(
        x = group, 
        y = n, 
        fill = group, 
        alpha = myc_focal
    )) + 
    geom_col() + 
    geom_label(
        aes(label = glue("{percent}%"), colour = group), 
        position = position_stack(vjust = 0.5), 
        show.legend = FALSE
    ) + 
    scale_fill_manual(
        values = colours$group_with_ebv
    ) + 
    scale_alpha_manual(
        values = c(
            "TRUE" = 1, "FALSE" = 0.7
        ), 
        name = "MYC Focal Break"
    ) + 
    scale_x_discrete(
        labels = str_replace(
            names(colours$group_with_ebv[c(2, 1, 5, 3)]), 
            "HGBCL-DH", 
            "HGBCL-\nDH"
        )
    ) +
    scale_colour_manual(
        values = c(
            "BL EBV-" = "black",
            "BL EBV+" = "black",
            "HGBCL-DH-BCL2" = "white",
            "DLBCL" = "black"
        )
    ) +
    theme_bw() + 
    theme(
        legend.position = "bottom"
    ) + 
    xlab("") + ylab("Count") +
    guides(fill = "none")
    
ggsave("results/myc_svs/myc_focal_break_by_group_with_ebv.pdf", height = 6, width = 4)


maf_in_myc <- maf_in_region %>% 
    filter(
        Chromosome == "chr8", 
        Start_Position >= 127735434, 
        End_Position <= 127741428
    ) %>% 
    count(Tumor_Sample_Barcode) %>% 
    full_join(distinct(select(maf_data, Tumor_Sample_Barcode))) %>% 
    mutate(n = replace_na(n, 0)) %>% 
    rename(
        sample_id = Tumor_Sample_Barcode, 
        count = n
    ) 

csr_shm <- myc_data %>%
    filter(partner == "IGH") %>% 
    mutate(IGH_partner_mechanism = ifelse(start_partner > 105861000, "SHM", "CSR")) %>% 
    select(
        patient_id, 
        biopsy_id, 
        IGH_partner_mechanism
    ) %>% 
    distinct()
colours_myc_with_csr <- c(
    "IGH-CSR" = unname(colours$myc_partner_group["IGH"]), 
    "IGH-SHM" = unname(colours$myc_partner_group["IGH"]), 
    colours$myc_partner_group[2:6]
)
colours_myc_with_csr["no break"] = "darkgrey"

myc_mut_counts <- sv_data_cat %>% 
    filter(ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2")) %>% 
    left_join(csr_shm) %>% 
    mutate(myc_partner_group = case_when(
        is.na(myc_partner_group) & myc_ba == "POS" ~ "not found", 
        is.na(myc_partner_group) & (myc_ba == "NEG") ~ "no break", 
        myc_partner == "TP63--IGH" ~ "IGH", 
        TRUE ~ myc_partner_group
    )) %>% 
    mutate(grp2 = factor(ifelse(
        myc_partner_group == "IGH", 
        glue("IGH-{IGH_partner_mechanism}"), 
        myc_partner_group
    ), levels = names(colours_myc_with_csr))) %>% 
    mutate(myc_partner_group = factor(
        myc_partner_group, levels = names(colours$myc_partner_group)
    )) %>% 
    filter(
        !is.na(myc_partner_group), 
        !is.na(grp2), 
        !(ICC_class == "BL" & myc_partner_group == "non-IG")
    ) %>% 
    left_join(maf_in_myc) 

counts_by_grp <- myc_mut_counts %>% 
    group_by(ICC_class, grp2) %>% 
    summarize(
        total_group = sum(!is.na(count)), 
        mutated = sum(count >= 1, na.rm = TRUE)
    ) %>% 
    mutate(percent_mut = round(mutated/total_group*100)) %>% 
    mutate(label = glue("N={mutated}\n{percent_mut}%"))

fish_test_ig_nonig <- myc_mut_counts %>% 
    filter(ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2")) %>% 
    mutate(myc_partner = case_when(
        myc_partner_group %in% c("IGH", "IGK/L") ~ "IG", 
        str_detect(myc_partner_group, "non-IG") ~ "non-IG" 
    )) %>% 
    filter(!is.na(myc_partner)) %>%
    mutate(mutated = ifelse(count >= 1, TRUE, FALSE))
    
fish_test_ig_nonig = broom::tidy(fisher.test(table(
    fish_test_ig_nonig$myc_partner, 
    fish_test_ig_nonig$mutated
)))

write_tsv(fish_test_ig_nonig, "results/myc_svs/myc_mutations_count_fish_test.tsv")

mut_test <- myc_mut_counts %>% 
    group_by(ICC_class) %>% 
    wilcox_test(count ~ grp2, ref.group = "IGH-CSR") %>% 
    filter(n2 > 4)

myc_mut_counts %>% 
    ggplot(aes(
        x = grp2, 
        y = count
    )) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_quasirandom(aes(
        colour = myc_partner_group
    )) + 
    geom_text(
        data = counts_by_grp, 
        aes(x = grp2, y = 57, label = label),
        inherit.aes = FALSE
    )+
    scale_colour_manual(
        values = c(
            colours$myc_partner_group[c(1:4, 6)], 
            "no break" = "darkgrey"
        ),
        name = "MYC Partner"
    ) + 
    stat_pvalue_manual(
        mut_test, 
        y = 43, 
        x = "group2"
    ) +
    facet_wrap(~ ICC_class, ncol = 1) + 
    xlab("") + ylab("MYC mutation count") +
    ylim(0,65) +
    theme_bw() + 
    theme(
        legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/myc_svs/myc_mutations_by_partner.pdf", height = 7, width = 7)

myc_mut_counts %>%
    mutate(mut_freq_per_100bp = count / 6) %>%
    group_by(ICC_class, grp2) %>%
    summarize(
        total_tumors = n(),
        min_mut_freq = min(mut_freq_per_100bp),
        max_mut_freq = max(mut_freq_per_100bp),
        avg_mut_freq = mean(mut_freq_per_100bp)
    ) %>% 
    write_tsv("results/myc_svs/myc_mut_freq_per_100bp.tsv")


sv_data_cat %>%
    filter(
        ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2"), 
        myc_partner == "IGH"
    ) %>% 
    mutate(myc_partner_igh_discrete = case_when(
        myc_partner_igh_discrete %in% c("Variable", "Emu", "IGHM") ~ "Variable/Emu/IGHM", 
        str_detect(myc_annotation_partner, "IGHJ") ~ "SHM", 
        !is.na(myc_partner_igh_discrete) ~ "IGHC"
    )) %>% 
    mutate(myc_partner_igh_discrete = factor(
        myc_partner_igh_discrete,
        levels = c("Variable/Emu/IGHM", "IGHC")
    )) %>% 
    filter(!is.na(myc_partner_igh_discrete)) %>%
    left_join(maf_in_myc) %>%
    ggplot(aes(
        x = myc_partner_igh_discrete,
        y = count
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(
        colour = myc_partner_igh_discrete
    )) +
    ggpubr::stat_compare_means(
        ref.group = "Variable/Emu/IGHM", 
        label = "p.signif", 
        label.y = 50
    ) +
    ggsci::scale_colour_d3(
        name = "MYC-IGH Partner"
    ) +
    facet_wrap(~ICC_class, ncol = 1) +
    xlab("") +
    ylab("MYC mutation count") +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave("results/myc_svs/myc_mutations_by_IGH_partner_binary.pdf", height = 6, width = 5)


# Write a maf of MYC mutations in DLBCL-MYC-neg and DLBCL-Ig and look for
# aSHM WRCY mutations
myc_no_break <- myc_mut_counts %>%
    filter(ICC_class == "DLBCL", grp2 == "no break") %>% 
    pull(sample_id)

myc_no_break <- maf_in_region %>%
    filter(
        Chromosome == "chr8",
        Start_Position >= 127735434,
        End_Position <= 127741428
    ) %>% 
    filter(Tumor_Sample_Barcode %in% myc_no_break) %>% 
    select(all_of(colnames(read_tsv("data/maf/region_mafs/FL_genome_capture.hg38.maf"))))
write_tsv(myc_no_break, "data/maf/region_mafs/DLBCL_MYC_no_break.hg38.maf")

myc_ig <- myc_mut_counts %>%
    filter(ICC_class == "DLBCL", str_detect(grp2, "IG[HKL]")) %>% 
    pull(sample_id)
    
myc_ig <- maf_in_region %>%
    filter(
        Chromosome == "chr8",
        Start_Position >= 127735434,
        End_Position <= 127741428
    ) %>% 
    filter(Tumor_Sample_Barcode %in% myc_ig)%>% 
    select(all_of(colnames(read_tsv("data/maf/region_mafs/FL_genome_capture.hg38.maf"))))
write_tsv(myc_ig, "data/maf/region_mafs/DLBCL_MYC_IG.hg38.maf")
