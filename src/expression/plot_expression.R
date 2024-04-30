source("src/libs.R")
library(tidyverse)
library(GAMBLR)
library(rstatix)
library(ggpubr)
library(glue)
library(readxl)
library(ggbeeswarm)
library(ComplexHeatmap)

md <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    filter(seq_type == "mrna") %>%
    mutate(dlbcl_call = ifelse(is.na(dlbcl_call), "ND", dlbcl_call)) %>%
    group_by(patient_id, biopsy_id) %>%
    slice_max(protocol == "Ribodepletion", n = 1, with_ties = FALSE) %>%
    ungroup()


sv_data <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv") %>%
    select(-matches("_ba$|_cn$|dlbcl_|dhitsig_"))

expression <- read_tsv("results/expression/myc_bcl2_bcl6_expr.tidy.tsv")

md_sv <- md %>%
    inner_join(sv_data) %>%
    mutate(myc_partner = case_when(
        myc_partner %in% c("BCL6", "PAX5") ~ "PAX5/ BCL6",
        str_detect(myc_partner, "IG[HKL]") ~ str_remove(myc_partner, ".*--"),
        !is.na(myc_partner) & ICC_class != "BL" ~ "Other non-IG",
        myc_bp_status == "FALSE_NEG" ~ "False Negative",
        TRUE ~ "None"
    )) %>%
    left_join(expression) %>%
    distinct() %>%
    filter(!(ICC_class == "BL" & myc_partner == "None"))

colours_partners <- c(
    "IGH" = "#87BC14",
    "IGK" = "#2E86AB",
    "IGL" = "#FC5130",
    "PAX5/ BCL6" = "#E40066",
    "Other non-IG" = "#C28CAE",
    "False Negative" = "#595959",
    "None" = "#BDBDBD"
)
median_bl_igh <- md_sv %>%
    filter(ICC_class == "BL" & myc_partner == "IGH" & hgnc_symbol == "MYC") %>%
    summarize(median = median(expression)) %>%
    pull(median)

median_dlbcl_mycneg <- md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(ICC_class == "DLBCL" & myc_partner == "None" & hgnc_symbol == "MYC") %>%
    summarize(median = median(expression)) %>%
    pull(median)

comp_means <- md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(
        ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol == "MYC"
    ) %>%
    mutate(partner_class = glue("{ICC_class}_{myc_partner}")) %>%
    rstatix::wilcox_test(expression ~ partner_class, p.adjust.method = "BH", conf.level = 0.9, ref.group = "DLBCL_None") %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "ns" & p.adj < 0.1, "*", p.adj.signif)) %>%
    separate(group2, into = c("ICC_class", "group2"), sep = "_")

write_tsv(comp_means, "results/expression/comp_means_myc.tsv")

comp_means_vs_IGH <- md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(
        ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol == "MYC"
    ) %>%
    mutate(myc_partner = factor(myc_partner, levels = names(colours_partners))) %>%
    group_by(ICC_class) %>%
    rstatix::wilcox_test(expression ~ myc_partner, p.adjust.method = "BH", conf.level = 0.9, ref.group = "IGH") %>%
    add_y_position() %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "ns" & p.adj < 0.1, "*", p.adj.signif))

md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(
        ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol == "MYC"
    ) %>%
    mutate(myc_partner = factor(myc_partner, levels = names(colours_partners))) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "DLBCL", "HGBCL-DH-BCL2")
        )
    ) %>%
    ggplot(aes(x = myc_partner, y = expression, colour = myc_partner)) +
    geom_hline(
        yintercept = median_bl_igh,
        colour = colours_partners["IGH"],
        lty = 2,
        lwd = 1
    ) +
    geom_hline(
        yintercept = median_dlbcl_mycneg,
        colour = colours_partners["None"],
        lty = 2,
        lwd = 1
    ) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.5) +
    ggpubr::stat_pvalue_manual(
        comp_means,
        x = "group2",
        y.position = 7,
        label = "p.adj.signif",
        vjust = 2
    ) +
    stat_pvalue_manual(
        comp_means_vs_IGH,
        label = "p.adj.signif",
        x = "group2",
        y.position = 15.75,
        vjust = 2
    ) +
    scale_colour_manual(values = colours_partners) +
    facet_grid(~ICC_class, scales = "free_x", space = "free") +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    ylim(6, 16) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/myc_expr_by_partner.pdf", height = 5, width = 7)



comp_means_dh_bcl6 <- md_sv %>%
    filter(
        ICC_class == "HGBCL-DH-BCL6" | (ICC_class == "DLBCL" & myc_partner == "None"),
        hgnc_symbol == "MYC"
    ) %>%
    mutate(myc_partner = factor(
        ifelse(myc_partner %in% c("PAX5/ BCL6", "Other non-IG"), "Non-IG", myc_partner),
        levels = c("IGH", "Non-IG", "False Negative", "None")
    )) %>%
    ungroup() %>%
    rstatix::wilcox_test(expression ~ myc_partner, p.adjust.method = "BH", conf.level = 0.9, ref.group = "None") %>%
    add_y_position() %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "ns" & p.adj < 0.1, "*", p.adj.signif)) %>%
    mutate(real_group2 = group2) %>%
    mutate(group2 = ifelse(group2 == "Non-IG", "Other non-IG", group2))

write_tsv(comp_means_dh_bcl6, "results/expression/comp_means_myc_DH-BCL6.tsv")

md_sv %>%
    filter(
        ICC_class == "HGBCL-DH-BCL6" | (ICC_class == "DLBCL" & myc_partner == "None"),
        hgnc_symbol == "MYC"
    ) %>%
    mutate(myc_partner = factor(myc_partner, levels = names(colours_partners))) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "DLBCL", "HGBCL-DH-BCL2")
        )
    ) %>%
    ggplot(aes(x = myc_partner, y = expression, colour = myc_partner)) +
    geom_hline(
        yintercept = median_bl_igh,
        colour = colours_partners["IGH"],
        lty = 2,
        lwd = 1
    ) +
    geom_hline(
        yintercept = median_dlbcl_mycneg,
        colour = colours_partners["None"],
        lty = 2,
        lwd = 1
    ) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.5) +
    ggpubr::stat_pvalue_manual(
        comp_means_dh_bcl6,
        x = "group2",
        y.position = 16,
        label = "p.adj.signif",
        vjust = 2
    ) +
    scale_colour_manual(values = colours_partners) +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/myc_expr_by_partner_dh_bcl6.pdf", height = 4, width = 4)


colours_partners2 <- colours_partners
names(colours_partners2)[7] <- "DLBCL MYC-R Neg"

md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(
        ICC_class == "HGBCL-DH-BCL2" |
            (ICC_class == "DLBCL" & myc_partner == "None"),
        hgnc_symbol == "MYC"
    ) %>%
    mutate(myc_partner = ifelse(ICC_class == "DLBCL", "DLBCL MYC-R Neg", myc_partner)) %>%
    mutate(myc_partner = factor(myc_partner, levels = names(colours_partners2))) %>%
    filter(!is.na(myc_partner)) %>%
    ggplot(aes(x = myc_partner, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 1, aes(colour = myc_partner)) +
    ggpubr::stat_compare_means(
        ref.group = "DLBCL MYC-R Neg",
        label.y = 16,
        label = "p.signif",
        vjust = 2
    ) +
    scale_colour_manual(values = colours_partners2) +
    cowplot::theme_cowplot() +
    xlab("") +
    ylab("Normalized Expression") +
    ylim(6, 16) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/myc_expr_by_partner_DH_vs_DLBCL_neg.pdf",
    height = 6, width = 5
)


comp_means_igh <- md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(
        ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol == "MYC",
        myc_partner == "IGH" | myc_bp_status == "TRUE_NEG"
    ) %>%
    mutate(igh_partner = case_when(
        str_detect(myc_annotation_partner, "[VDJ]") ~ "Variable",
        myc_annotation_partner == "Emu" ~ "Emu/IGHM",
        myc_bp_status == "TRUE_NEG" ~ "MYC-R Neg",
        TRUE ~ "CSR"
    )) %>%
    filter(!(igh_partner %in% c("Variable", "Emu/IGHM") & ICC_class == "HGBCL-DH-BCL2")) %>%
    mutate(partner_class = glue("{ICC_class}_{igh_partner}")) %>%
    rstatix::wilcox_test(expression ~ partner_class, p.adjust.method = "BH", conf.level = 0.9, ref.group = "DLBCL_MYC-R Neg") %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "ns" & p.adj < 0.1, "*", p.adj.signif)) %>%
    separate(group2, into = c("ICC_class", "group2"), sep = "_") %>%
    mutate(group2 = factor(group2, levels = c("Variable", "Emu/IGHM", "CSR", "MYC-R Neg")))


md_sv %>%
    filter(!(dlbcl_call == "ABC" & ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2"))) %>%
    filter(
        ICC_class %in% c("BL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol == "MYC",
        myc_partner == "IGH" | myc_bp_status == "TRUE_NEG"
    ) %>%
    mutate(igh_partner = case_when(
        str_detect(myc_annotation_partner, "[VDJ]") ~ "Variable",
        myc_annotation_partner == "Emu" ~ "Emu/IGHM",
        myc_bp_status == "TRUE_NEG" ~ "MYC-R Neg",
        TRUE ~ "CSR"
    )) %>%
    group_by(igh_partner, ICC_class) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(igh_partner = factor(igh_partner, levels = c("Variable", "Emu/IGHM", "CSR", "MYC-R Neg"))) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "HGBCL-DH-BCL2", "DLBCL")
        )
    ) %>%
    ggplot(aes(x = igh_partner, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 1, aes(colour = igh_partner)) +
    facet_grid(~ICC_class, scales = "free_x", space = "free") +
    ggpubr::stat_pvalue_manual(
        comp_means_igh,
        x = "group2",
        y.position = 15.75,
        label = "p.adj.signif",
        vjust = 2
    ) +
    ggsci::scale_color_d3() +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    ylim(6, 16) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/myc_expr_by_igh_partner.pdf", height = 4, width = 6)


comp_means_bcl2 <- md_sv %>%
    filter(
        ICC_class %in% c("FL", "HGBCL-DH-BCL2") |
            (ICC_class == "DLBCL" & dlbcl_call == "GCB"),
        hgnc_symbol == "BCL2"
    ) %>%
    mutate(bcl2_partner = case_when(
        bcl2_partner %in% c("IGH", "IGK", "IGL") |
            bcl2_ba == "POS" |
            ICC_class == "HGBCL-DH-BCL2" ~ "BCL2-R",
        TRUE ~ "None"
    )) %>%
    mutate(partner_class = glue("{ICC_class}_{bcl2_partner}")) %>%
    rstatix::wilcox_test(expression ~ partner_class, p.adjust.method = "BH", conf.level = 0.9, ref.group = "DLBCL_None") %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "ns" & p.adj < 0.1, "*", p.adj.signif)) %>%
    separate(group2, into = c("ICC_class", "group2"), sep = "_")

write_tsv(comp_means_bcl2, "results/expression/comp_means_bcl2.tsv")

colours_bcl2_partners <- c(
    "BCL2-R" = unname(colours_partners["IGH"]),
    "None" = unname(colours_partners["None"])
)

md_sv %>%
    filter(
        ICC_class %in% c("FL", "HGBCL-DH-BCL2") |
            (ICC_class == "DLBCL" & dlbcl_call == "GCB"),
        hgnc_symbol == "BCL2"
    ) %>%
    mutate(bcl2_partner = case_when(
        bcl2_partner %in% c("IGH", "IGK") |
            bcl2_ba == "POS" |
            ICC_class == "HGBCL-DH-BCL2" ~ "BCL2-R",
        TRUE ~ "None"
    )) %>%
    ggplot(aes(x = bcl2_partner, y = expression)) +
    # geom_hline(
    #     yintercept = median_bl_igh,
    #     colour = colours_bcl2_partners["IGH"],
    #     lty = 2,
    #     lwd = 1
    # ) +
    # geom_hline(
    #     yintercept = median_dlbcl_mycneg,
    #     colour = colours_partners["None"],
    #     lty = 2,
    #     lwd = 1
    # ) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 1, aes(colour = bcl2_partner)) +
    ggpubr::stat_pvalue_manual(
        comp_means_bcl2,
        x = "group2",
        y.position = 15.75,
        label = "p.adj.signif",
        vjust = 2
    ) +
    scale_colour_manual(values = colours_bcl2_partners) +
    facet_grid(~ICC_class, scales = "free_x", space = "free") +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    ylim(6, 16) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/bcl2_expr_by_igh_partner.pdf", height = 4, width = 5)

colours_bcl6 <- colours_partners[
    c("IGH", "PAX5/ BCL6", "Other non-IG", "False Negative", "None")
]
names(colours_bcl6) <- c(
    "IG",
    "MYC",
    "Other non-IG",
    "False Negative",
    "None"
)


comp_means_bcl6 <- md_sv %>%
    filter(
        ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2", "HGBCL-DH-BCL6"),
        hgnc_symbol == "BCL6"
    ) %>%
    mutate(bcl6_partner = case_when(
        bcl6_partner %in% c("IGH", "IGL", "IGK") ~ "IG",
        bcl6_partner == "MYC" ~ "MYC",
        !is.na(bcl6_partner) ~ "Other non-IG",
        bcl6_ba == "POS" & is.na(bcl6_partner) ~ "False Negative",
        TRUE ~ "None"
    )) %>%
    mutate(partner_class = glue("{ICC_class}_{bcl6_partner}")) %>%
    rstatix::wilcox_test(expression ~ partner_class, p.adjust.method = "BH", conf.level = 0.9, ref.group = "DLBCL_None") %>%
    mutate(p.adj.signif = ifelse(p.adj.signif == "ns" & p.adj < 0.1, "*", p.adj.signif)) %>%
    separate(group2, into = c("ICC_class", "group2"), sep = "_")
write_tsv(comp_means_bcl6, "results/expression/comp_means_bcl6.tsv")

md_sv %>%
    filter(
        ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2", "HGBCL-DH-BCL6"),
        hgnc_symbol == "BCL6"
    ) %>%
    mutate(bcl6_partner = case_when(
        bcl6_partner %in% c("IGH", "IGL", "IGK") ~ "IG",
        bcl6_partner == "MYC" ~ "MYC",
        !is.na(bcl6_partner) ~ "Other non-IG",
        bcl6_ba == "POS" & is.na(bcl6_partner) ~ "False Negative",
        TRUE ~ "None"
    )) %>%
    mutate(bcl6_partner = factor(
        bcl6_partner,
        levels = names(colours_bcl6)
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("DLBCL", "HGBCL-DH-BCL2", "HGBCL-DH-BCL6")
        )
    ) %>%
    ggplot(aes(x = bcl6_partner, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 1, aes(colour = bcl6_partner)) +
    scale_colour_manual(values = colours_bcl6) +
    stat_pvalue_manual(
        comp_means_bcl6,
        y.position = 14.5,
        x = "group2"
    ) +
    facet_grid(~ICC_class, scales = "free_x", space = "free") +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    ylim(6, 15) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/bcl6_expr_by_igh_partner.pdf", height = 4, width = 6)


md_sv %>%
    filter(
        ((ICC_class == "DLBCL") | (ICC_class == "HGBCL-DH-BCL2" & dlbcl_call == "GCB")),
        hgnc_symbol == "BCL6",
        dlbcl_call %in% c("ABC", "GCB")
    ) %>%
    mutate(bcl6_partner = case_when(
        bcl6_partner %in% c("IGH", "IGL", "IGK") ~ "IG",
        bcl6_partner == "MYC" ~ "MYC",
        !is.na(bcl6_partner) ~ "Other non-IG",
        bcl6_ba == "POS" & is.na(bcl6_partner) ~ "False Negative",
        TRUE ~ "None"
    )) %>%
    mutate(bcl6_partner = factor(
        bcl6_partner,
        levels = names(colours_bcl6)
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("DLBCL", "HGBCL-DH-BCL2", "HGBCL-DH-BCL6")
        )
    ) %>%
    ggplot(aes(x = bcl6_partner, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 1, aes(colour = bcl6_partner)) +
    scale_colour_manual(values = colours_bcl6) +
    facet_grid(dlbcl_call ~ ICC_class) +
    stat_compare_means(
        ref.group = "None",
        method = "wilcox",
        label = "p.signif",
        label.y = 15,
        hide.ns = TRUE,
        size = 5,
        method.args = list(alternative = "less", p.adjust.method = "BH")
    ) +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    ylim(6, 16) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

ggsave("results/expression/bcl6_expr_by_igh_partner_coo.pdf", height = 6, width = 6)

apex_test <- md_sv %>%
    filter(
        ICC_class %in% c("BL", "FL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol %in% c("AICDA", "UNG", "POLH", "STAT6")
    ) %>%
    mutate(ICC_class = ifelse(
        ICC_class == "DLBCL",
        glue("{dlbcl_call}-DLBCL"),
        ICC_class
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "FL", "GCB-DLBCL", "HGBCL-DH-BCL2", "ABC-DLBCL")
        )
    ) %>%
    filter(!is.na(ICC_class)) %>%
    group_by(hgnc_symbol) %>%
    pairwise_wilcox_test(
        expression ~ ICC_class,
        p.adjust.method = "BH"
    ) %>%
    add_y_position()


md_sv %>%
    filter(
        ICC_class %in% c("BL", "FL", "DLBCL", "HGBCL-DH-BCL2"),
        hgnc_symbol %in% c("AICDA", "UNG", "POLH", "STAT6")
    ) %>%
    mutate(ICC_class = ifelse(
        ICC_class == "DLBCL",
        glue("{dlbcl_call}-DLBCL"),
        ICC_class
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "FL", "GCB-DLBCL", "HGBCL-DH-BCL2", "ABC-DLBCL")
        )
    ) %>%
    filter(!is.na(ICC_class)) %>%
    ggplot(aes(
        x = ICC_class,
        y = expression
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = ICC_class), size = 0.5) +
    facet_wrap(~hgnc_symbol) +
    scale_colour_manual(
        values = colours$group[c("BL", "FL", "GCB-DLBCL", "HGBCL-DH-BCL2", "ABC-DLBCL")],
        name = ""
    ) +
    stat_pvalue_manual(
        apex_test,
        hide.ns = TRUE,
        tip.length = 0
    ) +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    theme(
        legend.position = "none",
        legend.justification = "left",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14)
    )

ggsave("results/expression/aid_related_expression.pdf", height = 8, width = 6)


md_sv %>%
    filter(
        ICC_class == "HGBCL-DH-BCL2",
        hgnc_symbol %in% c("APEX1", "APEX2", "AICDA", "UNG", "POLH", "STAT6")
    ) %>%
    filter(!is.na(myc_partner_group)) %>%
    mutate(myc_partner_group = factor(
        myc_partner_group,
        levels = names(colours$myc_partner_group)[c(1:4)]
    )) %>%
    ggplot(aes(
        x = myc_partner_group,
        y = expression
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = myc_partner_group), size = 0.5) +
    facet_wrap(~hgnc_symbol) +
    stat_compare_means(
        ref.group = "non-IG",
        label.y = 15,
        label = "p.signif",
        vjust = 2
    ) +
    scale_colour_manual(
        values = colours$myc_partner_group[1:4],
        name = "MYC Partner"
    ) +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave("results/expression/aid_related_dh_by_myc_partner.pdf", height = 6, width = 6)

colours_igh <- colours$ighc[names(colours$ighc) %in% unique(md_sv$hgnc_symbol)]

md_sv %>%
    filter(
        ICC_class %in% c("BL", "FL", "DLBCL", "HGBCL-DH-BCL2"),
        str_detect(hgnc_symbol, "^IGH")
    ) %>%
    mutate(ICC_class = ifelse(
        ICC_class == "DLBCL",
        glue("{dlbcl_call}-DLBCL"),
        ICC_class
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")
        )
    ) %>%
    mutate(hgnc_symbol = factor(hgnc_symbol, levels = names(colours_igh))) %>%
    filter(!is.na(ICC_class)) %>%
    ggplot(aes(
        x = ICC_class,
        y = expression
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = ICC_class), size = 0.5) +
    facet_wrap(~hgnc_symbol) +
    scale_colour_manual(
        values = colours$group[c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")],
        name = ""
    ) +
    stat_compare_means(
        ref.group = "HGBCL-DH-BCL2",
        label.y = 24,
        label = "p.signif",
        vjust = 2
    ) +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave("results/expression/IGH_expression_by_ICC_class.pdf", height = 7, width = 7)

igh_matrix <- md_sv %>%
    filter(
        ICC_class %in% c("BL", "FL", "DLBCL", "HGBCL-DH-BCL2"),
        str_detect(hgnc_symbol, "^IGH")
    ) %>%
    mutate(ICC_class = ifelse(
        ICC_class == "DLBCL",
        glue("{dlbcl_call}-DLBCL"),
        ICC_class
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")
        )
    ) %>%
    filter(!is.na(ICC_class)) %>%
    select(sample_id, hgnc_symbol, expression) %>%
    pivot_wider(
        names_from = "sample_id",
        values_from = "expression"
    ) %>%
    column_to_rownames("hgnc_symbol")

igh_matrix_scaled <- t(apply(igh_matrix, 1, scale))
colnames(igh_matrix_scaled) <- colnames(igh_matrix)
igh_matrix_scaled <- igh_matrix_scaled[rev(rownames(igh_matrix_scaled)), ]

csr <- read_tsv("data/ig_rearrangements/mixcr_all_loci.tsv") %>%
    select(patient_id, biopsy_id, IGH_C_Gene) %>%
    mutate(CSR = case_when(
        IGH_C_Gene %in% c("IGHM", "IGHD") ~ "No",
        str_detect(IGH_C_Gene, "IGH") ~ "Yes"
    )) %>%
    select(-IGH_C_Gene)


md_heatmap <- md_sv %>%
    left_join(csr) %>%
    select(ICC_class, sample_id, dlbcl_call, CSR) %>%
    distinct() %>%
    filter(
        ICC_class %in% c("BL", "FL", "DLBCL", "HGBCL-DH-BCL2")
    ) %>%
    mutate(ICC_class = ifelse(
        ICC_class == "DLBCL",
        glue("{dlbcl_call}-DLBCL"),
        ICC_class
    )) %>%
    mutate(
        ICC_class = factor(
            ICC_class,
            levels = c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")
        )
    ) %>%
    filter(!is.na(ICC_class)) %>%
    column_to_rownames("sample_id")

md_heatmap <- md_heatmap[colnames(igh_matrix_scaled), ]
identical(colnames(igh_matrix_scaled), rownames(md_heatmap))

ha <- HeatmapAnnotation(
    " " = md_heatmap$ICC_class,
    "CSR" = md_heatmap$CSR,
    col = list(
        " " = colours$group[c("BL", "FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")],
        "CSR" = c("Yes" = "#c41230", "No" = "#115284")
    ),
    annotation_legend_param = list(
        "CSR" = list(nrow = 2, title_position = "topleft")
    ),
    show_legend = c(FALSE, TRUE)
)

heatmap_ighc <- Heatmap(
    igh_matrix_scaled,
    name = " ",
    row_labels = rownames(igh_matrix_scaled),
    column_order = colnames(
        igh_matrix_scaled[, order(igh_matrix_scaled[1, ], decreasing = T)]
    ),
    row_names_side = "left",
    top_annotation = ha,
    column_split = md_heatmap$ICC_class,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    heatmap_legend_param = list(
        title = "Normalized Expression",
        legend_height = unit(5, "cm"),
        title_position = "leftcenter-rot"
    ),
    column_title_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10)
)

pdf("results/expression/IGH_expression_by_ICC_class_heatmap.pdf", height = 4, width = 8)
draw(
    heatmap_ighc,
    annotation_legend_side = "right"
)
dev.off()

md_sv %>%
    filter(
        ICC_class == "HGBCL-DH-BCL2",
        str_detect(hgnc_symbol, "^IGH")
    ) %>%
    filter(!is.na(myc_partner_group)) %>%
    mutate(myc_partner_group = factor(
        myc_partner_group,
        levels = names(colours$myc_partner_group)[c(1:4)]
    )) %>%
    mutate(hgnc_symbol = factor(hgnc_symbol, levels = names(colours_igh))) %>%
    filter(!is.na(ICC_class)) %>%
    ggplot(aes(
        x = myc_partner_group,
        y = expression
    )) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = myc_partner_group), size = 0.5) +
    facet_wrap(~hgnc_symbol) +
    scale_colour_manual(
        values = colours$myc_partner_group[1:4],
        name = "MYC Partner"
    ) +
    stat_compare_means(
        ref.group = "IGH",
        label.y = 22,
        label = "p.signif",
        vjust = 2
    ) +
    theme_bw() +
    xlab("") +
    ylab("Normalized Expression") +
    theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave("results/expression/IGH_expression_dh_by_myc_partner.pdf", height = 7, width = 7)
