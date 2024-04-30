source("src/libs.R")
library(tidyverse)
library(glue)
library(GAMBLR)
library(ggalluvial)
library(ggbeeswarm)


md <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv")

mixcr_IGH_filt <- read_tsv("data/ig_rearrangements/mixcr_IGH_filtered.tsv") %>%
    left_join(md) %>%
    mutate(group = case_when(
        ICC_class == "DLBCL" &
            dlbcl_call %in% c("ABC", "GCB", "UNCLASS") ~ glue("{dlbcl_call}-{ICC_class}"),
        TRUE ~ ICC_class
    ))
p1 <- mixcr_IGH_filt %>%
    filter(!is.na(Productive)) %>%
    ggplot(aes(x = log2(totalCloneCount), y = totalCloneFraction)) +
    geom_point() +
    facet_grid(Productive ~ C_Gene)

group_levels <- c(
    "BL",
    "FL",
    "HGBCL-DH-BCL2",
    "GCB-DLBCL",
    "ABC-DLBCL"
)

mixcr_IGH_filt %>%
    filter(
        group %in% group_levels,
        V_Gene != "NULL"
    ) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    group_by(group) %>%
    count(V_Gene) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ungroup() %>%
    filter(Percent > 2) %>%
    ggplot(aes(x = V_Gene, y = Percent, fill = group)) +
    geom_col() +
    geom_text(aes(x = V_Gene, y = Percent + 2, label = ifelse(n > 2, n, ""))) +
    facet_grid(~group) +
    coord_flip() +
    scale_fill_manual(values = colours$group[group_levels]) +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(
    "results/ig_rearrangement/mixcrTopVAllele_usage.pdf",
    height = 6,
    width = 10
)
mixcr_IGH_filt %>%
    filter(!is.na(C_Gene)) %>%
    filter(group %in% c(group_levels)) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    mutate(C_Gene = factor(C_Gene, levels = rev(names(colours$ighc)))) %>%
    ggplot(aes(x = C_Gene, fill = C_Gene)) +
    geom_bar() +
    facet_grid(Productive ~ group) +
    coord_flip() +
    xlab("") +
    scale_fill_manual(values = colours$ighc) +
    theme_bw() +
    theme(legend.position = "none")

mixcr_IGH_filt %>%
    filter(!is.na(C_Gene)) %>%
    filter(group %in% c(group_levels)) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    group_by(group) %>%
    count(C_Gene) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(C_Gene = ifelse(C_Gene == "NULL", "IGHE", C_Gene)) %>%
    mutate(Percent = ifelse(C_Gene == "IGHE", 0, Percent)) %>%
    mutate(n = ifelse(C_Gene == "IGHE", NA, n)) %>%
    mutate(C_Gene = factor(C_Gene, levels = rev(names(colours$ighc)))) %>%
    ggplot(aes(x = C_Gene, y = Percent, fill = C_Gene)) +
    geom_col() +
    geom_text(aes(x = C_Gene, y = Percent + 10, label = n)) +
    facet_grid(~group) +
    coord_flip() +
    scale_fill_manual(values = colours$ighc) +
    ylim(0, 100) +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none")
ggsave(
    "results/ig_rearrangement/mixcr_bestCGene_usage.pdf",
    height = 4,
    width = 10
)

mixcr_IGH_CSR <- mixcr_IGH_filt %>%
    filter(group %in% group_levels) %>%
    mutate(CSR = case_when(
        C_Gene %in% c("IGHM", "IGHD") ~ "FALSE",
        str_detect(C_Gene, "IGH") ~ "TRUE"
    )) %>%
    filter(!is.na(CSR))

mixcr_IGH_pairwise <- table(
    mixcr_IGH_CSR$group,
    mixcr_IGH_CSR$CSR
)

mixcr_IGH_pairwise <- pairwise_fisher_test(
    mixcr_IGH_pairwise,
    detailed = TRUE, p.adjust.method = "BH"
) %>%
    mutate(y.position = c(
        210, # BL ABC ns
        100, # FL ABC
        140, # GCB ABC ns
        110, # ABC DH
        180, # BL FL
        230, # BL GCB ns
        170, # BL DH
        140, # FL GCB ns
        160, # FL DH ns
        120 # GCB DH
    )) %>%
    mutate(CSR = TRUE, group = group1)

mixcr_IGH_CSR %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    group_by(group) %>%
    count(CSR) %>%
    mutate(percent = paste0(round(n / sum(n) * 100), "%")) %>%
    ggplot(aes(
        x = group, y = n, alpha = CSR, fill = group
    )) +
    geom_col(colour = "black") +
    geom_label(aes(
        label = percent,
        colour = group
    ), position = position_stack(vjust = 0.5), show.legend = FALSE) +
    scale_fill_manual(values = colours$group) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.7)) +
    scale_colour_manual(
        values = c(
            "BL" = "black",
            "FL" = "black",
            "HGBCL-DH-BCL2" = "white",
            "ABC-DLBCL" = "black",
            "GCB-DLBCL" = "black"
        )
    ) +
    stat_pvalue_manual(
        mixcr_IGH_pairwise,
        hide.ns = TRUE,
        tip.length = 0,
        show.legend = FALSE
    ) +
    scale_x_discrete(labels = str_replace(group_levels, "-D", "-\nD")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(
        fill = "none"
    ) +
    xlab("") +
    ylab("Number of Tumors")

ggsave("results/ig_rearrangement/mixcr_CSR_proportions.pdf", height = 6, width = 4)

mixcr_IGH_test <- mixcr_IGH_filt %>%
    filter(group %in% group_levels) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    mutate(IGH = factor(
        ifelse(C_Gene != "NULL", "Expressed", "NULL"),
        levels = c("Expressed", "NULL")
    )) %>%
    mutate(ffpe = factor(
        ifelse(preservation == "FFPE", "FFPE", "Frozen"),
        levels = c("Frozen", "FFPE")
    )) %>%
    mutate(myc_partner_igh_binary = factor(
        myc_partner_igh_binary,
        levels = c("non-IGH", "IGH")
    )) %>%
    select(sample_id, IGH, group, ffpe, cloneCount, myc_partner_igh_binary)

mixcr_IGH_test %>%
    count(ffpe, IGH) %>%
    group_by(ffpe) %>%
    mutate(percent = round(n / sum(n) * 100)) %>%
    ungroup() %>%
    ggplot(aes(x = ffpe, y = percent, fill = IGH)) +
    geom_col() +
    xlab("") +
    ylab("Percent")
ggsave("results/ig_rearrangement/mixcr_glm_IGH_null_vs_ffpe.pdf", height = 4, width = 4)

mixcr_IGH_test %>%
    count(ffpe, group, IGH) %>%
    group_by(ffpe, group) %>%
    mutate(percent = round(n / sum(n) * 100)) %>%
    ungroup() %>%
    ggplot(aes(x = ffpe, y = n, fill = IGH)) +
    geom_col(position = "dodge") +
    xlab("") +
    ylab("Percent") +
    coord_flip() +
    facet_wrap(~group, ncol = 1)
ggsave("results/ig_rearrangement/mixcr_glm_IGH_null_vs_ffpe_group.pdf", height = 6, width = 5)

broom::tidy(glm(
    IGH ~ ffpe + group,
    data = mixcr_IGH_test,
    family = "binomial"
)) %>%
    write_tsv("results/ig_rearrangement/mixcr_glm_IGH_null_vs_ffpe_group.tsv")

mixcr_IGH_test %>%
    ggplot(aes(x = ffpe, y = log2(cloneCount), colour = IGH)) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.8) +
    xlab("")
ggsave(
    "results/ig_rearrangement/mixcr_HGBCL-DH-BCL2_CSR_vs_MYC_partner.pdf",
    height = 5, width = 7
)

broom::tidy(glm(
    IGH ~ myc_partner_igh_binary + ffpe,
    data = filter(
        mixcr_IGH_test,
        group == "HGBCL-DH-BCL2"
    ),
    family = "binomial"
)) %>%
    write_tsv("results/ig_rearrangement/mixcr_glm_IGH_null_vs_ffpe_MYC_partner.tsv")

mixcr_igh_dh <- mixcr_IGH_filt %>%
    filter(group %in% c("HGBCL-DH-BCL2")) %>%
    filter(!is.na(myc_partner_igh_binary))

mixcr_igh_dh %>%
    group_by(myc_partner_igh_binary) %>%
    count(C_Gene) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(C_Gene = ifelse(C_Gene == "NULL", "IGHE", C_Gene)) %>%
    mutate(Percent = ifelse(C_Gene == "IGHE", 0, Percent)) %>%
    mutate(n = ifelse(C_Gene == "IGHE", NA, n)) %>%
    mutate(C_Gene = factor(C_Gene, levels = rev(names(colours$ighc)))) %>%
    ggplot(aes(x = C_Gene, y = Percent, fill = C_Gene)) +
    geom_col() +
    geom_text(aes(x = C_Gene, y = Percent + 2, label = n)) +
    facet_grid(~myc_partner_igh_binary) +
    coord_flip() +
    scale_fill_manual(values = colours$ighc) +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none")
ggsave("results/ig_rearrangement/mixcr_DH_BCL2_MYC_IG_non-IG.pdf", height = 5, width = 5)

table(
    MYC_Partner = mixcr_igh_dh$myc_partner_igh_binary,
    Class_Switched = mixcr_igh_dh$CSR
)
#            Class_Switched
# MYC_Partner FALSE TRUE
#     IGH        12    8
#     non-IGH    15   50

broom::tidy(
    fisher.test(
        table(
            MYC_Partner = mixcr_igh_dh$myc_partner_igh_binary,
            Class_Switched = mixcr_igh_dh$CSR
        )
    )
) %>%
    write_tsv("results/ig_rearrangement/DH_CSR_vs_MYC_partner_IGH_non-IGH.fish_test.tsv")

#   estimate p.value conf.low conf.high
#      <dbl>   <dbl>    <dbl>     <dbl>
# 1     4.89 0.00482     1.52      16.7

library(ggsankey)

mixcr_igh_dh %>%
    filter(myc_partner == "IGH", !is.na(myc_partner_igh_discrete)) %>%
    select(sample_id, C_Gene, myc_partner_igh_discrete) %>%
    filter(C_Gene != "NULL") %>%
    pivot_longer(
        -sample_id,
        names_to = "event",
        values_to = "value"
    ) %>%
    mutate(value = factor(value, levels = names(colours$ighc))) %>%
    mutate(event = ifelse(
        event == "C_Gene",
        "Constant Gene\nUsed",
        "MYC Partner"
    )) %>%
    ggplot(aes(
        x = event,
        stratum = value,
        alluvium = sample_id,
        fill = value
    )) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback") +
    geom_stratum() +
    scale_fill_manual(values = colours$ighc, name = "") +
    xlab("") +
    guides(alpha = "none") +
    theme_bw() +
    coord_flip()

ggsave("results/ig_rearrangement/mixcr_C_Gene_vs_MYC_Partner.pdf",
    height = 5, width = 8
)

mixcr_igh_dh %>%
    filter(
        myc_partner == "IGH",
        bcl2_partner == "IGH",
        !is.na(myc_partner_igh_discrete),
        C_Gene != "NULL"
    ) %>%
    select(
        sample_id,
        "Expressed\nConstant Gene" = C_Gene,
        "MYC Partner" = myc_partner_igh_discrete
    ) %>%
    make_long(`Expressed\nConstant Gene`, `MYC Partner`) %>%
    mutate(node = factor(node, levels = rev(names(colours$ighc)))) %>%
    ggplot(aes(
        x = x,
        next_x = next_x,
        node = node,
        next_node = next_node,
        fill = node,
        label = node
    )) +
    geom_sankey(flow.alpha = 0.75) +
    geom_sankey_label(show.legend = FALSE, colour = "white") +
    scale_fill_manual(values = colours$ighc, name = "IGH Region") +
    theme_sankey(base_size = 16) +
    xlab("")

ggsave(
    "results/ig_rearrangement/mixcr_C_Gene_vs_MYC_Partner_sankey.pdf",
    height = 5, width = 8
)
