source("src/libs.R")
library(tidyverse)
library(glue)
library(GAMBLR)
library(ggalluvial)
library(ggbeeswarm)
library(rstatix)
library(ggrepel)


md <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>% 
    filter(seq_type == "mrna") %>% 
    group_by(patient_id, biopsy_id) %>% 
    slice_max(preservation != "FFPE", n=1, with_ties = FALSE) %>%
    ungroup() %>% 
    select(patient_id, biopsy_id, preservation) %>% 
    left_join(read_tsv("data/metadata/breakpoint_capture_biopsies.tsv"))

mixcr_summary <- read_tsv("data/ig_rearrangements/mixcr_all_loci.tsv") %>% 
    left_join(md) %>% 
    mutate(group = ifelse(
        ICC_class == "DLBCL", 
        glue("{dlbcl_call}-DLBCL"), 
        ICC_class
    ))

flow <- read_tsv("data/flow_data/kappa_lambda_flow_categorized.tsv")

mixcr_flow <- inner_join(mixcr_summary, flow) %>% 
    mutate(mixcr_predicted = case_when(
        IGH_Productive == "Yes" & IGK_Productive == "Yes" & (IGK_totalCloneCount > IGL_totalCloneCount) ~ "Lambda-Kappa+", 
        IGH_Productive == "Yes" & IGL_Productive == "Yes" & (IGL_totalCloneCount > IGK_totalCloneCount) ~ "Lambda+Kappa-",
        IGH_Productive == "Yes" & IGK_Productive == "Yes" & IGL_Productive == "Yes" ~ "Lambda+Kappa+", 
        IGH_Productive == "Yes" & IGK_Productive == "Yes" ~ "Lambda-Kappa+", 
        IGH_Productive == "Yes" & IGL_Productive == "Yes" ~ "Lambda+Kappa-", 
        TRUE ~ "Lambda-Kappa-"
    )) 
    
write_tsv(mixcr_flow, "data/ig_rearrangements/mixcr_vs_flow.tsv")

mixcr_flow %>% 
    ggplot(aes(
        x = surface_IG, 
        y = percent
    )) +
    geom_quasirandom() + 
    facet_wrap(~ mixcr_predicted) + 
    coord_flip() + 
    ylim(0, 100) + 
    ylab("Percent of CD19+/CD20+ CD3- population") + 
    xlab("Surface Ig Category") + 
    theme_bw()
    
ggsave("results/ig_rearrangement/flow_sIg_mixcr.pdf", height = 5, width = 6)

table(sIG = mixcr_flow$surface_IG, mixcr = mixcr_flow$mixcr_predicted)

mixcr_flow %>% 
    group_by(mixcr_predicted) %>% 
    count(surface_IG) %>% 
    mutate(percent = round(n/sum(n) * 100)) %>% 
    ungroup() %>% 
    mutate(label = ifelse(mixcr_predicted == surface_IG, paste0(percent, "%"), "")) %>% 
    ggplot(aes(
        x = mixcr_predicted, 
        y = n,
        label = label, 
        fill = surface_IG
    )) + 
    geom_col() + 
    geom_text(position = position_stack(vjust = 0.5)) + 
    ggsci::scale_fill_d3() +
    coord_flip() + 
    xlab("MiXCR Predicted Light Chain") + 
    ylab("Count") + 
    theme_bw() 
    
ggsave("results/ig_rearrangement/flow_sIg_mixcr_counts.pdf", height = 3, width = 8)

mixcr_flow %>%
    filter(group %in% c("HGBCL-DH-BCL2", "ABC-DLBCL", "GCB-DLBCL", "FL")) %>% 
    group_by(mixcr_predicted, group) %>%
    count(surface_IG) %>%
    mutate(percent = round(n / sum(n) * 100)) %>%
    ungroup() %>%
    mutate(label = ifelse(mixcr_predicted == surface_IG, paste0(percent, "%"), "")) %>%
    mutate(group = factor(group, levels = c("FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL"))) %>% 
    ggplot(aes(
        x = mixcr_predicted,
        y = n,
        label = label,
        fill = surface_IG
    )) +
    geom_col() +
    geom_text(position = position_stack(vjust = 0.5)) + 
    ggsci::scale_fill_d3() +
    facet_wrap(~ group) +
    coord_flip() +
    xlab("MiXCR Predicted Light Chain") +
    ylab("Count") +
    theme_bw()

ggsave("results/ig_rearrangement/flow_sIg_mixcr_counts.pdf", height = 3, width = 8)


mixcr_flow %>%
    group_by(mixcr_predicted, preservation) %>%
    count(surface_IG) %>%
    mutate(percent = round(n / sum(n) * 100)) %>%
    ungroup() %>%
    mutate(label = ifelse(mixcr_predicted == surface_IG, paste0(percent, "%"), "")) %>%
    ggplot(aes(
        x = mixcr_predicted,
        y = n,
        label = label,
        fill = surface_IG
    )) +
    geom_col() +
    geom_text(position = position_stack(vjust = 0.5)) + 
    ggsci::scale_fill_d3() +
    facet_wrap(~preservation) +
    coord_flip() +
    xlab("MiXCR Predicted Light Chain") +
    ylab("Count") +
    theme_bw()

ggsave("results/ig_rearrangement/flow_sIg_mixcr_preservation.pdf", height = 3, width = 8)


test_data <- mixcr_flow %>% 
    filter(group %in% c("GCB-DLBCL", "ABC-DLBCL", "FL", "HGBCL-DH-BCL2")) %>% 
    mutate(sIg_NULL = ifelse(surface_IG == "Lambda-Kappa-", TRUE, FALSE)) %>% 
    mutate(group = factor(group, levels = c("FL", "HGBCL-DH-BCL2", "GCB-DLBCL", "ABC-DLBCL")))


test_data %>% 
    group_by(group) %>% 
    count(surface_IG) %>% 
    mutate(percent = str_c(round(n/sum(n) * 100), "%")) %>% 
    ungroup() %>% 
    ggplot(aes(
        x = surface_IG, 
        y = n,
        fill = surface_IG
    )) + 
    geom_col() + 
    geom_text(aes(label = percent), nudge_y = 3) +
    facet_wrap(~ group) + 
    coord_flip() + 
    theme_bw() + 
    theme(legend.position = "none") + 
    ggsci::scale_fill_d3() +
    xlab("") + ylab("Count")
    
ggsave("results/ig_rearrangement/flow_sIg_by_group.pdf", height = 3, width = 8)
    
fish_tab <- table(test_data$sIg_NULL, test_data$group)

fishers_res <- pairwise_fisher_test(fish_tab, p.adjust.method = "BH", detailed = TRUE) %>% 
    mutate(across(
        c(estimate, conf.low, conf.high), 
        log
    )) %>%
    arrange(p.adj) %>% 
    mutate(y.position = 40) %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.1 & p.adj > 0.01, "*", p.adj.signif))


test_data %>% 
    group_by(group) %>% 
    count(sIg_NULL) %>% 
    mutate(percent = round(n/sum(n) * 100)) %>% 
    # filter(sIg_NULL == TRUE) %>% 
    ggplot(aes(x = group, y = n, label = ifelse(sIg_NULL == TRUE, str_c(percent, "%"), ""))) + 
    geom_col(aes(fill = group, alpha = sIg_NULL)) + 
    geom_text(nudge_y = 1) +
    stat_pvalue_manual(
        fishers_res, 
        step.increase = .05, 
        hide.ns = TRUE
    ) + 
    xlab("") + ylab("Count") + 
    theme_bw() + 
    scale_fill_manual(
        values = colours$group[levels(test_data$group)], 
        guide = "none",
    ) + 
    scale_alpha_manual(
        values = c("TRUE" = 1, "FALSE" = 0.6), 
        name = "Surface Ig NULL"
    ) +
    scale_x_discrete(
        labels = str_replace(
            levels(test_data$group), 
            "DH-BCL2", "DH\n-BCL2"
        )
    ) +
    theme(legend.position = c(0.8, 0.8)) 
    
ggsave("results/ig_rearrangement/sIg_NULL_barplot.pdf", height = 5, width = 4)
