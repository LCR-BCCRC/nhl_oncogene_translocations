source("src/libs.R")
library(tidyverse)
library(GAMBLR)
library(data.table)
library(cowplot)
library(caret)
library(pROC)

# Load metadata. Select one representative sample.
md <- read_tsv("data/metadata/breakpoint_capture_best.tsv") 

# Load SV data
bcl2_data_raw <- read_tsv("data/sv_data/bcl2_annotated_breaks_unique.tsv")
bcl2_data <- bcl2_data_raw %>%
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

# Quantify mutations in BCL2 TSS

bcl2_tss_muts <- maf_in_region %>%
    filter(Region == "BCL2", Start_Position > 63312793) %>%
    count(sample_id = Tumor_Sample_Barcode) %>% 
    rename(TSS = n)

# BCL2 exon1: chr18:63318046-63318668
# Pad by 200 bp (generous)
bcl2_exon1_muts <- maf_in_region %>% 
    filter(
        Region == "BCL2", 
        Start_Position > 63317846, 
        End_Position < 63318868
    ) %>%
    count(sample_id = Tumor_Sample_Barcode) %>% 
    rename(exon1 = n)
    
bcl2_tss_muts <- bcl2_tss_muts %>% 
    left_join(bcl2_exon1_muts) %>% 
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
    
write_tsv(bcl2_tss_muts, "data/maf/bcl2_tss_mut_counts.tsv")

these_levels <- c("FL", "DLBCL", "HGBCL-DH-BCL2")
these_counts <- c("TSS", "exon1")

bcl2_table <- md %>%
    select(sample_id, patient_id, biopsy_id, bcl2_ba, ICC_class, matches("ihc|_call"), seq_type) %>%
    filter(ICC_class %in% these_levels) %>%
    mutate(ICC_class = factor(ICC_class, levels = these_levels)) %>% 
    left_join(select(sv_data_cat, patient_id, biopsy_id, bcl2_bp_status, partner = bcl2_partner)) %>% 
    left_join(bcl2_tss_muts) %>% 
    mutate(across(all_of(these_counts), ~ replace_na(.x, 0))) %>% 
    mutate(bcl2_ba = ifelse(bcl2_ba %in% c("FAIL", NA), "ND", bcl2_ba)) %>% 
    mutate(bcl2_ba = factor(bcl2_ba, levels = c("NEG", "POS", "ND"))) %>% 
    mutate(Rearrangement = ifelse(!is.na(partner), "POS", "NEG")) %>% 
    pivot_longer(
        c(TSS, exon1), 
        names_to = "which_count", 
        values_to = "n"
    )
    
mut_counts <- bcl2_table %>%
    filter(which_count == "TSS") %>%
    group_by(ICC_class, Rearrangement, bcl2_ba) %>%
    mutate(num_over_threshold = sum(n > 3)) %>%
    mutate(num_in_group = n()) %>% 
    mutate(pct_over_threshold = round(num_over_threshold / num_in_group * 100, digits = 1)) %>%
    ungroup() %>%
    mutate(label = glue("N={num_over_threshold}\n{pct_over_threshold}%")) %>% 
    select(ICC_class, Rearrangement, bcl2_ba, num_in_group, num_over_threshold, pct_over_threshold, label) %>%
    distinct() %>% 
    arrange(ICC_class, bcl2_ba, Rearrangement)
    
write_tsv(select(mut_counts, -label), "results/bcl2_svs/bcl2_sv_shm_counts.tsv")

bcl2_table %>% 
    filter(which_count == "TSS") %>% 
    ggplot(aes(
        x = Rearrangement, 
        y = n, 
        colour = ICC_class
    )) + 
    geom_hline(yintercept = 3.5) + 
    ggbeeswarm::geom_quasirandom(size = 0.5) +
    geom_text(
        data = mut_counts, 
        aes(x = Rearrangement, label = label, y = 120), 
        inherit.aes = FALSE
    ) +
    facet_grid(ICC_class ~ bcl2_ba) + 
    scale_colour_manual(values = colours$ICC_class[these_levels]) + 
    ylab("BCL2 TSS Mutation Count") + ylim(0, 140) +
    theme_bw() +
    theme(legend.position = "none")
    
# ggsave("results/bcl2_svs/bcl2_sv_shm.pdf", height = 6, width = 5)


plot_roc <- function(value, column, table){
    # value = "exon1"
    # column = "which_count"
    # table = bcl2_table
    rocobj <- roc(
        table[table[[column]] == value, ]$Rearrangement,
        table[table[[column]] == value, ]$n
    )

    roccoords <- coords(rocobj, "best", best.method = "youden") %>%
        mutate(across(everything(), ~ round(.x, digits = 2)))

    rocauc <- round(auc(rocobj), 2)

    ggroc(rocobj) +
        geom_point(
            data = roccoords,
            aes(
                x = specificity,
                y = sensitivity
            ),
            colour = "firebrick"
        ) +
        ggrepel::geom_text_repel(
            data = roccoords,
            aes(
                x = specificity,
                y = sensitivity,
                label = glue("{threshold} ({specificity}, {sensitivity})")
            ),
            nudge_x = .15,
            nudge_y = -.03,
            colour = "firebrick"
        ) +
        ggtitle(paste0("ROC Curve ", value, " (AUC = ", rocauc, ")")) +
        theme_bw()

    ggsave(paste0("results/bcl2_svs/bcl2_mut_vs_sv_roc_", value, ".pdf"), height = 4, width = 4)
}

lapply(
    c("exon1", "TSS"), 
    plot_roc,
    column = "which_count", 
    table = bcl2_table
)

# Split by lymphoma entity and ROC again for TSS mutations
lapply(
    unique(bcl2_table$ICC_class), 
    plot_roc,
    column = "ICC_class", 
    table = filter(bcl2_table, which_count == "TSS")
)

# Split by seq_type
lapply(
    c("genome", "capture"),
    plot_roc,
    column = "seq_type",
    table = filter(bcl2_table, which_count == "TSS")
)

# Using the Youden optimal cutoff of 3.5 mutations,
# compare sensitivity and specificity of FISH vs mutations
# across ICC_Class or seq_type

bcl2_table_TSS <- bcl2_table %>%
    filter(which_count == "TSS") %>%
    mutate(mutation_threshold = ifelse(n > 3, "POS", "NEG")) %>%
    mutate(across(
        all_of(c("Rearrangement", "mutation_threshold", "bcl2_ba")),
        ~ factor(.x, levels = c("POS", "NEG"))
    )) 

mut_vs_R <- broom::tidy(caret::confusionMatrix(
    table(
        Rearrangement = bcl2_table_TSS$Rearrangement,
        `Mutations > 3` = bcl2_table_TSS$mutation_threshold
    )
))



write_tsv(mut_vs_R, "results/bcl2_svs/bcl2_R_vs_mut_accuracy_TSS.tsv")


get_accuracy <- function(value, column, table){
    # value = "DLBCL"
    # column = "ICC_class"
    # table = bcl2_table[bcl2_table$which_count == "TSS", ]
    subset_table <- table[table[[column]] == value, ] %>%
        mutate(mutation_threshold = ifelse(n > 3, "POS", "NEG")) %>%
        select(Rearrangement, mutation_threshold, bcl2_ba) %>% 
        mutate(across(
            everything(),
            ~ factor(.x, levels = c("POS", "NEG"))
        ))%>%
        mutate(
            FISH_correct = bcl2_ba == Rearrangement,
            muts_correct = mutation_threshold == Rearrangement
        )
    mut_vs_R <- broom::tidy(caret::confusionMatrix(
        data = subset_table$mutation_threshold,
        reference = subset_table$Rearrangement
    )) %>% 
    mutate(
        group = value,
        test = "mutations_vs_rearrangement"
    )%>% 
        filter(term == "accuracy")
    fish_vs_R <- broom::tidy(caret::confusionMatrix(
        data = subset_table$bcl2_ba,
        reference = subset_table$Rearrangement
    )) %>%     
        mutate(
            group = value,
            test = "FISH_vs_rearrangement"
        ) %>% 
        filter(term == "accuracy")
        
    mcnemar_result <- broom::tidy(mcnemar.test(
        table(
            FISH = subset_table$FISH_correct,
            muts = subset_table$muts_correct
        )
    )) %>%
        select(p.value) %>%
            mutate(
                term = "mcnemar",
                group = value,
                test = "FISH_vs_mutation_accuracy"
            )
    return(bind_rows(mut_vs_R, fish_vs_R, mcnemar_result))
}

region_accuracy <- get_accuracy(
    value = "TSS", 
    column = "which_count", 
    table = bcl2_table
)

ICC_accuracy <- bind_rows(lapply(
    unique(bcl2_table$ICC_class),
    get_accuracy,
    column = "ICC_class",
    table = bcl2_table[bcl2_table$which_count == "TSS", ]
)) 

seq_accuracy <- bind_rows(lapply(
    unique(bcl2_table$seq_type),
    get_accuracy,
    column = "seq_type",
    table = bcl2_table[bcl2_table$which_count == "TSS", ]
)) 

bind_rows(
    region_accuracy,
    ICC_accuracy,
    seq_accuracy
) %>%
    mutate(across(matches("estimate|conf"), ~ round(.x, digits = 3))) %>% 
    mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", as.character(round(p.value, digits = 4)))) %>% 
    mutate(accuracy = glue("{estimate} ({conf.low}-{conf.high})")) %>% 
    write_tsv("results/bcl2_svs/bcl2_R_vs_mut_accuracy_subsets.tsv")

# Repeat for >=1 mutation for exon1

bcl2_table_ex1 <- bcl2_table %>%
    filter(which_count == "exon1") %>%
    mutate(mutation_threshold = ifelse(n > 0, "POS", "NEG")) %>%
    mutate(across(
        all_of(c("Rearrangement", "mutation_threshold", "bcl2_ba")),
        ~ factor(.x, levels = c("POS", "NEG"))
    ))

mut_vs_R <- broom::tidy(caret::confusionMatrix(
    table(
        Rearrangement = bcl2_table_ex1$Rearrangement,
        `Exon 1 Mutations > 1` = bcl2_table_ex1$mutation_threshold
    )
))

write_tsv(mut_vs_R, "results/bcl2_svs/bcl2_R_vs_mut_accuracy_exon1.tsv")

# Compare FISH to WGS/capture status
fish_vs_R <- broom::tidy(caret::confusionMatrix(
    table(
        Rearrangement = bcl2_table_TSS$Rearrangement,
        FISH = bcl2_table_TSS$bcl2_ba, 
        exclude = NA
    )
))

write_tsv(fish_vs_R, "results/bcl2_svs/bcl2_R_vs_FISH_accuracy.tsv")


tss_tx <- bcl2_data %>% 
    filter(ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2", "FL")) %>% 
    mutate(TSS_break = start_target > 63300000) %>% 
    count(TSS_break, partner, bp_status)
    
write_tsv(tss_tx, "results/bcl2_svs/bcl2_tx_location_partner_status.tsv")

bcl2_table %>% 
    filter(which_count == "TSS") %>% 
    group_by(ICC_class) %>% 
    count(bcl2_bp_status) %>% 
    filter(ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2", "FL")) %>% 
    filter(bcl2_bp_status != "UNKNOWN") %>% 
    group_by(ICC_class) %>% 
    mutate(total = sum(n)) %>% 
    pivot_wider(
        names_from = "bcl2_bp_status", 
        values_from = "n", 
        values_fill = 0
    ) %>% 
    write_tsv("results/bcl2_svs/bcl2_bp_status_by_ICC_class.tsv")


# point of reference: Of 226 BLs, 9 have a single mutation in BCL2 TSS, 217 have 0!

bcl2_table %>%
    filter(ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2")) %>% 
    filter(
        !is.na(bcl2_ihc), 
        dlbcl_call == "GCB", 
        bcl2_ba %in% c("NEG", "POS"), 
        which_count == "TSS"
    ) %>% 
    group_by(bcl2_ba, Rearrangement) %>% 
    count(bcl2_ihc) %>% 
    mutate(percent = round(n/sum(n) * 100)) %>% 
    ungroup() %>% 
    ggplot(aes(
        x = Rearrangement, 
        y = percent,
        label = glue("N={n}")
    )) +
    geom_col(aes(fill = bcl2_ihc)) + 
    ggrepel::geom_label_repel(
        aes(colour = bcl2_ihc), 
        position = position_stack(vjust = 0.5),
        direction = "y",
        vjust = 0.7, 
        show.legend = FALSE
    ) +
    facet_grid(~ bcl2_ba) + 
    ggsci::scale_fill_d3(name = "BCL2 IHC") +
    ggsci::scale_colour_d3() +
    theme_bw() + 
    ylab("") +
    theme(legend.position = "bottom")

# ggsave("results/bcl2_svs/bcl2_ihc_vs_r.pdf", height = 4, width = 4)

bcl2_table %>% 
    filter(ICC_class == "HGBCL-DH-BCL2", which_count == "TSS") %>% 
    count(bcl2_ihc)

bcl2_table %>%
    filter(ICC_class %in% c("DLBCL", "HGBCL-DH-BCL2")) %>%
    filter(
        !is.na(bcl2_ihc),
        dlbcl_call == "GCB",
        bcl2_ba %in% c("NEG", "POS"), 
        which_count == "TSS"
    ) %>%
    ggplot(aes(
        x = Rearrangement,
        y = n,
        colour = bcl2_ihc
    )) +
    geom_boxplot() +
    geom_quasirandom(dodge.width = 0.8) +
    facet_grid(~bcl2_ba) +
    ggsci::scale_colour_d3(name = "BCL2 IHC") +
    theme_bw() +
    ylab("BCL2 TSS Mutations (count)") +
    theme(legend.position = "bottom")

ggsave("results/bcl2_svs/bcl2_ihc_vs_mut_count.pdf", height = 4, width = 4)
