source("src/libs.R")

library(tidyverse)
library(GAMBLR)

md <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv")

# Could modify to include HGBCL-NOS.
fct_levels <- c(
    "BL",
    "FL",
    "DLBCL",
    "HGBCL-DH-BCL2",
    "HGBCL-DH-BCL6"
)

md_filt <- md %>%
    filter(ICC_class %in% fct_levels) %>%
    mutate(seq_type = case_when(
        capture & genome ~ "both",
        capture ~ "capture",
        genome ~ "genome"
    )) %>%
    mutate(seq_type = factor(
        seq_type,
        levels = c("capture", "genome", "both")
    ))


md_count <- md_filt %>%
    count(ICC_class) %>%
    mutate(ICC_class = factor(ICC_class, levels = fct_levels))

md_filt %>%
    count(ICC_class, seq_type) %>%
    mutate(ICC_class = factor(ICC_class, levels = fct_levels)) %>%
    ggplot(aes(
        x = reorder(ICC_class, desc(ICC_class)),
        y = n,
        alpha = seq_type,
        fill = ICC_class
    )) +
    geom_col(colour = "black") +
    scale_fill_manual(values = colours$ICC_class[fct_levels], name = "ICC Class") +
    scale_alpha_manual(values = c(
        "both" = 1,
        "genome" = 0.8,
        "capture" = 0.5
    ), name = "") +
    geom_text(
        data = md_count,
        aes(
            x = reorder(ICC_class, desc(ICC_class)),
            y = n + 20,
            label = n
        ), inherit.aes = FALSE
    ) +
    coord_flip() +
    ylab("Count") +
    xlab("") +
    cowplot::theme_cowplot() +
    guides(fill = "none") +
    theme(legend.position = c(0.7, 0.15))

ggsave("results/cohort/seq_type.pdf", height = 4, width = 6)

# Plot SV status

md_sv_status <- md_filt %>%
    select(
        patient_id,
        biopsy_id,
        seq_type,
        ICC_class,
        matches("bp_status")
    ) %>%
    pivot_longer(
        matches("bp_status"),
        names_to = "which_bp",
        values_to = "bp_status"
    ) %>%
    mutate(which_bp = str_remove(which_bp, "_bp_status")) %>%
    filter(bp_status %in% c("TRUE_POS", "FALSE_NEG")) %>%
    mutate(which_bp = factor(
        toupper(which_bp),
        levels = c("MYC", "BCL2", "BCL6")
    )) %>%
    filter(
        !(ICC_class == "BL" & which_bp %in% c("BCL2", "BCL6")),
        !(ICC_class == "FL" & which_bp %in% c("MYC", "BCL6"))
    )

overall_recall <- md_sv_status %>%
    group_by(which_bp) %>%
    count(status = bp_status) %>%
    mutate(percent = round(n / sum(n) * 100, digits = 0)) %>%
    ungroup() %>%
    mutate(ICC_class = "All")

md_sv_status %>%
    group_by(ICC_class, which_bp) %>%
    count(status = bp_status) %>%
    mutate(percent = round(n / sum(n) * 100, digits = 0)) %>%
    ungroup() %>%
    bind_rows(overall_recall) %>%
    ggplot(aes(
        x = ICC_class,
        y = n,
        fill = ICC_class,
        label = str_c(percent, "%"),
        group = status
    )) +
    geom_col(aes(alpha = status), colour = "black") +
    ggrepel::geom_label_repel(
        aes(fill = status, colour = status),
        position = position_stack(vjust = 0.5),
        direction = "y",
        vjust = 0.7
    ) +
    facet_grid(~which_bp, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(
        colours$ICC_class[fct_levels],
        "TRUE_POS" = "#595959",
        "FALSE_NEG" = "#BDBDBD"
    )) +
    scale_colour_manual(
        values = c("FALSE_NEG" = "black", "TRUE_POS" = "white")
    ) +
    scale_alpha_manual(
        values = c(
            "FALSE_NEG" = 0.4,
            "TRUE_POS" = 1
        ),
        name = "",
        labels = c("FISH Pos/Seq Neg", "FISH Pos/Seq Pos")
    ) +
    xlab("") +
    ylab("Count") +
    guides(fill = "none", colour = "none") +
    cowplot::theme_cowplot() +
    theme(
        legend.position = c(0.7, 0.9),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ggtitle("FISH Rearrangements Detected by Sequencing")

ggsave("results/cohort/sv_recall.pdf", height = 6, width = 9)


md_sv_status %>%
    group_by(seq_type, which_bp) %>%
    count(status = bp_status) %>%
    mutate(percent = round(n / sum(n) * 100, digits = 0)) %>%
    ungroup() %>%
    ggplot(aes(
        x = seq_type,
        y = n,
        fill = seq_type,
        label = str_c(percent, "%"),
        group = status
    )) +
    geom_col(aes(alpha = status), colour = "black") +
    ggrepel::geom_label_repel(
        aes(fill = status, colour = status),
        position = position_stack(vjust = 0.5),
        direction = "y",
        vjust = 0.7
    ) +
    facet_wrap(~which_bp) +
    scale_fill_manual(values = c(
        # colours$ICC_class[fct_levels],
        "TRUE_POS" = "#595959",
        "FALSE_NEG" = "#BDBDBD"
    )) +
    scale_colour_manual(
        values = c("FALSE_NEG" = "black", "TRUE_POS" = "white")
    ) +
    scale_alpha_manual(
        values = c(
            "FALSE_NEG" = 0.4,
            "TRUE_POS" = 1
        ),
        name = "",
        labels = c("FISH Pos/Seq Neg", "FISH Pos/Seq Pos")
    ) +
    xlab("") +
    ylab("Count") +
    guides(fill = "none", colour = "none") +
    cowplot::theme_cowplot() +
    theme(
        legend.position = c(0.7, 0.9),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ggtitle("FISH Rearrangements Detected by Sequencing")

ggsave("results/cohort/sv_recall_seq_type.pdf", height = 6, width = 9)

# Test if recall is significantly different for genome vs capture
seq_type_fish_test <- md_sv_status %>%
    filter(
        seq_type %in% c("capture", "genome")
    ) %>%
    group_by(which_bp) %>%
    summarise(table = list(
        table(bp_status, as.character(seq_type))
    )) %>%
    mutate(
        test = map(table, fisher.test),
        tidy = map(test, tidy)
    ) %>%
    unnest(tidy)

write_tsv(seq_type_fish_test, "results/cohort/sv_recall_by_seq_type_test.tsv")

# Get preservation
preservation <- read_tsv("data/metadata/breakpoint_capture_best.tsv") %>%
    filter(seq_type %in% c("genome", "capture")) %>%
    select(patient_id, biopsy_id, preservation) %>%
    mutate(preservation = ifelse(preservation == "FFPE", "FFPE", "Frozen"))

seq_type_glm <- md_sv_status %>%
    filter(
        seq_type %in% c("capture", "genome")
    ) %>%
    mutate(
        bp_status = bp_status == "TRUE_POS"
    ) %>%
    left_join(preservation) %>%
    group_by(which_bp) %>%
    nest() %>%
    mutate(
        test = map(data, ~ glm(bp_status ~ ICC_class + seq_type + preservation, data = .x, family = "binomial")),
        tidy = map(test, tidy)
    ) %>%
    unnest(tidy)

write_tsv(seq_type_glm, "results/cohort/sv_recall_by_seq_type_ICC_class_glm.tsv")

# Get coverage

coverage_files <- c(
    dir(
        "/projects/dscott_prj/CCSRI_1500/capture/results/QC/04-mergedMetrics/",
        pattern = "[.]hs_metrics.txt",
        full.names = TRUE
    ),
    "/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/gambl/qc-1.0/99-outputs/genome.qc_metrics.tsv"
)

get_coverage <- function(file) {
    # file <- coverage_files[1]
    cov <- read_tsv(file) %>%
        select(
            matches("sampleID|UID"),
            matches("MeanCorrectedCoverage|MEAN_TARGET_COVERAGE")
        )
    colnames(cov) <- c("sample_id", "coverage")
    return(cov)
}

coverage <- bind_rows(lapply(coverage_files, get_coverage)) %>%
    left_join(read_tsv("data/metadata/breakpoint_capture_best.tsv")) %>%
    filter(!is.na(patient_id)) %>%
    select(patient_id, biopsy_id, coverage) %>%
    group_by(patient_id, biopsy_id) %>%
    slice_max(coverage) %>%
    ungroup()

coverage_test <- md_sv_status %>%
    left_join(coverage) %>%
    group_by(which_bp) %>%
    rstatix::wilcox_test(
        coverage ~ bp_status,
        alternative = "less",
        ref.group = "TRUE_POS"
    )

write_tsv(coverage, "data/metadata/biopsy_best_coverage.tsv")
write_tsv(coverage_test, "data/metadata/coverage_bp_status_test.tsv")

md_filt %>%
    group_by(ICC_class, genome) %>%
    filter(!is.na(bcl2_ba)) %>%
    count(bcl2_ba) %>%
    mutate(percent = round(n / sum(n) * 100))

md_to_pub <- md_filt %>%
    select(
        -who_class,
        -morphology,
        -pathology,
        -cpr_complete,
        -specs_id,
        -matches("discrete|binary|group")
    ) %>%
    mutate(transformed = ifelse(transformed == "YES", "YES", NA))

md_to_pub %>%
    count(ICC_class, transformed) %>%
    group_by(ICC_class) %>%
    mutate(percent = round(n / sum(n) * 100)) %>%
    write_tsv("data/metadata/transformation_summary.tsv")

seq_md_to_pub <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    select(
        patient_id,
        biopsy_id,
        sample_id,
        ICC_class,
        tissue_status,
        seq_type,
        protocol,
        preservation,
        capture_version,
        publication,
        library_id
    ) %>%
    mutate(protocol = ifelse(seq_type == "mrna", protocol, NA)) %>%
    inner_join(select(md_to_pub, patient_id, biopsy_id))

seq_md_normals <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    filter(tissue_status == "normal") %>%
    select(patient_id, seq_type, normal_sample_id = sample_id) %>%
    distinct(patient_id, seq_type, .keep_all = TRUE)

seq_md_to_pub <- seq_md_to_pub %>%
    left_join(seq_md_normals)

write_tsv(seq_md_to_pub, "data/metadata/breakpoint_capture_seq_for_pub.tsv")

seq_md_to_pub %>%
    count(ICC_class, seq_type, tissue_status, publication) %>%
    write_tsv("data/metadata/sequencing_summary.tsv")

md %>%
    mutate(myc_r = case_when(
        myc_ba == "POS" ~ "FISH",
        !is.na(myc_partner) & myc_ba == "NEG" ~ "CRYPTIC",
        !is.na(myc_partner) ~ "SEQ",
        TRUE ~ "NEG"
    )) %>%
    mutate(bcl2_r = case_when(
        bcl2_ba == "POS" ~ "FISH",
        !is.na(bcl2_partner) & bcl2_ba == "NEG" ~ "CRYPTIC",
        !is.na(bcl2_partner) ~ "SEQ",
        TRUE ~ "NEG"
    )) %>%
    mutate(bcl6_r = case_when(
        bcl6_ba == "POS" ~ "FISH",
        !is.na(bcl6_partner) & bcl6_ba == "NEG" ~ "CRYPTIC",
        !is.na(bcl6_partner) ~ "SEQ",
        TRUE ~ "NEG"
    )) %>%
    select(
        patient_id,
        biopsy_id,
        ICC_class,
        matches("_r$|_ba")
    ) %>%
    filter(
        ICC_class == "HGBCL-DH-BCL2" & (!bcl2_ba %in% c("POS") | !myc_ba %in% c("POS")) |
            ICC_class == "HGBCL-DH-BCL6" & (!bcl6_ba %in% c("POS") | !myc_ba %in% c("POS"))
    ) %>%
    write_tsv(
        "data/metadata/rearrangements_found_by_sequencing.tsv"
    )
