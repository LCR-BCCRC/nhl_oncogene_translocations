source("src/libs.R")
library(tidyverse)
library(ComplexHeatmap)
library(GAMBLR)
library(data.table)
library(cowplot)

source("src/snvs/functions/calc_mut_freq_window.R")
source("src/snvs/functions/get_mutation_frequency_bin_matrix.r")
source("src/snvs/functions/process_regions.R")
source("src/snvs/functions/heatmap_mutation_frequency_bin.R")


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
    ))

# Load MAF data
maf_data <- read_tsv("data/maf/genome_capture.WRCY.hg38.maf")
maf <- maf_data


panel_regions <- fread(
    "data/region_data/capture_TE99028370--hg38.bed",
    col.names = c("chrom", "start", "end", "name")
)

gambl_svs_raw <- get_combined_sv(
    sample_ids = sv_data_cat$sample_id, 
    projection = "hg38"
)

regions <- list(
    IGH = "chr14:105550000-107043718",
    IGHV = "chr14:105861682-107043718"
)

regions <- lapply(regions, GAMBLR:::region_to_chunks)

gambl_svs_igh <- gambl_svs_raw %>%
    filter(CHROM_A == "chr14" & CHROM_B == "chr14") %>% 
    filter(
        START_A > regions$IGH$start, 
        START_B > regions$IGH$start
    )

capture_sv_files <- dir(
    "../CCSRI_1500/capture/results/annotate_sv_vcf/02-bedpe/", 
    pattern = ".bedpe$", 
    full.names = TRUE, 
    recursive = TRUE
)

colspec <- cols_condense(spec(read_tsv(
    "../CCSRI_1500/capture/results/annotate_sv_vcf/02-bedpe//capture_TE99028370--hg38/CLC02034--None--no_normal.annotated.bedpe"
)))

capture_svs_raw <- lapply(capture_sv_files, function(x){
    read_tsv(x, col_types = colspec)
})

capture_svs_igh <- bind_rows(capture_svs_raw) %>%
    filter(CHROM_A == "chr14" & CHROM_B == "chr14") %>%
    filter(
        START_A > regions$IGH$start,
        START_B > regions$IGH$start
    )
    

svs <- bind_rows(gambl_svs_igh, capture_svs_igh) 

# # Are SVs in the VDJ region enriched for RSS signal sequences per https://www.itb.cnr.it/rss/analyze.html?

# # Create a bed file with 50 bp padding around SV sites in VDJ region
# vdj_svs <- svs %>% 
#     dplyr::filter(
#         START_A > as.numeric(regions$IGHV$start) & 
#         END_A > as.numeric(regions$IGHV$start) &
#         START_B > as.numeric(regions$IGHV$start) & 
#         END_B > as.numeric(regions$IGHV$start)
#     ) %>% 
#     mutate(across(matches("START_"), ~ .x - 50)) %>% 
#     mutate(across(matches("END_"), ~ .x + 50)) %>% 
#     mutate(sv_name = glue("{tumour_sample_id}_{START_A}_{START_B}")) %>% 
#     distinct(sv_name, .keep_all = TRUE)

# # Write each half of the bedpe to file. Use bedtools getfasta to create fasta file of sequences.     
# write_tsv(
#     select(vdj_svs, CHROM_A, START_A, END_A, sv_name), 
#     "results/ig_rearrangement/vdj/vdj_1.bed", 
#     col_names = FALSE
#     )
    
# write_tsv(
#     select(vdj_svs, CHROM_B, START_B, END_B, sv_name), 
#     "results/ig_rearrangement/vdj/vdj_2.bed", 
#     col_names = FALSE
#     )
    
# # Load RSS predictions from "RSS Site". 

# load_rss <- function(filepath){
#     df <- read_tsv(
#         filepath, 
#         skip = 2, 
#         col_names = c(
#             "sequence_name", 
#             "start", 
#             "end",
#             "RSS_sequence", 
#             "strand", 
#             "RIC_score", 
#             "RIC_pass_fail"
#         )
#     )
#     spacer_length <- nchar(df$RSS_sequence)[1] - 16
#     return(
#         mutate(df, spacer_length = spacer_length)
#     )
# }

# rss_files_A <- dir(
#     "results/ig_rearrangement/vdj/", 
#     pattern = "vdj_1.*bp.txt", 
#     full.names = TRUE
# )

# rss_A <- bind_rows(lapply(rss_files_A, load_rss)) %>% 
#     filter(RIC_pass_fail == "PASS") %>% 
#     rename_with(~ str_c(.x, "_A"), -matches("name")) %>% 
#     select(sequence_name, RIC_pass_fail_A, spacer_length_A) %>% 
#     distinct()
    
# rss_files_B <- dir(
#     "results/ig_rearrangement/vdj/",
#     pattern = "vdj_2.*bp.txt",
#     full.names = TRUE
# )

# rss_B <- bind_rows(lapply(rss_files_B, load_rss)) %>%
#     filter(RIC_pass_fail == "PASS") %>%
#     rename_with(~ str_c(.x, "_B"), -matches("name")) %>%
#     select(sequence_name, RIC_pass_fail_B, spacer_length_B) %>%
#     distinct()

# # Combine A and B half of annotated bedpe   
# rss_combined <- full_join(rss_A, rss_B, relationship = "many-to-many") %>% 
#     mutate(
#         meets_12_23 = case_when(
#             spacer_length_A != spacer_length_B ~ "unequal", 
#             spacer_length_A == spacer_length_B ~ "equal"
#         )
#     ) %>% 
#     filter(!is.na(meets_12_23)) %>% 
#     select(-matches("RIC"))

# # Combine RSS annotations with bedpe
# vdj_svs_rss <- vdj_svs %>% 
#     mutate(across(matches("START_"), ~ .x + 50)) %>% 
#     mutate(across(matches("END_"), ~ .x - 50)) %>% 
#     left_join(rss_combined, by = c("sv_name" = "sequence_name")) %>% 
#     group_by(sv_name) %>% 
#     slice_max(!is.na(meets_12_23), n=1, with_ties = TRUE)%>% 
#     slice_max(meets_12_23 == "unequal", n=1, with_ties = FALSE) %>% 
#     ungroup()

# write_tsv(vdj_svs_rss, "results/ig_rearrangement/vdj/vdj_svs_rss_annotated_pairs.tsv")    


##### CSR Rearrangements #####
svs %>% 
    filter((ANNOTATION_B %in% c("IGHM", "Emu") & str_detect(ANNOTATION_A, "IGH[AGEM]|Emu")) |
        (ANNOTATION_A %in% c("IGHM", "Emu") & str_detect(ANNOTATION_B, "IGH[AGEM]|Emu"))) %>% 
    group_by(tumour_sample_id) %>% 
    summarize(num_csr_events = n()) %>% 
    ungroup() %>% 
    right_join(sv_data_cat, by = c("tumour_sample_id" = "sample_id")) %>% 
    filter(!is.na(ICC_class)) %>% 
    mutate(num_csr_events = replace_na(num_csr_events, 0)) %>% 
    mutate(CSR_event = ifelse(num_csr_events > 0, "CSR", "No CSR")) %>% 
    select(sample_id = tumour_sample_id, num_csr_events, CSR_event) %>% 
    write_tsv("data/ig_rearrangements/svar_master_evidence_of_csr.tsv")
    # count(ICC_class, CSR_event = num_csr_events > 0)

#chr14:105855000 - 105865861

svs_emu <- svs %>% 
    filter(ANNOTATION_B %in% c("Emu", "IGHM")) %>% 
    filter(START_B > 105855000) %>% 
    select(
        sample_id = tumour_sample_id, 
        position = START_B
    ) %>% 
    left_join(sv_data_cat) %>%
    filter(!is.na(ICC_class)) %>%
    mutate(type = "CSR")

maf_emu <- maf_data %>% 
    filter(
        Chromosome == "chr14", 
        Start_Position > 105855000, 
        End_Position < 105865861
    ) %>% 
    mutate(VAF = t_alt_count / t_depth) %>% 
    select(
        sample_id = Tumor_Sample_Barcode, 
        position = Start_Position
    ) %>% 
    left_join(sv_data_cat) %>% 
    filter(!is.na(ICC_class)) %>% 
    mutate(type = "SHM")

ggplot() + 
    geom_histogram(
        data = maf_emu,
        aes(x = position, fill = ICC_class), 
        alpha = 0.5
    ) + 
    geom_histogram(
        data = svs_emu, 
        aes(x = position, colour = ICC_class), 
        fill = NA
    ) + 
    scale_fill_manual(values = colours$ICC_class[1:5]) +
    scale_colour_manual(values = colours$ICC_class[1:5]) + 
    facet_grid(ICC_class ~ 1) + 
    theme(
        strip.background.x = element_blank(), 
        strip.text.x = element_blank()
    )

myc_svs <- read_tsv(
    "data/sv_data/myc_annotated_breaks_unique.tsv"
) %>% 
    filter(
        chrom_partner == "chr14", 
        start_partner > 105855000, 
        start_partner < 105865861
    ) %>% 
    select(patient_id, biopsy_id, position = start_partner) %>%
    left_join(sv_data_cat) %>% 
    mutate(type = "MYC")

bind_rows(maf_emu, svs_emu) %>% 
    bind_rows(myc_svs) %>% 
    mutate(ICC_class = case_when(
        ICC_class == "BL" & BL_EBV == "EBV-positive" ~ "BL EBV+", 
        ICC_class == "BL" & BL_EBV == "EBV-negative" ~ "BL EBV-", 
        TRUE ~ ICC_class
    )) %>% 
    mutate(BL_EBV = ifelse(!str_detect(ICC_class, "BL"), "EBV-negative", BL_EBV)) %>% 
    mutate(ICC_class = str_remove(ICC_class, "HGBCL-")) %>% 
    mutate(ICC_class = factor(ICC_class, levels = c(
        "BL EBV+", 
        "BL EBV-",
        "FL", 
        "DH-BCL2", 
        "DLBCL"
    ))) %>% 
    filter(!ICC_class %in% c("HGBCL-DH-BCL6", NA)) %>% 
    ggplot(aes(x = position, fill = ICC_class, alpha = BL_EBV)) + 
    geom_histogram(binwidth = 100, position = "dodge") + 
    geom_vline(xintercept = 105861000, colour = "firebrick") + 
    scale_fill_manual(values = c(
        "BL EBV-" = unname(colours$group["BL"]),
        "BL EBV+" = unname(colours$group["BL"]),
        "DH-BCL2" = unname(colours$ICC_class["HGBCL-DH-BCL2"]),
        colours$ICC_class[c(
            "FL", 
            "DLBCL"
    )])) +
    scale_alpha_manual(values = c("EBV-negative" = 1, "EBV-positive" = 0.7)) +
    # facet_grid(rows = vars(ICC_class, type), scales = "free_y") +
    facet_grid(type ~ ICC_class, scales = "free_y") + 
    xlab("chr8 Position (Mb)") + ylab("") + 
    theme_bw() +
    theme(
        # strip.background.x = element_blank(), 
        # strip.text.x = element_blank(), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) 

ggsave("results/shm_heatmaps/Emu_muts_vs_CSR_svs.pdf", height = 5, width = 12)

csr_boundary <- 105861000

count_csr <- read_tsv(
    "data/sv_data/myc_annotated_breaks_unique.tsv"
) %>% 
    filter(partner == "IGH") %>% 
    select(patient_id, biopsy_id, position = start_partner) %>% 
    mutate(myc_csr_shm = ifelse(
        position > csr_boundary, 
        "SHM", 
        "CSR"
    )) %>%
    left_join(sv_data_cat) %>%
        mutate(ICC_class = case_when(
            ICC_class == "BL" & BL_EBV == "EBV-positive" ~ "BL EBV+",
            ICC_class == "BL" & BL_EBV == "EBV-negative" ~ "BL EBV-",
            TRUE ~ ICC_class
        )) %>%
        mutate(BL_EBV = ifelse(!str_detect(ICC_class, "BL"), "EBV-negative", BL_EBV)) %>%
        mutate(ICC_class = str_remove(ICC_class, "HGBCL-")) %>%
        mutate(ICC_class = factor(ICC_class, levels = c(
            "BL EBV+",
            "BL EBV-",
            "DH-BCL2",
            "DLBCL"
        ))) %>%
        filter(!ICC_class %in% c("HGBCL-DH-BCL6", NA, "FL")) 
        
fish_tab <- table(count_csr$ICC_class, count_csr$myc_csr_shm)

fish_test <- rstatix::pairwise_fisher_test(fish_tab) %>% 
    mutate(y.position = c(
        92, 
        96, 
        100, 
        80,
        84, 
        60 
    )) %>% 
    mutate(ICC_class = group1)

count_csr %>% 
    count(ICC_class, myc_csr_shm) %>% 
    group_by(ICC_class) %>% 
    mutate(percent = paste0(round(n/sum(n) * 100), "%")) %>% 
    ungroup() %>% 
    ggplot(aes(
        x = ICC_class, 
        y = n, 
        fill = ICC_class, 
        alpha = myc_csr_shm
    )) + 
    geom_col() +
    geom_label(
        aes(label = percent, colour = ICC_class), 
        position = position_stack(vjust = 0.8), 
        show.legend = FALSE
    ) +
    scale_fill_manual(values = c(
        "BL EBV-" = unname(colours$group["BL"]),
        "BL EBV+" = "#B398C6",
        "DH-BCL2" = unname(colours$ICC_class["HGBCL-DH-BCL2"]),
        colours$ICC_class[c(
            "FL", 
            "DLBCL"
    )])) +
    scale_colour_manual(
        values = c(
            "BL EBV-" = "black",
            "BL EBV+" = "black",
            "DH-BCL2" = "white",
            "DLBCL" = "black"
        )
    ) +
    scale_alpha_manual(
        values = c("SHM" = 1, "CSR" = 0.7), 
        name = "MYC-R Mechanism"
    ) + 
    ggpubr::stat_pvalue_manual(
        data = fish_test, 
        inherit.aes = FALSE, 
        hide.ns = TRUE
    ) +
    guides(fill = "none", colour = "none") +
    theme_bw() + 
    theme(legend.position = "bottom") +
    xlab("") + ylab("Count")

ggsave("results/shm_heatmaps/CSR_vs_SHM_MYCR.pdf", height = 5, width = 3.5)
