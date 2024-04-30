source("src/libs.R")
library(tidyverse)
library(glue)
library(GAMBLR)


# Load full project metadata
metadata <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    filter(seq_type == "mrna") %>%
    mutate(group = case_when(
        ICC_class == "DLBCL" ~ glue("{dlbcl_call} DLBCL"),
        TRUE ~ ICC_class
    ))

# Load and tidy MYC partner data
myc_data <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv") 

# Paths to mixcr outputs
gambl_root <- "/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/gambl/"
llmpp_root <- "/projects/dscott_prj/CCSRI_1500/transcriptomes/results/"
mixcr_path <- "mixcr-1.2/99-outputs/txt/mrna/mixcr.{sample_id}.clonotypes.{chain}.igblast.txt"

colspec <- cols(
    .default = col_character(),
    cloneId = col_double(),
    cloneCount = col_double(),
    cloneFraction = col_double(),
    nSeqFR1 = col_logical(),
    minQualFR1 = col_logical(),
    nSeqCDR1 = col_logical(),
    minQualCDR1 = col_logical(),
    nSeqFR2 = col_logical(),
    minQualFR2 = col_logical(),
    minQualCDR2 = col_double(),
    nSeqFR3 = col_logical(),
    minQualFR3 = col_logical(),
    minQualCDR3 = col_double(),
    minQualFR4 = col_double(),
    aaSeqFR1 = col_logical(),
    aaSeqCDR1 = col_logical(),
    aaSeqFR2 = col_logical(),
    aaSeqFR3 = col_logical(),
    numMissingRegions = col_double()
)

read_mixcr <- function(id, path) {
    mixcr <- read_tsv(path, col_types = colspec) %>%
        mutate(sample_id = id)
    if (dim(mixcr)[1] == 0) {
        mixcr <- add_row(mixcr, tibble_row(sample_id = id))
    }
    return(mixcr)
}

load_mixcr <- function(sample_id, chain) {
    this_mixcr <- glue(mixcr_path)
    gambl_mixcr <- paste0(gambl_root, this_mixcr)
    llmpp_mixcr <- paste0(llmpp_root, this_mixcr)
    if (file.exists(gambl_mixcr)) {
        mixcr <- read_mixcr(sample_id, gambl_mixcr)
        return(mixcr)
    } else if (file.exists(llmpp_mixcr)) {
        mixcr <- read_mixcr(sample_id, llmpp_mixcr)
        return(mixcr)
    } else {
        message(glue("No data available for {sample_id}"))
    }
}

mixcr_IGH <- lapply(metadata$sample_id, load_mixcr, chain = "IGH")
mixcr_IGH <- bind_rows(mixcr_IGH)

mixcr_IGH_filt <- mixcr_IGH %>%
    select(
        sample_id,
        Productive,
        mutatedStatus,
        matches("clone"),
        sequenceIdentity = igblastnTopSequenceIdentity,
        matches("[VJ]Allele"),
        matches("Gene"),
        matches("Score"),
        -matches("^Seq", ignore.case = FALSE)
    ) %>%
    mutate(
        V_Gene = str_remove(str_remove(igblastnTopVAllele, "[*].*"), ",.*"),
        V_Gene = case_when(
            !str_detect(V_Gene, "IGH") ~ bestVGene,
            TRUE ~ V_Gene
        )
    ) %>%
    mutate(
        J_Gene = str_remove(str_remove(igblastnJAllele, "[*].*"), ",.*"),
        J_Gene = case_when(
            !str_detect(J_Gene, "IGH") ~ bestJGene,
            TRUE ~ J_Gene
        )
    ) %>%
    mutate(C_Gene = bestCGene) %>% 
    group_by(sample_id, V_Gene, J_Gene, C_Gene) %>% 
    mutate(cloneId = paste0(cloneId, collapse = ",")) %>% 
    mutate(totalCloneCount = sum(cloneCount)) %>% 
    mutate(totalCloneFraction = sum(cloneFraction)) %>% 
    mutate(meanSequenceIdentity = mean(as.numeric(sequenceIdentity), na.rm = TRUE)) %>% 
    ungroup() %>% 
    filter(
        totalCloneCount >= 100,
        totalCloneFraction > 0.1
    ) %>%
    group_by(sample_id) %>%
    slice_max(Productive == "Yes", with_ties = TRUE) %>%
    slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(
        sample_id,
        Productive,
        mutatedStatus,
        matches("totalClone"),
        matches("[sS]equenceIdentity"),
        matches("clone"),
        matches("_Gene$"),
        everything()
    ) %>%
    distinct() %>%
    full_join(distinct(select(mixcr_IGH, sample_id))) %>%
    left_join(select(
        metadata, 
        sample_id, 
        patient_id, 
        biopsy_id, 
        preservation
    )) %>%
    rowwise() %>%
    mutate(
        across(
            matches("_Gene"),
            ~ ifelse(is.na(.x) | Productive == "N/A", "NULL", .x)
        )
    ) %>%
    mutate(Productive = str_replace(Productive, "N/A", "No")) %>%
    mutate(CSR = case_when(
        C_Gene == "IGHM" ~ FALSE,
        C_Gene != "NULL" ~ TRUE
    )) %>%
    ungroup() %>%
    distinct(sample_id, .keep_all = TRUE) %>% 
    group_by(patient_id, biopsy_id) %>% 
    slice_max(!is.na(totalCloneFraction), n=1, with_ties = FALSE) %>% 
    ungroup() %>%
    mutate(locus = "IGH")

write_tsv(mixcr_IGH_filt, "data/ig_rearrangements/mixcr_IGH_filtered.tsv")


mixcr_IGK <- lapply(metadata$sample_id, load_mixcr, chain = "IGK")
mixcr_IGK <- bind_rows(mixcr_IGK)

mixcr_IGK_filt <- mixcr_IGK %>%
    select(
        sample_id,
        Productive,
        mutatedStatus,
        matches("clone"),
        matches("[VJ]Allele"),
        matches("Gene"),
        matches("Score"),
        -matches("Seq")
    ) %>%
    mutate(
        V_Gene = str_remove(str_remove(igblastnTopVAllele, "[*].*"), ",.*"),
        V_Gene = case_when(
            !str_detect(V_Gene, "IGK") | is.na(V_Gene) ~ bestVGene,
            TRUE ~ V_Gene
        )
    ) %>%
    mutate(
        J_Gene = str_remove(str_remove(igblastnJAllele, "[*].*"), ",.*"),
        J_Gene = case_when(
            !str_detect(J_Gene, "IGK") | is.na(J_Gene) ~ bestJGene,
            TRUE ~ J_Gene
        )
    ) %>%
    mutate(C_Gene = bestCGene) %>%
    group_by(sample_id, V_Gene, J_Gene, C_Gene) %>%
    mutate(cloneId = paste0(cloneId, collapse = ",")) %>%
    mutate(totalCloneCount = sum(cloneCount)) %>%
    mutate(totalCloneFraction = sum(cloneFraction)) %>%
    ungroup() %>%
    filter(
        totalCloneCount >= 100,
        totalCloneFraction > 0.1
    ) %>%
    group_by(sample_id) %>%
    slice_max(Productive == "Yes", with_ties = TRUE) %>% 
    slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(
        sample_id,
        Productive,
        mutatedStatus,
        matches("totalClone"),
        matches("clone"),
        matches("_Gene$"),
        everything()
    ) %>%
    distinct() %>%
    full_join(distinct(select(mixcr_IGK, sample_id))) %>%
    left_join(select(
        metadata, 
        sample_id, 
        patient_id, 
        biopsy_id, 
        preservation
    )) %>%
    rowwise() %>%
    mutate(
        across(
            matches("_Gene"),
            ~ ifelse(is.na(.x) | Productive == "N/A", "NULL", .x)
        )
    ) %>%
    mutate(Productive = str_replace(Productive, "N/A", "No")) %>%
    ungroup() %>%
    distinct(sample_id, .keep_all = TRUE) %>% 
    group_by(patient_id, biopsy_id) %>% 
    slice_max(!is.na(totalCloneFraction), n=1, with_ties = FALSE) %>% 
    ungroup() %>%
    mutate(locus = "IGK")

write_tsv(mixcr_IGK_filt, "data/ig_rearrangements/mixcr_IGK_filtered.tsv")


mixcr_IGL <- lapply(metadata$sample_id, load_mixcr, chain = "IGL")
mixcr_IGL <- bind_rows(mixcr_IGL)

mixcr_IGL_filt <- mixcr_IGL %>%
    select(
        sample_id,
        Productive,
        mutatedStatus,
        matches("clone"),
        matches("[VJ]Allele"),
        matches("Gene"),
        matches("Score"),
        -matches("Seq")
    ) %>%
    mutate(
        V_Gene = str_remove(str_remove(igblastnTopVAllele, "[*].*"), ",.*"),
        V_Gene = case_when(
            !str_detect(V_Gene, "IGL") | is.na(V_Gene) ~ bestVGene,
            TRUE ~ V_Gene
        )
    ) %>%
    mutate(
        J_Gene = str_remove(str_remove(igblastnJAllele, "[*].*"), ",.*"),
        J_Gene = case_when(
            !str_detect(J_Gene, "IGL") | is.na(J_Gene) ~ bestJGene,
            TRUE ~ J_Gene
        )
    ) %>%
    mutate(C_Gene = bestCGene) %>%
    group_by(sample_id, V_Gene, J_Gene, C_Gene) %>%
    mutate(cloneId = paste0(cloneId, collapse = ",")) %>%
    mutate(totalCloneCount = sum(cloneCount)) %>%
    mutate(totalCloneFraction = sum(cloneFraction)) %>%
    ungroup() %>%
    filter(
        totalCloneCount >= 100,
        totalCloneFraction > 0.1
    ) %>%
    group_by(sample_id) %>%
    slice_max(Productive == "Yes", with_ties = TRUE) %>%
    slice_max(cloneCount, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(
        sample_id,
        Productive,
        mutatedStatus,
        matches("totalClone"),
        matches("clone"),
        matches("_Gene$"),
        everything()
    ) %>%
    distinct() %>%
    full_join(distinct(select(mixcr_IGL, sample_id))) %>%
    left_join(select(
        metadata, 
        sample_id, 
        patient_id, 
        biopsy_id, 
        preservation
    )) %>%
    rowwise() %>%
    mutate(
        across(
            matches("_Gene"),
            ~ ifelse(is.na(.x) | Productive == "N/A", "NULL", .x)
        )
    ) %>%
    mutate(Productive = str_replace(Productive, "N/A", "No")) %>%
    ungroup() %>%
    distinct(sample_id, .keep_all = TRUE) %>% 
    group_by(patient_id, biopsy_id) %>% 
    slice_max(!is.na(totalCloneFraction), n=1, with_ties = FALSE) %>% 
    ungroup() %>% 
    mutate(locus  = "IGL")

write_tsv(mixcr_IGL_filt, "data/ig_rearrangements/mixcr_IGL_filtered.tsv")

# Combine light and heavy chain data

mixcr_all <- mixcr_IGH_filt %>% 
    bind_rows(mixcr_IGK_filt) %>% 
    bind_rows(mixcr_IGL_filt) %>% 
    select(
        patient_id, 
        biopsy_id, 
        locus,
        Productive, 
        mutatedStatus, 
        totalCloneCount, 
        totalCloneFraction, 
        V_Gene, 
        J_Gene, 
        C_Gene
    ) %>% 
    pivot_wider(
        names_from = "locus", 
        values_from = Productive:C_Gene, 
        names_glue = "{locus}_{.value}"
    )

write_tsv(mixcr_all, "data/ig_rearrangements/mixcr_all_loci.tsv")

