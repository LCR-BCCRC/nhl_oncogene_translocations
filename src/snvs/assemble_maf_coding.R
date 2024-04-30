library(tidyverse)
library(GAMBLR)
library(data.table)
library(readxl)

source("src/metadata/swap_biopsy_id.R")
source("src/libs.R")

# Load project metdata
md <- read_tsv("data/metadata/breakpoint_capture_md.tsv")

# Subset to one of genome or capture per biopsy
# prioritizing genome

md <- md %>%
    filter(
        seq_type %in% c("genome", "capture"),
        tissue_status == "tumour"
    ) %>%
    group_by(patient_id, biopsy_id) %>%
    slice_max(seq_type == "genome", n=1, with_ties = FALSE) %>%
    ungroup() %>%
    select(patient_id, biopsy_id, ICC_class, seq_type)

md %>% count(ICC_class, seq_type)
#    ICC_class     seq_type     n
#  1 BL            capture      7
#  2 BL            genome     219
#  3 COMFL         genome       8
#  4 DLBCL         capture     86
#  5 DLBCL         genome     182
#  6 FL            capture     13
#  7 FL            genome     130
#  8 HGBCL-DH-BCL2 capture    190
#  9 HGBCL-DH-BCL2 genome      31
# 10 HGBCL-DH-BCL6 capture     34
# 11 HGBCL-DH-BCL6 genome       5
# 12 HGBCL-NOS     capture     23
# 13 HGBCL-NOS     genome       8
# 14 PBL           capture     63
# 15 PBL           genome       5

# Obtain available exome data

md_exomes <- read_tsv("../CCSRI_1500/exomes/data/metadata/LLMPP_exomes_samples.tsv") %>%
    filter(tissue_status == "tumour") %>%
    swap_biopsy_id() %>%
    mutate(biopsy_id = case_when(
        str_detect(sample_id, "PBL") ~ patient_id,
        str_detect(seq_sample_id, "^EX") & !is.na(llmpp_specs_id) ~ llmpp_specs_id,
        TRUE ~ biopsy_id
    )) %>%
    select(patient_id, biopsy_id, Tumor_Sample_Barcode = sample_id)

md_exomes <- md %>%
    filter(seq_type == "capture") %>%
    left_join(md_exomes) %>%
    filter(!is.na(Tumor_Sample_Barcode))

# Load regions bed file
# ARID1A chr1:27,001,002-27,130,121
# KMT2D chr12:49,412,758-49,449,107
# MEF2B chr19:19,256,376-19,281,098
# TNFRSF14 chr1:2,487,805-2,495,267
# EZH2 chr7:148,504,464-148,581,441
# ATP6V1B2 chr8:20,054,704-20,079,207
# CREBBP chr16:3,775,056-3,930,121
# STAT6 chr12:57,489,187-57,505,196
regions <- data.frame(
    "chrom" = c("chr1", "chr12", "chr19", "chr1", "chr7", "chr8", "chr16", "chr12"),
    "start" = c(27000000, 49412000, 19250000, 2487000, 148504000, 20054000, 3775000, 57489187),
    "end" = c(27130000, 49450000, 19282000, 2496000, 148582000, 20080000, 3931000, 57505196),
    "name" = c("ARID1A", "KMT2D", "MEF2B", "TNFRSF14", "EZH2*", "ATP6V1B2", "CREBBP*", "STAT6")
)

# Load all GAMBL SSM data for the genes of interest

gambl_md <- get_gambl_metadata(seq_type_filter = "genome") %>%
    filter(patient_id %in% md$patient_id) %>%
    swap_biopsy_id() %>%
    select(patient_id, biopsy_id, Tumor_Sample_Barcode = sample_id)

genome_md <- gambl_md %>%
    inner_join(filter(md, seq_type == "genome")) %>%
    distinct()

genome_maf <- get_ssm_by_regions(
    projection = "grch37",
    regions_bed = regions,
    streamlined = FALSE,
    basic_columns = TRUE
) %>%
    filter(Tumor_Sample_Barcode %in% genome_md$Tumor_Sample_Barcode) %>%
    filter(Variant_Classification %in% GAMBLR:::coding_vc) %>%
    mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "EZH2" & str_detect(HGVSp_Short, "Y646") ~ "EZH2*",
        Hugo_Symbol == "CREBBP" &
            Start_Position > 3785000 &
            End_Position < 3791000 &
            Variant_Classification == "Missense_Mutation" ~ "CREBBP*",
        TRUE ~ Hugo_Symbol
    )) %>% 
    select(1:45, -Entrez_Gene_Id)


# Load exome maf
exome_maf <- read_tsv(
    "/projects/dscott_prj/CCSRI_1500/exomes/results/slms3_with_unmatchedNormal/vcf2maf_noncanonical/99-outputs/deblacklisted/exomes.maf",
    col_types = spec(genome_maf)
)


# Subset exome maf to genes of interest
exome_maf_filt <- exome_maf %>%
    filter(Hugo_Symbol %in% str_remove(regions$name, "[*]")) %>%
    filter(Variant_Classification %in% GAMBLR:::coding_vc) %>%
    filter(Tumor_Sample_Barcode %in% md_exomes$Tumor_Sample_Barcode) %>%
    mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "EZH2" & str_detect(HGVSp_Short, "Y646") ~ "EZH2*",
        Hugo_Symbol == "CREBBP" &
            Start_Position > 3785000 &
            End_Position < 3791000 &
            Variant_Classification == "Missense_Mutation" ~ "CREBBP*",
        TRUE ~ Hugo_Symbol
    )) %>% 
    select(1:45, -Entrez_Gene_Id)


combined_maf <- bind_rows(exome_maf_filt, genome_maf)

write_tsv(combined_maf, "data/maf/genome_capture.ARID1A_KMT2D.maf")

combined_md <- bind_rows(md_exomes, genome_md)


fish_test <- prettyForestPlot(
    maf = combined_maf,
    metadata = combined_md,
    genes = regions$name,
    rm_na_samples = FALSE,
    comparison_column = "ICC_class",
    comparison_values = c("HGBCL-DH-BCL2", "FL"),
    custom_colours = colours$group[c("FL", "HGBCL-DH-BCL2")]
)

fish_test$arranged

# Load and tidy MYC partner data
myc_data <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv") 

md_myc <- combined_md %>%
    left_join(myc_data) %>%
    filter(!is.na(myc_partner_group)) %>%
    filter(ICC_class == "HGBCL-DH-BCL2") %>%
    distinct(patient_id, biopsy_id, .keep_all = TRUE)

fish_test <- prettyForestPlot(
    maf = combined_maf,
    metadata = md_myc,
    genes = c(regions$name, "CREBBP"),
    comparison_column = "myc_partner_igh_binary",
    comparison_values = c("IGH", "non-IGH")
)

fish_test$arranged

fish_test$fisher

ggsave("results/ig_rearrangement/DH_MYC_partner_vs_ARID1A_KMT2D.pdf")
