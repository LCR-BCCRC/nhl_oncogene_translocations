source("src/libs.R")
library(tidyverse)
library(GAMBLR)
library(data.table)

# Load project metdata with best sequencing sample
md <- read_tsv("data/metadata/breakpoint_capture_best.tsv")

md %>% count(ICC_class, seq_type)
#    ICC_class     seq_type     n
#  1 BL            capture      7
#  2 BL            genome     219
#  3 COMFL         genome       8
#  4 DLBCL         capture     83
#  5 DLBCL         genome     180
#  6 FL            capture     13
#  7 FL            genome     130
#  8 HGBCL-DH-BCL2 capture    193
#  9 HGBCL-DH-BCL2 genome      34
# 10 HGBCL-DH-BCL6 capture     34
# 11 HGBCL-DH-BCL6 genome       4
# 12 HGBCL-NOS     capture     23
# 13 HGBCL-NOS     genome       8
# 14 PBL           capture     63
# 15 PBL           genome       5

# Load regions bed file
regions <- read_tsv(
    "data/region_data/region_bed_simple.hg38.bed",
    col_names = c("chrom", "start", "end", "name")
)

# Create data.table from regions bed file for foverlaps
regions_bed <- data.table(regions)
colnames(regions_bed) <- c("Chromosome", "Start_Position", "End_Position", "name")
setkey(regions_bed, Chromosome, Start_Position, End_Position)

# Load all GAMBL SSM data for the breakpoint capture space
genome_maf <- get_ssm_by_regions(
    projection = "hg38",
    regions_bed = regions,
    streamlined = FALSE,
    basic_columns = TRUE
) %>%
    filter(Tumor_Sample_Barcode %in% md$sample_id) %>%
    left_join(regions_bed, by = "Chromosome", suffix = c("", "_region")) %>%
    filter(Start_Position >= Start_Position_region, End_Position < End_Position_region) %>% 
    select(-matches("_region"))

# Load capture maf
capture_maf <- read_tsv(
    "/projects/dscott_prj/CCSRI_1500/capture/results/slms_3-1.0_vcf2maf-1.3/level_3/capture--hg38/breakpoint_capture_slms3_merged.maf",
    col_types = spec(genome_maf)
)

# Subset capture maf to variants in capture space
capture_maf <- foverlaps(data.table(capture_maf), regions_bed) %>%
    filter(!is.na(Start_Position)) %>%
    select(-Start_Position, -End_Position) %>%
    rename_with(~ str_remove(.x, "i[.]")) %>%
    select(all_of(colnames(genome_maf))) %>%
    mutate(Tumor_Sample_Barcode = str_remove(Tumor_Sample_Barcode, "--.*")) %>%
    filter(Tumor_Sample_Barcode %in% md[md$seq_type == "capture", ]$sample_id) %>% 
    mutate(Entrez_Gene_Id = as.character(Entrez_Gene_Id))

combined_maf <- bind_rows(capture_maf, genome_maf)

write_tsv(combined_maf, "data/maf/genome_capture.hg38.maf")

# Write individual ICC_class mafs for motif checking

for(class in unique(md$ICC_class)){
    ids_this_class <- filter(md, ICC_class == class) %>% pull(sample_id)
    combined_maf %>%
        filter(Tumor_Sample_Barcode %in% ids_this_class) %>% 
        write_tsv(paste0("data/maf/region_mafs/", class, "_genome_capture.hg38.maf"))
}

# To annotate with WRCY motif overlap:
# src/snvs/CheckMotifMutBias.sh 