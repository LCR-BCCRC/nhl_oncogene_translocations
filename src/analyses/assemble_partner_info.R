## DEPRECATED ##
## Instead use src/svs/assemble_partner_info_func.R as part of src/metadata/create_metadata.R ##

library(tidyverse)
library(readxl)
source("src/metadata/swap_biopsy_id.R")

# Load tidy metadata table
md <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv")
# Obtain biopsy IDs for WCM samples
wcm <- md %>%
    filter(str_detect(patient_id, "WCM")) %>%
    select(patient_id, biopsy_id_wcm = biopsy_id, sample_id = capture)

# Load SV file and tidy up columns
llmpp_myc_raw <- read_tsv(
    "data/raw_sv_data/myc_annotated_unique_with_bcl2_bcl6_partners.tsv"
)
llmpp_myc <- llmpp_myc_raw %>%
    # Anonymize biopsy ID
    swap_biopsy_id() %>%
    left_join(wcm) %>%
    mutate(biopsy_id = case_when(
        !is.na(biopsy_id_wcm) ~ biopsy_id_wcm,
        str_detect(sample_id, "BLGSP-71-30") ~ str_sub(sample_id, 1, 20),
        str_detect(sample_id, "PBL") ~ sample_id,
        str_detect(sample_id, "^SPECS") ~ sample_id,
        TRUE ~ biopsy_id
    )) %>%
    select(
        matches("sample_id"),
        patient_id,
        biopsy_id,
        seq_type,
        matches("target|partner|status"),
        vaf, dp,
        manta_name, gridss_name
    )

# Join SV data with metadata

sv_md <- left_join(md, llmpp_myc) %>%
    mutate(myc_ba = ifelse(myc_ba == "FAIL", NA, myc_ba)) %>%
    mutate(myc_bp_status = case_when(
        !is.na(myc_partner) & myc_ba == "POS" ~ "TRUE_POS",
        is.na(myc_partner) & myc_ba == "NEG" ~ "TRUE_NEG",
        !is.na(myc_partner) & myc_ba == "NEG" ~ "CRYPTIC_FOUND",
        is.na(myc_partner) & myc_ba == "POS" ~ "FALSE_NEG",
        TRUE ~ "UNKNOWN"
    ))

sv_md %>%
    filter(!is.na(myc_bp_status)) %>%
    count(!is.na(myc_partner), myc_ba, myc_bp_status) %>%
    rename(has_myc_partner = `!is.na(myc_partner)`)
# A tibble: 6 x 4
#   has_myc_partner myc_ba myc_bp_status     n
#   <lgl>           <chr>  <chr>         <int>
# 1 FALSE           NEG    TRUE_NEG        293
# 2 FALSE           POS    FALSE_NEG        97
# 3 FALSE           NA     UNKNOWN          83
# 4 TRUE            NEG    CRYPTIC_FOUND    34
# 5 TRUE            POS    TRUE_POS        395
# 6 TRUE            NA     UNKNOWN         102

# Some cases have "UNKNOWN" for bp_status even though
# they have FISH and genome/capture. Will fix after filling
# out missingness.

# Is this cryptic count accurate?
# Capture samples have been checked - genomes need to be reviewed.
# Most genome bams are in grch37 but these breakpoints have been
# lifted to hg38. Need to get original breakpoint coordinates to check.

gsc_dir <- "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/data/genome_bams/"
# Identify breakpoints by the unique IDs output by SV callers
cryptic_check <- sv_md %>%
    filter(myc_bp_status == "CRYPTIC_FOUND") %>%
    filter(seq_type == "genome", !str_detect(sample_id, "BLGSP")) %>%
    select(sample_id = genome, gridss_name, manta_name)
# Obtain grch37 bedpe
svs_grch37 <- GAMBLR::get_combined_sv(sample_ids = cryptic_check$sample_id, oncogenes = "MYC")
# Join updated coordinates and write out in IGV-friendly format
cryptic_check %>%
    left_join(mutate(svs_grch37, sample_id = tumour_sample_id)) %>%
    rowwise() %>%
    mutate(path = dir(
        gsc_dir,
        pattern = sample_id,
        full.names = TRUE
    )[1]) %>%
    mutate(
        r1 = str_c(CHROM_A, ":", START_A - 100, "-", END_A + 100),
        r2 = str_c(CHROM_B, ":", START_B - 100, "-", END_B + 100)
    ) %>%
    select(path, r1, r2) # %>%
# write_tsv("genomes_for_igv.tsv")

sv_md %>%
    filter(sample_id %in% c(
        "CCS_2042",
        "CCS_1884",
        "SPECS_4953",
        "SPECS_5358",
        "CLC07777"
    )) %>%
    rowwise() %>%
    mutate(path = dir(
        "/projects/dscott_prj/CCSRI_1500/bam_archive/capture/merged_bams/",
        pattern = sample_id,
        full.names = TRUE,
        recursive = TRUE
    )[1]) %>%
    ungroup() %>%
    mutate(
        r1 = str_c(chrom_target, ":", start_target - 100, "-", end_target + 100),
        r2 = str_c(chrom_partner, ":", start_partner - 100, "-", end_partner + 100)
    ) %>%
    select(path, r1, r2) # %>%
# write_tsv("capture_for_igv.tsv")


# After checking these cryptic genome rearrangements, the following
# criteria will be used to determine if SV is truly cryptic to
# BA FISH:
# DLBCL Trio genome (ends in _tumor) -> artifact
# Intrachromosomal -> doesn't meet definition

sv_md_final <- sv_md %>%
    mutate(myc_bp_status = case_when(
        chrom_partner == "chr8" & myc_ba == "NEG" ~ "TRUE_NEG",
        str_detect(genome, "_tumor") & myc_ba == "NEG" ~ "TRUE_NEG",
        TRUE ~ myc_bp_status
    )) %>%
    mutate(across(matches("target|partner"), ~ ifelse(myc_bp_status == "TRUE_NEG", NA, .x))) %>%
    select(-manta_name, -gridss_name)

write_tsv(sv_md_final, "data/sv_data/myc_annotated_breaks_unique_bcl2_bcl6.tsv")

sv_md_final %>%
    filter(!is.na(myc_bp_status)) %>%
    count(!is.na(myc_partner), myc_ba, myc_bp_status) %>%
    rename(has_myc_partner = `!is.na(myc_partner)`)
# A tibble: 6 x 4
#   has_myc_partner myc_ba myc_bp_status     n
# 1 FALSE           NEG    TRUE_NEG        305
# 2 FALSE           POS    FALSE_NEG        97
# 3 FALSE           NA     UNKNOWN          83
# 4 TRUE            NEG    CRYPTIC_FOUND    22
# 5 TRUE            POS    TRUE_POS        395
# 6 TRUE            NA     UNKNOWN         102

sv_md_final %>%
    count(ICC_class, myc_bp_status) %>%
    pivot_wider(
        names_from = myc_bp_status,
        values_from = n,
        values_fill = 0
    )
#   ICC_class     CRYPTIC_FOUND FALSE_NEG TRUE_NEG TRUE_POS UNKNOWN
# 1 BL                        6        41        3       63     113
# 2 COMFL                     0         0        8        0       0
# 3 DLBCL                     5        13      188       51      11
# 4 FL                        1         6       72       12      52
# 5 HGBCL-DH-BCL2             5        27        0      189       0
# 6 HGBCL-DH-BCL6             0         5        0       34       0
# 7 HGBCL-NOS                 1         2       14       13       1
# 8 PBL                       1         3       23       33       8

# Which BLs have apparently cryptic MYC translocations?


sv_md_final %>%
    count(ICC_class, bcl2_bp_status) %>%
    pivot_wider(
        names_from = bcl2_bp_status,
        values_from = n,
        values_fill = 0
    )

#   ICC_class     TRUE_NEG UNKNOWN  `NA` CRYPTIC_FOUND FALSE_NEG TRUE_POS
# 1 BL                  10     172    44             0         0        0
# 2 COMFL                0       8     0             0         0        0
# 3 DLBCL              140      81     0             3         2       42
# 4 FL                   8      70     0             5         2       58
# 5 HGBCL-DH-BCL2        2       0     0             3         9      207
# 6 HGBCL-DH-BCL6       33       4     0             2         0        0
# 7 HGBCL-NOS           25       2     0             1         1        2
# 8 PBL                  8      60     0             0         0        0

sv_md_final %>%
    count(ICC_class, bcl6_bp_status) %>%
    pivot_wider(
        names_from = bcl6_bp_status,
        values_from = n,
        values_fill = 0
    )
