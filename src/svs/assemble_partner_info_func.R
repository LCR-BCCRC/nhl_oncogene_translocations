library(tidyverse)
library(readxl)
source("src/metadata/swap_biopsy_id.R")
source("src/libs.R")

get_partners <- function(target, metadata){
    
        md <- metadata %>% 
        select(
            patient_id, 
            biopsy_id, 
            capture, 
            genome, 
            ICC_class, 
            ba  = str_c(target, "_ba")
        ) 
        
    # Obtain biopsy IDs for WCM samples
    wcm <- md %>%
        filter(str_detect(patient_id, "WCM")) %>%
        select(patient_id, biopsy_id_wcm = biopsy_id, sample_id = capture)

    # Load SV file and tidy up columns
    
    breakpoint_files <- list()
    breakpoint_files$myc <- "data/raw_sv_data/myc_annotated_unique_with_bcl2_bcl6_partners.tsv"
    breakpoint_files$bcl2 <- "data/raw_sv_data/bcl2_annotated_breaks_unique.tsv"
    breakpoint_files$bcl6 <- "data/raw_sv_data/bcl6_annotated_breaks_unique.tsv"
    
    llmpp_raw <- read_tsv(breakpoint_files[[target]])

    llmpp_bp <- llmpp_raw %>%
        filter(
            !(annotation_target == "BCL2" & annotation_partner == "IGLL4P"), 
            !(annotation_target == "BCL2" & who_class %in% c("BL"))
        ) %>% 
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
            patient_id,
            biopsy_id,
            matches("target|partner|status"),
            vaf, dp,
            manta_name, gridss_name
        ) %>% 
        mutate(
            across(matches("_name"), 
            ~ ifelse(is.na(.x), FALSE, TRUE)
        )) %>% 
        rename_with(~ str_replace(.x, "_name", "_called")) %>% 
        rename_with(~ str_remove(.x, paste0(target, "_"))) %>% 
        select(-matches("bcl[26]_")) %>% 
        group_by(patient_id, biopsy_id) %>%
        slice_max(!is.na(chrom_target)) %>% 
        ungroup()

    # Join SV data with metadata

    sv_md <- left_join(md, llmpp_bp) %>%
        mutate(ba = ifelse(ba == "FAIL", NA, ba)) %>%
        mutate(bp_status = case_when(
            !is.na(partner) & ba == "POS" ~ "TRUE_POS",
            is.na(partner) & ba == "NEG" ~ "TRUE_NEG",
            !is.na(partner) & ba == "NEG" ~ "CRYPTIC_FOUND",
            is.na(partner) & ba == "POS" ~ "FALSE_NEG",
            TRUE ~ "UNKNOWN"
        )) %>% 
        distinct(patient_id, biopsy_id, .keep_all = TRUE) %>% 
        mutate(across(
            c(partner, annotation_partner), 
            ~ ifelse(.x == "LRMP", "IRAG2", .x)
        )) %>% 
        mutate(across(
            matches("target|partner|vaf|dp"), 
            ~ ifelse(is.na(partner), NA, .x)
        )) %>% 
        mutate(across(
            matches("called"), 
            ~ ifelse(is.na(partner), FALSE, .x)
        ))

    if(target == "myc"){
        
        sv_md <- sv_md %>%
            mutate(bp_status = case_when(
                # Intrachromosomal SVs on background of FISH neg are true negative
                chrom_partner == "chr8" & ba == "NEG" ~ "TRUE_NEG",
                # Cryptic SVs in FFPE tumours from LSARP_Trios were always artifacts
                str_detect(genome, "_tumor") & ba == "NEG" ~ "TRUE_NEG",
                patient_id == "SPECS_5334" ~ "TRUE_NEG",
                TRUE ~ bp_status
            )) %>%
            mutate(
                across(matches("target|partner|vaf|dp"), 
                ~ ifelse(bp_status == "TRUE_NEG", NA, .x))
            )  %>% 
            mutate(
                across(matches("manta|gridss"), 
                ~ ifelse(bp_status == "TRUE_NEG", FALSE, .x))
            )
    }
    
    sv_md <- sv_md %>% 
        mutate(across(
            matches("capture|genome"), 
            ~ ifelse(is.na(.x), "FALSE", "TRUE")
        )) 
    
    write_tsv(sv_md, glue("data/sv_data/{target}_annotated_breaks_unique.tsv"))
    
    sv_simplified <- sv_md %>% 
        select(patient_id, biopsy_id, ba, bp_status, partner, annotation_partner) %>%
        mutate(partner_group = case_when(
            partner %in% c("BCL6", "PAX5") ~ "recurrent non-IG",
            partner %in% c("IGK", "IGL") ~ "IGK/L",
            partner %in% c("IGH", "MYC") ~ partner,
            !is.na(partner) ~ "non-IG"
        )) %>%
        mutate(partner_igh_binary = case_when(
            partner == "IGH" ~ "IGH",
            !is.na(partner) ~ "non-IGH"
        )) %>%
        mutate(partner_igh_discrete = case_when(
            str_detect(annotation_partner, "IGH[VDJ][:digit:]+-") ~ "Variable",
            annotation_partner %in% names(colours$ighc) ~ annotation_partner
        )) %>% 
        rename_with(~ glue("{target}_{.x}"), matches("partner|ba|status"))
    
    return(sv_simplified)

}
