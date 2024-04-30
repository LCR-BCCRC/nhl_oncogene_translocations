# Load raw untidy FlowJo gated results tables and create tidy/anonymized tables

source("src/libs.R")
library(tidyverse)
library(readxl)
source("src/metadata/swap_biopsy_id.R")

# Load original 60 samples all gated with the same method
flow_results_raw <- read_csv("data/flow_data/flowjo_table_raw.csv")

flow_key <- read_tsv("data/flow_data/flow_id_key.tsv")

flow_key <- swap_biopsy_id(flow_key)

flow_tidy1 <- flow_results_raw %>% 
    mutate(ID = str_remove_all(ID, " [Tt]ube.*")) %>% 
    left_join(flow_key, by = c("ID" = "Filename_trim")) %>% 
    select(
        patient_id, 
        biopsy_id, 
        flow_id = ID, 
        parent = `45side_19pos_3neg_freq_of_parent_percent`, 
        scatter = `45side_freq_of_parent_percent`,
        matches("Lambda")
    ) %>% 
    mutate(gate = "19+3-", scatter_method = "45side") %>% 
    distinct(patient_id, biopsy_id, .keep_all = TRUE)

# Load older flow data with variable assays and gating strategies

flow_raw2 <- read_tsv("data/flow_data/flowjo_table_raw2.tsv")

flow_tidy2 <- flow_raw2 %>% 
    mutate(ID = str_remove_all(`...1`, " [Tt]ube.*")) %>% 
    left_join(
        mutate(flow_key, Filename_trim = str_remove(Filename_trim, "^B")), 
        by = c("ID" = "Filename_trim")
    ) %>% 
    select(-`...1`) %>% 
    pivot_longer(
        -c(patient_id, biopsy_id, scatter, ID), 
        names_to = "gate", 
        values_to = "percent"
    ) %>% 
    mutate(gate = str_remove_all(gate, "[.].*")) %>% 
    separate(gate, into = c("gate", "quadrant"), sep = ": ", fill = "right") %>% 
    mutate(quadrant = replace_na(quadrant, "parent")) %>% 
    distinct() %>% 
    mutate(gate = str_remove(gate, "/Q[1234]")) %>% 
    pivot_wider(
        names_from = "quadrant", 
        values_from = "percent"
    ) %>% 
    filter(!is.na(parent)) %>% 
    group_by(patient_id, biopsy_id) %>% 
    slice_max(!is.na(`polyLambda- , polyKappa+`), n=1, with_ties = TRUE) %>% 
    slice_max(parent, n=1, with_ties = FALSE) %>% 
    ungroup() %>% 
    mutate(scatter_method = "scatter") %>% 
    rename_all(~str_remove_all(.x, "poly| , ")) %>% 
    rename(flow_id = ID)
    
flow_manual <- read_tsv("data/flow_data/25-Sep-2023-flow-manual-gate.tsv")

censor <- flow_manual %>% 
    filter(
        
    )

flow_manual <- read_tsv("data/flow_data/25-Sep-2023-flow-manual-gate.tsv")

censor <- flow_manual %>% 
    filter(final_sIg %in% c("censor", "indeterminant")) %>% 
    pull(flow_id)

flow_results <- flow_tidy2 %>% 
    bind_rows(flow_tidy1) %>% 
    # Remove cases censored after manual review for low B cell count
    filter(!flow_id %in% censor) %>% 
    filter(parent >= 25) %>% 
    select( 
        patient_id, 
        biopsy_id, 
        flow_id,
        scatter_method, 
        scatter_percent = scatter, 
        parent_gate = gate, 
        parent_percent = parent, 
        matches("Lambda")
    )

write_tsv(select(flow_results, -flow_id), "data/flow_data/kappa_lambda_flow.tsv")


flow_results_cat <- flow_results %>% 
    pivot_longer(
        matches("Lambda"),
        names_to = "surface_IG",
        values_to = "percent"
    ) %>%
    group_by(biopsy_id) %>%
    slice_max(percent) %>%
    ungroup() %>% 
    left_join(flow_manual) %>% 
    mutate(surface_IG = case_when(
        final_sIg == "lambda" ~ "Lambda+Kappa-", 
        final_sIg == "kappa" ~ "Lambda-Kappa+", 
        final_sIg == "null" ~ "Lambda-Kappa-", 
        TRUE ~ surface_IG
    )) %>% 
    select(-final_sIg, -flow_id)

write_tsv(flow_results_cat, "data/flow_data/kappa_lambda_flow_categorized.tsv")
