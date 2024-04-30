swap_biopsy_id <- function(df, biopsy_col = "biopsy_id", patient_col = "patient_id") {
    biopsy_key <- "/projects/dscott_prj/CCSRI_1500/shared/metadata/LCR_database_tables/Internal_Biopsies.xlsx"
    biopsy_key <- read_excel(biopsy_key) %>%
        select(Patient_ID, Biopsy_ID = Source_ID, Biopsy_res)%>%
        distinct(Patient_ID, Biopsy_ID, .keep_all = TRUE)
    join_key <- c("Patient_ID", "Biopsy_ID")
    names(join_key) <- c(patient_col, biopsy_col)
    df2 <- df %>%
        left_join(biopsy_key, by = join_key) %>%
        mutate(., !!biopsy_col := ifelse(!is.na(Biopsy_res), Biopsy_res, get(biopsy_col))) %>%
        select(all_of(colnames(df)))
}
