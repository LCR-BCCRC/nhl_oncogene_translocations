library(tidyverse)

md <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    filter(seq_type == "mrna")
    

mat_tidy <- read_tsv("results/expression/myc_bcl2_bcl6_expr.tidy.tsv")

multi <- md %>%
    group_by(biopsy_id) %>%
    filter(n() > 1) %>%
    ungroup()  %>% 
    left_join(mat_tidy) %>% 
    dplyr::select(
        biopsy_id,
        hgnc_symbol, 
        expression, 
        preservation, 
        protocol
    )  

multi %>%
    mutate(preservation = ifelse(preservation == "FFPE", preservation, "Frozen")) %>% 
    distinct(
        biopsy_id,
        hgnc_symbol, 
        preservation, 
        protocol, 
        .keep_all = TRUE
    ) %>% 
    group_by(biopsy_id, hgnc_symbol, preservation) %>% 
    slice_head(n=1) %>% 
    ungroup() %>% 
    pivot_wider(
        names_from = preservation, 
        values_from = c(expression, protocol), 
        names_glue = "{preservation}_{.value}"
    ) %>% 
    drop_na() %>% 
    mutate(same_protocol = case_when(
        Frozen_protocol == "PolyA" & Frozen_protocol == FFPE_protocol ~ "both_PolyA", 
        Frozen_protocol == "Ribodepletion" & Frozen_protocol == FFPE_protocol ~ "both_Ribodepletion", 
        TRUE ~ "different_protocol"
    )) %>% 
    ggplot(aes(
        x = Frozen_expression, 
        y = FFPE_expression, 
        colour = same_protocol
    )) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    ggpubr::stat_cor() +
    facet_wrap(~ hgnc_symbol, scales = "free")
    
ggsave("tmp/PolyA_vs_Ribodepletion_expression.pdf", height = 20, width = 23)


multi %>%
    mutate(preservation = ifelse(preservation == "FFPE", preservation, "Frozen")) %>%
    distinct(
        biopsy_id,
        hgnc_symbol,
        protocol,
        .keep_all = TRUE
    ) %>% 
    group_by(biopsy_id, hgnc_symbol, protocol) %>% 
    slice_head(n=1) %>% 
    ungroup() %>%
    pivot_wider(
        names_from = protocol,
        values_from = c(expression, preservation),
        names_glue = "{protocol}_{.value}"
    ) %>%
    drop_na() %>%
    # filter(hgnc_symbol == "AICDA") %>%
    # dplyr::count(PolyA_preservation,Ribodepletion_preservation)
    mutate(same_preservation = case_when(
        PolyA_preservation == "FFPE" & Ribodepletion_preservation == PolyA_preservation ~ "both_FFPE",
        PolyA_preservation == "Frozen" & Ribodepletion_preservation == PolyA_preservation ~ "both_Frozen",
        TRUE ~ "different_preservation"
    )) %>% 
    ggplot(aes(
        x = PolyA_expression,
        y = Ribodepletion_expression,
        colour = same_preservation
    )) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    ggpubr::stat_cor() +
    facet_wrap(~hgnc_symbol, scales = "free") 

ggsave("tmp/FFPE_vs_Frozen_expression.pdf", , height = 20, width = 23)
