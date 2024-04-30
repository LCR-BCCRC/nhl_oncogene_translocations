library(biomaRt)
library(tidyverse)

# Make GRCh37 MART object
grch37 = useEnsembl(biomart="ensembl", GRCh = 37, dataset = "hsapiens_gene_ensembl")

# Create a biotype filter to extract all IG genes
filter_values <- searchFilterValues(grch37, "transcript_biotype", "IG")
filter_values <- filter_values[str_detect(filter_values, "_gene")]


IG <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcript_biotype", "chromosome_name", "start_position", "end_position"), 
                        filters = "transcript_biotype", values = filter_values, mart = grch37)

# Remove all pseudogenes
IG <- IG %>% 
  filter(chromosome_name == "14" | 
           chromosome_name == "2" | 
           chromosome_name == "22")

# Combine 
IGH <- IG %>% 
  filter(chromosome_name == "14") %>% 
  mutate(transcript_biotype = ifelse(transcript_biotype == "IG_C_gene", "IGH C", "IGH V")) %>%
  group_by(transcript_biotype) %>% 
  summarise(start = min(start_position), 
            end = max(end_position)) %>% 
  rename(name = transcript_biotype) %>% 
  mutate(chrom = "chr14") %>% 
  dplyr::select(chrom, start, end, name)

IGHC <- IG %>% 
  filter(chromosome_name == "14" & 
           transcript_biotype == "IG_C_gene") %>% 
  mutate(chromosome_name = "chr14") %>%
  dplyr::select(chrom = chromosome_name, 
         start = start_position, 
         end = end_position, 
         name = external_gene_name)

IGH <- bind_rows(IGH, IGHC)

IGK <- IG %>% 
  filter(chromosome_name == "2") %>% 
  mutate(transcript_biotype = ifelse(transcript_biotype == "IG_C_gene", "IGL C", "IGL V")) %>% 
  group_by(transcript_biotype) %>% 
  summarise(start = min(start_position), 
            end = max(end_position)) %>% 
  mutate(chromosome_name = "chr2") %>%
  dplyr::select(chrom = chromosome_name, 
                start = start, 
                end = end, 
                name = transcript_biotype)

IGL <- IG %>% 
  filter(chromosome_name == "22") %>% 
  mutate(transcript_biotype = ifelse(transcript_biotype == "IG_C_gene", "IGK C", "IGK V")) %>%
  group_by(transcript_biotype) %>% 
  summarise(start = min(start_position), 
            end = max(end_position)) %>% 
  mutate(chromosome_name = "chr22") %>%
  dplyr::select(chrom = chromosome_name, 
                start = start, 
                end = end, 
                name = transcript_biotype)


# Get coding genes for ZCCHC7, LRMP, SOCS1, and BCL6 partner regions

genes <- c("MYC", "BCL2", "BCL6", "CASC8", "PVT1", "RFTN1", "ZCCHC7", "PAX5", "LPP", "CASC1", "LRMP", "SOCS1")

gene_list <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcript_biotype", "chromosome_name", "start_position", "end_position"), 
      filters = "external_gene_name", values = genes, mart = grch37)

gene_bed <- gene_list %>% 
  distinct(external_gene_name, .keep_all = TRUE) %>%
  dplyr::select(chrom = chromosome_name, 
         start = start_position, 
         end = end_position, 
         name = external_gene_name) %>% 
  mutate(chrom = str_c("chr", chrom))

bed_all <- bind_rows(gene_bed, IGH, IGK, IGL)


write_tsv(bed_all, "data/transloc_genes.tsv")
