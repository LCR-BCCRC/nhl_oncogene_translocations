library(tidyverse)
library(circlize)
library(data.table)

source("src/R/circlize/prepare_circos_HGBL.R")
source("src/R/circlize/plot_circos_HGBL.R")

colours <- c("chr8" = "grey", 
             "chr2" = "#2E86AB", 
             "chr22" = "#FC5130", 
             "chr14" = "#87BC14")

breakpoints <- read_tsv("data/assembled_breakpoints_with_PMBCL_manual.txt")
breakpoints_IG_MYC_HGBL <- breakpoints %>% 
  filter(method == "CAPTURE" & 
           target == "MYC" & 
           PATH == "HGBL-HG" &
           str_detect(chr_partner, "14|2$|22")) 

genes <- read_tsv("data/transloc_genes.tsv")

genes_ig_myc <- genes %>% 
  filter(str_detect(name, "MYC|CASC8|PVT1|IGH|IGK|IGL")) %>% 
  filter(name != "IGH C") %>% 
  mutate(chrom = ifelse(!str_detect(chrom, "chr"), str_c("chr", chrom), chrom)) %>% 
  mutate(name = ifelse(name == "IGHG1", "IGH C", name), 
         name = ifelse(str_detect(name, "IGH[A-Z]"), "", name), 
         name = ifelse(name == "IGK C", "C", name), 
         name = ifelse(name == "IGL C", "C", name)) %>%
  as.data.table()


circos.clear()
pdf("HGBL_HG_IG_MYC.pdf", height = 8.5, width = 8.5)
circos <- prepare_circos(breakpoints_IG_MYC_HGBL, colours)
plot_circos(circos, genes_ig_myc)
dev.off()

breakpoints_IG_MYC_DLBCL <- breakpoints %>% 
  filter(method == "Chong et al." & 
           target == "MYC" & 
           chr_target == "8" &
           str_detect(chr_partner, "14|^2$|22")) %>% 
  distinct(Patient_ID, .keep_all = TRUE)

circos.clear()
pdf("HGBL_DLBCL_IG_MYC.pdf", height = 8.5, width = 8.5)
circos <- prepare_circos(breakpoints_IG_MYC_DLBCL, colours)
plot_circos(circos, genes_ig_myc)
dev.off()


breakpoints_IG_PBL <- breakpoints %>% 
  filter(method != "Chong et al." & 
           target == "MYC" & 
           chr_target == "8" & 
           PATH == "PBL" & 
           str_detect(chr_partner, "14|^2$|22"))

circos.clear()
pdf("PBL_IG_MYC.pdf", height = 8.5, width = 8.5)
circos <- prepare_circos(breakpoints_IG_PBL, colours)
plot_circos(circos, genes_ig_myc)
dev.off()


