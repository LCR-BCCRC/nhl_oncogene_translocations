library(tidyverse)
library(ggbio)
library(GenomicRanges)
library(AnnotationHub)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(cowplot)

# Load capspace bed file
capspace <- read_tsv(
    "data/region_data/capture_TE99028370--hg38.bed", 
    col_names = c("chr", "start", "end", "name")
) %>% 
mutate(chr = str_remove(chr, "chr")) 

# Collapse regions to one per chromosome and add some padding
regions <- capspace %>% 
    group_by(name) %>% 
    mutate(
        start = min(start), 
        end = max(end)
    ) %>% 
    distinct() %>% 
    mutate(start = ifelse(name == "ZCCHC7", start - 50000, start))

chroms <- unique(regions$chr)

# Convert regions and capspace to granges
region_gr <- makeGRangesFromDataFrame(
    regions, 
    keep.extra.columns = TRUE
)

capspace_gr <- makeGRangesFromDataFrame(
    capspace, 
    keep.extra.columns = TRUE
)

# Load FISH probe coordinates
fish <- read_tsv("data/region_data/fish_probes.hg38.bed") %>%
    separate(name, into = c("name", "colour"), sep = "-") %>% 
    dplyr::filter(name %in% c("Vysis", "Dako", "Metasystems")) %>% 
    mutate(region = case_when(
        chrom == "chr8" ~ "MYC", 
        chrom == "chr3" ~ "BCL6", 
        chrom == "chr18" ~ "BCL2"
    ))

fish_gr <- makeGRangesFromDataFrame(
    fish, 
    keep.extra.columns = TRUE
)

# Create EnsDb from AnnotationHub
ah <- AnnotationHub()
## Get EnsDb for homo sapiens, Ensembl release 105
edb <- ah[["AH98047"]]

# Specify biotypes for plotting
biotypes <- c(
    "IG_C_gene",
    "IG_J_gene",
    "IG_D_gene",
    "IG_V_gene",
    "lncRNA",
    "protein_coding"
)

# Sanity check transcript filters
Tx <- transcripts(
    edb, 
    filter = list(
        GRangesFilter(capspace_gr), 
        AnnotationFilter(~ tx_is_canonical == TRUE), 
        GeneBiotypeFilter(biotypes)
    ), 
    columns = c(listColumns(edb, "tx"), "gene_name", "gene_biotype")
)

# Set desired gene_names to include
protein_coding <- Tx[Tx$gene_biotype == "protein_coding"]$gene_name
protein_coding <- protein_coding[!protein_coding == ""]

lncRNA <- c("PVT1", "CASC11", "CASC8", "CCDC29", "CASC19")

ig_genes <- Tx[str_detect(Tx$gene_biotype, "IG_.*_gene")]$gene_name

gene_names <- c(protein_coding, lncRNA, ig_genes, "PAX5")


plot_regions <- function(region){
    labels <- data.frame(
        capspace_gr[capspace_gr$name == region]
    ) %>% 
    rowwise() %>%
        mutate(x = mean(c(start, end)), y = 0) %>%
        ungroup() %>%
        mutate(label = ifelse(
            width > 500000,
            paste0(round(width / 1e6, digits = 1), " MB"),
            paste0(round(width / 1e3), " kB")
        ))
    
    gene_track <- autoplot(
        edb,
        AnnotationFilterList(
            GRangesFilter(region_gr[region_gr$name == region]),
            AnnotationFilter(~ tx_is_canonical == TRUE),
            GeneNameFilter(gene_names), 
            GeneBiotypeFilter(biotypes)
        ),
        names.expr = "gene_name",
        stat = "identity",
        label.color = "black", 
        color = "darkblue", 
        fill = "darkblue", 
        rect.height = 0.2
    ) +
        geom_rect(
            capspace_gr[capspace_gr$name == region], 
            alpha = 0.2, 
            rect.height = 1, 
            fill = "firebrick" 
        ) + 
        geom_text(data = labels, aes(
            x = x, 
            y = y,
            label = label
        )) +
        theme_cowplot() +
        ggtitle(paste0(region, " (chr", seqnames(region_gr[region_gr$name == region]), ")"))
        
    if(region %in% fish_gr$region){
        fish_region <- fish_gr[fish_gr$region == region]
        palette <- fish_region$colour
        names(palette) <- fish_region$colour
        gene_track <- gene_track +
            geom_rect(
                fish_region,
                rect.height = 0.05,
                aes(
                    fill = colour, 
                    colour = colour,
                    label = name
                )
            ) + 
            scale_fill_manual(
                values = palette
            ) +
            scale_colour_manual(
                values = palette
            ) +
                theme(legend.position = "none")
    }
    
    return(gene_track@ggplot)

}

gene_plots <- lapply(regions$name, plot_regions)

pdf("results/cohort/capspace.pdf", height = 20, width = 10)
cowplot::plot_grid(plotlist = gene_plots, ncol = 1)
dev.off()
