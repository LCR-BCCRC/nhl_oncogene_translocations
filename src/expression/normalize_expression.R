library(tidyverse)
library(DESeq2)
library(limma)
library(BiocParallel)
library(biomaRt)

# Set up parallelization

register(MulticoreParam(64))

# Load the metadata and set paths to matrices

md <- read_tsv("data/metadata/breakpoint_capture_md.tsv") %>%
    filter(seq_type == "mrna") %>%
    distinct() %>%
    filter(preservation != "sorted_cells")

gambl_matrix <-
    "/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/gambl/salmon-1.0/99-outputs/quant_to_hg38_matrix/mrna/salmon.genes.counts.tsv"

llmpp_matrix <-
    "/projects/dscott_prj/CCSRI_1500/transcriptomes/results/salmon-1.1/99-outputs/quant_to_hg38_matrix/mrna/salmon.genes.counts.tsv"


# Filter the matrix to samples meeting the minimum read threshold
# of 2M reads
matrix <- lapply(c(gambl_matrix, llmpp_matrix), function(x) {
    read.delim(
        x,
        row.names = 1,
        check.names = FALSE,
        sep = "\t"
    )
})
matrix <- bind_cols(matrix)
matrix <- matrix[, colnames(matrix) %in% md$sample_id]
matrix_filt <- matrix[, colSums(matrix) > 2e6]
dropped <- colnames(matrix)[!colnames(matrix) %in% colnames(matrix_filt)]
dropped <- colSums(matrix[, dropped, drop = FALSE])

# Manipulate metadata to match matrix
meta <- md %>%
    dplyr::select(sample_id, ICC_class, preservation, protocol) %>%
    column_to_rownames("sample_id")

meta <- meta[rownames(meta) %in% colnames(matrix_filt), ]

# Run DESeq2

matrix_filt <- as.matrix(matrix_filt[, colnames(matrix_filt) %in% rownames(meta)])
meta <- meta[colnames(matrix_filt), ]
identical(colnames(matrix_filt), rownames(meta))

dds <- DESeqDataSetFromMatrix(matrix_filt,
    colData = meta,
    design = ~1
)

dds


dds <- DESeq(dds,
    parallel = TRUE
)

# Perform variance stabilzing transformation
vsd <- vst(dds, blind = FALSE)

# Plot PCA for each uncorrected batch

bio_var <- "ICC_class"
final_batches <- c("preservation", "protocol")

for (batch in c(final_batches, bio_var)) {
    plot_path <- file.path("results/expression", paste0("PCA_uncorrected_", batch, ".pdf"))
    pdf(plot_path, height = 10, width = 10)
    print(plotPCA(vsd, intgroup = batch))
    dev.off()
}

mat <- assay(vsd)

# Perform batch effect correction on each batch

vsd_corr <- vsd

batch_matrix <- model.matrix(
    as.formula(paste("~", paste(final_batches, collapse = " + "))),
    meta
)
design_matrix <- model.matrix(
    as.formula(paste("~", paste(bio_var, collapse = " + "))),
    meta
)

mat <- removeBatchEffect(
    mat,
    covariates = batch_matrix,
    design = design_matrix
)


assay(vsd_corr) <- mat

# Plot PCA for each corrected batch plus cohort
for (batch in c(final_batches, bio_var)) {
    plot_path <- file.path("results/expression", paste0("PCA_corrected_", batch, ".pdf"))
    pdf(plot_path, height = 10, width = 10)
    print(plotPCA(vsd_corr, intgroup = batch))
    dev.off()
}

# Obtain human-readable gene IDs from biomaRt
mat <- data.frame(mat, check.names = FALSE)

mat_renamed <- mat
mat_renamed$gene_id <- rownames(mat)

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gencode_to_symbol <- data.frame(
    gene_id = rownames(mat),
    ensembl_gene_id = gsub("[.].*", "", rownames(mat)),
    stringsAsFactors = FALSE
)

gene_list <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    values = (gencode_to_symbol$ensembl_gene_id), mart = ensembl
)

gencode_to_symbol <- merge(
    gencode_to_symbol,
    gene_list
)

gencode_to_symbol <- gencode_to_symbol[!duplicated(gencode_to_symbol$gene_id), ]

mat_renamed <- merge(
    gencode_to_symbol,
    mat_renamed,
    by = "gene_id"
)

write.table(
    mat_renamed,
    "results/expression/vst_matrix.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)

# Also write an uncorrected vst matrix for PCA analysis
mat_uncor <- data.frame(assay(vsd), check.names = FALSE)
mat_uncor$gene_id <- rownames(mat_uncor)

mat_uncor <- merge(
    gencode_to_symbol,
    mat_uncor,
    by = "gene_id"
)

write.table(
    mat_uncor,
    "results/expression/uncorrected_vst_matrix.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)

# Select key genes for a pared down expression matrix

# mat_renamed <- fread("results/expression/vst_matrix.tsv")
d90_dzsig <- read_tsv("data/region_data/dlbcl90_dzsig_nanostring_genes.tsv")
mat_min <- mat_renamed[
    mat_renamed$hgnc_symbol %in% c(
        "MYC", 
        "BCL2", 
        "BCL6", 
        "APEX1", 
        "APEX2", 
        "POLH", 
        "AICDA", 
        "UNG", 
        "TET2", 
        "STAT6", 
        "KMT2D", 
        "ARID1A", 
        d90_dzsig$hgnc_symbol
    ) |
    str_detect(mat_renamed$hgnc_symbol, "IGH[MDE]$|IGH[GA][:digit:]$"), 
]

mat_tidy <- mat_min %>%
    pivot_longer(
        -c(gene_id, ensembl_gene_id, hgnc_symbol),
        names_to = "sample_id",
        values_to = "expression"
    )

write_tsv(mat_tidy, "results/expression/myc_bcl2_bcl6_expr.tidy.tsv")
