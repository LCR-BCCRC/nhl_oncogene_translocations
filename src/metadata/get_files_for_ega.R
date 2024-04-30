library(tidyverse)
library(glue)
library(GAMBLR)

# Columns in the samples table (* Required):
# alias*
# title
# description
# biological_sex*
# subject_id*
# phenotype*
# biosample_id
# case_control
# organism_part
# cell_line

# Columns in the runs table:
# sample	file1	file2
# This table can be customized, and other fields can be added
# file1 and file2 are only needed for fastq
# Otherwise the can just be called file

seq_md <- read_tsv("data/metadata/breakpoint_capture_seq_for_pub.tsv")

sex_md <- read_tsv("data/metadata/breakpoint_capture_biopsies.tsv") %>%
    select(patient_id, sex) %>%
    drop_na() %>%
    distinct()


# Get info for unpublished genomes
genomes <- seq_md %>%
    filter(is.na(publication)) %>%
    filter(seq_type == "genome")

genomes_gambl <- get_gambl_metadata(seq_type_filter = "genome") %>%
    filter(sample_id %in% genomes$sample_id) %>%
    select(sample_id, data_path)

genomes_tidy <- genomes %>%
    select(sample_id, patient_id, ICC_class) %>%
    left_join(genomes_gambl) %>%
    left_join(sex_md) %>%
    mutate(link_name = glue("{sample_id}.genome.grch37.cram")) %>%
    mutate(file = glue("/{link_name}")) %>%
    mutate(alias = sample_id) %>%
    select(
        alias = sample_id,
        biological_sex = sex,
        phenotype = ICC_class,
        subject_id = patient_id,
        data_path,
        link_name,
        file
    )

genomes_tidy %>%
    select(file_name = link_name, data_path) %>%
    write_tsv("data/EGA/genomes_files.tsv")

genomes_tidy %>%
    select(alias, biological_sex, subject_id, phenotype) %>%
    write_tsv("data/EGA/genomes_samples.tsv")

genomes_tidy %>%
    select(
        sample = alias,
        file
    ) %>%
    write_tsv("data/EGA/genomes_runs.tsv")


# Upload all capture data since published Chong data doesn't seem to exist in EGA
captures <- seq_md %>%
    filter(seq_type == "capture")

file_pattern <- "/projects/dscott_prj/CCSRI_1500/bam_archive/capture/merged_bams/{capture_version}--hg38/{sample_id}.merged.cram"

captures_tidy <- captures %>%
    select(sample_id, patient_id, ICC_class, capture_version) %>%
    left_join(sex_md) %>%
    mutate(link_name = glue("{sample_id}.{capture_version}.hg38.cram")) %>%
    mutate(data_path = glue(file_pattern)) %>%
    mutate(available = file.exists(data_path)) %>%
    mutate(file = glue("/{link_name}")) %>%
    mutate(alias = sample_id) %>%
    select(
        alias = sample_id,
        biological_sex = sex,
        phenotype = ICC_class,
        subject_id = patient_id,
        data_path,
        link_name,
        file,
        capture_version
    )

captures_tidy %>%
    select(file_name = link_name, data_path) %>%
    write_tsv("data/EGA/capture_files.tsv")

captures_tidy %>%
    select(alias, biological_sex, subject_id, phenotype) %>%
    write_tsv("data/EGA/capture_samples.tsv")

captures_tidy %>%
    select(
        sample = alias,
        file
    ) %>%
    write_tsv("data/EGA/capture_runs.tsv")

# Upload unpublished RNAseq
rnaseqs <- seq_md %>%
    filter(seq_type == "mrna") %>%
    filter(is.na(publication)) %>%
    filter(protocol == "Ribodepletion")

file_rnaseq <- "/projects/dscott_prj/CCSRI_1500/transcriptomes/data/gsc_bam/{sample_id}.bam"
file_gambl <- "/projects/rmorin/projects/gambl-repos/gambl-lhilton/data/mrna_bams/{sample_id}.hg38.bam"

rnaseqs_tidy <- rnaseqs %>%
    select(sample_id, patient_id, ICC_class, capture_version) %>%
    left_join(sex_md) %>%
    mutate(link_name = glue("{sample_id}.mrna.bam")) %>%
    mutate(data_path = glue(file_rnaseq)) %>%
    mutate(available = file.exists(data_path)) %>%
    mutate(data_path = ifelse(
        !available,
        glue(file_gambl),
        data_path
    )) %>%
    mutate(available = file.exists(data_path)) %>%
    mutate(link_name = glue("{sample_id}.mrna.bam")) %>%
    mutate(file1 = glue("/{sample_id}.mrna.R1.fastq.gz")) %>%
    mutate(file2 = glue("/{sample_id}.mrna.R2.fastq.gz")) %>%
    mutate(alias = sample_id) %>%
    select(
        alias = sample_id,
        biological_sex = sex,
        phenotype = ICC_class,
        subject_id = patient_id,
        data_path,
        link_name,
        file1,
        file2
    )

rnaseqs_tidy %>%
    select(file_name = link_name, data_path) %>%
    write_tsv("data/EGA/rnaseq_files_symlinks.tsv")

rnaseqs_tidy %>%
    select(file_name = file1) %>%
    bind_rows(select(rnaseqs_tidy, file_name = file2)) %>%
    arrange(file_name) %>%
    write_tsv("data/EGA/rnaseq_fastqs.tsv")

rnaseqs_tidy %>%
    select(alias, biological_sex, subject_id, phenotype) %>%
    write_tsv("data/EGA/rnaseqs_samples.tsv")

rnaseqs_tidy %>%
    select(
        sample = alias,
        file1,
        file2
    ) %>%
    write_tsv("data/EGA/rnaseqs_runs.tsv")
