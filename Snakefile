import pandas as pd
import oncopipe as op

genes = ["myc", "bcl2", "bcl6"]
ICC = ["BL", "DLBCL", "HGBCL-DH-BCL2", "FL", "HGBCL-DH-BCL6"]
best = "data/metadata/breakpoint_capture_best.tsv"
biopsies = "data/metadata/breakpoint_capture_biopsies.tsv"
samples = "data/metadata/breakpoint_capture_md.tsv"
svs = expand("data/sv_data/{gene}_annotated_breaks_unique.tsv", gene = genes)
expression = "results/expression/myc_bcl2_bcl6_expr.tidy.tsv"


# Cohort summary
rule plot_case_counts: 
    input: 
        sv = svs, 
        biopsies = biopsies
    output: 
        seq_type = "results/cohort/seq_type.pdf", 
        recall = "results/cohort/sv_recall.pdf", 
        recall_seq_type = "results/cohort/sv_recall_seq_type.pdf"
    params: 
        script = "src/metadata/plot_case_counts.R"
    script: 
        "{params.script}"

# MAF data from genomes and breakpoint capture        
rule generate_maf: 
    input: 
        capture_maf = "/projects/dscott_prj/CCSRI_1500/capture/results/slms_3-1.0_vcf2maf-1.3/level_3/capture--hg38/breakpoint_capture_slms3_merged.maf"
    output: 
        maf = "data/maf/genome_capture.hg38.maf"
    params: 
        script = "src/snvs/assemble_maf.R"
    script: 
        "{params.script}"

rule wrcy_maf: 
    input: 
        maf = rules.generate_maf.output.maf
    output: 
        wrcy_maf = "data/maf/genome_capture.WRCY.hg38.maf"
    params: 
        script = "src/snvs/CheckMotifMutBias.sh"
    shell: 
        "{params.script}"

# Circos plots
rule circos: 
    input: 
        best = best, 
        biopsies = biopsies,
        svs = svs, 
        maf = rules.wrcy_maf.output.wrcy_maf
    output: 
        expand("results/sv_circos/{ICC}.all.pdf", ICC = ICC)
    params: 
        script = "src/circos/generate_circos_plots.R"
    script: 
        "{params.script}"
        
# Cryptic BCL2 and mutation data

rule cryptic_bcl2: 
    input: 
        best = best, 
        biopsies = biopsies, 
        maf_data = str(rules.wrcy_maf.output.wrcy_maf),
        svs = svs[2],  
        regions = "data/region_data/capture_TE99028370--hg38.bed"
    output: 
        "data/maf/bcl2_tss_mut_counts.tsv", 
        "results/bcl2_svs/bcl2_R_vs_mut_accuracy_TSS.tsv", 
        "results/bcl2_svs/bcl2_R_vs_FISH_accuracy.tsv", 
        "results/bcl2_svs/bcl2_tx_location_partner_status.tsv", 
        "results/bcl2_svs/bcl2_bp_status_by_ICC_class.tsv", 
        "results/bcl2_svs/bcl2_ihc_vs_mut_count.pdf"
    params: 
        script = "src/svs/cryptic_BCL2_svs.R", 
    script: 
        "{params.script}"
        
# Additional MYC analyses

rule myc_svs: 
    input: 
        best = best, 
        svs = svs[1], 
        biopsies = biopsies, 
        maf = rules.wrcy_maf.output.wrcy_maf, 
        regions = "data/region_data/capture_TE99028370--hg38.bed"
    output: 
        "results/myc_svs/myc_break_by_partner_hist.pdf", 
        "results/myc_svs/myc_break_by_IGH_mechanism.pdf", 
        "results/myc_svs/myc_mutations_by_partner.pdf", 
        "results/myc_svs/myc_mutations_by_IGH_partner_binary.pdf"
    params: 
        script = "src/svs/myc_svs.R"
    script: 
        "{params.script}"

# Gene expression        
rule normalize_expression: 
    output: 
        expression = expression
    params: 
        script = "src/expression/normalize_expression.R"
    script: 
        "{params.script}"
        
rule plot_expression: 
    input: 
        samples = samples, 
        biopsies = biopsies, 
        expression = expression
    output: 
        "results/expression/myc_expr_by_partner.pdf", 
        "results/expression/myc_expr_by_partner_DH_vs_DLBCL_neg.pdf", 
        "results/expression/myc_expr_by_igh_partner.pdf", 
        "results/expression/bcl2_expr_by_igh_partner.pdf", 
        "results/expression/bcl6_expr_by_igh_partner.pdf", 
        "results/expression/bcl6_expr_by_igh_partner_coo.pdf", 
        "results/expression/aid_related_expression.pdf", 
        "results/expression/aid_related_dh_by_myc_partner.pdf", 
        "results/expression/IGH_expression_by_ICC_class.pdf", 
        "results/expression/IGH_expression_dh_by_myc_partner.pdf"
    params: 
        script = "src/expression/plot_expression.R", 
    script: 
        "{params.script}"

# aSHM and CSR       
rule ashm_heatmaps_IGH: 
    input: 
        best = best, 
        biopsies = biopsies,
        svs = svs, 
        maf = rules.wrcy_maf.output.wrcy_maf
    output: 
        "results/shm_heatmaps/IGHC_shm_heatmap.pdf", 
        "results/shm_heatmaps/count_shm_ighc.pdf", 
        "results/shm_heatmaps/count_shm_ighc_dh_bcl2_igh.pdf", 
    params: 
        script = "src/snvs/generate_ashm_heatmaps.R"
    script: 
        "{params.script}"
        
rule ashm_heatmaps_all: 
    input: 
        best = best, 
        maf = rules.generate_maf.output.maf, 
        regions = "data/region_data/regions_for_mutsig.tsv"
    output: 
        "results/shm_heatmaps/mutation_counts_per_region.pdf", 
        "results/shm_heatmaps/mutation_counts_per_region_dzsig.pdf"
    params: 
        script = "src/snvs/mutation_count_per_ICC.R"
    script: 
        "{params.script}"
        
rule igh_csr_vs_shm: 
    input: 
        best = best, 
        biopsies = biopsies, 
        maf = rules.wrcy_maf.output.wrcy_maf, 
        svs = svs, 
        panel_regions = "data/region_data/capture_TE99028370--hg38.bed"
    output: 
        "results/shm_heatmaps/Emu_muts_vs_CSR_svs.pdf"
    params: 
        script = "src/svs/all_igh_svs.R"
    script: 
        "{params.script}"
        

# MiXCR and flow cytometry        
rule process_mixcr: 
    input: 
        samples = samples, 
        biopsies = biopsies
    output: 
        IGH = "data/ig_rearrangements/mixcr_IGH_filtered.tsv", 
        IGK = "data/ig_rearrangements/mixcr_IGK_filtered.tsv", 
        IGL = "data/ig_rearrangements/mixcr_IGL_filtered.tsv", 
        ALL = "data/ig_rearrangements/mixcr_all_loci.tsv"
    params: 
        script = "src/ig_rearrangement/process_mixcr.R"
    script: 
        "{params.script}"
        
rule mixcr_vs_flow: 
    input: 
        samples = samples, 
        mixcr = str(rules.process_mixcr.output.ALL), 
        flow = "data/flow_data/kappa_lambda_flow_categorized.tsv"
    output: 
        "data/ig_rearrangements/mixcr_vs_flow.tsv", 
        "results/ig_rearrangement/flow_sIg_mixcr.pdf", 
        "results/ig_rearrangement/flow_sIg_mixcr_counts.pdf", 
        "results/ig_rearrangement/flow_sIg_mixcr_preservation.pdf", 
        "results/ig_rearrangement/flow_sIg_by_group.pdf", 
        "results/ig_rearrangement/sIg_NULL_barplot.pdf"
    params: 
        script = "src/ig_rearrangement/mixcr_vs_flow.R"
    script: 
        "{params.script}"
        
rule plot_mixcr: 
    input: 
        biopsies = biopsies, 
        IGH = str(rules.process_mixcr.output.IGH)
    output: 
        "results/ig_rearrangement/mixcr_glm_IGH_null_vs_ffpe_MYC_partner.tsv", 
        "results/ig_rearrangement/DH_CSR_vs_MYC_partner_IGH_non-IGH.fish_test.tsv", 
        "results/ig_rearrangement/mixcr_glm_IGH_null_vs_ffpe_group.tsv", 
        "results/ig_rearrangement/mixcr_DH_BCL2_MYC_IG_non-IG.pdf", 
        "results/ig_rearrangement/mixcr_C_Gene_vs_MYC_Partner.pdf", 
        "results/ig_rearrangement/mixcr_C_Gene_vs_MYC_Partner_sankey.pdf"
    params: 
        script = "src/ig_rearrangement/plot_mixcr.R"
    script: 
        "{params.script}"

        

rule all: 
    input: 
        # Figures to summarize case counts and SV recall
        str(rules.plot_case_counts.output.seq_type), 
        str(rules.plot_case_counts.output.recall),   
        str(rules.plot_case_counts.output.recall_seq_type),
        # Circos plots
        rules.circos.output,
        # Expression plots
        rules.plot_expression.output, 
        # MiXCR and flow
        rules.process_mixcr.output, 
        rules.mixcr_vs_flow.output, 
        rules.plot_mixcr.output, 
        # aSHM heatmaps and frequency plots
        rules.ashm_heatmaps_IGH.output, 
        rules.ashm_heatmaps_all.output, 
        # Additional SV analyses 
        rules.cryptic_bcl2.output,
        rules.myc_svs.output
        
        