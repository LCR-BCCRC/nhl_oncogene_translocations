# Oncogene rearrangements in mature B-cell lymphomas

This repository contains all code used in the analysis of *MYC*, *BCL2*, and *BCL6* rearrangements from targeted capture, whole genome, and RNAseq data of Burkitt lymphoma (BL), diffuse large B-cell lymphoma (DLBCL), follicular lymphoma (FL), and high-grade B-cell lymphoma, double hit with *MYC* and *BCL2* or *BCL6* rearrangements (HGBCL-DH-*BCL2*/HGBCL-DH-*BCL6*). All R packages used in the analysis can be installed reproducibly with `renv`, using R version 4.1.\*. Outputs can be regenerated with the included Snakefile:

```         
snakemake -j5 all
```

## Sequencing data analysis

Sequencing files were processed outside of this repository using modules from [LCR-Modules](https://github.com/LCR-BCCRC/lcr-modules) as follows:

| Sequencing Type | Output                             | Module(s)               |
|------------------|--------------------------------|-----------------------|
| genome, capture | Somatic variant calls              | slms_3-1.0, vcf2maf-1.3 |
| genome, capture | Somatic structural variants        | svar_master-1.0         |
| genome, capture | Quality control (coverage metrics) | qc-1.0                  |
| RNAseq          | Immunoglobulin rearrangements      | mixcr-1.2               |
| RNAseq          | Expression                         | salmon-1.1              |