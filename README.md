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

## Citing this data
All data and experiments are described in our paper: 
Hilton, L.K., Collinge, B., Ben-Neriah, S., Alduaij, W., Shaalan, H., Weng, A.P., Cruz, M., Slack, G.W., Farinha, P., Miyata-Takata, T., Boyle, M., Meissner, B., Cook, J.R., Ondrejka, S.L., Ott, G., Rosenwald, A., Campo, E., Amador, C., Greiner, T.C., Raess, P.W., Song, J.Y., Inghirami, G., Jaffe, E.S., Weisenburger, D.D., Chan, W.C., Beiske, K., Fu, K., Delabie, J., Pittaluga, S., Iqbal, J., Wright, G., Sehn, L.H., Savage, K.J., Mungall, A.J., Feldman, A.L., Staudt, L.M., Steidl, C., Rimsza, L.M., Morin, R.D., Scott, D.W., 2024. Motive and opportunity: MYC rearrangements in high-grade B-cell lymphoma with MYC and BCL2 rearrangements (an LLMPP study). Blood 144, 525–540. https://doi.org/10.1182/blood.2024024251
