src/snvs/CheckMotifMutBias.py \
    -m data/maf/genome_capture.hg38.maf \
    -r /projects/rmorin/reference/lcr-modules-references/genomes/hg38/genome_fasta/genome.fa \
    -s "WRCY" \
    -i 3 \
    -o data/maf/genome_capture.WRCY.hg38.tsv \
    --annotate_maf data/maf/genome_capture.WRCY.hg38.maf \
    -t data/region_data/region_bed_simple.hg38.bed
    
# for maf in data/maf/region_mafs/*_genome_capture.hg38.maf; do 

#     out=${maf/hg38.maf/WRCY.hg38.tsv}
#     echo "$out"
    
#     src/snvs/CheckMotifMutBias.py \
#         -m $maf \
#         -r /projects/rmorin/reference/lcr-modules-references/genomes/hg38/genome_fasta/genome.fa \
#         -s "WRCY" \
#         -i 3 \
#         -o $out \
#         -t data/region_data/region_bed_simple.hg38.bed
            
# done