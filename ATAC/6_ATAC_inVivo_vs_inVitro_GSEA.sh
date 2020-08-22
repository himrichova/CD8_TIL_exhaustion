#########################################################################################
# Compare in vivo and in vitro ATAC peaks (using GSEA)
# 
# Ranking of all in vivo consensus peaks compared to diff.accessible peaks in vitro
#########################################################################################

# 1) Generate a bed file and rnk file for all in VIVO peaks


awk -F "\t" '{if(NR>1){print $1 "\t" $3}}' ATAC_inVivo_DESeq2.Tumor_high_vs_Tumor_low.txt|sort -k2 -gr > inVivo_DESeq2_ATACseq_PEAKS.log2FC.rnk


paste <(cat inVivo_DESeq2_ATACseq_PEAKS.log2FC.rnk|cut -f1 |tr ":" "\t" |tr "-" "\t") <(cat  inVivo_DESeq2_ATACseq_PEAKS.log2FC.rnk) > inVivo_DESeq2_ATACseq_PEAKS.bed




##############################################################################
# 2) Assign signif. diff. accessible peaks in VITRO to ATAC peaks in VIVO


# assign diff.accessible peaks in vitro peaks to in vivo peaks 
# (per each diff.access. in vitro peak, get a region ID "chr:start-stop" of the closest in vivo peak)

for signif_bed_inVitro in ATAC_inVitro_High.bed ATAC_inVitro_Low.bed
do
closestBed -a ${signif_bed_inVitro} -b inVivo_DESeq2_ATACseq_PEAKS.bed |awk -F "\t" '{print $8 }' |sort -u  > ${signif_bed_inVitro%.bed}.assigned_to_inVivo_PEAKS.txt
done


paste ATAC_inVitro_High.assigned_to_inVivo_PEAKS.txt ATAC_inVitro_Low.assigned_to_inVivo_PEAKS.txt > signif_lfc05_p05_ATAC_inVITRO_assigned_to_inVivo_PEAKS.gmx


##############################################################################
# 3) Run GSEA

i=inVivo_DESeq2_ATACseq_PEAKS.log2FC.rnk

java -cp /home/himrichova/software/gsea/gsea-3.0.jar -Xmx3200m -Xms1600m xtools.gsea.GseaPreranked -rpt_label ${i%.rnk}_ATAC_in_VITRO_assigned_PEAKS -gmx signif_lfc05_p05_ATAC_inVITRO_assigned_to_inVivo_PEAKS.gmx -rnk ${i} -set_min 10 -set_max 5000 -plot_top_x 40 -nperm 10000


