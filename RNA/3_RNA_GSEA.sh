#########################################
# GSEA analysis
#########################################

# Generate a .rnk file of genes (mapped to human symbols) and log2FC values
awk -F "\t" '{if(NR>1 && $3!="NA"){print $9 "\t" $3}}' DESeq2_RNAseq_CooksT_FiltT.Tumor_high_vs_Tumor_low.mart_annot.txt > DESeq2_RNAseq_CooksT_FiltT.Tumor_high_vs_Tumor_low.mart_annot.rnk

# Run GSEA
i=DESeq2_RNAseq_CooksT_FiltT.Tumor_high_vs_Tumor_low.mart_annot.rnk
java -cp /home/himrichova/software/gsea/gsea-3.0.jar -Xmx200G xtools.gsea.GseaPreranked -rpt_label gsea_${i%.rnk} -gmx exhaustion_related_genesets.gmt -rnk ${i} -set_min 20 -set_max 2500 -plot_top_x 40 -nperm 1000000 

