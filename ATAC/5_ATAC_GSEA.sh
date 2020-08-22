######################################################
# Assign peaks to the closest genes and GSEA
######################################################

# Specify ATAC DESEq output file

DESeq_results=ATAC_inVivo_DESeq2.Tumor_high_vs_Tumor_low.txt    # in vivo peaks
#DESeq_results=ATAC_inVitro_DESeq2.Tumor_high_vs_Tumor_low.txt  # in vitro peaks


# Create a bed file including the DESeq results:

paste <(cat ${DESeq_results}|cut -f1 |tr ":" "\t" |tr "-" "\t") <(cat  ${DESeq_results}) |awk -F "\t" '{if(NR>1){print $1 "\t" $2 "\t" $3 "\t" $6 }}' > ${DESeq_results%.txt}.log2FC.bed


# Assign peaks to genes (mm10) and rank genes according to log2FC values

closestBed -a ${DESeq_results%.txt}.log2FC.bed -b /home/himrichova/1_Projects/resources/genomes/mm10/mm10_genes.bed |awk -F "\t" '{print $8 "\t" $4}' |sort -u |sort -k2 -gr > genes_${DESeq_results%.txt}.log2FC.rnk




DESeq2_rnk=genes_${DESeq_results%.txt}.log2FC.rnk

# Remove NA values and selet max/min score per mm10 gene
cat ${DESeq2_rnk}|grep -v -w NA| sort -k1 -V | groupBy -g 1 -c 2,2 -o max,min > ${DESeq2_rnk%.txt}.max_min.txt



# Map mouse genes to human symbols and add columns with absolute values of max/min log2FC values
R --vanilla --args ${DESeq2_rnk%.txt}.max_min.txt mouse_to_human Mouse_gene_name < /data/groups/lab_winter/himrichova/resources/mapping_genes_between_human_and_mouse/map_genes_between_mouse_and_human.R 



rnk_human_genes=${DESeq2_rnk%.txt}.max_min.txt.ToHuman.txt

# Per each peak, select the log2FC value with the highest log2FC from the ATAC DESeq analysis:

awk -F "\t" '{if(NR>1){print $0}}' ${rnk_human_genes} |awk -F "\t" '{if($6>$7){print $1 "\t" $2 "\t" $6 "\t" $0};if($6<=$7){print $1 "\t" $3 "\t" $7 "\t" $0}}' |sort -k8 -V |groupBy -g 8 -c 2,2 -o max,min > ${rnk_human_genes%.ToHuman.txt}.ToHuman.Max_Min.txt

R --vanilla --args ${rnk_human_genes%.ToHuman.txt}.ToHuman.Max_Min.txt < /data/groups/lab_winter/himrichova/resources/mapping_genes_between_human_and_mouse/abs_value.R


awk -F "\t" '{if(NR>1){print $0}}' ${rnk_human_genes%.ToHuman.txt}.ToHuman.Max_Min.txt.abs.txt|awk -F "\t" '{if($4>$5){print $1 "\t" $2};if($4<=$5){print $1 "\t" $3}}'  |sort -u|sort -k2 -gr > ${rnk_human_genes%.ToHuman.txt}.absMax_absMin.ToHuman.absMax_abxMin.human.rnk


# Run GSEA

i=${rnk_human_genes%.ToHuman.txt}.absMax_absMin.ToHuman.absMax_abxMin.human.rnk
java -cp /home/himrichova/software/gsea/gsea-3.0.jar -Xmx200G xtools.gsea.GseaPreranked -rpt_label ${i%.rnk}_gsea -gmx exhaustion_related_genesets.gmt -rnk ${i} -set_min 20 -set_max 2500 -plot_top_x 40 -nperm 1000000 



