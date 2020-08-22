#########################################################################################################################
# Processing output from the icisTarget TF-binding motif enrichment analysis
########################################################################################################################
module load multigrep/2018

cd icistarget_analysis_results

mkdir predicted_targets


# Path to the motif database annotation file used by iRegulon (Janky et al., 2014) and i-cisTarget (Imrichova et al., 2015). TFs identified as "Unlikely to be sequence specific TF" or "ssDNA/RNA binding" were filetered out according to (Lambert et al., 2018).

motifs_tbl=/home/himrichova/1_Projects/Motifs/mm9/motifs-v8-nr.mgi-m0.001-o0.0.true_TFs.tbl


# Specify the file of RNA DESeq results (ranking of gene names based on the log2FC values)
RNA_DESeq_results=DESeq2_RNAseq.Tumor_high_vs_Tumor_low.mart_annot.geneNames.rnk



# ctx_results_dir: Specify directory names where icistarget results (directories) were downloaded

# - icisTarget results for differentially methylated regions based on the following thresholds: adj.p-value <= 0.05 & |meanquotlog2| >= 0.5 
# - icisTarget results for differentially accessible regions in vivo (ATAC-seq) based on the following thresholds: adj.p-value <= 0.05 & |log2FC| >= 0.5 
# - icisTarget results for differentially accessible regions in vitro (ATAC-seq) based on the following thresholds: p-value <= 0.05 & |log2FC| >= 0.5 

for ctx_results_dir in ctx_ATACvivo_Tumor_Low_log2FC05_adjP05 ctx_ATACvivo_Tumor_High_log2FC05_adjP05 ctx_ATACvitro_Tumor_Low_log2FC05_P05 ctx_ATACvitro_Tumor_High_log2FC05_P05 ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05 ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05 ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05 ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05
do

echo ${ctx_results_dir}

for cluster_id in $(cat ${ctx_results_dir}/clusters.tbl|cut -f1|sort -u)
do

multigrep.sh -f 1 -w -g <(cat ${ctx_results_dir}/clusters.tbl ${motifs_tbl} |multigrep.sh -f 1 -w -p ${cluster_id} |cut -f3|sort -u) |cut -f6 |sort -u > predicted_targets/${ctx_results_dir}_possible_TFs_cluster_${cluster_id}.txt

multigrep.sh -f 3 -w -g <(cat ${ctx_results_dir}/clusters.tbl ${ctx_results_dir}/statistics.tbl|multigrep.sh -f 1 -w -p ${cluster_id} |cut -f3|sort -u) |cut -f3,8 > predicted_targets/${ctx_results_dir}_NES_scores_cluster_${cluster_id}.txt


multigrep.sh -f 3 -w -g <(cat ${ctx_results_dir}/clusters.tbl ${ctx_results_dir}/statistics.tbl |multigrep.sh -f 1 -w -p ${cluster_id} |cut -f3|sort -u) |cut -f9|tr ";" "\n"|sort -u > predicted_targets/${ctx_results_dir}_target_regions_cluster_${cluster_id}.txt


multigrep.sh -f 1 -w -g <(cat predicted_targets/${ctx_results_dir}_possible_TFs_cluster_${cluster_id}.txt) ${RNA_DESeq_results} > predicted_targets/${ctx_results_dir}_any_log2FC_of_possible_TFs_cluster_${cluster_id}.txt


paste <(echo cluster_${cluster_id}) <(cat predicted_targets/${ctx_results_dir}_possible_TFs_cluster_${cluster_id}.txt|wc -l) <(cut -f2 predicted_targets/${ctx_results_dir}_NES_scores_cluster_${cluster_id}.txt |sort -k1 -gr |head -1) <(cut -f2 predicted_targets/${ctx_results_dir}_NES_scores_cluster_${cluster_id}.txt|wc -l) <(cat predicted_targets/${ctx_results_dir}_target_regions_cluster_${cluster_id}.txt |wc -l) <(cat predicted_targets/${ctx_results_dir}_any_log2FC_of_possible_TFs_cluster_${cluster_id}.txt|cut -f2 |sort -k1 -gr |head -1)  <(cat predicted_targets/${ctx_results_dir}_any_log2FC_of_possible_TFs_cluster_${cluster_id}.txt|cut -f2|sort -k1 -gr |tail -1) <(cat predicted_targets/${ctx_results_dir}_any_log2FC_of_possible_TFs_cluster_${cluster_id}.txt|sort -k2 -gr |head -1|cut -f1) <(cat predicted_targets/${ctx_results_dir}_any_log2FC_of_possible_TFs_cluster_${cluster_id}.txt|sort -k2 -gr |tail -1|cut -f1) <(cat predicted_targets/${ctx_results_dir}_any_log2FC_of_possible_TFs_cluster_${cluster_id}.txt |sort -k2 -gr |cut -f1 |tr "\n" ",") <(cat predicted_targets/${ctx_results_dir}_possible_TFs_cluster_${cluster_id}.txt |sort -k1 -V  |tr "\n" ",") >> true_TFs_${ctx_results_dir}.MotifCluster_numberTFs_topNES_numberMotif_numberTGs_highestExpression_lowestExpression_highTF_lowTF_anyExpressedTF_anyTF.txt


done 




# Filter out cluster motifs having no expressed TF assigned:
awk -F "\t" '{if($6!=""){print $0}}' true_TFs_${ctx_results_dir}.MotifCluster_numberTFs_topNES_numberMotif_numberTGs_highestExpression_lowestExpression_highTF_lowTF_anyExpressedTF_anyTF.txt > TFs_true_TFs_${ctx_results_dir}.MotifCluster_numberTFs_topNES_numberMotif_numberTGs_highestExpression_lowestExpression_highTF_lowTF_anyExpressedTF_anyTF.txt


# Generate a list of all possible expressed TFs (a file with 1 column = TF names)
cut -f10 true_TFs_${ctx_results_dir}.MotifCluster_numberTFs_topNES_numberMotif_numberTGs_highestExpression_lowestExpression_highTF_lowTF_anyExpressedTF_anyTF.txt|tr "," "\n"|sort -u > all_possible_TFs_${ctx_results_dir}.txt


# Generate a table without non-expressed TFs:
awk -F "\t" '{if($10!=""){print $0}}' true_TFs_${ctx_results_dir}.MotifCluster_numberTFs_topNES_numberMotif_numberTGs_highestExpression_lowestExpression_highTF_lowTF_anyExpressedTF_anyTF.txt |sort -k5 -gr |awk -F "\t" '{print $1 "\t" $2 "\t" $4 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }' > complete_list_true_TFs_${ctx_results_dir}.txt

done









#####################
# ATAC in vivo
#####################
# Generate a list of predicted TFs specific for MDR/MG_Hi or MDR/MG_Lo ATAC peaks in vivo

grep -w -E -f all_possible_TFs_ctx_ATACvivo_Tumor_High_log2FC05_adjP05.txt -v all_possible_TFs_ctx_ATACvivo_Tumor_Low_log2FC05_adjP05.txt > TFs_ATACvivo_tumor_Low_specific.txt
grep -w -E -f all_possible_TFs_ctx_ATACvivo_Tumor_Low_log2FC05_adjP05.txt -v all_possible_TFs_ctx_ATACvivo_Tumor_High_log2FC05_adjP05.txt > TFs_ATACvivo_tumor_High_specific.txt



# TFs binding peaks accessible in MDR/MG_Hi
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_ATACvivo_Tumor_High_log2FC05_adjP05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_ATACvivo_Tumor_High_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_ATACvivo_Tumor_High_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_ATACvivo_tumor_High_specific.txt |tr "\n" ",")
done |grep -w -E -f TFs_ATACvivo_tumor_High_specific.txt > complete_list_true_TFs_ctx_ATACvivo_Tumor_High_log2FC05_adjP05.specific_TFs.txt
 


# TFs binding peaks accessible in MDR/MG_Low
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_ATACvivo_Tumor_Low_log2FC05_adjP05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_ATACvivo_Tumor_Low_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_ATACvivo_Tumor_Low_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_ATACvivo_tumor_Low_specific.txt |tr "\n" ",")
done |grep -w -E -f TFs_ATACvivo_tumor_Low_specific.txt > complete_list_true_TFs_ctx_ATACvivo_Tumor_Low_log2FC05_adjP05.specific_TFs.txt




#####################
# ATAC in vitro
#####################
# Generate a list of predicted TFs specific for MDR/MG_Hi or MDR/MG_Lo ATAC peaks in vitro

grep -w -E -f all_possible_TFs_ctx_ATACvitro_Tumor_High_log2FC05_adjP05.txt -v all_possible_TFs_ctx_ATACvitro_Tumor_Low_log2FC05_adjP05.txt > TFs_ATACvitro_tumor_Low_specific.txt
grep -w -E -f all_possible_TFs_ctx_ATACvitro_Tumor_Low_log2FC05_adjP05.txt -v all_possible_TFs_ctx_ATACvitro_Tumor_High_log2FC05_adjP05.txt > TFs_ATACvitro_tumor_High_specific.txt



# TFs binding peaks accessible in MDR/MG_Hi
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_ATACvitro_Tumor_High_log2FC05_adjP05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_ATACvitro_Tumor_High_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_ATACvitro_Tumor_High_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_ATACvitro_tumor_High_specific.txt |tr "\n" ",")
done |grep -w -E -f TFs_ATACvitro_tumor_High_specific.txt > complete_list_true_TFs_ctx_ATACvitro_Tumor_High_log2FC05_adjP05.specific_TFs.txt
 


# TFs binding peaks accessible in MDR/MG_Low
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_ATACvitro_Tumor_Low_log2FC05_adjP05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_ATACvitro_Tumor_Low_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_ATACvitro_Tumor_Low_log2FC05_adjP05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_ATACvitro_tumor_Low_specific.txt |tr "\n" ",")
done |grep -w -E -f TFs_ATACvitro_tumor_Low_specific.txt > complete_list_true_TFs_ctx_ATACvitro_Tumor_Low_log2FC05_adjP05.specific_TFs.txt






#####################
# WGBS CpGs
#####################
# Generate a list of predicted TFs specific for hypo-methylated CpGs in MDR/MG_Hi or MDR/MG_Lo samples

grep -w -E -f all_possible_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt -v all_possible_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt > TFs_tumor_Low_specific.CPG.txt

grep -w -E -f all_possible_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt -v all_possible_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt > TFs_tumor_High_specific.CPG.txt



# WGBS CPGs hypomethyl in MDR/MG_Hi
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_tumor_High_specific.CPG.txt |tr "\n" ",")
done |grep -w -E -f TFs_tumor_High_specific.CPG.txt > complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.specific_TFs.txt



# WGBS CPGs hypomethyl in MDR/MG_Lo
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_tumor_Low_specific.CPG.txt |tr "\n" ",")
done |grep -w -E -f TFs_tumor_Low_specific.CPG.txt > complete_list_true_TFs_ctx_WGBS_CPG_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.specific_TFs.txt



#####################
# WGBS Promoters
#####################
# Generate a list of predicted TFs specific for hypo-methylated promoters in MDR/MG_Hi or MDR/MG_Lo samples

grep -w -E -f all_possible_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt -v all_possible_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt > TFs_tumor_Low_specific.PROMOTERS.txt

grep -w -E -f all_possible_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt -v all_possible_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt > TFs_tumor_High_specific.PROMOTERS.txt



# WGBS Promoters hypomethyl in MDR/MG_Hi
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_tumor_High_specific.PROMOTERS.txt |tr "\n" ",")
done |grep -w -E -f TFs_tumor_High_specific.PROMOTERS.txt > complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_High_meanquotlog2_05_FDR_05.specific_TFs.txt



# WGBS Promoters hypomethyl in MDR/MG_Low
for R in $(seq 1 $(cat complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt|wc -l))
do
paste <(cat complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $0}}' ) <(cat complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.txt| awk -F "\t" -v row=${R} '{if(NR==row){print $10}}' |tr "," "\n"|grep -w -E -f  TFs_tumor_Low_specific.PROMOTERS.txt |tr "\n" ",")
done |grep -w -E -f TFs_tumor_Low_specific.PROMOTERS.txt > complete_list_true_TFs_ctx_WGBS_PROMOTERS_hypomethyl_in_Tumor_Low_meanquotlog2_05_FDR_05.specific_TFs.txt




######################################################################
######################################################################
######################################################################