################################################################################################
# Define differentially accessible peaks in vivo and in vitro
################################################################################################


# Thresholds applied on in vivo ATAC DESeq results: adjusted P-value ≤ 0.05 & |log2FC| ≥ 0.5

inVivo_results=ATAC_inVivo_DESeq2.Tumor_high_vs_Tumor_low.txt

awk -F "\t" '{if(NR>1 && $3>=0.5 && $7 <= 0.05){print $1}}' ${i} > ATAC_inVivo_High.txt
awk -F "\t" '{if(NR>1 && $3<=-0.5 && $7 <= 0.05){print $1}}' ${i} > ATAC_inVivo_Low.txt

# generate files with coordinates of significantly differentially accessible peaks in vivo

awk -F "\t" '{if(NR>1 && $3>=0.5 && $7 <= 0.05){print $0}}' ${i} | cut -f1 | tr ":" "\t" | tr "-" "\t" > ATAC_inVivo_High.bed
awk -F "\t" '{if(NR>1 && $3<=-0.5 && $7 <= 0.05){print $0}}' ${i} | cut -f1 | tr ":" "\t" | tr "-" "\t" > ATAC_inVivo_Low.bed




# Thresholds applied on in vitro ATAC DESeq results: P-value ≤ 0.05 & |log2FC| ≥ 0.5

inVitro_results=ATAC_inVitro_DESeq2.Tumor_high_vs_Tumor_low.txt

awk -F "\t" '{if(NR>1 && $3>=0.5 && $6 <= 0.05){print $1}}' ${i} > ATAC_inVitro_High.txt
awk -F "\t" '{if(NR>1 && $3<=-0.5 && $6 <= 0.05){print $1}}' ${i} > ATAC_inVitro_Low.txt

# generate files with coordinates of significantly differentially accessible peaks in vitro

awk -F "\t" '{if(NR>1 && $3>=0.5 && $6 <= 0.05){print $0}}' ${i} | cut -f1 | tr ":" "\t" | tr "-" "\t" > ATAC_inVitro_High.bed
awk -F "\t" '{if(NR>1 && $3<=-0.5 && $6 <= 0.05){print $0}}' ${i} | cut -f1 | tr ":" "\t" | tr "-" "\t" > ATAC_inVitro_Low.bed



# Convert genome coordinates of differentially accessible regions from mm10 to mm9
# (liftovered coordinates are used as input to icisTarget)

for i in ATAC*.bed
do
echo $i
liftOver -minMatch=0.95 ${i} /home/himrichova/1_Projects/resources/liftover_over_chain/mm10ToMm9.over.chain ${i%.bed}.mm10ToMm9.bed ${i%.bed}.mm10ToMm9.log
done


