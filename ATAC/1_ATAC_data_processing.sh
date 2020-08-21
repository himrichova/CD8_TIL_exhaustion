##################################################################
# ATAC-seq data processing 
##################################################################

module load bedtools/2.26.0



# Specify the location of the raw data files (raw bam files that are available in the GEO)
local_raw_data_location=

# Specify the directory where the output data will be generated
mapping_results_dir=

# Specify path to a file with sample names that are supposed to be processed 
sample_names=

# Example of the samples_names.txt
#Sample_Name     Condition
#Tumor_high_1    Tumor_high
#Tumor_high_2    Tumor_high
#Tumor_low_1     Tumor_low



# Create subdirectories for the results
mkdir -p ${mapping_results_dir}
mkdir -p ${mapping_results_dir}/fastq
mkdir -p ${mapping_results_dir}/fastqc
mkdir -p ${mapping_results_dir}/bam
mkdir -p ${mapping_results_dir}/peaks
mkdir -p ${mapping_results_dir}/peaks/log10p_9_peak_files
mkdir -p ${mapping_results_dir}/../data_analysis/
mkdir -p ${mapping_results_dir}/../data_analysis/coverage_merged_macs2_peaks



cd ${mapping_results_dir}

# Adapt "ATAC_sample_name_XX.sh" script: add sample name and paths to the local raw data and mapping results directories:

for sample_name in $( awk -F "\t" '{if(NR>1){print $1}}' $sample_names)
do
echo $sample_name
sed 's/sample_name_XX/'${sample_name}'/g' ATAC_sample_name_XX.sh | sed 's+SPECIFY_local_raw_data_location+'${local_raw_data_location}'+g' | sed 's+SPECIFY_sample_names+'${sample_names}'+g' | sed 's+SPECIFY_mapping_results_dir+'${mapping_results_dir}'+g' > ATAC_${sample_name}.sh
done




# Submit jobs:

for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
echo $sample_name
sbatch ATAC_${sample_name}.sh
done



############################################################################################################
# Generate a consensus peak list from either all in vivo or in vitro samples (MDR/MG high and MDR/MG low)

# move to the directory with log10p_9_peak_files (either in vivo or in vitro)
cd ${mapping_results_dir}/peaks/log10p_9_peak_files/

# merge ATAC peaks from all in vivo samples or all in vitro samples
bedtools merge -i <(cat *.blacklists_filtered.bed|sort -k1,1V -k2,2n -k3,3n) | sort -k1,1V -k2,2n -k3,3n |awk -F "\t" '{print $0 "\t" $1":"$2"-"$3}' > ATAC_hi_lo_merged_log10p_9.mm10.bed




# Quantify the chromatin accessibility of each consensus peak in each sample -> generate a count files, counting the number of reads from the filtered ATAC BAM file that overlap each region from the consensus peak list 

OUTPUT_DIR=${mapping_results_dir}/../data_analysis/coverage_merged_macs2_peaks
REGIONS=${mapping_results_dir}/peaks/log10p_9_peak_files/ATAC_hi_lo_merged_log10p_9.mm10.bed
MAPPING_RESULTS_BAM=${mapping_results_dir}/bam

mkdir ${OUTPUT_DIR}

################################################################
# Coverage from Bed 
# (bamToBed -> extended 200bp -> coverage)
##############################################################

for B in ${MAPPING_RESULTS_BAM}/*mm10.filtered.bam
do
echo "${B}"

bedtools bamtobed -i ${B} | bedtools slop -g /data/groups/lab_bock/shared/resources/genomes/mm10/mm10.chromSizes -b 200 > ${OUTPUT_DIR}/$(basename ${B%.bam}).extended200bp.bed

if [ ! -f ${OUTPUT_DIR}/coverage-$(basename ${B%.bam}).extended200bp.bed ]
then
echo "coverage-$(basename ${B%.bam}).extended200bp.bed does NOT exist yet."
coverageBed -a ${OUTPUT_DIR}/$(basename ${B%.bam}).extended200bp.bed -b ${REGIONS} > ${OUTPUT_DIR}/coverage-$(basename ${B%.bam}).extended200bp.bed
else
echo "$(basename ${B%.bam}).extended200bp.bed exists."
fi
done



################################
# Generate the read count matrix


cd ${mapping_results_dir}/../data_analysis/coverage_merged_macs2_peaks/DESeq2_High_Low_samples

rm read_counts_matrix.tmp

first_sample_name=$(awk -F "\t" '{if(NR==2){print $1}}' ${sample_names})
cat ../coverage-${first_sample_name}.trimmed.bowtie2.mm10.filtered.extended200bp.bed |cut -f1-4 > read_counts_matrix.tmp



for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
echo $sample_name
paste <(cat read_counts_matrix.tmp) <(cat ../coverage-${sample_name}.trimmed.bowtie2.mm10.filtered.extended200bp.bed  |cut -f5 ) > read_counts_matrix.2.tmp
mv read_counts_matrix.2.tmp read_counts_matrix.tmp
done


paste <(echo chr) <(echo start) <(echo stop) <(echo region_ID) <(awk -F "\t" '{if(NR>1){print $1}}' sample_names_final.txt |tr "\n" "\t") > colnames.txt


cat colnames.txt read_counts_matrix.tmp > read_counts_matrix.txt

cat read_counts_matrix.txt |cut -f4- > read_counts_matrix.reg_IDs.txt
