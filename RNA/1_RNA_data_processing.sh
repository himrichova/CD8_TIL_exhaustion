##################################################################
# RNA data processing
##################################################################

# Specify the location of a file with sample names that are supposed to be processed 
sample_names=

# Specify the location of the raw data files (raw bam files that are available in the GEO)
local_raw_data_location=

# Specify the directory where the output data will be generated
mapping_results_dir=




# Create subdirectories for the results
mkdir -p ${mapping_results_dir}
mkdir -p ${mapping_results_dir}/fastq
mkdir -p ${mapping_results_dir}/fastqc
mkdir -p ${mapping_results_dir}/star_mapping
mkdir -p ${mapping_results_dir}/../data_analysis/
mkdir -p ${mapping_results_dir}/../data_analysis/star_gene_counts


cd ${mapping_results_dir}

# Adapt "RNA_sample_name_XX.sh" script: add sample name and paths to the local raw data and mapping results directories:

for sample_name in $( awk -F "\t" '{if(NR>1){print $1}}' $sample_names)
do
echo $sample_name
sed 's/sample_name_XX/'${sample_name}'/g' RNA_sample_name_XX.sh | sed 's+SPECIFY_local_raw_data_location+'${local_raw_data_location}'+g' | sed 's+SPECIFY_sample_names+'${sample_names}'+g' | sed 's+SPECIFY_mapping_results_dir+'${mapping_results_dir}'+g' > RNA_${sample_name}.sh
done



# Submit jobs:

for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
echo $sample_name
sbatch RNA_${sample_name}.sub
done





##################################
# Generate a read count matrix
##################################

cd ${mapping_results_dir}/../data_analysis/star_gene_counts

for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
echo $sample_name
awk -F "\t" '{if(NR>1) {print $1 "\t" $2}}' ${mapping_results_dir}/star_mapping/${sample_name}/Gene_counts.tsv > ${sample_name}_read_counts_exons.txt
done





first_sample_name=$(awk -F "\t" '{if(NR==2){print $1}}' ${sample_names})

cat <(echo gene_name) <(awk -F "\t" '{if(NR>0){print $1}}' ${first_sample_name}_read_counts_exons.txt) > count_matrix.txt


for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
paste <(cat count_matrix.txt) <(cat <(echo ${sample_name}) <(awk -F "\t" '{if(NR>0){print $2}}' ${sample_name}_read_counts_exons.txt )) >> count_matrix.tmp
mv count_matrix.tmp count_matrix.txt
done




