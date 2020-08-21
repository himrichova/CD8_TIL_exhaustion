#!/bin/bash
#SBATCH --job-name='RNA_sample_name_XX'
#SBATCH --output='RNA_sample_name_XX.log'
#SBATCH --time='12:00:00'
#SBATCH --mem-per-cpu=50G
#SBATCH --partition='shortq'
#SBATCH -m block
/bin/hostname

module load star/2.5.2b
module load FastQC
module load htslib/1.3.1
module load samtools/1.3.1

genome=mm10
sample_names=SPECIFY_sample_names
local_raw_data_location=SPECIFY_local_raw_data_location
mapping_results_dir=SPECIFY_mapping_results_dir

sample_name=sample_name_XX


echo $sample_name


# Generate fastq file
# Target to produce: ${mapping_results_dir}/fastq/${sample_name}_R1.fastq
fastq_out=$(echo ${mapping_results_dir}/fastq/${sample_name}_R1.fastq)

`samtools view ${local_raw_data_location}/${sample_name}.bam | awk -v fastq_out=${fastq_out} '{ print "@"$1"\n"$10"\n+\n"$11 > fastq_out; }'`

# Target to produce fastqc report
`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_R1.fastq`

# Trimming
# Target to produce: `{out_dir}/fastq/${sample_name}_R1_trimmed.fastq`

`/cm/shared/apps/java/jdk/1.7.0_80/bin/java -Xmx60000m -jar /cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar SE -phred33 -threads 2 ${mapping_results_dir}/fastq/${sample_name}_R1.fastq ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq HEADCROP:14 ILLUMINACLIP:/data/groups/lab_bock/shared/resources/adapters/epignome_adapters_2_add.fa:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:36`


# Target to produce fastqc report
`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq`




# STAR mapping
# Target to produce: `${mapping_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam`

`STAR --runThreadN 4 --genomeDir /data/groups/lab_bock/shared/resources/genomes/mm10/indexed_STAR/ --readFilesIn ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq  --outFileNamePrefix ${mapping_results_dir}/star_mapping/${sample_name}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --outFilterMismatchNoverReadLmax 0.04 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate`

samtools index ${mapping_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam

# Target to produce bigwig file: ${mapping_results_dir}/hisat2/${sample_name}.bw
`bamCoverage -b ${mapping_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam -o ${mapping_results_dir}/hisat2/${sample_name}.bw --normalizeUsing RPKM`

# generate a table with read counts per transcript
mkdir ${mapping_results_dir}/star_mapping/${sample_name}
/cm/shared/apps/R/3.5.1/bin/Rscript /home/himrichova/scripts/epigen_pipelines/pipelines_trim7_bowtieLocal/pipelines/tools/rna_Overlaps.R ${mapping_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.ensembl.gtf ${mapping_results_dir}/star_mapping/${sample_name} single

awk -F "\t" '{if(NR>1) {print $1 "\t" $2}}' ${mapping_results_dir}/star_mapping/${sample_name}/Gene_counts.tsv > ${mapping_results_dir}/star_mapping/${sample_name}/${sample_name}_read_counts_exons.txt

####### end of the job
