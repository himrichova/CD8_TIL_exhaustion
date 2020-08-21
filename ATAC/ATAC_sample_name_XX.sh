#!/bin/bash
#SBATCH --job-name='ATAC_sample_name_XX'
#SBATCH --output='ATAC_sample_name_XX.log'
#SBATCH --time='12:00:00'
#SBATCH --partition='shortq'
#SBATCH -m block
/bin/hostname


module load bowtie/2.2.4
module load sambamba/0.5.5
module load htslib/1.3.1
module load samtools/1.3.1
module load bedtools/2.26.0
module load multigrep/2018
module load MACS/2.1.0
#macs2 --version
#macs2 2.1.1.20160309




genome=mm10
local_raw_data_location=SPECIFY_local_raw_data_location
sample_names=SPECIFY_sample_names
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
`/cm/shared/apps/java/jdk/1.7.0_80/bin/java -Xmx60000m -jar /cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar SE -phred33 -threads 2 ${mapping_results_dir}/fastq/${sample_name}_R1.fastq ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq HEADCROP:13 ILLUMINACLIP:/data/groups/lab_bock/shared/resources/adapters/epignome_adapters_2_add.fa:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:18`

# Target to produce fastqc report
`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq`


# Generate aligned filtered bam file and index
input_fastq=${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq

bowtie2 --very-sensitive --no-discordant -p 8 -x /data/groups/lab_bock/shared/resources/genomes/${genome}/indexed_bowtie2/${genome} --met-file ${mapping_results_dir}/bam/${sample_name}.aln_metrics.${genome}.txt ${input_fastq}  2> ${mapping_results_dir}/bam/${sample_name}.aln_rates.${genome}.txt | samtools view -S -b - | samtools sort -o ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.bam -


sambamba markdup -t 8 -r --compression-level=0 ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.bam ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.nodups.nofilter.bam 2> ${mapping_results_dir}/bam/${sample_name}.dups_metrics.${genome}.txt

sambamba view -t 8 -f bam --valid -F "not unmapped and not (secondary_alignment or supplementary) and mapping_quality >= 30" ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.nodups.nofilter.bam |sambamba sort -t 8 /dev/stdin -o ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.bam

samtools index ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.bam
samtools index ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.bam

sambamba flagstat ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.bam > ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.flagstat
sambamba flagstat ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.bam > ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.flagstat




# Generate bigWig
bedtools bamtobed -i ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.bam | bedtools slop -i stdin -g /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes -s -l 0 -r 130 | /home/himrichova/1_Projects/0_scripts/fix_bedfile_genome_boundaries.py ${genome} | genomeCoverageBed -bg -g /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes -i stdin > ${mapping_results_dir}/bam/${sample_name}_${genome}.cov

awk 'NR==FNR{sum+= $4; next}{ $4 = ($4 / sum) * 1000000; print}' ${mapping_results_dir}/bam/${sample_name}_${genome}.cov ${mapping_results_dir}/bam/${sample_name}_${genome}.cov| sort -k1,1 -k2,2n > ${mapping_results_dir}/bam/${sample_name}_${genome}.normalized.cov

bedGraphToBigWig ${mapping_results_dir}/bam/${sample_name}_${genome}.normalized.cov /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes ${mapping_results_dir}/bam/${sample_name}_${genome}.bigWig

chmod 755 ${mapping_results_dir}/bam/${sample_name}_${genome}.bigWig




# Peak calling
macs2 callpeak -t ${mapping_results_dir}/bam/${sample_name}.trimmed.bowtie2.${genome}.filtered.bam --nomodel --p 1e-09 --extsize 147 -g mm -n ${sample_name}_ATAC_${genome} --outdir ${mapping_results_dir}/peaks/log10p_9_peak_files

intersectBed -a ${mapping_results_dir}/peaks/log10p_9_peak_files/${sample_name}_ATAC_${genome}_peaks.narrowPeak -b /data/groups/lab_bock/shared/resources/regions/blacklists/blacklists_v2/mm10-blacklist.v2.bed -v | multigrep.sh -f 1 -w -g ~/1_Projects/annotation/mouse/mm10/chromosomes_mm10.txt > ${mapping_results_dir}/peaks/log10p_9_peak_files/${sample_name}_ATAC_${genome}.narrowPeak.blacklists_filtered.bed


