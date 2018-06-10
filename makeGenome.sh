#!/bin/bash

#----QSUB INFO----#
##choose queue
####PBS -q
##list - node are nodes: ppn are cpus per node: walltime=walltime
#PBS -l nodes=1:ppn=8,mem=100gb,walltime=600:00:00
##job name
#PBS -N MakeGenome

#----LOAD MODULES----#
module load fastqc
module load trimmomatic/0.36
module load samtools/1.3.1
module load python/2.7.12
module load bwa/0.7.15
module load java/1.8.0_91 
module load picard/2.4.1
module load gatk/3.6
module load bedtools/2 
module load vcftools/v0.1.14-14
module load snpeff/4.3p
module load bcftools
module load perl/5.26.0
module load xz/5.2.2
module load python/2.7.12
module load java/1.8.0_91
module load bcftools
module load bedtools/2
module load vcftools/v0.1.14-14
module load htslib

#----ENVIRONMENT----#
ref=""	#name of reference genome without file extension
wd=""	#path to directory
cd ${wd}

species="" #specific identifier of data without extensions

#----------------------------#
#----BEGIN MAKING GENOMES----#
#----------------------------#

#!!PREP REFERENCE - Only do once!!
###Indexes reference fasta for use by bwa, samtools, picard
###Remove lines if fasta is already indexed by these packages
bwa index ${ref}.fa
samtools faidx ${ref}.fa
java -Xms2g -Xmx14g -jar /tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=${ref}.fa O=${ref}.dict

#MAP - bwa mem
###Store raw reads in directory named 1_raw-reads
mkdir 2_mapping
bwa mem -M -t 8 -R "@RG\tID:${species}\tLB:lib\tPL:ILLUMINA\tSM:${species}" -v 2 ${ref}.fa ${wd}/1_raw-reads/${species}_1.fastq ${wd}/1_raw-reads/${species}_2.fastq > ${wd}/2_mapping/${species}.sam

#SORT & generate mapping STATS
mkdir Stats
samtools view -Sbu ${wd}/2_mapping/${species}.sam | samtools sort -@ 8 -o ${wd}/2_mapping/${species}.bam
samtools index ${wd}/2_mapping/${species}.bam
touch ${wd}/Stats/${species}.stats.txt
touch ${wd}/Stats/${species}.stats.txt
samtools flagstat ${wd}/2_mapping/${species}.bam > ${wd}/Stats/${species}.stats.txt
samtools depth ${wd}/2_mapping/${species}.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${wd}/Stats/${species}.depth.txt

###REALIGN INDELS
mkdir 4_indels
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${wd}/${ref}.fa \
	-I ${wd}/2_mapping/${species}.bam \
	-o ${wd}/4_indels/${species}.intervals
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ${wd}/${ref}.fa \
	-I ${wd}/2_mapping/${species}.bam \
	-targetIntervals ${wd}/4_indels/${species}.intervals \
	-o ${wd}/4_indels/${species}.indel.bam
samtools index ${wd}/4_indels/${species}.indel.bam

#GATK: CALL VARIANTS
mkdir 5_variants
###retrieve raw vcf
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R ${ref}.fa \
	-I ${wd}/4_indels/${species}.indel.bam \
	-o ${wd}/5_variants/${species}.GATK.rawvar.vcf

###split into raw files for indels and snps
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${ref}.fa \
	-V ${wd}/5_variants/${species}.GATK.rawvar.vcf \
	-selectType SNP \
	-o ${wd}/5_variants/${species}.GATK.rawsnp.vcf
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R ${ref}.fa \
        -V ${wd}/5_variants/${species}.GATK.rawvar.vcf \
        -selectType INDEL \
        -o ${wd}/5_variants/${species}.GATK.rawindel.vcf
	
###filter based on GATK recommended stats
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R ${ref}.fa \
	-V ${wd}/5_variants/${species}.GATK.rawsnp.vcf \
	--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
	--filterName "basic_snp_filter" \
	-o ${wd}/5_variants/${species}.GATK.filtersnp.vcf
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
       -T VariantFiltration \
        -R ${ref}.fa \
        -V ${wd}/5_variants/${species}.GATK.rawindel.vcf \
	--filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' \
        --filterName "basic_indel_filter" \
        -o ${wd}/5_variants/${species}.GATK.filterindel.vcf

###CREATE NEW FASTA
snp="5_variants/${species}.GATK.filtersnp.vcf"        #vcf of SNPs
indel="5_variants/${species}.GATK.filterindel.vcf"    #vcf of indels
mkdir 6_new-fasta

###zip and index file of SNP variants
bgzip ${wd}/${snp}
bcftools index ${wd}/${snp}.gz

###insert SNPs into original reference
bcftools consensus -f ${wd}/${ref}.fa ${wd}/${snp}.gz -o ${wd}/6_new-fasta/${species}.snp.fa

###mask indels/gaps onto new reference
bedtools maskfasta -fi ${wd}/6_new-fasta/${species}.snp.fa -bed ${wd}/${indel} -fo ${wd}/6_new-fasta/${species}.fa
