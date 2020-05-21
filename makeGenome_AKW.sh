#!/bin/bash

#----QSUB INFO----#
##choose queue
####PBS -q
##list - node are nodes: ppn are cpus per node: walltime=walltime
#PBS -l nodes=1:ppn=10,mem=100gb,walltime=80:00:00:00
##job name
#PBS -N MakeGenome
##number of array elements
#PBS -t 34
##run on schwartz nodes to prevent preemptions
#PBS -W x=FLAGS:ADVRES:tss0019_lab

#----LOAD MODULES----#
module load fastqc/11.5
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
module load bcftools/1.3.1
module load perl/5.26.0
module load xz/5.2.2
module load python/2.7.12
module load htslib

#----ENVIRONMENT----#
wd="/scratch/adc0032/Sceloporus"	#path to directory
datdir="${wd}/SRRs" 	#path to SRRs
jobfile="${wd}/Scelops" # list of SRR#s to set up jobs for
ref="${wd}/Sund1"	#name of reference genome without file extension
meta="${wd}/Scelop_meta.csv" #matrix to associate names with SRR# and other metadata


##Setting up array; start 34 jobs, assigning the task ID to the amount of jobs requested in the qsub -t; job number, pull SRR with that same line number in job file for analysis
itemsToProcess=34
taskID=$PBS_ARRAYID
startLineNumber=$(($taskID * $itemsToProcess))
endLineNumber=$(( $startLineNumber + $itemsToProcess ))
startLineNumber=$(( $startLineNumber + 1))

for line in $taskID
do

cd ${wd}
sm=$( head -n $line $jobfile | tail -n 1) # pull SRR from file associated with PBS job ID
species=`grep "${sm}" ${meta} |awk -F, '{print $1}'` #specific identifier of data without extensions; getting this information from the metadata file
SRR=`grep "${sm}" ${meta} |awk -F, '{print $2}'` #assigning SRR # to a variable Using this second source to verify my list in my log file (code removed for space)

#setting up directory, organized by species, with subfolders for holding rawread data
mkdir ${species}
fdir="${wd}/${species}" # assigning species directory to variable for easy reference
mkdir ${fdir}/1_raw-reads

#put raw data in species director subfolder
cp ${datdir}/${SRR}_1.fastq.gz ${fdir}/1_raw-reads/${species}_1.fastq.gz &&
cp ${datdir}/${SRR}_2.fastq.gz ${fdir}/1_raw-reads/${species}_2.fastq.gz

#run a quality check on the raw reads using fastqc
zcat ${fdir}/1_raw-reads/${species}_1.fastq.gz|fastqc  --outdir=${fdir}/1_raw-reads stdin
 
#----------------------------#
#----BEGIN MAKING GENOMES----#
#----------------------------#

#!!PREP REFERENCE - Only do once!!
###Indexes reference fasta for use by bwa, samtools, picard
###Remove lines if fasta is already indexed by these packages
#bwa index ${ref}
#samtools faidx ${ref}
#java -Xms2g -Xmx14g -jar /tools/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=${ref} O=${ref}.dict

#MAP - bwa mem
###Store raw reads in directory named 1_raw-reads
mkdir ${fdir}/2_mapping
bwa mem -M -t 8 -R "@RG\tID:${species}\tLB:lib\tPL:ILLUMINA\tSM:${species}" -v 2 ${ref}.fasta ${fdir}/1_raw-reads/${species}_1.fastq.gz ${fdir}/1_raw-reads/${species}_2.fastq.gz > ${fdir}/2_mapping/${species}.sam

#SORT & generate mapping STATS
mkdir ${fdir}/Stats
samtools view -Sbu ${fdir}/2_mapping/${species}.sam | samtools sort -@ 8 -o ${fdir}/2_mapping/${species}.bam
samtools index ${fdir}/2_mapping/${species}.bam
touch ${fdir}/Stats/${species}.stats.txt
touch ${fdir}/Stats/${species}.depth.txt
samtools flagstat ${fdir}/2_mapping/${species}.bam > ${fdir}/Stats/${species}.stats.txt
samtools depth ${fdir}/2_mapping/${species}.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${fdir}/Stats/${species}.depth.txt

###REALIGN INDELS
mkdir ${fdir}/4_indels
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${ref}.fasta \
	-I ${fdir}/2_mapping/${species}.bam \
	-o ${fdir}/4_indels/${species}.intervals
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ${ref}.fasta \
	-I ${fdir}/2_mapping/${species}.bam \
	-targetIntervals ${fdir}/4_indels/${species}.intervals \
	-o ${fdir}/4_indels/${species}.indel.bam
samtools index ${fdir}/4_indels/${species}.indel.bam

#GATK: CALL VARIANTS
mkdir ${fdir}/5_variants
###retrieve raw vcf
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R ${ref}.fasta \
	-I ${fdir}/4_indels/${species}.indel.bam \
	-o ${fdir}/5_variants/${species}.GATK.rawvar.vcf

###split into raw files for indels and snps
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${ref}.fasta \
	-V ${fdir}/5_variants/${species}.GATK.rawvar.vcf \
	-selectType SNP \
	-o ${fdir}/5_variants/${species}.GATK.rawsnp.vcf
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R ${ref}.fasta \
        -V ${fdir}/5_variants/${species}.GATK.rawvar.vcf \
        -selectType INDEL \
        -o ${fdir}/5_variants/${species}.GATK.rawindel.vcf

###create bed file with zero coverage coordinates; used to mask no-coverage sites in final species fasta
bedtools genomecov -ibam ${fdir}/4_indels/${species}.indel.bam -bga |awk '$4 == 0' > ${fdir}/5_variants/${species}.zero.bed

###filter based on GATK recommended stats
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R ${ref}.fasta \
	-V ${fdir}/5_variants/${species}.GATK.rawsnp.vcf \
	--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
	--filterName "basic_snp_filter" \
	-o ${fdir}/5_variants/${species}.GATK.filtersnp.vcf
java -Xms2g -Xmx14g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
       -T VariantFiltration \
        -R ${ref}.fasta \
        -V ${fdir}/5_variants/${species}.GATK.rawindel.vcf \
	--filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' \
        --filterName "basic_indel_filter" \
        -o ${fdir}/5_variants/${species}.GATK.filterindel.vcf

###CREATE NEW FASTA
snp="5_variants/${species}.GATK.filtersnp.vcf"         #vcf of SNPs
indel="5_variants/${species}.GATK.filterindel.vcf"     #vcf of indels
zero="/5_variants/${species}.zero.bed"            #bed of zero-coverage sites

mkdir ${fdir}/6_new-fasta

###zip and index file of SNP variants
bgzip ${fdir}/${snp}
bcftools index ${fdir}/${snp}.gz

###insert SNPs into original reference
bcftools consensus -f ${ref}.fasta ${fdir}/${snp}.gz -o ${fdir}/6_new-fasta/${species}.snp.fa

###mask indels & no-coverage regions for new species genome
bedtools maskfasta -fi ${fdir}/6_new-fasta/${species}.snp.fa -bed ${fdir}/${indel} -fo ${fdir}/6_new-fasta/${species}.in.fa
bedtools maskfasta -fi ${fdir}/6_new-fasta/${species}.in.fa -bed $zero -fo ${fdir}/6_new-fasta/${species}.fa

done
