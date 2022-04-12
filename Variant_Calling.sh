#!/bin/bash
module load gatk-4.2.2.0-gcc-8.5.0-7fb72tl
~/miniconda3/envs/Trimming/bin/jar
~/miniconda3/envs/Trimming/bin/java
conda activate /etc/ace-data/home/smthande/miniconda3/envs/Trimming
#======================================================================================================
#
#                                    Using HaplotypeCaller for varaint calling
#
#======================================================================================================

mkdir GVC
for sample in `cat sample_id.txt`
do
gatk --java-options "-Xmx4g" HaplotypeCaller  \                           #calling the tool for variant calling
 -R ../Reference/hg38.fa \                                                #path to the reference genome
  -I BQRS/${sample}_calibrated.bam \                                      #path to the already produced calibrated file from data preprocessing steps
  --output GVC/${sample}.g.vcf.gz \                                       #output which is a genomic vcf
  -ERC GVCF                                                               #specifying that it should be a gvcf
 done
#combing the files(joint variant calling) into one file. It produces <NON REF> in the alternate allele position
gatk CombineGVCFs \                                                      #Combining tool
  -R ../Reference/hg38.fa \                                              #reference pathway
   --variant GVC/sample_1.g.vcf.gz \                                     #sample 1
  --variant GVC/sample_2.g.vcf.gz \                                      #sample 2
  --variant GVC/sample_3.g.vcf.gz\                                       #sample 3
    --variant GVC/sample_n.g.vcf.gz                                      #any number of samples can be added
  --output Combined_samples.g.vcf.gz                                     #the output of the combined vcf files
#Genotyping produces the normal vcf  file that will be used for further analysis
gatk --java-options "-Xmx4g" GenotypeGVCFs \                             #GenotypeGVCF is used for convertion
   -R ../Reference/hg38.fa \                                             #Reference files
  -V  Combined_samples.g.vcf.gz \                                        #File to be converted
 -O Combined_samples.vcf                                                 #The output file
