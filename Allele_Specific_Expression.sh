#!/bin/bash
mkdir SITES
for sample in `cat sample_id.txt`
do
bcftools view GVC/${sample}.g.vcf -m2 -M2 -v snps > SITES/${sample}_bi_sites.vcf #using the files produced in Haplotypecaller to filter biallelic sites only
bgzip SITES/${sample}_bi_sites.vcf                                               #zipping the biallelic file
bcftools index SITES/${sample}_bi_sites.vcf.gz                                   #Indexing the vcf file
tabix SITES/${sample}_bi_sites.vcf.gz                                            #indexing still
done
#Conducting the allele specific expression analysis

mkdir Allelic_Imbalance
for sample in `cat sample_id.txt`
gatk ASEReadCounter -I BQSR/${sample}_calibrated.bam \                           #Alignent file that is calibrated
                -V SITES/${sample}_bi_sites.vcf.gz \                             #input file which is the bam file
                   -R ../Reference/hg38.fa \                                     #pathway to the reference file
                --output-format CSV \                                            #output format is the csv file
                    -O Allelic_Imbalance/${sample}_AI.csv                        #Output files
done
