#!/bin/bash
#===========================================================================================================================================
#
#                                  Annotating the VCF files with SnpSift, snpEff and annovar
#
#============================================================================================================================================


#Annotating with snpsift and snpEff
mkdir NEW_SIFT
java -jar /etc/ace-data/home/smthande/miniconda3/share/snpsift-4.3.1t-3/SnpSift.jar annotate  \         #the path to the SnpSift jar file 
../Reference/dbsnp_138.hg38.vcf\                                              #specifying the database to use for the annotation of the variants
 Combined_samples.vcf \                                                       #file to annotate(from HaplotypeCaller)
 > Combined_samples.sift.vcf                                                  #The output file


java -jar /etc/ace-data/home/smthande/miniconda3/share/snpeff-5.0-1/snpEff.jar ann -v -o vcf -c  #path to the jar file and the annotation options
../snpEff.config hg19 \                                                                          #path to the config.file
Combined_samples_sift.vcf \                                                                      #File to annotate
 > Combined_samples_snpeff.vcf                                                                   #Output file 

#Annoatting with annovar
mkdir Annovar_annotation
perl table_annovar.pl /etc/ace-data/home/smthande/Master_Thesis/annovar/Combined_samples.vcf.gz humandb/ -buildver hg19 -out Annovar_annotation/Combined_samples.vcf.gz\
 -remove -protocol \
refGene,exac03,gwasCatalog,AFR.sites.2014_10,ALL.sites.2014_10,AMR.sites.2014_10,EUR.sites.2014_10,SAS.sites.2014_10,tfbsConsSites,refGeneWithVer,gnomad211_genome -operation g,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput 

#It uses perl language , withe the various databases that are accesible through the annovar documantation site.
