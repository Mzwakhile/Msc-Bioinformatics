#!/bin/bash  
#SBATCH --job-name={}                                 #specifying the job name i.e Quality_Control
# SBATCH --nodes=1                                    #specifying the number of nodes to use for a particular process
#SBATCH --error=Calib.err                             # If you encounter an error , the error should be written here
#SBATCH --output=Calib.log                            # this is the file showing the output

#loading modules to be used for the tasks
#=================================================================================================
#
#                                       Loading Modules and conda environment activation
#
#=================================================================================================
module load "anaconda3-2019.10-gcc-9.3.0-7gh72na"
eval "$(conda shell.bash hook)"
module load "hisat2-2.1.0-gcc-9.3.0-ld56fsj"
module load "samtools-1.10-gcc-9.3.0-egzhq7w"
module load "vcftools-0.1.14-gcc-9.3.0-vdqtx"
conda activate /etc/ace-data/home/smthande/miniconda3/envs/Project      # path of the environment
conda activate /etc/ace-data/home/smthande/miniconda3/envs/Trimming      # path to the environment
module load gatk-4.2.2.0-gcc-8.5.0-7fb72tl
#=================================================================================================
#                                                                                                #
#                                      Samples download and the metadata file                    #
#                                                                                                #
#================================================================================================#
##Required in the working directory
# 1. sample_id.txt: contains a list of samples                         #Contains the ID of the sequences to be used in a for loop for example

# 2. Sample sequences (R1 and R2)                                      #The foward and reverse read
#==================================================================================================
#                                                                                                 #
#                       Data Preprocessing                                                        #
#                                                                                                 #
#==================================================================================================

echo -e "\n TB Data Preprocessing... \n"

mkdir -p QUALITY                                    #create directory for the fastqc output

Quality Check using 'FastQC'
for sample in `cat sample_id.txt`
do
       fastqc ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz --outdir QUALITY
done
Aggregate FastQC results using 'MultiQC'
multiqc QUALITY -o QUALITY

Quality Trimming with 'Trim Galore'
mkdir -p GALORE   #Directory for the output of Trim Galore

for sample in `cat sample_id.txt`
do
        trim_galore -j 10 --paired ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz -q 25 --length 20 --fastqc -o GALORE
done
# -j 10 specifies the number of cores to be used for trimming
#--paired runs a paired-end validation on both trimmed R1 and R2 FastQ files to removes entire read pairs if at least one of the two sequences became shorter than the $# -q 25 to remove portions of the reads with phred score less than 25 at the 3' end
# --length 20 to filter trimmed reads less than 20 bp


Quality Check the trimmed Reads with 'Multiqc'
multiqc GALORE -o GALORE                                        #An aggregate of the quality control for all the files

#====================================================================================================
#                                                                                                   #
#                             HAVE TO DO HISAT ALIGNMENT                                             #
#                                                                                                   #
#====================================================================================================

#Defining the number of threads to be used
threads=20
#Alignment of sequences to the reference genome using 'HISAT'

echo -e "\n Now Running Alignement using Hisat2... \n"

mkdir -p ref
# Download the human reference genome
wget -c ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/hg38.fa.gz -O ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gunzip  ref/hg38.fa.gz


HISAT2 indexing to build the reference genome index
hisat2-build -p $threads hisat2/hg38.fa hisat2/ref.idx

#HISAT2 Alignment
mkdir hisat2
for sample in `cat sample_id.txt`
do
      hisat2 -p $threads -x /etc/ace-data/home/smthande/Master_Thesis/Reference/HISAT2/ref.idx  \  #path to the Reference 
              -1 GALORE/${sample}_R1_001_val_1.fq.gz  -2 GALORE/${sample}_R2_001_val_2.fq.gz  \    # the foward and reverse trimmed  reads
              -S   hisat2/${sample}_hisat2.sam                                                     # output of alignement

  samtools view -Sb hisat2/${sample}_hisat2.sam  | samtools sort  > hisat2/${sample}_hisat2_sorted.bam  #compress the sam files to binary format
     samtools index hisat2/${sample}_hisat2_sorted.bam                                                   #Indexing the compressed files

     rm  hisat2/${sample}_hisat2.sam #remove the .sam files for storage purposes
done
#-x >The basename of the index for the reference genome.
#-1 > the foward reads of the samples
#-2 > the reverse reads of all the samples
#-S >File to write SAM alignments to, which by default its written to "stdout" or "standard out"

#=====================================================================================================
#########################################ADDING READ GROUPS TO THE BAM FLES###########################
~/miniconda3/envs/Trimming/bin/java                                                #java pathfile
mkdir -p Read_Group                                                                #Direcory for the output
for sample in `cat sample_id.txt`
do
gatk AddOrReplaceReadGroups \                                                       # GATK tool for adding read groups
     I=hisat2/${sample}_hisat2_sorted.bam \                                        #Input files
       O=Read_Group/${sample}_hisat2_read_group.bam \                               #Output directory for the files
      RGID=4 \                                                                      #ReadGroup ID
       RGLB=lib1 \                                                                  #ReadGroup Library
       RGPL=ILLUMINA \                                                              #ReadGroup Platform
       RGPU=unit1 \                                                                 #ReadGroup Process Units
       RGSM=2
done

#======================================================BASE CALIBRATION==========================
mkdir Calibration
for sample in `cat sample_id.txt`
do
gatk BaseRecalibrator \                                                             #GATK tool calling
   -I Read_Group/${sample}_hisat2_read_group.bam \                                  #Input files to be calibrated
   -R ../Reference/hg38.fa \                                                        # specifying the human genome reference
   --known-sites /etc/ace-data/home/smthande/Master_Thesis/Reference/1000G_omni2.5.hg38.vcf \  #path to the known site file to calibrate with
   --known-sites /etc/ace-data/home/smthande/Master_Thesis/Reference/dbsnp_138.hg38.vcf\       #path to another known site to calibrate with
   -O Calibration/${sample}_recal_data.table                                               # Calibration table 
done

mkdir BQSR
for sample in `cat BWC_id.txt`
do
gatk ApplyBQSR \                                                              #Tool for  Applying scores
  -R ../Reference/hg38.fa \                                                   #Reference genome path
  -I Read_Group/${sample}_hisat2_read_group.bam \                             #Input file
  --bqsr-recal-file Calibration/${sample}_recal_data.table \                  #calibration table to be used
   -O ${sample}_calibrated.bam                                                #output file which is calibrated bam file for downstream analysis
done


mkdir Bam_Stats
for sample in `cat UG.txt`
do
bamtools stats -in BQSR/${sample}_calibrated.bam >  Bam_Stats/${sample}_calibrated.bam.txt #obtaining the alignment statistics like number of reads mapped, singletons etc
done
