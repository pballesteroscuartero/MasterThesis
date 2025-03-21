#!/bin/bash

#This code takes the RNAseq data in fastq format and performs basic preprocessing steps such as alignment, duplicate removal, and counting of reads.
#The steps follow GATK best practices for RNAseq data processing.


#Commands:
#working_dir_root: Directory where the generated data will be stored
#data_folder: Name of the folder within working_dir_root where the input files are stored
#toolds_dir_root: directory containing picard and gatk installed. It should also contain reference genome, STAR index of the genome, and dbSNP

working_dir_root=$1
data_folder=$2
tools_dir_root=$3
ref_genome=$tools_dir_root/genome_38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
genome_index=$tools_dir_root/genome_38/STAR_index
gtf=$tools_dir_root/genome_38/Homo_sapiens.GRCh38.111.gtf
dbSNP=$tools_dir_root/genome_38/dbSNP_homosapiens_renamed.vcf.gz

mkdir -p $working_dir_root/STAR_aligned
mkdir -p $working_dir_root/finals
mkdir -p $working_dir_root/counts
mkdir -p  $working_dir_root/MuTect


# set up output filenames and locations

for file in "$working_dir_root/$data_folder"/*
do 
    if [ -f "$file" ]
    then 
        echo "Processing file: $file"

        IFS='/' read -ra parts <<< "$file"
        file_name="${parts[-1]}"

        mkdir -p $working_dir_root/STAR_aligned/$file_name
        mkdir -p $working_dir_root/finals/$file_name
        mkdir -p $working_dir_root/counts/$file_name
        mkdir -p $working_dir_root/MuTect/$file_name
                
        align_out=$working_dir_root/STAR_aligned/$file_name/${file_name}_
        addorreplace_input_bam=$working_dir_root/STAR_aligned/$file_name/${file_name}_Aligned.sortedByCoord.out.bam
        duplicates_input_bam=$working_dir_root/STAR_aligned/$file_name/${file_name}_AddInfo.bam
        split_input_bam=$working_dir_root/STAR_aligned/$file_name/${file_name}_RemoveDuplicates.bam
        realign_input_bam=$working_dir_root/STAR_aligned/$file_name/${file_name}_Splitted.bam
        realign_table=$working_dir_root/STAR_aligned/$file_name/${file_name}_SJ.out.tab
        counts_input_bam=$working_dir_root/STAR_aligned/$file_name/${file_name}_BQSR.bam
        metrics_file=$working_dir_root/STAR_aligned/$file_name/${file_name}_metrics.txt
        counts=$working_dir_root/counts/$file_name/${file_name}_featurecounts.txt
        final=$working_dir_root/finals/$file_name/${file_name}_final.bam


        STAR --runThreadN 4 --genomeDir $genome_index --readFilesIn $file --outFileNamePrefix $align_out \
        --twopassMode Basic --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within --outSAMattributes NH HI NM MD AS --readFilesCommand zcat --quantMode TranscriptomeSAM
                
        java -jar $tools_dir_root/picard.jar AddOrReplaceReadGroups -I $addorreplace_input_bam -O $duplicates_input_bam -SO coordinate -RGID id -RGLB library -RGPL platform -RGPU machine -RGSM $file_name
        # Remove duplicates with Picard
        java -jar $tools_dir_root/picard.jar MarkDuplicates -INPUT $duplicates_input_bam -OUTPUT $split_input_bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -METRICS_FILE $metrics_file
        
        # SplitNCigarsReads for exon segment splitting of reads
        java -jar $tools_dir_root/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar  SplitNCigarReads -R $ref_genome -I $split_input_bam -O $realign_input_bam --tmp-dir $tools_dir_root/temp_files

        # Base realignment with gatk
        java -jar $tools_dir_root/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar BaseRecalibrator -I $realign_input_bam --known-sites $dbSNP -R $ref_genome -O $realign_table
        
        #Apply BQSR
        java -jar $tools_dir_root/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar ApplyBQSR -R $ref_genome -I $realign_input_bam --bqsr-recal-file $realign_table -O $counts_input_bam
        
        # Count mapped reads
        featureCounts -T 4 -s 2 -a $gtf -o $counts $counts_input_bam

        java -jar $tools_dir_root/picard.jar ReorderSam -I $counts_input_bam -O $final -R $ref_genome -CREATE_INDEX TRUE  -SEQUENCE_DICTIONARY $tools_dir_root/genome_38/Homo_sapiens.GRCh38.dna.primary_assembly.dict


        echo "File $fq processed!"

    else
        echo "File not found: $file"

    fi 
done
