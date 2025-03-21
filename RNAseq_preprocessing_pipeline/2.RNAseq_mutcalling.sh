#!/bin/bash

#This code takes an aligned BAM file and finds the somatic mutations identified in the patients, and the normalized gene expressions of each feature

#Some paths
working_dir_root=$1
tools_dir_root=$2
ref_genome=$tools_dir_root/genome_38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=$tools_dir_root/genome_38/Homo_sapiens.GRCh38.111.gtf
data_folder=$working_dir_root/finals
rsem_out=$working_dir_root/ref_rsem


#It takes time so do only once
rsem-prepare-reference --gtf $tools_dir_root/genome_38/Homo_sapiens.GRCh38.111.gtf --star \
 $tools_dir_root/genome_38/Homo_sapiens.GRCh38.dna.primary_assembly.fa $rsem_out/ref_rsem

for file in $data_folder/*
do
    
    echo "Processing file: $file"
    IFS='/' read -ra parts <<< "$file"
    file_name="${parts[-1]}"

    input_file=$file/${file_name}_final.bam
    if [ -f ${input_file} ]
    then
        mkdir -p $working_dir_root/MuTect/$file_name

        mutect_out=$working_dir_root/MuTect/$file_name/${file_name}.vcf
        mutect_filter=$working_dir_root/MuTect/$file_name/${file_name}_filtered.vcf
        mutect_final=$working_dir_root/MuTect/$file_name/${file_name}_final.vcf
        rsem_input=$working_dir_root/STAR_aligned/$file_name/${file_name}_Aligned.toTranscriptome.out.bam

        java -jar $tools_dir_root/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 -R $ref_genome -I $input_file \
         --germline-resource $tools_dir_root/genome_38/germline_homosapiens_renamed.vcf.gz --panel-of-normals $tools_dir_root/genome_38/pon_homosapiens_renamed.vcf.gz \
          --dont-use-soft-clipped-bases -O $mutect_out

        java -jar $tools_dir_root/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar FilterMutectCalls -V $mutect_out -R $ref_genome -O $mutect_filter
        grep -E "PASS|##|#" $mutect_filter > $mutect_final

        rsem-calculate-expression --alignments $rsem_input $rsem_out/ref_rsem ${file_name}_rsem

    fi
done





