#!/bin/bash

#User should have PERL5LIB running

working_dir_root=$1
tools_dir_root=$2
input_vep_folder=$working_dir_root/MuTect

for file in "$input_vep_folder"/*
do 
    echo "Processing file: $file"

    IFS='/' read -ra parts <<< "$file"
    file_name="${parts[-1]}"

    mkdir -p $working_dir_root/VEP/$file_name
    vep_in=$working_dir_root/MuTect/$file_name/${file_name}_final.vcf
    vep_out_raw=$working_dir_root/VEP/$file_name/${file_name}_vep_out_raw.txt
    vep_out_filtered=$working_dir_root/VEP/$file_name/${file_name}_vep_out_filtered.txt

    $tools_dir_root/ensembl-vep/vep -i $vep_in -o $vep_out_raw --species homo_sapiens --cache --dir_cache $tools_dir_root --uniprot --protein --symbol --transcript_version --tsl --offline --tab --sift p --polyphen p --hgvsg --check_existing
    
    $tools_dir_root/ensembl-vep/filter_vep -i $vep_out_raw -o $vep_out_filtered -filter "(IMPACT is HIGH) or (IMPACT is MODERATE)" --force_overwrite

    echo "VEP filtered"
done
