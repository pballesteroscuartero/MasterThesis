library(dplyr)
library(stringr)
library(rstatix)
library(ggplot2)
library(tibble)
library(tidyr)
library(ggpubr)
library(discover)

output_path "../outputs"
filter_vep = function(input_path) {
  folders = list.dirs(input_path)[-1]
  colnames = c("Uploaded_variation",	"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",	"cDNA_position",	"CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"IMPACT",	"DISTANCE",	"STRAND",	"FLAGS",	"SYMBOL",	"SYMBOL_SOURCE",	"HGNC_ID",	"TSL",	"ENSP",	"SWISSPROT",	"TREMBL",	"UNIPARC",	"UNIPROT_ISOFORM",	"SIFT",	"PolyPhen",	"HGVSg",	"CLIN_SIG",	"SOMATIC",	"PHENO")
  vep_combined = data.frame()
  genes_remove = c()
  
  #Concat all the VEP mutations so that we can later filter them out with the COSMIC Cancer mutation census
  for (i in folders){
    file_name = tail(str_split(i, "/")[[1]], n=1)
    sample_id = str_split(file_name, "_")[[1]][1]
    vep_path = paste0(i, "/", file_name, "_vep_out_filtered.txt")
    vep_file = read.table(vep_path, sep = "\t", comment.char = "#")
    colnames(vep_file)=colnames
    vep_file$chromosome = unlist(lapply(vep_file$Location, function(x) str_split(x, ":")[[1]][1]))
    vep_file$sample = sample_id
    vep_remove = vep_file %>%
      subset(chromosome == "Y" | chromosome == "MT")
    genes_remove = append(genes_remove, unique(vep_remove$SYMBOL))
    vep_filt = vep_file %>%
      subset(!SYMBOL %in% genes_remove)
    vep_combined = rbind(vep_combined, vep_filt)
  }
  
  #Filter out mutants with no known existing variants
  vep_census = vep_combined %>%
    subset(Existing_variation != "-") #We remove 129765 variants that were not reported in literature
  
  #Keep only variants known in the large intestine
  vep_census_li = vep_census %>%
    subset(HGVSg %in% census_mut_filt$HGVSG) #We remove 91414 variants that are not listed as large intestine specific
  
  return(vep_census_li)
}

create_mutation_matrix = function(input_df, mutation_type, patient_threshold) {
  #Mutation must be snv or deletion
  mutation_matrix = data.frame()
  
  for (i in unique(input_df$sample)){
    vep_file = input_df %>%
      subset(sample == i) %>%
      subset(mut_type == mutation_type)
    mut_gene = unique(vep_file$SYMBOL)
    vect = rep(1, length(mut_gene)) %>%
      t() %>%
      data.frame()
    colnames(vect)=mut_gene
    rownames(vect) = i
    mutation_matrix = dplyr::bind_rows(mutation_matrix, vect)
  }
  
  mutation_matrix[is.na(mutation_matrix)] = 0
  mutations_filt = mutation_matrix_snv[, colSums(mutation_matrix) >= patient_threshold]
  
  return(mutations_filt)
  
  
}

census_mut = read.table("../Cosmic_MutantCensus_Tsv_v99_GRCh38_largeintestine.tar", sep = "\t", header = T)
census_mut_filt = census_mut %>%
  subset(MUTATION_SOMATIC_STATUS == "Confirmed somatic variant" | MUTATION_SOMATIC_STATUS == "Reported in another cancer sample as somatic")
mut_genes = unique(census_mut$GENE_SYMBOL)

vep_census_li = filter_vep("Y:/Data/Masters/Pilar/CRC/VEP")
vep_snv_indel = vep_census_li %>%
  dplyr::distinct(Uploaded_variation, sample, .keep_all = TRUE)%>%
  mutate(amino_original = lapply(Amino_acids, function(x) str_split(x, "/")[[1]][1])) %>%
  mutate(amino_change = lapply(Amino_acids, function(x) str_split(x, "/")[[1]][2])) %>%
  mutate(mut_type = ifelse(amino_change == "*", "deletion", "snv")) %>%
  data.frame() 

patient_threshold = 3
mutation_matrix_snv = create_mutation_matrix(vep_snv_indel, "snv", patient_threshold)
mutation_matrix_del = create_mutation_matrix(vep_snv_indel, "deletion", patient_threshold)

saveRDS(mutations_filt_snv, paste0(output_path, "binary_mutations_",patient_threshold, "_or_more_patients_snv.rds"))
saveRDS(mutations_filt_del, paste0(output_path, "binary_mutations_",patient_threshold, "_or_more_patients_del.rds"))


