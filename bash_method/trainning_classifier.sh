#!/usr/bin/env bash

# Immediately stop on errors. keep track of the commands as they get executed. 
set -ue
# dont show commands
set +x


# This bash file runs analysis for fungal ITS sequencing data generated in Ion S platform
# through qiime2 using dada2

#####################################################################
############################### help message ########################
#####################################################################
help()
{
  echo -en "This shell script runs taxonomy classification using qiime2\n through fit-classifier-naive-bayes\n"
  echo -en  "\n"
  echo -en "Intermediate files will be stored in tmp folder\n"
  echo -en "Sintax: bash trainning_classifier.sh [-help|i|p|r]\n" 
  echo -en  "\n"
  echo -en "options:\n"
  echo -en "--help    show help message\n"
  echo -en "-f    Full path to folder containing reference files for sequences and taxonomy\n"
  echo -en "-s    FeatureData[Sequence] - fasta file file\n"
  echo -en "-t    FeatureData[Taxonomy] - file with taxonomy assignment\n"
  echo -en "-o    Output filename. The filename must have the extension .qza\n"
  echo -en  "\n"
} 

#####################################################################
########## define parameters ########################################
#####################################################################
while getopts ":h:f:s:t:o:" flag
do
  case "${flag}" in
      h) help;;
      f) folder=$OPTARG;;
      s) sequences=$OPTARG;;
      t) taxonomy=$OPTARG;;
      o) output_file=$OPTARG;;
  esac
done


echo "Folder: $folder";
echo "FeatureData[Sequence]: $sequences";
echo "FeatureData[Taxonomy]: $taxonomy";
echo "output file: $output_file";



#####################################################################
########## Main program #############################################
#####################################################################
mkdir -p "$folder"/output_folder
mkdir -p "$folder"/tmp

# train classifier
echo "Import seq reference file"
qiime tools import \
--type FeatureData[Sequence] \
--input-path "$folder"/"$sequences" \
--output-path "$folder"/tmp/reference_sequences.qza

echo "Import taxonomy reference file"
qiime tools import \
--type FeatureData[Taxonomy] \
--input-path "$folder"/"$taxonomy" \
--output-path "$folder"/tmp/reference_taxonomy.qza \
--input-format HeaderlessTSVTaxonomyFormat

echo "Start trainning classifier"
time qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads "$folder"/tmp/reference_sequences.qza \
--i-reference-taxonomy "$folder"/tmp/reference_taxonomy.qza \
--o-classifier "$folder"/output_folder/"$output_file"
