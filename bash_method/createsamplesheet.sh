#!/usr/bin/env bash

# Immediately stop on errors. keep track of the commands as they get executed. 
set -ue
# dont show commands
set +x




#####################################################################
############################### help message ########################
#####################################################################
help()
{
  echo -en "This shell script creates samplesheet.csv \n"
  echo -en  "\n"
  echo -en "Sintax: bash createsamplesheet.sh [-h|i] \n" 
  echo -en "options:\n"
  echo -en "--h    show help message\n"
  echo -en "--i    folder with fastq files - full_path/ \n"
  echo -en  "\n"
} 


#####################################################################
########## define parameters ########################################
#####################################################################
while getopts ":h:i:" flag
do
  case "${flag}" in
      h) help;;
      i) fastq_folder=$OPTARG;;
  esac
done



echo "sample" > "$fastq_folder"/sample-id.txt
echo "r1" > "$fastq_folder"/r1.txt
echo "r2" > "$fastq_folder"/r2.txt
for i in $(ls "$fastq_folder"/*fastq); do echo "$i" | awk -F _bp_ '{print $2}' | awk -F . '{print $1}' >> "$fastq_folder"/sample-id.txt; done
#for i in $(ls "$fastq_folder"/*fastq); do echo "$i" | awk -F . '{print $1}' >> "$fastq_folder"/output_"$outdir"/temp/sample-id.txt; done
find $(pwd)/"$fastq_folder"/*fastq >> "$fastq_folder"/r1.txt
paste -d , "$fastq_folder"/sample-id.txt "$fastq_folder"/r1.txt "$fastq_folder"/r2.txt > "$fastq_folder"/samplesheet.csv
