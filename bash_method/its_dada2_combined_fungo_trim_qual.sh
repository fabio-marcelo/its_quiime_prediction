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
  echo -en "This shell script runs taxonomy classification using qiime2 \n"
  echo -en  "\n"
  echo -en "Sintax: bash its_dada2_iontorrent.sh [-h|i|p|r] \n" 
  echo -en "options:\n"
  echo -en "-h    show help message\n"
  echo -en "-i    folder with fastq files - full_path/ \n"
  echo -en "-p    primer sequence \n"
  echo -en "-r    full_path/to/trainned_classifier.qza \n"
  echo -en "-a    full_path/to/reference_reads.fasta \n"
  echo -en "-t    full_path/to/tax_file.txt \n"
  echo -en "-o    output folder name \n"
  echo -en  "\n"
} 


#####################################################################
########## define parameters ########################################
#####################################################################
while getopts ":h:i:p:r:a:t:o:" flag
do
  case "${flag}" in
      h) help;;
      i) fastq_folder=$OPTARG;;
      p) primer_seq=$OPTARG;;
      r) trainned_classifier=$OPTARG;;
      a) ref_reads=$OPTARG;;
      t) tax_file=$OPTARG;;
      o) outdir=$OPTARG;;
  esac
done


echo "set variables"
echo "fastq_folder: $fastq_folder";
echo "primer_seq: $primer_seq";
echo "trainned classifier: $trainned_classifier";
echo "reference reads: $ref_reads";
echo "tax file: $tax_file";
echo "output dir name: $outdir";
echo -e "\n"



#####################################################################
########## Main program #############################################
#####################################################################
# cria pasta necessarias
echo -e "set infraestructure \n"
echo -e "\n"
mkdir -p "$fastq_folder"/output_"$outdir"
mkdir -p "$fastq_folder"/output_"$outdir"/temp
mkdir -p "$fastq_folder"/output_"$outdir"/log

# setup file
echo -e "write setup file \n"
echo -e "\n"
echo "fastq_folder: $fastq_folder" >> "$fastq_folder"/output_"$outdir"/temp/setup_file.txt
echo "primer_seq: $primer_seq" >> "$fastq_folder"/output_"$outdir"/temp/setup_file.txt
echo "trainned classifier: $trainned_classifier" "$fastq_folder"/output_"$outdir"/temp/setup_file.txt
echo "reference reads: $ref_reads" >> "$fastq_folder"/output_"$outdir"/temp/setup_file.txt
echo "tax file: $tax_file" >> "$fastq_folder"/output_"$outdir"/temp/setup_file.txt
echo "output dir name: $outdir" >> "$fastq_folder"/output_"$outdir"/temp/setup_file.txt



# create manifest file
echo -e "create manifest-file \n"
echo -e "\n"
echo "sample-id" > "$fastq_folder"/output_"$outdir"/temp/sample-id.txt
echo "absolute-filepath" > "$fastq_folder"/output_"$outdir"/temp/filepath.txt
for i in $(ls "$fastq_folder"/*fastq); do echo "$i" | awk -F _bp_ '{print $2}' | awk -F . '{print $1}' >> "$fastq_folder"/output_"$outdir"/temp/sample-id.txt; done
#for i in $(ls "$fastq_folder"/*fastq); do echo "$i" | awk -F . '{print $1}' >> "$fastq_folder"/output_"$outdir"/temp/sample-id.txt; done
find $(pwd)/"$fastq_folder"/*fastq >> "$fastq_folder"/output_"$outdir"/temp/filepath.txt
paste "$fastq_folder"/output_"$outdir"/temp/sample-id.txt "$fastq_folder"/output_"$outdir"/temp/filepath.txt > "$fastq_folder"/output_"$outdir"/temp/manifest-file.tsv
echo -e "manifest-file created \n"


# create metadata file
echo -e "create metadata-file \n"
echo "SampleID" > "$fastq_folder"/output_"$outdir"/temp/SampleID.txt
echo "BarcodeSequence" > "$fastq_folder"/output_"$outdir"/temp/BarcodeSequence.txt
echo "LinkerPrimerSequence" > "$fastq_folder"/output_"$outdir"/temp/LinkerPrimerSequence.txt
echo "Path" > "$fastq_folder"/output_"$outdir"/temp/Path.txt
for i in $(ls "$fastq_folder"/*fastq); do echo "$i" | awk -F _bp_ '{print $2}' | awk -F . '{print $1}' >> "$fastq_folder"/output_"$outdir"/temp/SampleID.txt; done
find $(pwd)/"$fastq_folder"/*fastq >> "$fastq_folder"/output_"$outdir"/temp/Path.txt
paste "$fastq_folder"/output_"$outdir"/temp/SampleID.txt "$fastq_folder"/output_"$outdir"/temp/BarcodeSequence.txt "$fastq_folder"/output_"$outdir"/temp/LinkerPrimerSequence.txt "$fastq_folder"/output_"$outdir"/temp/Path.txt > "$fastq_folder"/output_"$outdir"/temp/metadata-file.txt
echo -e "metadata-file created \n"
echo -e "\n"


# import fastq files
echo "import fastq files"
qiime tools import --type 'SampleData[SequencesWithQuality]' \
--input-path "$fastq_folder"/output_"$outdir"/temp/manifest-file.tsv \
--output-path "$fastq_folder"/output_"$outdir"/temp/fastq_imported.qza \
--input-format SingleEndFastqManifestPhred33V2 2> "$fastq_folder"/output_"$outdir"/log/import-files_erro.txt
echo "finishing import fastq files"
echo -e "\n"


# visualize fastq files imported
echo "visualize fastq files imported"
qiime demux summarize \
--i-data "$fastq_folder"/output_"$outdir"/temp/fastq_imported.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/view-inspec_import.qzv 2> "$fastq_folder"/output_"$outdir"/log/visualize-imported_erro.txt
echo "open qzv file in https://view.qiime2.org/"
echo -e "\n"

# trimm primers
echo "start trimming"
qiime cutadapt trim-single \
--i-demultiplexed-sequences "$fastq_folder"/output_"$outdir"/temp/fastq_imported.qza \
--p-cores $(grep -c processor /proc/cpuinfo) \
--p-front "$primer_seq" \
--p-quality-cutoff-5end 20 \
--p-quality-cutoff-3end 20 \
--p-error-rate 0.2 \
--p-minimum-length 80 \
--o-trimmed-sequences "$fastq_folder"/output_"$outdir"/trimmed-seqs.qza \
--verbose 2> "$fastq_folder"/output_"$outdir"/log/trimming_erro.txt
echo "finishing trimming"
echo -e "\n"

qiime demux summarize \
--i-data "$fastq_folder"/output_"$outdir"/trimmed-seqs.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/view-trimmed_sequences.qzv


# denoising
## com a inclusao de --p-quality-cutoff-5end e --p-quality-cutoff-3end, 
## "--p-trim-left 6 \" foi removido.
echo "start denoising"
qiime dada2 denoise-single \
--p-n-threads $(grep -c processor /proc/cpuinfo) \
--i-demultiplexed-seqs "$fastq_folder"/output_"$outdir"/trimmed-seqs.qza \
--p-max-ee 2 \
--p-trunc-len 0 \
--p-pooling-method 'independent' \
--p-chimera-method 'consensus' \
--o-representative-sequences "$fastq_folder"/output_"$outdir"/representative-seqs.qza \
--o-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
--o-denoising-stats "$fastq_folder"/output_"$outdir"/denoise-stats.qza 2> "$fastq_folder"/output_"$outdir"/log/denoising_erro.txt
echo "finishing denoising"
echo -e "\n"

# visualize denoising
echo "visualize denoising"
qiime metadata tabulate \
--m-input-file "$fastq_folder"/output_"$outdir"/denoise-stats.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/view-inspect_denoise-stats.qzv 2> "$fastq_folder"/output_"$outdir"/log/visualize-denoised_erro.txt

# Denoise stats to tsv
echo "set denoise stats tsv file"
qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/denoise-stats.qza \
--output-path "$fastq_folder"/output_"$outdir"/denoise-stats_tsv
echo "open qzv file in https://view.qiime2.org/"
echo -e "\n"

# Feature table summarize
qiime feature-table summarize \
--i-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/table-denoised_summary.qzv


# Generate rarefaction curves
## --p-max-depth     you can only rarefy to as many features are present in your most populated sample.
#qiime diversity alpha-rarefaction \
#--i-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
#--p-max-depth 36298 \
#--p-steps 20 \
#--m-metadata-file "$fastq_folder"/output_"$outdir"/temp/metadata-file.txt \
#--o-visualization "$fastq_folder"/output_"$outdir"/rarefaction_plot.qzv



## New
#qiime feature-table summarize \
#--i-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
#--m-sample-metadata-file metadata.tsv \
#--o-visualization analysis/visualisations/16s_table.qzv \
#--verbose



echo "Import seq reference file"
qiime tools import \
--type FeatureData[Sequence] \
--input-path "$ref_reads" \
--output-path "$fastq_folder"/output_"$outdir"/temp/reference_sequences.qza 2> "$fastq_folder"/output_"$outdir"/log/import_seq_file_erro.txt
echo "seq reference file imported"
echo -e "\n"

echo "Import taxonomy reference file"
qiime tools import \
--type FeatureData[Taxonomy] \
--input-path "$tax_file" \
--output-path "$fastq_folder"/output_"$outdir"/temp/reference_taxonomy.qza \
--input-format HeaderlessTSVTaxonomyFormat 2> "$fastq_folder"/output_"$outdir"/log/import_tax_file_erro.txt
echo "seq reference file imported"
echo -e "\n"

# taxonomy classification
echo "Starting Taxonomic identification - vsearch"
qiime feature-classifier classify-consensus-vsearch \
--i-query "$fastq_folder"/output_"$outdir"/representative-seqs.qza \
--i-reference-reads "$fastq_folder"/output_"$outdir"/temp/reference_sequences.qza \
--i-reference-taxonomy "$fastq_folder"/output_"$outdir"/temp/reference_taxonomy.qza \
--p-perc-identity 0.99 \
--p-threads 18 \
--o-classification "$fastq_folder"/output_"$outdir"/vsearch-taxonomyITS.qza \
--output-dir "$fastq_folder"/output_"$outdir"/temp/vsearch 2> "$fastq_folder"/output_"$outdir"/log/vsearch_log.txt
echo -e "\n"

qiime metadata tabulate \
--m-input-file "$fastq_folder"/output_"$outdir"/vsearch-taxonomyITS.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/view-vsearch-taxonomyITS.qzv

echo "Starting Taxonomic identification - blast"
qiime feature-classifier classify-consensus-blast \
--p-perc-identity 0.99 \
--i-query "$fastq_folder"/output_"$outdir"/representative-seqs.qza \
--i-reference-reads "$fastq_folder"/output_"$outdir"/temp/reference_sequences.qza \
--i-reference-taxonomy "$fastq_folder"/output_"$outdir"/temp/reference_taxonomy.qza \
--o-classification "$fastq_folder"/output_"$outdir"/blast-taxonomyITS.qza \
--verbose \
--o-search-results "$fastq_folder"/output_"$outdir"/temp/blast 2> "$fastq_folder"/output_"$outdir"/log/blast_log.txt
echo -e "\n"

qiime metadata tabulate \
--m-input-file "$fastq_folder"/output_"$outdir"/blast-taxonomyITS.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/view-blast-taxonomyITS.qzv

#echo "Start trainning classifier"
#time qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads "$fastq_folder"/output_"$outdir"/temp/reference_sequences.qza \
#--i-reference-taxonomy "$fastq_folder"/output_"$outdir"/temp/reference_taxonomy.qza \
#--o-classifier "$fastq_folder"/output_"$outdir"/temp/trainned_classifier_qiime_release_s_29.11.2022.qza  

echo "Starting Taxonomic identification - sklearn"
qiime feature-classifier classify-sklearn \
--p-n-jobs 1 \
--p-confidence 0.99 \
--i-classifier "$trainned_classifier" \
--p-reads-per-batch 10000 \
--i-reads "$fastq_folder"/output_"$outdir"/representative-seqs.qza \
--verbose \
--o-classification "$fastq_folder"/output_"$outdir"/sklearn-taxonomyITS.qza 2> "$fastq_folder"/output_"$outdir"/log/sklearn_log.txt
echo -e "\n"

qiime metadata tabulate \
--m-input-file "$fastq_folder"/output_"$outdir"/sklearn-taxonomyITS.qza \
--o-visualization "$fastq_folder"/output_"$outdir"/view-sklearn-taxonomyITS.qzv
echo "open qzv file in https://view.qiime2.org/"
echo "Finishing Taxonomic identification"
echo -e "\n"

# make table
## export taxa
echo "Starting export"
qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/table-denoised.qza \
--output-path "$fastq_folder"/output_"$outdir"/feature-table

qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/vsearch-taxonomyITS.qza \
--output-path "$fastq_folder"/output_"$outdir"/vsearch-tabela-formato-tabular

qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/blast-taxonomyITS.qza \
--output-path "$fastq_folder"/output_"$outdir"/blast-tabela-formato-tabular

qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/sklearn-taxonomyITS.qza \
--output-path "$fastq_folder"/output_"$outdir"/sklearn-tabela-formato-tabular

## export taxonomy
qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/vsearch-taxonomyITS.qza \
--output-path "$fastq_folder"/output_"$outdir"/vsearch-taxonomy

qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/blast-taxonomyITS.qza \
--output-path "$fastq_folder"/output_"$outdir"/blast-taxonomy

qiime tools export \
--input-path "$fastq_folder"/output_"$outdir"/sklearn-taxonomyITS.qza \
--output-path "$fastq_folder"/output_"$outdir"/sklearn-taxonomy

## replace header
cp "$fastq_folder"/output_"$outdir"/vsearch-taxonomy/*.tsv "$fastq_folder"/output_"$outdir"/vsearch-taxonomy/vsearch-taxonomy_header.tsv
sed -i 's/Feature ID/#otu-id/g' "$fastq_folder"/output_"$outdir"/vsearch-taxonomy/vsearch-taxonomy_header.tsv
sed -i 's/Taxon/taxonomy/g' "$fastq_folder"/output_"$outdir"/vsearch-taxonomy/vsearch-taxonomy_header.tsv

cp "$fastq_folder"/output_"$outdir"/blast-taxonomy/*.tsv "$fastq_folder"/output_"$outdir"/blast-taxonomy/blast-taxonomy_header.tsv
sed -i 's/Feature ID/#otu-id/g' "$fastq_folder"/output_"$outdir"/blast-taxonomy/blast-taxonomy_header.tsv
sed -i 's/Taxon/taxonomy/g' "$fastq_folder"/output_"$outdir"/blast-taxonomy/blast-taxonomy_header.tsv

cp "$fastq_folder"/output_"$outdir"/sklearn-taxonomy/*.tsv "$fastq_folder"/output_"$outdir"/sklearn-taxonomy/sklearn-taxonomy_header.tsv
sed -i 's/Feature ID/#otu-id/g' "$fastq_folder"/output_"$outdir"/sklearn-taxonomy/sklearn-taxonomy_header.tsv
sed -i 's/Taxon/taxonomy/g' "$fastq_folder"/output_"$outdir"/sklearn-taxonomy/sklearn-taxonomy_header.tsv

## add metadata
biom add-metadata \
--input-fp "$fastq_folder"/output_"$outdir"/feature-table/feature-table.biom \
--observation-metadata-fp "$fastq_folder"/output_"$outdir"/vsearch-taxonomy/vsearch-taxonomy_header.tsv \
--output-fp "$fastq_folder"/output_"$outdir"/vsearch-biom-with-taxonomy.biom

biom add-metadata \
--input-fp "$fastq_folder"/output_"$outdir"/feature-table/feature-table.biom \
--observation-metadata-fp "$fastq_folder"/output_"$outdir"/blast-taxonomy/blast-taxonomy_header.tsv \
--output-fp "$fastq_folder"/output_"$outdir"/blast-biom-with-taxonomy.biom

biom add-metadata \
--input-fp "$fastq_folder"/output_"$outdir"/feature-table/feature-table.biom \
--observation-metadata-fp "$fastq_folder"/output_"$outdir"/sklearn-taxonomy/sklearn-taxonomy_header.tsv \
--output-fp "$fastq_folder"/output_"$outdir"/sklearn-biom-with-taxonomy.biom

## convert biom to text
biom convert \
--input-fp "$fastq_folder"/output_"$outdir"/vsearch-biom-with-taxonomy.biom \
--output-fp "$fastq_folder"/output_"$outdir"/vsearch-biom-with-taxonomy.tsv \
--to-tsv \
--observation-metadata-fp "$fastq_folder"/output_"$outdir"/vsearch-taxonomy/vsearch-taxonomy_header.tsv \
--header-key taxonomy

biom convert \
--input-fp "$fastq_folder"/output_"$outdir"/blast-biom-with-taxonomy.biom \
--output-fp "$fastq_folder"/output_"$outdir"/blast-biom-with-taxonomy.tsv \
--to-tsv \
--observation-metadata-fp "$fastq_folder"/output_"$outdir"/blast-taxonomy/blast-taxonomy_header.tsv \
--header-key taxonomy

biom convert \
--input-fp "$fastq_folder"/output_"$outdir"/sklearn-biom-with-taxonomy.biom \
--output-fp "$fastq_folder"/output_"$outdir"/sklearn-biom-with-taxonomy.tsv \
--to-tsv \
--observation-metadata-fp "$fastq_folder"/output_"$outdir"/sklearn-taxonomy/sklearn-taxonomy_header.tsv \
--header-key taxonomy

## print columns except column one
cat "$fastq_folder"/output_"$outdir"/blast-biom-with-taxonomy.tsv | cut --complement -f1 > "$fastq_folder"/output_"$outdir"/blast-final_output.csv
cat "$fastq_folder"/output_"$outdir"/vsearch-biom-with-taxonomy.tsv | cut --complement -f1 > "$fastq_folder"/output_"$outdir"/vsearch-final_output.csv
cat "$fastq_folder"/output_"$outdir"/sklearn-biom-with-taxonomy.tsv | cut --complement -f1 > "$fastq_folder"/output_"$outdir"/sklearn-final_output.csv
echo "Finishing export"

## create a barplot
echo -e "create a barplot \n"
qiime taxa barplot \
--i-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
--i-taxonomy "$fastq_folder"/output_"$outdir"/vsearch-taxonomyITS.qza \
--m-metadata-file "$fastq_folder"/output_"$outdir"/temp/metadata-file.txt \
--o-visualization "$fastq_folder"/output_"$outdir"/vsearch-taxa-bar-plots-filtered.qzv 2> "$fastq_folder"/output_"$outdir"/log/create-barplot_erro.txt
echo -e "\n"

qiime taxa barplot \
--i-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
--i-taxonomy "$fastq_folder"/output_"$outdir"/blast-taxonomyITS.qza \
--m-metadata-file "$fastq_folder"/output_"$outdir"/temp/metadata-file.txt \
--o-visualization "$fastq_folder"/output_"$outdir"/blast-taxa-bar-plots-filtered.qzv 2> "$fastq_folder"/output_"$outdir"/log/create-barplot_erro.txt

qiime taxa barplot \
--i-table "$fastq_folder"/output_"$outdir"/table-denoised.qza \
--i-taxonomy "$fastq_folder"/output_"$outdir"/sklearn-taxonomyITS.qza \
--m-metadata-file "$fastq_folder"/output_"$outdir"/temp/metadata-file.txt \
--o-visualization "$fastq_folder"/output_"$outdir"/sklearn-taxa-bar-plots-filtered.qzv 2> "$fastq_folder"/output_"$outdir"/log/create-barplot_erro.txt
echo -e "barplot created \n"
echo -e "open qzv file in https://view.qiime2.org/ \n"
echo -e "Pipeline finished"


