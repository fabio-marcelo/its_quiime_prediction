 #! /usr/bin/env nextflow
 nextflow.enable.dsl=2

 // Define help message
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:

        nextflow run main.nf --primer_seq 'primer_sequence' --fastq_folder '/path/to/my/fastq_files' \
          --trainned_classifier 'full_path/and_name/to/my/trainned_classifier' \
          --ref_reads 'full_path/to/reference_reads.fasta' \
          --tax_file 'full_path/to/tax_file.txt' \
          --outdir 'output-folder-name' \
          --threads "15"

        If using default parameters, you inly need to provide "fastq_folder" and a csv file (named samplesheet.csv)
        containing the sample names in the first column "named "sample", reads1 in the second column named "r1" and
        a third column empty column "r2".
  
  
  

        Mandatory arguments:
         --fastq_folder                 Folder with fastq files (full path required)
         --primer_seq                   Primer sequence ("AACTCCG") [must be surrounded with quotes]
                                        [Default: "GGAAGTAAAAGTCGTAACAAGG"]
         --trainned_classifier          full_path/to/trainned_classifier.qza (full path required)
                                        [Default: trainned_qiime-2023.7_ver9_99_s_all_29.11.2022_dev.qza.gz]
         --ref_reads                    Full path to reference reads file (.fasta) (full path required)
                                        [Default: sh_refs_qiime_ver9_99_s_all_29.11.2022_dev.fasta.gz]
         --tax_file                     Full path to taxonomy file (.txt) (full path required)
                                        [Default: sh_taxonomy_qiime_ver9_99_s_all_29.11.2022_dev.txt.gz]
         --outdir                       Output folder name

       Optional arguments:
        --threads                      Number of CPUs to use [Default: 15]

          qiime cutadapt trim-single arguments
          --p_quality_cutoff_5end        Quality cut off in 5' end [Default: 20]
          --p_quality_cutoff_3end        Quality cut off in 3' end [Default: 20]
          --p_error_rate                 Maximum allowed error rate [Default: 0.2]
          --p_minimum_length             Discard reads shorter than specified value. Note, the 
                                         cutadapt default of 0 has been overridden, because
                                         that value produces empty sequence records [Default: 80]
          
          qiime dada2 denoise-single arguments
          --p_max_ee                     Reads with number of expected errors higher than
                                         this value will be discarded [Default: 2] 
          --p_trunc_len                  Position at which sequences should be truncated due
                                         to decrease in quality. This truncates the 3' end of
                                         the of the input sequences, which will be the bases
                                         that were sequenced in the last cycles. Reads that
                                         are shorter than this value will be discarded. If 0
                                         is provided, no truncation or length filtering will
                                         be performed [Default: 0] 
        
        """
        .stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


process SAMPLESHEET{
    publishDir "${params.outdir}", mode:'copy', pattern: "samplesheet.csv"

    input:
        path(folder) 

    output:
        path ("samplesheet.csv"), emit: samplesheet
        path ("sample.txt")
        path ("r1.txt")
        path ("r2.txt")

    script:
    """
    touch sample.txt
    echo "sample" >> sample.txt
    touch r1.txt
    echo "r1" >> r1.txt
    touch r2.txt
    echo "r2" >> r2.txt
    touch samplesheet.csv
    for i in "\$(ls ${params.fastq_folder}/*fastq)"; do echo "\$i" | awk -F _280_ '{print \$2}' | awk -F . '{print \$1}' >> sample.txt; done
    find "\$(pwd)"/${params.fastq_folder}/*fastq >> r1.txt
    paste -d , sample.txt r1.txt r2.txt > samplesheet.csv
    """
    
}


process FASTQC {
    publishDir "${params.outdir}/fastqc_dir", mode:'copy'
    tag "$meta.id"

    input:
        tuple val(meta), path(reads)

    output:
      //file("*.{html,zip}")
      path("*.html")
      path("*.zip"), emit: zip

    script:
    
    // The meta map holds the single/paired information that can now easily be parsed in the module in a standardized way
    // cat ${reads} | fastqc stdin:${meta.id} forces the filename to be the string which fills meta.id
    if(meta.single_end){

        """
        cat ${reads} | fastqc stdin:"${meta.id}"
        """

    } else {

        """
        cat ${reads}[0] | fastqc stdin:${meta.id}
        cat ${reads}[1] | fastqc stdin:${meta.id}
        """
        
    }

}

// Manifest file
  // create manifest file to import data into qiime2
process MANIFEST {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path folder
    
  output:
    path ("manifest-file.tsv", emit: manFile)
    path ("sample-id.txt")
    path ("filepath.txt")

  script:
  """
  touch sample-id.txt
  touch filepath.txt
  echo "sample-id" >> sample-id.txt
  echo "absolute-filepath" >> filepath.txt
  for i in "\$(ls ${params.fastq_folder}/*fastq)"; do echo "\$i" | awk -F _280_ '{print \$2}' | awk -F . '{print \$1}' >> sample-id.txt; done
  find "\$(pwd)"/${params.fastq_folder}/*fastq >> filepath.txt
  paste sample-id.txt filepath.txt > manifest-file.tsv
  """
}

// Metadata file
  // create manifest file to import data into qiime2
process METADATA_FILE {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path folder3
    
  output:
    path ("metadata-file.txt"), emit: metadataFile
    path ("SampleID.txt")
    path ("BarcodeSequence.txt")
    path ("LinkerPrimerSequence.txt")
    path ("Path.txt")

  script:
  """
  touch SampleID.txt
  touch BarcodeSequence.txt
  touch LinkerPrimerSequence.txt
  touch Path.txt
  echo "SampleID" >> SampleID.txt
  echo "BarcodeSequence" >> BarcodeSequence.txt
  echo "LinkerPrimerSequence" >> LinkerPrimerSequence.txt
  echo "Path" >> Path.txt
  for i in "\$(ls ${params.fastq_folder}/*fastq)"; do echo "\$i" | awk -F _280_ '{print \$2}' | awk -F . '{print \$1}' >> SampleID.txt; done
  find "\$(pwd)"/${params.fastq_folder}/*fastq >> Path.txt
  paste SampleID.txt BarcodeSequence.txt LinkerPrimerSequence.txt Path.txt > metadata-file.txt
  """
}

// Setup file
  // outputs a txt file with parameters used in the run
process SETUP {
  publishDir "${params.outdir}", mode:'copy'

  output:
  file "setup_file.txt"

  script:
  """
  touch setup_file.txt
  echo "fastq_folder: ${params.fastq_folder}" >> setup_file.txt
  echo "primer_seq: ${params.primer_seq}" >> setup_file.txt
  echo "trainned classifier: ${params.trainned_classifier}" >> setup_file.txt
  echo "reference reads: ${params.ref_reads}" >> setup_file.txt
  echo "tax file: ${params.tax_file}" >> setup_file.txt
  echo "output dir name: ${params.outdir}" >> setup_file.txt
  echo "threads: ${params.threads}" >> setup_file.txt
  echo "p_quality_cutoff_5end: ${params.p_quality_cutoff_5end}" >> setup_file.txt
  echo "p_quality_cutoff_3end: ${params.p_quality_cutoff_3end}" >> setup_file.txt
  echo "p_error_rate: ${params.p_error_rate}" >> setup_file.txt
  #echo "p_minimum_length: ${params.p_minimum_length}" >> setup_file.txt
  """
}

// Import data from manifest file
process IMPORTDATA {
  publishDir "${params.outdir}", mode:'copy'
    
  input:
  file manFile

  output:
  path ("fastq_imported.qza"), emit: fastq_imported

  script:
  """
  qiime tools import --type 'SampleData[SequencesWithQuality]' \
  --input-path ${manFile} \
  --output-path fastq_imported.qza \
  --input-format SingleEndFastqManifestPhred33V2
  """
}

// visualize fastq imported files
process VIEW_IMPORTED {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "fastq_imported.qza"

  output:
  file "view-inspec_import.qzv"

  script:
  """
  qiime demux summarize --i-data fastq_imported.qza \
  --o-visualization view-inspec_import.qzv
  """
}

// trimming
  // trim and filter using cutadapt
params.p_minimum_length = params.p_minimum_length ?: [20, 80]  
process TRIM {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  tuple path(file), val(k)

  output:
  path ("*minLen_trimmed-seqs.qza"), emit: trimmed_seqs
  //path ("*minLen_trimmed-seqs.qza"), emit: trimmed_seqs

  script:
  """
  qiime cutadapt trim-single --i-demultiplexed-sequences $file \
  --p-cores "${params.threads}" \
  --p-front "${params.primer_seq}" \
  --p-quality-cutoff-5end "${params.p_quality_cutoff_5end}" \
  --p-quality-cutoff-3end "${params.p_quality_cutoff_3end}" \
  --p-error-rate "${params.p_error_rate}" \
  --p-minimum-length "${k}" \
  --o-trimmed-sequences "${k}minLen_trimmed-seqs.qza"
  """
}

// inspect trimmed
process VIEW_TRIMMED {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  path(trimmed_seqs)

  output:
  //file "*_view-trimmed_sequences.qzv"
  path ("*_view-trimmed_sequences.qzv"), emit: view_trimmed_seqs

  script:
  def prefix  = files(trimmed_seqs).first().name - '_trimmed-seqs.qza'
  """
  qiime demux summarize \
  --i-data ${prefix}_trimmed-seqs.qza \
  --o-visualization ${prefix}_view-trimmed_sequences.qzv
  """
}

// denoising
  // denoise reads using dada2
process DENOISE {
  publishDir "${params.outdir}", mode:'copy'
  tag "$trimmed_seqs"
  
  input:
  file(trimmed_seqs)

  output:
  path ("*representative-seqs.qza"), emit: rep_seqs
  path ("*table-denoised.qza"), emit: table_denoised
  path ("*denoise-stats.qza"), emit: denoise_stats
  // tuple file("representative-seqs.qza"), file("table-denoised.qza"), file("denoise-stats.qza")


  script:
  def prefix  = files(trimmed_seqs).first().name - '_trimmed-seqs.qza'
  """
  qiime dada2 denoise-single \
  --p-n-threads "${params.threads}" \
  --i-demultiplexed-seqs ${prefix}_trimmed-seqs.qza \
  --p-max-ee "${params.p_max_ee}" \
  --p-trunc-len "${params.p_trunc_len}" \
  --p-pooling-method 'independent' \
  --p-chimera-method 'consensus' \
  --o-representative-sequences ${prefix}_representative-seqs.qza \
  --o-table ${prefix}_table-denoised.qza \
  --o-denoising-stats ${prefix}_denoise-stats.qza
  """
}

// view denoising output
process VIEW_DENOISE {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file(denoise_stats)

  output:
  file ("*view-inspect_denoise-stats.qzv")
  
  script:
  def prefix  = files(denoise_stats).first().name - '_denoise-stats.qza'
  """
  qiime metadata tabulate \
  --m-input-file ${prefix}_denoise-stats.qza \
  --o-visualization ${prefix}_view-inspect_denoise-stats.qzv
  """
}

// Denoise stats to tsv
process DENOISE_TSV {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  path(denoise_stats)

  output:
  file "*denoise-stats_tsv"
  
  script:
  def prefix  = files(denoise_stats).first().name - '_denoise-stats.qza'
  """
  qiime tools export \
  --input-path ${prefix}_denoise-stats.qza \
  --output-path ${prefix}_denoise-stats_tsv
  """
}

// Denoise stats to qzv
process DENOISE_SUMMARY {
  publishDir "${params.outdir}", mode:'copy'
  tag "$table_denoised"
  
  input:
  file(table_denoised)

  output:
  file "*table-denoised_summary.qzv"
  
  script:
  def prefix  = files(table_denoised).first().name - '_table-denoised.qza'
  """
  qiime feature-table summarize \
  --i-table ${prefix}_table-denoised.qza \
  --o-visualization ${prefix}_table-denoised_summary.qzv
  """
}

// Import seq reference file
process DECOMPRESS_REF {
    input:
    file(import_ref_reads_ch)

    output:
    path ("*"), emit: ref_uncompressed
    
    script:
    """
    gzip --decompress --force ${import_ref_reads_ch}
    """
}

process IMPORT_REF {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  path(ref_uncompressed)

  output:
  path ("reference_sequences.qza"), emit: ref_seqs
  
  script:
  """
  qiime tools import \
  --type FeatureData[Sequence] \
  --input-path ${ref_uncompressed} \
  --output-path reference_sequences.qza
  """
}

// Import taxonomy reference file
process DECOMPRESS_TAX {
    input:
    path(import_tax_ch)

    output:
    path ("*"), emit: tax_uncompressed
    
    script:
    """
    gzip --decompress --force ${import_tax_ch}
    """
}

process IMPORT_TAX {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  path(tax_uncompressed)

  output:
  path ("reference_taxonomy.qza"), emit: ref_tax
  
  script:
  """
  qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path ${tax_uncompressed} \
  --output-path reference_taxonomy.qza \
  --input-format HeaderlessTSVTaxonomyFormat
  """
}

// taxonomy classification
process VSEARCH {
  tag "$rep_seqs"
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  path(rep_seqs)
  file "reference_sequences.qza"
  file "reference_taxonomy.qza"

  output:
  path ("*vsearch-taxonomyITS.qza"), emit: vsearch_taxonomyITS
  path ("*view-vsearch-taxonomyITS.qzv"), emit: view_vsearch_taxonomyITS
  path ("*vsearch"), emit: vsearch
  
  script:
  def prefix  = files(rep_seqs).first().name - '_representative-seqs.qza'
  """
  qiime feature-classifier classify-consensus-vsearch \
  --i-query ${prefix}_representative-seqs.qza \
  --i-reference-reads reference_sequences.qza \
  --i-reference-taxonomy reference_taxonomy.qza \
  --p-perc-identity "${params.vsearch_p_perc_identity}" \
  --p-threads "${params.threads}" \
  --o-classification ${prefix}_vsearch-taxonomyITS.qza \
  --output-dir ${prefix}_vsearch

  qiime metadata tabulate \
  --m-input-file ${prefix}_vsearch-taxonomyITS.qza \
  --o-visualization ${prefix}_view-vsearch-taxonomyITS.qzv
  """
}

process BLAST {
  publishDir "${params.outdir}", mode:'copy'
  tag "$rep_seqs"
  
  input:
  file(rep_seqs)
  file "reference_sequences.qza"
  file "reference_taxonomy.qza"

  output:
  path ("*blast-taxonomyITS.qza"), emit: blast_taxonomyITS
  path ("*view-blast-taxonomyITS.qzv"), emit: view_blast_taxonomyITS
  path ("*blast"), emit: blast
  
  script:
  def prefix  = files(rep_seqs).first().name - '_representative-seqs.qza'
  """
  qiime feature-classifier classify-consensus-blast \
  --i-query ${prefix}_representative-seqs.qza \
  --i-reference-reads reference_sequences.qza \
  --i-reference-taxonomy reference_taxonomy.qza \
  --p-perc-identity "${params.blast_p_perc_identity}" \
  --o-classification ${prefix}_blast-taxonomyITS.qza \
  --output-dir ${prefix}_blast

  qiime metadata tabulate \
  --m-input-file ${prefix}_blast-taxonomyITS.qza \
  --o-visualization ${prefix}_view-blast-taxonomyITS.qzv
  """
}

process DECOMPRESS_QZA {
    input:
    path(trainned_ch)

    output:
    path ("*.qza"), emit: qza_uncompressed
    
    script:
    """
    gzip --decompress --force ${trainned_ch}
    """
}

process SKLEARN{
  publishDir "${params.outdir}", mode:'copy'
  tag "$rep_seqs"
  
  input:
  file qza_uncompressed
  file(rep_seqs)
  
  output:
  path "*sklearn-taxonomyITS.qza", emit: sklearn_taxonomyITS
  path "*view-sklearn-taxonomyITS.qzv", emit: view_sklearn_taxonomyITS
  
  script:
  def prefix  = files(rep_seqs).first().name - '_representative-seqs.qza'
  """
  qiime feature-classifier classify-sklearn \
  --p-n-jobs -1 \
  --p-confidence "${params.sklearn_p_confidence}" \
  --i-classifier "${qza_uncompressed}" \
  --p-reads-per-batch "${params.sklearn_p_reads_per_batch}" \
  --i-reads ${prefix}_representative-seqs.qza \
  --o-classification ${prefix}_sklearn-taxonomyITS.qza

  qiime metadata tabulate \
  --m-input-file ${prefix}_sklearn-taxonomyITS.qza \
  --o-visualization ${prefix}_view-sklearn-taxonomyITS.qzv
  """
}

// export biom table
process EXPORT_BIOM1 {
  publishDir "${params.outdir}", mode:'copy'
  tag "$table_denoised"

  input:
    file(table_denoised)

  output:
    path ("*feature-table/feature-table.biom"), emit: feature_table_biom20
    path ("*feature-table.biom", emit: file_feature_table_biom20)

  script:
  def prefix  = files(table_denoised).first().name - '_table-denoised.qza'
  """
  qiime tools export \\
    --input-path ${prefix}_table-denoised.qza \\
    --output-path ${prefix}_feature-table
  
  cat ${prefix}_feature-table/feature-table.biom > ${prefix}_feature-table.biom
  """  
}

process EXPORT_BIOM2 {
  publishDir "${params.outdir}", mode:'copy'
  tag "$y"

  input:
    file(y)

  output:
    path ("*feature-table/feature-table.biom"), emit: feature_table_biom80
    path ("*feature-table.biom", emit: file_feature_table_biom80)

  script:
  def prefix  = files(y).first().name - '_table-denoised.qza'
  """
  qiime tools export \\
    --input-path ${prefix}_table-denoised.qza \\
    --output-path ${prefix}_feature-table
  
  cat ${prefix}_feature-table/feature-table.biom > ${prefix}_feature-table.biom
  """  
}

// export taxonomy
process TAXONOMY {
  publishDir "${params.outdir}", mode:'copy'

  input:    
    file(x)    
    //file(vsearch_taxonomyITS)
    //file(blast_taxonomyITS)

  output:
    path ("*-taxonomy"), emit: taxonomy_tsv
    //path ("*vsearch-taxonomy"), emit: file_vsearch_taxonomy_tsv
    //path ("*blast-taxonomy"), emit: file_blast_taxonomy_tsv

  script:  
  def prefix  = files(x).first().name - '-taxonomyITS.qza'  
  """
  qiime tools export \
  --input-path $x \
  --output-path ${prefix}-taxonomy

  
  """
}

// replace header
process REPLACE_HEADER1 {
  tag "${tax_ch_20}"
  publishDir "${params.outdir}", mode:'copy'
  
  input:
    path (tax_ch_20)
    
  output:
    path ("*_taxonomy.tsv"), emit: header_tsv_20
     
  script:
  def prefix  = files(tax_ch_20).first().name - '-taxonomy' 
  """
  sed 's/Feature ID/#out-id/g' ${prefix}-taxonomy/taxonomy.tsv > ${prefix}_taxonomy.tsv
  sed -i 's/Taxon/taxonomy/g' ${prefix}_taxonomy.tsv

  """  
}

process REPLACE_HEADER2 {
  tag "${tax_ch_80}"
  publishDir "${params.outdir}", mode:'copy'
  
  input:
    path (tax_ch_80)
    
  output:
    path ("*_taxonomy.tsv"), emit: header_tsv_80
     
  script:
  def prefix  = files(tax_ch_80).first().name - '-taxonomy' 
  """
  sed 's/Feature ID/#out-id/g' ${prefix}-taxonomy/taxonomy.tsv > ${prefix}_taxonomy.tsv
  sed -i 's/Taxon/taxonomy/g' ${prefix}_taxonomy.tsv
  """  
}


// add metadata
process METADATA20 {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path (file_feature_table_biom20)
    path (header_tsv_20)
    //path (sklearn_header_tsv)
    //path (vsearch_header_tsv)

  output:
    //path ("20minLen_blast_biom_with_taxonomy.biom"), emit: blast_biom_with_taxonomy20
    //path ("20minLen_sklearn_biom_with_taxonomy.biom"), emit: sklearn_biom_with_taxonomy20
    //path ("20minLen_vsearch_biom_with_taxonomy.biom"), emit: vsearch_biom_with_taxonomy20
    path ("*biom_with_taxonomy.biom"), emit:biom20

  script:
  def prefix  = files(header_tsv_20).first().name - '_taxonomy.tsv'
  """
  biom add-metadata \
  --input-fp 20minLen_feature-table.biom \
  --observation-metadata-fp  ${prefix}_taxonomy.tsv \
  --output-fp ${prefix}_biom_with_taxonomy.biom
  
  """
}

process METADATA80 {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path (file_feature_table_biom80)
    path (header_tsv_80)
    //path (sklearn_header_tsv)
    //path (vsearch_header_tsv)

  output:
    //path ("20minLen_blast_biom_with_taxonomy.biom"), emit: blast_biom_with_taxonomy20
    //path ("20minLen_sklearn_biom_with_taxonomy.biom"), emit: sklearn_biom_with_taxonomy20
    //path ("20minLen_vsearch_biom_with_taxonomy.biom"), emit: vsearch_biom_with_taxonomy20
    path ("*biom_with_taxonomy.biom"), emit:biom80

  script:
  def prefix  = files(header_tsv_80).first().name - '_taxonomy.tsv'
  """
  biom add-metadata \
  --input-fp 80minLen_feature-table.biom \
  --observation-metadata-fp  ${prefix}_taxonomy.tsv \
  --output-fp ${prefix}_biom_with_taxonomy.biom
  """
}


// Convert biom to text
process CONVERT_BIOM20 {
  publishDir "${params.outdir}", mode:'copy'

  input:
    /* path (blast_biom_with_taxonomy)
    path (sklearn_biom_with_taxonomy)
    path (vsearch_biom_with_taxonomy)
    path (vsearch_header_tsv)
    path (blast_header_tsv)
    path (sklearn_header_tsv) */
    path (convert_ch_20)
    path (header_ch_20)
  
  output:
    path ("*_biom_with_taxonomy.tsv"), emit: biom_with_taxonomy_tsv_20
    
  script:
    def prefix  = files(header_ch_20).first().name - '_taxonomy.tsv'
  """
  biom convert \
  --input-fp ${prefix}_biom_with_taxonomy.biom \
  --output-fp ${prefix}_biom_with_taxonomy.tsv \
  --to-tsv \
  --observation-metadata-fp ${prefix}_taxonomy.tsv \
  --header-key taxonomy
  """
}
process CONVERT_BIOM80 {
  publishDir "${params.outdir}", mode:'copy'

  input:
    /* path (blast_biom_with_taxonomy)
    path (sklearn_biom_with_taxonomy)
    path (vsearch_biom_with_taxonomy)
    path (vsearch_header_tsv)
    path (blast_header_tsv)
    path (sklearn_header_tsv) */
    path (convert_ch_80)
    path (header_ch_80)
  
  output:
    path ("*_biom_with_taxonomy.tsv"), emit: biom_with_taxonomy_tsv_80
    
  script:
    def prefix  = files(header_ch_80).first().name - '_taxonomy.tsv'
  """
  biom convert \
  --input-fp ${prefix}_biom_with_taxonomy.biom \
  --output-fp ${prefix}_biom_with_taxonomy.tsv \
  --to-tsv \
  --observation-metadata-fp ${prefix}_taxonomy.tsv \
  --header-key taxonomy
  """
}

// keep all columns but the first one
process ADJUST_COLS {
  publishDir "${params.outdir}", mode:'copy'

  input:
    /* path (vsearch_biom_with_taxonomy_tsv)
    path (blast_biom_with_taxonomy_tsv)
    path (sklearn_biom_with_taxonomy_tsv) */
    path (biom_tsv_files)
  
  output:
    path ("*_final_output.csv"), emit: final_output_csv
    //path ("vsearch_final_output.csv"), emit: vsearch_final_output_csv
    //path ("sklearn_final_output.csv"), emit: sklearn_final_output_csv

  script:
    def prefix  = files(biom_tsv_files).first().name - '_biom_with_taxonomy.tsv'
  """
  cat $biom_tsv_files | cut --complement -f1 > ${prefix}_final_output.csv
  """

}


// format table using a custom python script
process PYTHON_TASK20 {
  publishDir "${params.outdir}/final_report", mode:'copy'

  input:
    path (outdir_ch)
    /* path (sklearn_final_output_csv)
    path (vsearch_final_output_csv)
    path (blast_final_output_csv) */
    path (python20)
    
      
  output:
    //stdout (emit: python_out)
    path ("*final_table.csv"), emit: python_out_20

  script:
  """
  python_its.py ${outdir_ch} 20minLen_final_table.csv 20minLen
  """
}

process PYTHON_TASK80 {
  publishDir "${params.outdir}/final_report", mode:'copy'

  input:
    path (outdir_ch)
    /* path (sklearn_final_output_csv)
    path (vsearch_final_output_csv)
    path (blast_final_output_csv) */
    path (python80)
    
      
  output:
    //stdout (emit: python_out)
    path ("*final_table.csv"), emit: python_out_80

  script:
  """
  python_its.py ${outdir_ch} 80minLen_final_table.csv 80minLen
  """
}

// format table using a custom python script
process REPORT1 {
  //stageInMode 'copy'
  //stageOutMode 'copy'
  //publishDir "${params.outdir}/final_report", mode:'copy'
  
  debug true
  input:
    path (python_out_20)
    path (python_out_80)
    path (outdir_ch2)
          
  output:
    path ("${params.outdir}/final_report/report_html.html"), emit: report
    //stdout
        
  script:
  """
  cp $projectDir/bin/report_html.qmd ${params.outdir}
  quarto render ${params.outdir}/report_html.qmd -P fastqc_dir:"fastqc_dir" -P output:"report_html.html" -P folder_full:
  mv ${params.outdir}/report_html.html ${params.outdir}/final_report
  rm -rf ${params.outdir}/report_html.qmd 
  """
}

process REPORT2 {
  //stageInMode 'copy'
  //stageOutMode 'copy'
  //publishDir "${params.outdir}/final_report", mode:'copy'
  
  debug true
  input:
    path (python_out_20)
    path (python_out_80)
    path (outdir_ch2)
          
  output:
    path ("${params.outdir}/final_report/report_pdf.pdf"), emit: report
    //stdout
        
  script:
  """
  cp $projectDir/bin/report_pdf.qmd ${params.outdir}
  quarto render ${params.outdir}/report_pdf.qmd -P fastqc_dir:"fastqc_dir" -P output:"report_pdf.pdf" -P folder_full:
  mv ${params.outdir}/report_pdf.pdf ${params.outdir}/final_report
  rm -rf ${params.outdir}/report_pdf.qmd 
  """
}
// // WORKFLOW // // 
workflow {
    // Define channels for paired-end and single-end files
    def paired_channel = Channel
        .fromFilePairs("${params.fastq_folder}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}")
        .map { id, reads ->
            def sample = id.tokenize("_R")[0]
            def meta = [sample: sample, single_end: false]
            [meta, [r1, r2]]
        }

    def single_channel = Channel
        .fromPath("${params.fastq_folder}/*.{fastq,fastq.gz,fq,fq.gz}")
        .filter { path -> 
            !path.name.contains("_R1") && !path.name.contains("_R2")
        }
        .map { path ->
            def sample = path.name.split("_280_")[1]
            def meta = [id: sample, single_end: true]
            [meta, [path]]
        }

    // Merge paired-end and single-end channels
        // If use concat oeprator, itens will be emitted in the specified order
    paired_channel
        .mix(single_channel)
        .set { my_channel }

    channel.fromPath(params.projectDir)
      .set { db_ch  }
    
    channel.fromPath(params.fastq_folder)
      .set { folder }

    channel.fromPath(params.fastq_folder)
      .set { folder2 }

    channel.fromPath(params.fastq_folder)
      .set { folder3 }

    channel.fromPath(params.ref_reads)
      .set { import_ref_reads_ch }

    channel.fromPath(params.tax_file)
      .set { import_tax_ch }

    channel.fromPath(params.trainned_classifier)
      .set { trainned_ch }

    channel.fromPath(params.outdir)
      .set { outdir_ch }

    channel.fromPath(params.outdir)
      .set { outdir_ch2 }

    channel.fromPath(params.projectDir)
      .set { db_ch }

// Create configuration file
    SAMPLESHEET(folder)
// Run fastqc
    FASTQC(my_channel)
    MANIFEST(folder2)
    METADATA_FILE(folder3)
    SETUP()
// Import fastq files using read_ch
    IMPORTDATA(MANIFEST.out.manFile.collectFile())
// visualize imported fastq files 
    VIEW_IMPORTED(IMPORTDATA.out.collectFile(name: 'fastq_imported.qza'))
// trimming
    //channel to import qza file
    imported_qza = IMPORTDATA.out.fastq_imported.collect()
    // channel to import minimum lengths
    length = channel.of(params.p_minimum_length)
    // combine channels
    imported_qza
      .flatMap {file -> params.p_minimum_length.collect {k -> 
          [file, k] 
        }
         }
      .set {combine}

    TRIM(combine) 
// view trimmed
    VIEW_TRIMMED(TRIM.out.trimmed_seqs.collect().flatten())
// denoise
    DENOISE(TRIM.out.trimmed_seqs.collect().flatten())
    VIEW_DENOISE(DENOISE.out.denoise_stats.collect().flatten())
    DENOISE_TSV(DENOISE.out.denoise_stats.collect().flatten())
    DENOISE_SUMMARY(DENOISE.out.table_denoised.collect().flatten())
// Import
    DECOMPRESS_REF(import_ref_reads_ch)
    DECOMPRESS_TAX(import_tax_ch)
    IMPORT_REF(DECOMPRESS_REF.out.ref_uncompressed.collect())
    IMPORT_TAX(DECOMPRESS_TAX.out.tax_uncompressed.collect())
// Taxonomy classification
    DECOMPRESS_QZA(trainned_ch)
    //DENOISE.out.rep_seqs.collect().flatten().view()
    denoise_ch = DENOISE.out.rep_seqs.collect()
    //denoise_ch.view()

    VSEARCH(DENOISE.out.rep_seqs.collect().flatten(), IMPORT_REF.out.ref_seqs.collect(), IMPORT_TAX.out.ref_tax.collect())
    BLAST(DENOISE.out.rep_seqs.collect().flatten(), IMPORT_REF.out.ref_seqs.collect(), IMPORT_TAX.out.ref_tax.collect())
    SKLEARN(DECOMPRESS_QZA.out.qza_uncompressed.collect(), DENOISE.out.rep_seqs.collect().flatten())
// Export
    /* EXPORT_BIOM1(DENOISE.out
      .table_denoised
      .collectFile(name: '20minLen_table-denoised.qza'))
    
    EXPORT_BIOM2(DENOISE.out
      .table_denoised
      .collectFile(name: '80minLen_table-denoised.qza')) */
    denoised_20 = DENOISE.out.table_denoised
      .collect()
      .flatten()  // Achata as tuplas em um único nível
      .filter { filePath -> filePath.toString().contains('20minLen') }

    denoised_80 = DENOISE.out.table_denoised
      .collect()
      .flatten()  // Achata as tuplas em um único nível
      .filter { filePath -> filePath.toString().contains('80minLen') }
    EXPORT_BIOM1(denoised_20)
    
    EXPORT_BIOM2(denoised_80)
    
    sklearn_out_ch = SKLEARN.out.sklearn_taxonomyITS.collect().flatten()
    vsearch_out_ch = VSEARCH.out.vsearch_taxonomyITS.collect().flatten()
    blast_out_ch = BLAST.out.blast_taxonomyITS.collect().flatten()

    sklearn_out_ch
      .mix(vsearch_out_ch,blast_out_ch)
      .set { taxonomy_input }
    
    //sklearn_out_ch.view()
    
    //TAXONOMY(SKLEARN.out.sklearn_taxonomyITS.collect().flatten(), VSEARCH.out.vsearch_taxonomyITS.collect().flatten(), BLAST.out.blast_taxonomyITS.collect().flatten())
    TAXONOMY(taxonomy_input)
    //##################################################################

    tax_ch_20 = TAXONOMY.out.taxonomy_tsv
      .collect()
      .flatten()  // Achata as tuplas em um único nível
      .filter { filePath -> filePath.toString().contains('20minLen') }  // Filtra os caminhos que contêm "20minLen"
      
    
    tax_ch_80 = TAXONOMY.out.taxonomy_tsv
      .collect()
      .flatten()  // Achata as tuplas em um único nível
      .filter { filePath -> filePath.toString().contains('80minLen') }  // Filtra os caminhos que contêm "20minLen"
        
    //tax_ch_20.view()
    REPLACE_HEADER1(tax_ch_20)
    REPLACE_HEADER2(tax_ch_80)
    
// Add metadata

    METADATA20(EXPORT_BIOM1.out.file_feature_table_biom20.collect(),
      REPLACE_HEADER1.out.header_tsv_20.collect().flatten())
    METADATA80(EXPORT_BIOM2.out.file_feature_table_biom80.collect(),
      REPLACE_HEADER2.out.header_tsv_80.collect().flatten())

// Convert biom
    convert_ch_20 = METADATA20.out.biom20
      .collect()
    header_ch_20 = REPLACE_HEADER1.out.header_tsv_20
      .collect()
      .flatten()

    CONVERT_BIOM20(convert_ch_20, 
      header_ch_20)

    convert_ch_80 = METADATA80.out.biom80
      .collect()
    header_ch_80 = REPLACE_HEADER2.out.header_tsv_80
      .collect()
      .flatten()
    
    CONVERT_BIOM80(convert_ch_80, 
      header_ch_80)
// Adjust cols
    CONVERT_BIOM20.out.biom_with_taxonomy_tsv_20
        .mix(CONVERT_BIOM80.out.biom_with_taxonomy_tsv_80)
        .set { biom_tsv_files }

    ADJUST_COLS(biom_tsv_files)
// adjust sheets
// Python script
    python80 = ADJUST_COLS.out.final_output_csv
      .collect()
      .flatten()  // Achata as tuplas em um único nível
      .filter { filePath -> filePath.toString().contains('80minLen') }
    python20 = ADJUST_COLS.out.final_output_csv
      .collect()
      .flatten()  // Achata as tuplas em um único nível
      .filter { filePath -> filePath.toString().contains('20minLen') } 
    //python.view()
    PYTHON_TASK80(outdir_ch, python80)
    PYTHON_TASK20(outdir_ch, python20)    
// render report
    REPORT1(PYTHON_TASK20.out.python_out_20.collect(),
      outdir_ch2, PYTHON_TASK80.out.python_out_80.collect())
    
    REPORT2(PYTHON_TASK20.out.python_out_20.collect(),
      outdir_ch2, PYTHON_TASK80.out.python_out_80.collect())
      
      }
