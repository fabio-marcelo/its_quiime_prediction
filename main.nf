 #! /usr/bin/env nextflow
 nextflow.enable.dsl=2

 // Define help message
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:

        nextflow run main.nf --primer_seq "primer_sequence" --fastq_folder /path/to/my/fastq_files 
        --trainned_classifier "full_path/and_name/to/my/trainned_classifier 
        --ref_reads full_path/to/reference_reads.fasta --tax_file full_path/to/tax_file.txt 
        --outdir output-folder-name --threads "15"

  
  
  

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
    publishDir "${params.outdir}", mode:'copy'

    input:
        path(folder) 

    output:
        path ("samplesheet.csv")
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
    for i in "\$(ls ${params.fastq_folder}/*fastq)"; do echo "\$i" | awk -F _bp_ '{print \$2}' | awk -F . '{print \$1}' >> sample.txt; done
    find "\$(pwd)"/${params.fastq_folder}/*fastq >> r1.txt
    paste -d , sample.txt r1.txt r2.txt > samplesheet.csv
    """
    
}

// Channel for the samplesheet
ch_samplesheet = Channel.fromPath("${params.fastq_folder}/samplesheet.csv")

// Parse it line by line
ch_reads = ch_samplesheet.splitCsv(header:true).map {

    // This is the read1 and read2 entry
    r1 = it['r1']
    r2 = it['r2']

    // Detect wiether single-end or paired-end
    is_singleEnd = r2.toString()=='' ? true : false
    
    // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
    meta = [id: it['sample'], single_end: is_singleEnd]
    
    // We return a nested map, the first entry is the meta map, the second one is the read(s)
    r2.toString()=='' ? [meta, [r1]] : [meta, [r1, r2]]

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
        cat ${reads} | fastqc stdin:${meta.id}
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
  for i in "\$(ls ${params.fastq_folder}/*fastq)"; do echo "\$i" | awk -F _bp_ '{print \$2}' | awk -F . '{print \$1}' >> sample-id.txt; done
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
  for i in "\$(ls ${params.fastq_folder}/*fastq)"; do echo "\$i" | awk -F _bp_ '{print \$2}' | awk -F . '{print \$1}' >> SampleID.txt; done
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
  echo "p_minimum_length: ${params.p_minimum_length}" >> setup_file.txt
  """
}

// Import data from manifest file
process IMPORTDATA {
  publishDir "${params.outdir}", mode:'copy'
    
  input:
  file manFile

  output:
  file "fastq_imported.qza"

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
process TRIM {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "fastq_imported.qza"

  output:
  file "trimmed-seqs.qza"

  script:
  """
  qiime cutadapt trim-single --i-demultiplexed-sequences fastq_imported.qza \
  --p-cores "${params.threads}" \
  --p-front "${params.primer_seq}" \
  --p-quality-cutoff-5end "${params.p_quality_cutoff_5end}" \
  --p-quality-cutoff-3end "${params.p_quality_cutoff_3end}" \
  --p-error-rate "${params.p_error_rate}" \
  --p-minimum-length "${params.p_minimum_length}" \
  --o-trimmed-sequences trimmed-seqs.qza
  """
}

// inspect trimmed
process VIEW_TRIMMED {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "trimmed-seqs.qza"

  output:
  file "view-trimmed_sequences.qzv"

  script:
  """
  qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization view-trimmed_sequences.qzv
  """
}

// denoising
  // denoise reads using dada2
process DENOISE {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "trimmed-seqs.qza"

  output:
  path ("representative-seqs.qza"), emit: rep_seqs
  path ("table-denoised.qza"), emit: table_denoised
  path ("denoise-stats.qza"), emit: denoise_stats
  // tuple file("representative-seqs.qza"), file("table-denoised.qza"), file("denoise-stats.qza")


  script:
  """
  qiime dada2 denoise-single \
  --p-n-threads "${params.threads}" \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-max-ee "${params.p_max_ee}" \
  --p-trunc-len "${params.p_trunc_len}" \
  --p-pooling-method 'independent' \
  --p-chimera-method 'consensus' \
  --o-representative-sequences representative-seqs.qza \
  --o-table table-denoised.qza \
  --o-denoising-stats denoise-stats.qza
  """
}

// view denoising output
process VIEW_DENOISE {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "denoise-stats.qza"

  output:
  file "view-inspect_denoise-stats.qzv"
  
  script:
  """
  qiime metadata tabulate \
  --m-input-file denoise-stats.qza \
  --o-visualization view-inspect_denoise-stats.qzv
  """
}

// Denoise stats to tsv
process DENOISE_TSV {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  path "denoise-stats.qza"

  output:
  file "denoise-stats_tsv"
  
  script:
  """
  qiime tools export \
  --input-path denoise-stats.qza \
  --output-path denoise-stats_tsv
  """
}

// Denoise stats to qzv
process DENOISE_SUMMARY {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "table-denoised.qza"

  output:
  file "table-denoised_summary.qzv"
  
  script:
  """
  qiime feature-table summarize \
  --i-table table-denoised.qza \
  --o-visualization table-denoised_summary.qzv
  """
}

// Import seq reference file
process DECOMPRESS_REF {
    input:
    file(import_ref_reads_ch)

    output:
    path ("*.fasta"), emit: ref_uncompressed
    
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
    path ("*.txt"), emit: tax_uncompressed
    
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
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "representative-seqs.qza"
  file "reference_sequences.qza"
  file "reference_taxonomy.qza"

  output:
  path ("vsearch-taxonomyITS.qza"), emit: vsearch_taxonomyITS
  path ("view-vsearch-taxonomyITS.qzv"), emit: view_vsearch_taxonomyITS
  path ("vsearch"), emit: vsearch
  
  script:
  """
  qiime feature-classifier classify-consensus-vsearch \
  --i-query representative-seqs.qza \
  --i-reference-reads reference_sequences.qza \
  --i-reference-taxonomy reference_taxonomy.qza \
  --p-perc-identity "${params.vsearch_p_perc_identity}" \
  --p-threads "${params.threads}" \
  --o-classification vsearch-taxonomyITS.qza \
  --output-dir vsearch

  qiime metadata tabulate \
  --m-input-file vsearch-taxonomyITS.qza \
  --o-visualization view-vsearch-taxonomyITS.qzv
  """
}

process BLAST {
  publishDir "${params.outdir}", mode:'copy'
  
  input:
  file "representative-seqs.qza"
  file "reference_sequences.qza"
  file "reference_taxonomy.qza"

  output:
  path ("blast-taxonomyITS.qza"), emit: blast_taxonomyITS
  path ("view-blast-taxonomyITS.qzv"), emit: view_blast_taxonomyITS
  
  script:
  """
  qiime feature-classifier classify-consensus-blast \
  --i-query representative-seqs.qza \
  --i-reference-reads reference_sequences.qza \
  --i-reference-taxonomy reference_taxonomy.qza \
  --p-perc-identity "${params.blast_p_perc_identity}" \
  --o-classification blast-taxonomyITS.qza \
  --o-search-results blast

  qiime metadata tabulate \
  --m-input-file blast-taxonomyITS.qza \
  --o-visualization view-blast-taxonomyITS.qzv
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
  
  input:
  file qza_uncompressed
  file "representative-seqs.qza"
  
  output:
  path "sklearn-taxonomyITS.qza", emit: sklearn_taxonomyITS
  path "view-sklearn-taxonomyITS.qzv", emit: view_sklearn_taxonomyITS
  
  script:
  """
  qiime feature-classifier classify-sklearn \
  --p-n-jobs -1 \
  --p-confidence "${params.sklearn_p_confidence}" \
  --i-classifier "${qza_uncompressed}" \
  --p-reads-per-batch "${params.sklearn_p_reads_per_batch}" \
  --i-reads representative-seqs.qza \
  --o-classification sklearn-taxonomyITS.qza

  qiime metadata tabulate \
  --m-input-file sklearn-taxonomyITS.qza \
  --o-visualization view-sklearn-taxonomyITS.qzv
  """
}

// export biom table
process EXPORT_BIOM {
  publishDir "${params.outdir}", mode:'copy'

  input:
    file "table-denoised.qza"

  output:
    path ("feature-table/feature-table.biom"), emit: feature_table_biom
    path ("feature-table.biom", emit: file_feature_table_biom)

  script:
  """
  qiime tools export \\
    --input-path table-denoised.qza \\
    --output-path feature-table
  
  cat feature-table/feature-table.biom > feature-table.biom
  """  
}

// export taxonomy
process TAXONOMY {
  publishDir "${params.outdir}", mode:'copy'

  input:    
    path (sklearn_taxonomyITS)    
    path (vsearch_taxonomyITS)
    path (blast_taxonomyITS)

  output:
    path ("sklearn-taxonomy"), emit: file_sklearn_taxonomy_tsv
    path ("vsearch-taxonomy"), emit: file_vsearch_taxonomy_tsv
    path ("blast-taxonomy"), emit: file_blast_taxonomy_tsv

  script:    
  """
  qiime tools export \
  --input-path sklearn-taxonomyITS.qza \
  --output-path sklearn-taxonomy

  qiime tools export \
  --input-path vsearch-taxonomyITS.qza \
  --output-path vsearch-taxonomy

  qiime tools export \
  --input-path blast-taxonomyITS.qza \
  --output-path blast-taxonomy
  """
}
// replace header
process REPLACE_HEADER1 {
  tag "${file_sklearn_taxonomy_tsv}"
  publishDir "${params.outdir}", mode:'copy'
  
  input:
    path (file_sklearn_taxonomy_tsv)
    
  output:
    path ("sklearn_taxonomy.tsv"), emit: sklearn_header_tsv
    
  script:
  """
  sed 's/Feature ID/#out-id/g' sklearn-taxonomy/taxonomy.tsv > sklearn_taxonomy.tsv
  sed -i 's/Taxon/taxonomy/g' sklearn_taxonomy.tsv
  """
}

process REPLACE_HEADER2 {
  tag "${file_vsearch_taxonomy_tsv}"
  publishDir "${params.outdir}", mode:'copy'
  
  input:
    path (file_vsearch_taxonomy_tsv)
  
  output:
    path ("vsearch_taxonomy.tsv"), emit: vsearch_header_tsv
  
  script:
  """
  sed 's/Feature ID/#out-id/g' vsearch-taxonomy/taxonomy.tsv > vsearch_taxonomy.tsv
  sed -i 's/Taxon/taxonomy/g' vsearch_taxonomy.tsv
  """
}

process REPLACE_HEADER3 {
  tag "${file_blast_taxonomy_tsv}"
  publishDir "${params.outdir}", mode:'copy'
  
  input:
    path (file_blast_taxonomy_tsv)
  
  output:
    path ("blast_taxonomy.tsv", emit: blast_header_tsv)
  
  script:
  """
  sed 's/Feature ID/#out-id/g' blast-taxonomy/taxonomy.tsv > blast_taxonomy.tsv
  sed -i 's/Taxon/taxonomy/g' blast_taxonomy.tsv
  """
}

// add metadata
process METADATA {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path (file_feature_table_biom)
    path (blast_header_tsv)
    path (sklearn_header_tsv)
    path (vsearch_header_tsv)

  output:
    path ("blast_biom_with_taxonomy.biom"), emit: blast_biom_with_taxonomy
    path ("sklearn_biom_with_taxonomy.biom"), emit: sklearn_biom_with_taxonomy
    path ("vsearch_biom_with_taxonomy.biom"), emit: vsearch_biom_with_taxonomy

  script:
  """
  biom add-metadata \
  --input-fp feature-table.biom \
  --observation-metadata-fp blast_taxonomy.tsv \
  --output-fp blast_biom_with_taxonomy.biom

  biom add-metadata \
  --input-fp feature-table.biom \
  --observation-metadata-fp sklearn_taxonomy.tsv \
  --output-fp sklearn_biom_with_taxonomy.biom

  biom add-metadata \
  --input-fp feature-table.biom \
  --observation-metadata-fp vsearch_taxonomy.tsv \
  --output-fp vsearch_biom_with_taxonomy.biom 
  """
}

// Convert biom to text
process CONVERT_BIOM {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path (blast_biom_with_taxonomy)
    path (sklearn_biom_with_taxonomy)
    path (vsearch_biom_with_taxonomy)
    path (vsearch_header_tsv)
    path (blast_header_tsv)
    path (sklearn_header_tsv)
  
  output:
    path ("vsearch_biom_with_taxonomy.tsv"), emit: vsearch_biom_with_taxonomy_tsv
    path ("sklearn_biom_with_taxonomy.tsv"), emit: sklearn_biom_with_taxonomy_tsv
    path ("blast_biom_with_taxonomy.tsv"), emit: blast_biom_with_taxonomy_tsv

  script:
  """
  biom convert \
  --input-fp vsearch_biom_with_taxonomy.biom \
  --output-fp vsearch_biom_with_taxonomy.tsv \
  --to-tsv \
  --observation-metadata-fp vsearch_taxonomy.tsv \
  --header-key taxonomy

  biom convert \
  --input-fp sklearn_biom_with_taxonomy.biom \
  --output-fp sklearn_biom_with_taxonomy.tsv \
  --to-tsv \
  --observation-metadata-fp sklearn_taxonomy.tsv \
  --header-key taxonomy

  biom convert \
  --input-fp blast_biom_with_taxonomy.biom \
  --output-fp blast_biom_with_taxonomy.tsv \
  --to-tsv \
  --observation-metadata-fp blast_taxonomy.tsv \
  --header-key taxonomy
  """
}

// keep all columns but the first one
process ADJUST_COLS {
  publishDir "${params.outdir}", mode:'copy'

  input:
    path (vsearch_biom_with_taxonomy_tsv)
    path (blast_biom_with_taxonomy_tsv)
    path (sklearn_biom_with_taxonomy_tsv)
  
  output:
    path ("blast_final_output.csv"), emit: blast_final_output_csv
    path ("vsearch_final_output.csv"), emit: vsearch_final_output_csv
    path ("sklearn_final_output.csv"), emit: sklearn_final_output_csv

  script:
  """
  cat $blast_biom_with_taxonomy_tsv | cut --complement -f1 > blast_final_output.csv
  cat $vsearch_biom_with_taxonomy_tsv | cut --complement -f1 > vsearch_final_output.csv
  cat $sklearn_biom_with_taxonomy_tsv | cut --complement -f1 > sklearn_final_output.csv
  """

}

// format table using a custom python script
process PYTHON_TASK {
  publishDir "${params.outdir}/final_report", mode:'copy'

  input:
    path (outdir_ch)
    path (sklearn_final_output_csv)
    path (vsearch_final_output_csv)
    path (blast_final_output_csv)
    
      
  output:
    //stdout (emit: python_out)
    path ("final_table.csv"), emit: python_out

  script:
  """
  python_its.py ${outdir_ch} final_table.csv
  """
}

// format table using a custom python script
process REPORT {
  //stageInMode 'copy'
  //stageOutMode 'copy'
  //publishDir "${params.outdir}/final_report", mode:'copy'
  
  debug true
  input:
    path (python_out)
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
// // WORKFLOW // // 
workflow {
    db_ch = channel.fromPath(params.projectDir)
    folder = channel.fromPath(params.fastq_folder)
    folder2 = channel.fromPath(params.fastq_folder)
    folder3 = channel.fromPath(params.fastq_folder)
    import_ref_reads_ch = channel.fromPath(params.ref_reads)        
    import_tax_ch = channel.fromPath(params.tax_file)
    trainned_ch = channel.fromPath(params.trainned_classifier)
    outdir_ch = channel.fromPath(params.outdir)
    outdir_ch2 = channel.fromPath(params.outdir)
    db_ch = channel.fromPath(params.projectDir)
    //PREPARE_DB(db_ch)
// Create configuration file
    SAMPLESHEET(folder)
// Run fastqc
    FASTQC(ch_reads)
    MANIFEST(folder2)
    METADATA_FILE(folder3)
    SETUP()
// Import fastq files using read_ch
    IMPORTDATA(MANIFEST.out.manFile.collectFile())
// visualize imported fastq files 
    VIEW_IMPORTED(IMPORTDATA.out.collectFile(name: 'fastq_imported.qza'))
// trimming
    TRIM(IMPORTDATA.out.collectFile(name: 'fastq_imported.qza')) 
// view trimmed
    VIEW_TRIMMED(TRIM.out.collectFile(name: 'trimmed-seqs.qza'))
// denoise
    DENOISE(TRIM.out.collectFile(name: 'trimmed-seqs.qza'))
    VIEW_DENOISE(DENOISE.out.denoise_stats.collect())
    DENOISE_TSV(DENOISE.out.denoise_stats.collect())
    DENOISE_SUMMARY(DENOISE.out.table_denoised.collect())
// Import
    DECOMPRESS_REF(import_ref_reads_ch)
    DECOMPRESS_TAX(import_tax_ch)
    IMPORT_REF(DECOMPRESS_REF.out.ref_uncompressed.collect())
    IMPORT_TAX(DECOMPRESS_TAX.out.tax_uncompressed.collect())
// Taxonomy classification
    DECOMPRESS_QZA(trainned_ch)
    VSEARCH(DENOISE.out.rep_seqs.collect(), IMPORT_REF.out.ref_seqs.collect(), IMPORT_TAX.out.ref_tax.collect())
    BLAST(DENOISE.out.rep_seqs.collect(), IMPORT_REF.out.ref_seqs.collect(), IMPORT_TAX.out.ref_tax.collect())
    SKLEARN(DECOMPRESS_QZA.out.qza_uncompressed.collect(), DENOISE.out.rep_seqs.collect())
// Export
    EXPORT_BIOM(DENOISE.out.table_denoised.collect())
    TAXONOMY(SKLEARN.out.sklearn_taxonomyITS.collect(), VSEARCH.out.vsearch_taxonomyITS.collect(), BLAST.out.blast_taxonomyITS.collect())
    REPLACE_HEADER1(TAXONOMY.out.file_sklearn_taxonomy_tsv.collect())
    REPLACE_HEADER2(TAXONOMY.out.file_vsearch_taxonomy_tsv.collect())
    REPLACE_HEADER3(TAXONOMY.out.file_blast_taxonomy_tsv.collect())
// Add metadata
    METADATA(REPLACE_HEADER3.out.blast_header_tsv.collect(), EXPORT_BIOM.out.feature_table_biom.collect(), REPLACE_HEADER1.out.sklearn_header_tsv.collect(), REPLACE_HEADER2.out.vsearch_header_tsv.collect())
// Convert biom
    CONVERT_BIOM(METADATA.out.vsearch_biom_with_taxonomy.collect(), METADATA.out.sklearn_biom_with_taxonomy.collect(), METADATA.out.blast_biom_with_taxonomy.collect(), REPLACE_HEADER3.out.blast_header_tsv.collect(),REPLACE_HEADER2.out.vsearch_header_tsv.collect(), REPLACE_HEADER1.out.sklearn_header_tsv.collect())
// Adjust cols
    ADJUST_COLS(CONVERT_BIOM.out.vsearch_biom_with_taxonomy_tsv.collect(), CONVERT_BIOM.out.blast_biom_with_taxonomy_tsv.collect(), CONVERT_BIOM.out.sklearn_biom_with_taxonomy_tsv.collect())
// adjust sheets
    PYTHON_TASK(outdir_ch, ADJUST_COLS.out.sklearn_final_output_csv.collect(), ADJUST_COLS.out.blast_final_output_csv.collect(), ADJUST_COLS.out.vsearch_final_output_csv.collect())    
// render report
    REPORT(PYTHON_TASK.out.collect(), outdir_ch2)
}