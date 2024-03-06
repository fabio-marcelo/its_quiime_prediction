# FUNGUS_ITS_TAXONOMY
This is a pipeline implemented in nextflow to identify fungus through ITS2 metabarcoding generated in S5 IonTorrent. The outputs (one csv and one html report) are stored in the `final_report` folder inside the output directory.



# Prerequisites
* OS Ubuntu
* Docker [tutorial](https://docs.docker.com/engine/install/ubuntu/)
* Nextflow [tutorial](https://www.nextflow.io/docs/latest/getstarted.html)
* Nextflow deals with image download and container run by itself.

# Database
The database used for fungi taxonomy and included in this repository was download from [UNITE QIIME release_29.11.2022 for Fungi](https://dx.doi.org/10.15156/BIO/2483915). It includes the trainned classifier for `quay.io/qiime2/core:2023.7`.
Files are stored in the `arquivos_db` folder in gz format. The pipeline itself decompress the files.

# Pipeline
The pipeline is executed through `quay.io/qiime2/core:2023.7` image and includes a python script for data preprocessing and a R script to generate the html report through `quarto`.

## Parameters
* Due to variability in the length of the sequences for ITS,we opted for not using the parameter `--p-trunc-len` = 0 [tutorial](https://benjjneb.github.io/dada2/ITS_workflow.html);
* As the pattern shown by Ion S5 fastq files quality the default for parameters `p_quality_cutoff_5end` and `p_quality_cutoff_3end` is 20;
* `p_error_rate = 0.2`
* `p_minimum_length = 80`
* `p_max_ee` = 2
* Percent identity for `qiime feature-classifier classify-consensus-vsearch = 0.99`
* Percent identity for `qiime feature-classifier classify-consensus-blast = 0.99`
* Percent identity for `qiime feature-classifier classify-sklearn = 0.99`
* Percent identity for `qiime feature-classifier classify-sklearn --p-reads-per-batch = 10000`

# Usage
## Help message
```bash
nextflow run github/its_pipeline/main.nf --help

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
```
## Run the pipeline
You need to provide at least the `--fastq_folder` parameter.

`nextflow run github/its_pipeline/main.nf --fastq_folder 'its_reads'`
