#!/bin/bash
## Swift Biosciences 16S Qiime 2 workflow
## Author Benli Chai & Sukhinder Sandhu 20190618
## Modified to run on Ixodes, AEM updated 20191023
## Make executable by running: chmod +x qiime2_workflow.sh
## run: ./qiime2_workflow.sh

source /home/lymelab/miniconda2/etc/profile.d/conda.sh
conda --version
conda activate qiime2-2019.7

####SET UP PATHS####
WD=/home/lymelab/lab_members/mann/swiftTest/ #working directory
RAW=/home/lymelab/lab_members/mann/swiftTest/raw #raw fastq files directory
CUTADAPT=/usr/bin/cutadapt #path to cutadapt
PRIMERS=/home/lymelab/lab_members/mann/swiftTest/primers_16S_V1-9_anchored.fasta #primer file formatted paired-end primer trimming
READLEN=130 #minimum read length after primer trimming (before PE merging; use 280 PE300)
CLASSIFIER_seq=/home/lymelab/lab_members/mann/swiftTest/silva_132_99_16S.qza #path to qiime formatted reference database (sequences)
CLASSIFIER_tax=/home/lymelab/lab_members/mann/swiftTest/consensus_taxonomy_7_levels.qza #path to qiime formatted reference database taxonomy
MAP=/home/lymelab/lab_members/mann/swiftTest/map.txt #metadata file, should follow qiime2 formatting parameters
QIIME=/home/lymelab/miniconda2/envs/qiime2-2019.7/bin/qiime #path to your qiime install (I'm not sure this is necessary but a fail safe if conda environment doens't activate)

####MAKE OUTPUT DIRECTORIES####
rm -r ${WD}/Qobj
mkdir ${WD}/Qobj #folder for all qiime2 data objects
rm -r ${WD}/Fastq
mkdir ${WD}/Fastq #folder for intermediate fastq files
rm -r ${WD}/Export
mkdir ${WD}/Export #folder to export qiime2 data object

####SET UP LOG FILE####
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
runlog='log'
log=$runlog.$current_time
echo "Starting time: "$current_time >> $log
echo "Primer file used: "$PRIMERS >> $log

####PRIMER TRIM####
for R1 in ${RAW}/*_R1_*fastq.gz; do #the input directory
    R2=${R1/_R1_/_R2_} #the path to the read 2 file
    echo $R1 >> $log
    echo $R2 >> $log
    echo

    basenameR1=${R1##*/}
    basenameR2=${R2##*/}
    prefixR1=${basenameR1%".fastq.gz"}
    prefixR2=${basenameR2%".fastq.gz"}
    trimmedR1=${WD}/Fastq/${prefixR1}"_primerTrimmed.fastq"
    trimmedR2=${WD}/Fastq/${prefixR2}"_primerTrimmed.fastq"
    untrimmedR1=${WD}/Fastq/${prefixR1}"_primerNotFound.fastq" #won't go to downstream analysis
    untrimmedR2=${WD}/Fastq/${prefixR2}"_primerNotFound.fastq" #won't go to downstream analysis
    prefix=${prefixR1%_R1*} #the sample name prefix of R1 and R2
    echo Processing ${prefix} ...

    #PE option to trim primers off the reads
    echo
    echo Trimming primers ...
    $CUTADAPT -e 0.10 -g file:$PRIMERS -G file:$PRIMERS \
              -o $trimmedR1 -p $trimmedR2  \
              --untrimmed-output $untrimmedR1 \
              --untrimmed-paired-output $untrimmedR2 \
              $R1 $R2 \
              --max-n 0 \
              --minimum-length $READLEN \
              >& log.${prefix}.cutadapt.txt

done

mkdir ${WD}/Fastq/noprimers
mv ${WD}/Fastq/*NotFound.fastq ${WD}/Fastq/noprimers

####MAKE MANIFEST FILE####
rm manifest.tsv
manifest=${WD}/manifest.tsv
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" >> $manifest
printf "${prefix}\t${trimmedR1}\t${trimmedR2}\n" >> $manifest

### Starting the ASV approach
     ## Use dada2 for removal of noise including chimeras and generation of ASVs 
     $QIIME tools import \
             --type 'SampleData[PairedEndSequencesWithQuality]' \
             --input-path ${WD}/manifest.tsv \
             --input-format PairedEndFastqManifestPhred33V2 \
             --output-path ${WD}/Qobj/trimmed_sequences.qza \
     
     echo DADA2 processing ...
     $QIIME dada2 denoise-paired \
             --i-demultiplexed-seqs ${WD}/Qobj/trimmed_sequences.qza \
             --o-table ${WD}/Qobj/dada2_ASVs \
             --p-chimera-method 'pooled' \
             # trim N bases from 5' end of forward reads (shouldn't be an issue after primer trimming)
             --p-trim-left-f 0 \
             # trim N bases from the 5' end of reverse reads (shouldn't be an issue after primer trimming)
             --p-trim-left-r 0 \
             # truncate forward reads to N length
             --p-trunc-len-f 0 \
             # truncate reverse reads to N length
             --p-trunc-len-r 0 \
             --o-representative-sequences ${WD}/Qobj/dada2_repset \
             --o-denoising-stats ${WD}/Qobj/dada2_stats \
             --verbose \
    
     echo Assigning taxonomy with blast and ${CLASSIFIER_seq}....
     ## Classification of ASVs from dada2 using blast
     $QIIME feature-classifier classify-consensus-blast \
             --i-reference-reads $CLASSIFIER_seq \
             --i-reference-taxonomy $CLASSIFIER_tax \
             --i-query ${WD}/Qobj/dada2_repset.qza \
             --p-strand both \
             --o-classification ${WD}/Qobj/dada2_rep_taxonomy

     ## Collapse ASV-based feature table by taxonomy
     $QIIME taxa collapse \
             --i-table ${WD}/Qobj/dada2_ASVs.qza \
             --i-taxonomy ${WD}/Qobj/dada2_rep_taxonomy.qza \
             --o-collapsed-table ${WD}/Qobj/dada2_ASVs_collapsed.qza \
             --p-level 6
     
     echo Exporting ASV table ...
     ## Export ASV-based feature table
     $QIIME tools export \
             --input-path ${WD}/Qobj/dada2_ASVs_collapsed.qza \
             --output-path ${WD}/Export/dada2

####DEBLUR####
     ## Quality filtering
     $QIIME quality-filter q-score \
             --i-demux ${WD}/Qobj/trimmed_sequences.qza \
             --o-filtered-sequences ${WD}/Qobj/qfilter_seq \
             --o-filter-stats ${WD}/Qobj/qfilter_stats 
    
    echo Deblur processing ...
    $QIIME deblur denoise-16S \
            --i-demultiplexed-seqs ${WD}/Qobj/qfilter_seq.qza \
            # trim reads to N length
            --p-trim-length 120 \
            --o-representative-sequences ${WD}/Qobj/rep-seqs-deblur.qza \
            --o-table ${WD}/Qobj/table-deblur.qza \
            --p-sample-stats \
            --o-stats ${WD}/Qobj/deblur-stats.qza \

     ## Classification of deblur ASVs from dada2 using blast
     echo Assigning taxonomy with blast and ${CLASSIFIER_seq} ...
     $QIIME feature-classifier classify-consensus-blast \
             --i-reference-reads $CLASSIFIER_seq \
             --i-reference-taxonomy $CLASSIFIER_tax \
             --i-query ${WD}/Qobj/rep-seqs-deblur.qza \
             --p-strand both \
             --o-classification ${WD}/Qobj/deblur_rep_taxonomy

     ## Collapse deblur-based feature table by taxonomy
     $QIIME taxa collapse \
             --i-table ${WD}/Qobj/table-deblur.qza \
             --i-taxonomy ${WD}/Qobj/deblur_rep_taxonomy.qza \
             --o-collapsed-table ${WD}/Qobj/deblur_ASVs_collapsed.qza \
             --p-level 6

     ## Export deblur-based feature table
     echo Exporting ASV table ...
     $QIIME tools export \
             --input-path ${WD}/Qobj/deblur_ASVs_collapsed.qza \
             --output-path ${WD}/Export/deblur

     ## Enter Export folder
     cd ${WD}/Export

     ## Convert feature tables from biom format to tsv format
     echo Converting BIOM to TSV format ...
	 biom convert \
             --to-tsv \
             --input-fp dada2/feature-table.biom \
             --output-fp dada2/feature-table_ASV.tsv

	 biom convert \
             --to-tsv \
             --input-fp deblur/feature-table.biom \
             --output-fp deblur/feature-table_ASV.tsv

####MAKE TAXONOMY BARCHART####
     echo Generating taxonomy barcharts ...
     $QIIME taxa barplot \
             --i-table ${WD}/Qobj/dada2_ASVs.qza \
             --i-taxonomy ${WD}/Qobj/dada2_rep_taxonomy.qza \
             --o-visualization ${WD}/Export/dada2/dada2_tax_barplot \
             --m-metadata-file $MAP

     $QIIME taxa barplot \
             --i-table ${WD}/Qobj/table-deblur.qza \
             --i-taxonomy ${WD}/Qobj/deblur_rep_taxonomy.qza \
             --o-visualization ${WD}/Export/deblur/deblur_tax_barplot \
             --m-metadata-file $MAP             
     echo Processing complete!
    
