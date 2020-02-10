#documentation of nf-core/rnaseq workflow: https://github.com/nf-core/rnaseq
#citation: https://www.biorxiv.org/content/10.1101/610741v1

##################
#Install nextflow
##################
#following three lines only need to be run if nextflow not already installed
# java -version # should be v8+
# curl -fsSL get.nextflow.io | bash
# sudo mv nextflow /usr/local/bin
nextflow pull nf-core/rnaseq # update pipeline

###################
#Envrionment setup
###################
EMAIL=allison.e.mann@gmail.com # change to your email, will send you overview of results
REF=Rnor_6.0 # change to the reference genome you want to map your reads to -- for full list: https://github.com/nf-core/rnaseq/blob/master/conf/igenomes.config
OUT=/home/lymelab/lab_members/mann/rna_seq # output directory
RAW=/home/lymelab/data/raw_data/RNA-seq/raw # path to your folder that holds your raw fastq.gz 
MAX_MEM=100GB # set to the maximum amount of memory you want to use - this will be machine specific
ETYPE=--email # by default will send result email, can set to --email_on_fail if you only want to be sent an email if the workflow is not successful
cd $RAW

#############################
#RNA seq processing pipeline
#############################
#using nextflow and nf-core rnaseq pipeline
#single end
nextflow run nf-core/rnaseq --reads '*.fastq.gz' --genome $REF -profile conda $ETYPE $EMAIL --singleEnd --max_memory $MAX_MEM &> "$(date +%Y-%m-%d_%H-%M-%S.txt)"
#paired end
# nextflow run nf-core/rnaseq --reads '*_R{1,2}*.fastq.gz' --genome $REF -profile conda $ETYPE $EMAIL --outdir $OUT --max_memory $MAX_MEM &> "$(date +%Y-%m-%d_%H-%M-%S.txt)"

#output documentation: https://github.com/nf-core/rnaseq/blob/master/docs/output.md



