#Enviornment setup
WD=/home/lymelab/lab_members/mann/ # CHANGE ME to your working directory
BCL=/home/lymelab/lab_members/mann/cellranger-tiny-bcl-1.2.0 # CHANGE ME to your bcl folder
FASTQS=test_demultiplex # CHANGE ME to the name of the output folder where your fastqs will be written
SAMPS=/home/lymelab/lab_members/mann/cellranger-tiny-bcl-samplesheet-1.2.0.csv # CHANGE ME to the location of your sample sheet (Illumina format)
REF=/home/lymelab/reference_databases/refdata-cellranger-GRCh38-and-mm10-3.1.0 # CHANGE ME to the path to your reference transcriptome
EXP=test_exp # CHANGE ME to the name of the output folder where your gene expression data wil be written
NAME=test_sample # CHANGE ME to the name of the sample being analyzed in cellranger count
cd $WD

#Demultiplex and quality check from BCL output
cellranger mkfastq --id=$FASTQS --run=$BCL --samplesheet=$SAMPS --qc
#Generate feature/gene expression counts
cellranger count --id=$EXP --transcriptome=$REF --fastqs=$FASTQS/outs/fastq_path --sample=$NAME --expect-cells=1000