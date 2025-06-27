#this compilation of commands can be used to download annotations used in the paper by Hewel and Wierczeiko et al.,
#ont oligos and bed files with modified sites annotation
aws s3 sync --no-sign-request s3://ont-open-data/rna-modbase-validation_2025.03/references references
#GRCh38 GENCODE primary assembly
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
#corresponding GENCODE gtf annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
#chrR genome hg19
wget https://github.com/vikramparalkar/rDNA-Mapping-Genomes/blob/main/Human_hg19-rDNA_genome_v1.0.tar.gz
