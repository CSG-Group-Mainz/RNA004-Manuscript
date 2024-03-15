library(bambu)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("The following arguments must be given: BAM file, FASTA file, GTF file")
} 

bam.paths = args[1]
fa.path = args[2]
gtf.path = args[3]
out.path = args[4]
bam.paths = strsplit(bam.paths, ",")[[1]]

all.out = paste0(out.path, "/all")
novel.out = paste0(out.path, "/novel")
fulllen.out = paste0(out.path, "/fulllen")

print(bam.paths)
print(fa.path)
print(gtf.path)
print(all.out)
print(novel.out)
print(fulllen.out)
print(paste0(out.path, "/se_object.rds"))

NDR = 1
quant = TRUE

annotations = prepareAnnotations(gtf.path)
se = bambu(reads = bam.paths, annotations = annotations, genome = fa.path, quant = quant, ncore = 16) # NDR = NDR, 
writeBambuOutput(se, path = all.out)

# se.novel = se[mcols(se)$novelTranscript,]
# writeBambuOutput(se.novel, path = novel.out)

# se.fulllen = se[assays(se)$fullLengthCounts >= 1,]
# writeBambuOutput(se.fulllen, path = fulllen.out)

saveRDS(se, paste0(out.path, "/se_object.rds"))





# test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
# fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")

# gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

# bambuAnnotations <- prepareAnnotations(gtf.file)

# se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)

