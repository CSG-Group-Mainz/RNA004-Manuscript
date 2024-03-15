library(DSS)
require(bsseq)

if (length(commandArgs(trailingOnly = TRUE)) != 8) {
  stop("Usage: Rscript myscript.R input_path1.1 input_path1.2 input_path2.1 input_path2.2 name1 name2 output_path")
}
input_sample1_1 = commandArgs(trailingOnly = TRUE)[1]
input_sample1_2 = commandArgs(trailingOnly = TRUE)[2]
input_sample2_1 = commandArgs(trailingOnly = TRUE)[3]
input_sample2_2 = commandArgs(trailingOnly = TRUE)[4]

name1 = commandArgs(trailingOnly = TRUE)[5]
name2 = commandArgs(trailingOnly = TRUE)[6]
out_prefix = commandArgs(trailingOnly = TRUE)[7]

name1.1 = paste0(name1, "_1")
name1.2 = paste0(name1, "_2")
name2.1 = paste0(name2, "_1")
name2.2 = paste0(name2, "_2")
output_dir = commandArgs(trailingOnly = TRUE)[8]

print(paste0("###### Starting differential methylation analysis between samples ", name1, " & ", name2, " with two replicates ######"))

print(paste0("Reading: ", input_sample1_1))
data_sample1.1 = read.table(input_sample1_1, header = TRUE)
print(paste0("Reading: ", input_sample1_2))
data_sample1.2 = read.table(input_sample1_2, header = TRUE)

print(paste0("Reading: ", input_sample2_1))
data_sample2.1 = read.table(input_sample2_1, header = TRUE)
print(paste0("Reading: ", input_sample2_2))
data_sample2.2 = read.table(input_sample2_2, header = TRUE)

print("Creating data object")
BSobj = makeBSseqData(list(data_sample1.1, data_sample1.2, data_sample2.1, data_sample2.2), c(name1.1, name1.2, name2.1, name2.2))

print("Performing DML test")
dmltest = DMLtest(BSobj, group1 = c(name1.1, name1.2), group2 = c(name2.1, name2.2), smoothing = TRUE, ncores = 16)

print("Calling DML")
outpath = file.path(output_dir, paste0(out_prefix, "_", name1, "_", name2, "_dml.tsv"))
dmls = callDML(dmltest, p.threshold=1)
print(paste0("Writing DML to ", outpath))
write.table(dmls, file = outpath, quote = FALSE, row.names = FALSE, sep="\t")

print("######       DONE        ######")
