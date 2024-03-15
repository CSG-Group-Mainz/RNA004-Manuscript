##### FEATURE COUNTS 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import re
import pandas as pd 
import matplotlib.pyplot as plt

counts_file="/raid/awiercze_analysis/RNA002_RNA004_paper/RNA004.GRCh38.featurecounts.tsv"
counts = pd.read_csv(counts_file, sep = "\t", header = 0, comment='#')
counts = counts.set_index("Geneid")
counts = counts.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
counts.columns = [x.split('/')[4] for x in counts.columns]


counts.columns = [re.sub("(_basecall.GRCh38.bam.*)$", "", x) for x in counts.columns]
counts.columns
samples = pd.DataFrame({"ID": counts.columns, "Condition" : ["HEK293", "HEK293", "Blood", "Blood", "UHRR", "UHRR"]}).set_index("ID")
counts_df = counts
metadata = samples
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]
counts_df = counts_df.T

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="Condition",
    refit_cooks=True,
    inference=inference,
)
dds.deseq2()

x = pd.DataFrame(dds.layers["normed_counts"]).T
x.index = counts_df.T.index
x.columns = counts_df.T.columns
gtf = "/home/awiercze/gencode.v43.annotation.gtf.gz"
gtf_file = pd.read_csv(gtf, comment="#", sep = "\t", header = None)
gtf_file_gene = gtf_file[gtf_file[2] == "gene"]
gtf_file_gene["ENSEMBL"] = gtf_file_gene[8].str.split('; ', expand=True)[0].str.replace("gene_id ", "").str.replace('"', '')

transcripts = ["METTL16","METTL3","METTL14","WTAP","ZC3H13","RBM15","KIAA1429","FTO","ALKBH5"]
transcripts_out = []
for i in transcripts:
    gtf_sub = gtf_file_gene[gtf_file_gene[8].str.contains(i)]
    print(gtf_sub["ENSEMBL"])
    transcripts_out.append(gtf_sub["ENSEMBL"])
    
    
transcrpts_tmp = [i.to_list() for i in transcripts_out]
count = []
for t in transcrpts_tmp:
    count += t
y = x[x.index.isin(count)]
y.to_csv("/raid/awiercze_analysis/RNA002_RNA004_paper/writer_eraser_counts.tsv", sep = "\t")
y_out = y.reset_index()
plot_table = y_out.melt(id_vars="Geneid")

fig, ax = plt.subplots(figsize=(15,8))
sns.barplot(plot_table, y="value", x="Geneid", hue="variable", ax=ax)
ax.tick_params(axis='x', labelrotation=90)
fig.tight_layout()
fig.savefig("/raid/awiercze_analysis/RNA002_RNA004_paper/writer_eraser_counts_plot.png")