#_____________________________________________________________________________________________
# RNA002 vs RNA004 SCRIPT: FEATURE COUNTS AND QC PLOTS
import re
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import itertools


colors_input = {"RNA002_UHRR" : "#005f73", 
"RNA002_HEK293" :  "#0a9396",
 "RNA002_blood" :"#94d2bd", 
"RNA004_blood" : "#ee9b00", 
"RNA004_HEK293" : "#ca6702", 
"RNA004_UHRR" : "#bb3e03",
"Illumina" : "#e9d8a6"}

# Load feature counts table for all 9 samples 
counts_file="RNA002_RNA004.GRCh38.featurecounts.tsv"
counts = pd.read_csv(counts_file, sep = "\t", header = 0, comment='#')
counts = counts.set_index("Geneid")
counts = counts.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
counts.columns = [x.split('/')[4] for x in counts.columns]

counts.columns = [re.sub("(_basecall.GRCh38.bam.*)$", "", x) for x in counts.columns]
counts.columns = [re.sub("_A", "_HEK293_A", x) for x in counts.columns ]
counts.columns  = [re.sub("_B", "_HEK293_B", x) for x in counts.columns ]
counts.columns  = [re.sub("_Direct", "_normal", x) for x in counts.columns ]

# Create condition metadata file for coloring
condition_Val  = [re.sub("_[A|B]", "", x) for x in counts.columns ]
condition_Val = [re.sub("(blood_.*)$", "blood", x) for x in condition_Val]
condition_Val = [re.sub("(_[1|2].*)$", "", x) for x in condition_Val]
samples = pd.DataFrame({"ID": counts.columns, "Condition" : condition_Val}).set_index("ID")
counts_df = counts
metadata = samples


metadata['colors'] = [colors_input[i] for i in metadata.Condition]

# Load feature counts table from Illumina data
counts_file="GSE47774_count_table.gencodev43.txt"
counts_ngs = pd.read_csv(counts_file, sep = "\t", header = 0, comment='#')
counts_ngs = counts_ngs.set_index("Geneid")
counts_ngs = counts_ngs.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
counts_ngs.columns = [x.split('/')[2] for x in counts_ngs.columns]
counts_ngs.columns = [re.sub("(.gencodev43.Aligned.sortedByCoord.out.bam.*)$", "", x) for x in counts_ngs.columns]

# Calculate mean counts over all samples
counts_ngs['mean'] = counts_ngs.mean(axis=1)
counts_df['Illumina'] = counts_ngs['mean']
condition_Val_new = condition_Val
condition_Val_new.append("Illumina")

samples2 = pd.DataFrame({"ID": counts_df.columns, "Condition" : condition_Val_new}).set_index("ID")
metadata = samples2

# Plot Gene Counts for all genes with coverage >= 10 
df_covered_genes = (counts_df >= 10).sum(axis=0)

gene_counts_10 = pd.DataFrame(df_covered_genes)
gene_counts_10 = gene_counts_10.reset_index()
gene_counts_10.columns = ["sample_type", "# genes covered > 10"]

gene_counts_10  = pd.merge(metadata.reset_index(), gene_counts_10, left_on='ID', right_on='sample_type', how = "left")

fig, ax = plt.subplots(figsize=(4.5,4))
sns.barplot(gene_counts_10, y="# genes covered > 10", x="sample_type", hue = "Condition", ax=ax, palette = colors_input, edgecolor='black')
ax.axhline(12000, ls='--', color = "black")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
fig.tight_layout()

hatches = ['','/', '', '', '/', '', '/', '', '/']
for i, bar in enumerate(ax.patches):
   
    if i < 9: 
       bar.set_hatch(hatches[i]) 
ax.legend([],[], frameon=False)       
ax.set(xlabel='', ylabel='Number of genes covered >= 10')
fig.savefig("PAPER_gene_counts_plot_cov10.svg")
df_covered_genes = (counts_df >= 1).sum(axis=0)

gene_counts_1 = pd.DataFrame(df_covered_genes)
gene_counts_1 = gene_counts_1.reset_index()
gene_counts_1.columns = ["sample_type", "# genes covered > 0"]
gene_counts_1  = pd.merge(metadata.reset_index(), gene_counts_1, left_on='ID', right_on='sample_type', how = "left")

fig, ax = plt.subplots(figsize=(4.5,4))
sns.barplot(gene_counts_1, y="# genes covered > 0", x="sample_type", hue = "Condition", ax=ax, palette = colors_input, edgecolor='black')
ax.axhline(12000, ls='--', color = "black")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
fig.tight_layout()
hatches = ['','/', '', '', '/', '', '/', '', '/']
for i, bar in enumerate(ax.patches):
   
    if i < 9: 
       bar.set_hatch(hatches[i]) 
ax.legend([],[], frameon=False)       
ax.set(xlabel='', ylabel='Number of genes covered > 0')
fig.savefig("PAPER_gene_counts_plot_cov1.svg")

#_____________________________________________________________________________________________
# QUALITY PLOTS
# Load stat output from NanoComp
nanocomp_stats = pd.read_csv("NanoComp_RNA002_RNA004/NanoComp-data.tsv.gz", sep = "\t")
nanocomp_stats_filt = nanocomp_stats[(nanocomp_stats['lengths'] >= 200) & (nanocomp_stats['lengths'] <= 150000)]
nanocomp_stats_filt = pd.merge(nanocomp_stats_filt, metadata.reset_index(), left_on = "dataset", right_on = "ID", how = "left")

# Plot average base call quality
fig, ax = plt.subplots(figsize=(5,5))
sns.violinplot(nanocomp_stats_filt, y="quals", x="dataset", palette = samples['colors'].to_dict(), inner = None)
#ax.tick_params(axis='x', labelrotation=90)
#ax.axhline(12000, ls='--', color = "black")
hatches = ['','//', '', '', '//', '', '//', '', '//']
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
fig.tight_layout()
for c,i in enumerate(ax.get_children()):
    if isinstance(i, mpl.collections.PolyCollection):
        i.set_hatch(hatches[c])
        
ax.legend([],[], frameon=False)
ax.set_ylim(0,33)   
ax.set(xlabel='', ylabel='Average base call quality score')
fig.savefig("PAPER_violin_plot_base_quality.svg")

# Plot percent identity
fig, ax = plt.subplots(figsize=(5,5))

sns.violinplot(nanocomp_stats_filt, y="percentIdentity", x="dataset", palette = samples['colors'].to_dict(), inner = None)
ax.set_ylim(75, 101)
hatches = ['','//', '', '', '//', '', '//', '', '//']
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
fig.tight_layout()
for c,i in enumerate(ax.get_children()):
    if isinstance(i, mpl.collections.PolyCollection):
        i.set_hatch(hatches[c])
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
fig.tight_layout()        
ax.legend([],[], frameon=False)       
ax.set(xlabel='', ylabel='Percent reference identity\n(GRCh38)')
fig.savefig("PAPER_violin_plot_percentIdentity.svg")

# Plot N50 read length
nanocomp_stats = pd.read_csv("NanoComp_RNA002_RNA004/NanoStats.txt", sep = "\t")
nanostats = nanocomp_stats.set_index("Metrics").transpose()
nanostats = nanostats.reset_index()
nanostats.n50 = pd.to_numeric(nanostats.n50)

fig, ax = plt.subplots(figsize=(4.5,4))
sns.barplot(nanostats, y="n50", x="index", hue = "index", ax=ax, palette = samples['colors'].to_dict(), edgecolor='black')
#ax.tick_params(axis='x', labelrotation=45)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
fig.tight_layout()

hatches = ['','/', '', '', '/', '', '/', '', '/']
for i, bar in enumerate(ax.patches):
   
    if i < 9: 
       bar.set_hatch(hatches[i]) 
ax.legend([],[], frameon=False)       
ax.set(xlabel='', ylabel='N50')
fig.savefig("PAPER_n50_bar_plot.svg")