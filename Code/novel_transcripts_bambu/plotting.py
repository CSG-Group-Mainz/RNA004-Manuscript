import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import upsetplot 
from matplotlib_venn import venn2

def plot_venn(data: pd.DataFrame, a: str, b: str, title: str = "Number of novel transcripts"):
    a_set = set(data.loc[data[a] > 0, "TXNAME"])
    b_set = set(data.loc[data[b] > 0, "TXNAME"])
    fig, ax = plt.subplots(figsize=(5,5))
    venn2([a_set, b_set], set_labels=[a, b], set_colors=["#d62728", "#2ca02c"], alpha=0.8, ax=ax)
    ax.set_title(title)
    return fig

def plot_upset(data: pd.DataFrame, a: str, b: str, title: str = "Number of novel transcripts"):
    data = data.copy()
    data[a] = data[a] > 0
    data[b] = data[b] > 0

    data = data.set_index([a,b])["TXNAME"]

    fig = plt.figure()
    upsetplot.plot(data, subset_size="count", fig=fig)
    fig.suptitle(title)
    return fig

def plot_scatter(data: pd.DataFrame, a: str, b: str, title: str = "Transcript counts"):
    fig, ax = plt.subplots()
    sns.scatterplot(data, x=a, y=b, ax=ax)
    ax.set_xlabel(f"Count {a}")
    ax.set_ylabel(f"Count {b}")
    ax.set_title(title)
    return fig


# 1-2 novel transcripts 

counts = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/1_2_RNA002_RNA004/1_2_novel/counts_transcript.txt", sep="\t").rename(columns={"1_RNA002_all_hac.GRCh38": "1_RNA002",
                                                                                                                        "1_2_RNA004_sup.GRCh38": "1_2_RNA004"})
p = plot_venn(counts, "1_RNA002", "1_2_RNA004")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/1_2_RNA002_004_venn_novel.svg")

p = plot_upset(counts, "1_RNA002", "1_2_RNA004")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/1_2_RNA002_004_upset_novel.svg")

p = plot_scatter(counts, "1_RNA002", "1_2_RNA004", "Novel transcript counts")
counts.loc[counts["1_2_RNA004"]>800, "text"] = counts.loc[counts["1_2_RNA004"]>800, "TXNAME"] + "\n" + counts.loc[counts["1_2_RNA004"]>800, "GENEID"]
ax = p.axes[0]
for i, r in counts.loc[counts["text"].notna()].iterrows():
    ax.text(x=r[2]+5, y=r[3]+50, s=r[4], fontsize=8)
ax.set_xlim(xmax=650)
ax.set_ylim(ymax=8700)

p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/1_2_RNA002_004_scatter_novel.svg")


# 1-2 full-length transcripts 

counts = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/1_2_RNA002_RNA004/1_2_fulllen/counts_transcript.txt", sep="\t").rename(columns={"1_RNA002_all_hac.GRCh38": "1_RNA002",
                                                                                                                        "1_2_RNA004_sup.GRCh38": "1_2_RNA004"})

p = plot_venn(counts, "1_RNA002", "1_2_RNA004", "Number of full-length transcripts")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/1_2_RNA002_004_venn_fulllen.svg")

p = plot_upset(counts, "1_RNA002", "1_2_RNA004", "Number of full-length transcripts")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/1_2_RNA002_004_upset_fulllen.svg")

p = plot_scatter(counts, "1_RNA002", "1_2_RNA004", "Full-length transcript counts")
counts.loc[(counts["1_RNA002"]>10000) | (counts["1_2_RNA004"]>100000), "text"] = counts.loc[(counts["1_RNA002"]>3000) | (counts["1_2_RNA004"]>33000), "TXNAME"] + "\n" + counts.loc[(counts["1_RNA002"]>3000) | (counts["1_2_RNA004"]>33000), "GENEID"]
ax = p.axes[0]
for i, r in counts.loc[counts["text"].notna()].iterrows():
    ax.text(x=r[2]+80, y=r[3]+100, s=r[4], fontsize=8)
ax.set_xlim(xmax=35000)
ax.set_ylim(ymax=300000)

p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/1_2_RNA002_004_scatter_fulllen.svg")


# A-B novel transcripts 

counts = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/A_B_RNA002_RNA004/A_B_novel/counts_transcript.txt", sep="\t").rename(columns={"A_B_RNA002_hac.GRCh38": "A_B_RNA002",
                                                                                                                        "A_B_RNA004_sup.GRCh38": "A_B_RNA004"})

p = plot_venn(counts, "A_B_RNA002", "A_B_RNA004")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/A_B_RNA002_004_venn_novel.svg")

p = plot_upset(counts, "A_B_RNA002", "A_B_RNA004")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/A_B_RNA002_004_upset_novel.svg")

p = plot_scatter(counts, "A_B_RNA002", "A_B_RNA004", "Novel transcript counts")
counts.loc[(counts["A_B_RNA004"]>5000) | (counts["A_B_RNA002"]>900), "text"] = counts.loc[(counts["A_B_RNA004"]>5000) | (counts["A_B_RNA002"]>900), "TXNAME"] + "\n" + counts.loc[(counts["A_B_RNA004"]>5000) | (counts["A_B_RNA002"]>900), "GENEID"]
ax = p.axes[0]
for i, r in counts.loc[counts["text"].notna()].iterrows():
    ax.text(x=r[2]+10, y=r[3], s=r[4], fontsize=8)
ax.set_xlim(xmax=14000)
ax.set_ylim(ymax=50000)

p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/A_B_RNA002_004_scatter_novel.svg")


# A-B full-length transcripts 

counts = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/bambu_results_gencode/A_B_RNA002_RNA004/A_B_fulllen/counts_transcript.txt", sep="\t").rename(columns={"A_B_RNA002_hac.GRCh38": "A_B_RNA002",
                                                                                                                        "A_B_RNA004_sup.GRCh38": "A_B_RNA004"})

p = plot_venn(counts, "A_B_RNA002", "A_B_RNA004", "Number of full-length transcripts")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/A_B_RNA002_004_venn_fulllen.svg")

p = plot_upset(counts, "A_B_RNA002", "A_B_RNA004", "Number of full-length transcripts")
p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/A_B_RNA002_004_upset_fulllen.svg")

p = plot_scatter(counts, "A_B_RNA002", "A_B_RNA004", "Full-length transcript counts")
counts.loc[(counts["A_B_RNA002"]>50000) | (counts["A_B_RNA004"]>200000), "text"] = counts.loc[(counts["A_B_RNA002"]>50000) | (counts["A_B_RNA004"]>200000), "TXNAME"] + "\n" + counts.loc[(counts["A_B_RNA002"]>50000) | (counts["A_B_RNA004"]>200000), "GENEID"]
ax = p.axes[0]
for i, r in counts.loc[counts["text"].notna()].iterrows():
    ax.text(x=r[2], y=r[3], s=r[4], fontsize=8)
ax.set_xlim(xmax=120000)
ax.set_ylim(ymax=500000)

p.savefig("/home/vincent/projects/charlotte_rna2_4/novel_transcripts/plots/A_B_RNA002_004_scatter_fulllen.svg")