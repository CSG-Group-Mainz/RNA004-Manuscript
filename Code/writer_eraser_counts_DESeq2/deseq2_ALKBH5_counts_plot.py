import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("/home/vincent/projects/rna2_4/featureCounts/writer_eraser_counts.tsv", sep="\t")
data = data.loc[data["Geneid"] == "ENSG00000091542.9"]
data = data.drop("Geneid", axis=1).T
data = data.reset_index().rename(columns={12: "Norm. Counts", "index": "Sample"})

fig, ax = plt.subplots()
sns.barplot(data, x="Sample", y="Norm. Counts", hue="Sample", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=-40, ha="left", fontsize=12)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
ax.set_ylabel("Norm. Count", fontsize=14)
ax.set_xlabel("Sample", fontsize=14)
ax.set_title("ALKBH5 - Normalized counts", fontsize=16);
fig.tight_layout()
fig.savefig("/home/vincent/projects/rna2_4/featureCounts/ALKBH5_counts.svg")

data = data.loc[data["Sample"].isin(["RNA004_A", "RNA004_B", "RNA004_UHRR_1", "RNA004_UHRR_2"])]

fig, ax = plt.subplots()
sns.barplot(data, x="Sample", y="Norm. Counts", hue="Sample")#, ax=ax, palette=["#1f77b4", "#ff7f0e", "#9467bd", "#8c564b"])
ax.set_xticklabels(ax.get_xticklabels(), rotation=-40, ha="left", fontsize=12)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
ax.set_ylabel("Norm. Count", fontsize=14)
ax.set_xlabel("Sample", fontsize=14)
ax.set_title("ALKBH5 - Normalized counts", fontsize=16);
fig.tight_layout()

fig.savefig("/home/vincent/projects/rna2_4/featureCounts/ALKBH5_counts_UHRR_HEK.svg")