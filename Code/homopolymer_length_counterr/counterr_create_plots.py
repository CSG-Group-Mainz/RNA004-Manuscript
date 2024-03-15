import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_true_len_by_base(true_len: int, df: pd.DataFrame, sample1: str, sample2: str) -> plt.Figure:
    df = df.loc[(df["Sample"].isin([sample1, sample2])) & (df["True length"]==true_len)]
    base_counts = {sample1: {}, sample2: {}}
    for s in [sample1, sample2]:
        for b in ["A", "C", "G", "U"]:
            base_counts[s][b] = int(df.loc[(df["Sample"]==s) & (df["Base"]==b), "Count"].sum())

    fig, ax = plt.subplots(figsize=(10,8))

    ax.axvline(true_len, ymin=0, ymax=1, linestyle="--", color="black")

    sns.lineplot(df, x="Homopolymer length", y="Rel. count", hue="Base", style="Sample", 
                 hue_order=["A", "C", "G", "U"], palette=["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"], linewidth=2, alpha=0.75, dashes={sample1: (2,1), sample2: (1, 0)},
                 ax=ax)
    
    legend = ax.legend(title_fontsize=12, fontsize=11)
    for i, b in zip([1,2,3,4], ["A", "C", "G", "U"]):
        legend.get_texts()[i].set_text(f"{b}, counts:\nRNA2: {base_counts[sample1][b]}\nRNA4: {base_counts[sample2][b]}")

    ax.set_xlabel("Homopolymer Length", fontsize=14)
    ax.set_ylabel("Relative Count", fontsize=14)
    # ax.set_xlim(true_len-3,true_len+3)
    ax.set_xlim(0,15)
    ax.set_title(f"Observed homopolymer lengths for true length {true_len}", fontsize=16)

    fig.tight_layout()
    return fig

def plot_true_len(true_len: int, df: pd.DataFrame, sample1: str, sample2: str) -> plt.Figure:
    df = df.loc[(df["Sample"].isin([sample1, sample2])) & (df["True length"]==true_len)]
    base_counts = {sample1: {}, sample2: {}}
    for s in [sample1, sample2]:
        for b in ["A", "C", "G", "U"]:
            base_counts[s][b] = int(df.loc[(df["Sample"]==s) & (df["Base"]==b), "Count"].sum())

    fig, ax = plt.subplots(figsize=(10,8))

    ax.axvline(true_len, ymin=0, ymax=1, linestyle="--", color="black")

    sns.lineplot(df, x="Homopolymer length", y="Rel. count", hue="Sample", 
                 palette="tab10", linewidth=3, alpha=0.75, ax=ax)
    
    ax.set_xlabel("Observed homopolymer Length", fontsize=14)
    ax.set_ylabel("Frequency", fontsize=14)
    # ax.set_xlim(true_len-3,true_len+3)
    ax.set_xlim(0,15)
    ax.set_title(f"Observed homopolymer lengths for true length {true_len}", fontsize=16)

    fig.tight_layout()
    return fig

def create_subplots_by_base(df: pd.DataFrame, sample1: str, sample2: str):

    fig, ax = plt.subplots(ncols=2, nrows=4, figsize=(15,20), sharey=True)

    for i, (r, c) in zip(range(3,11), [(0,0), (1,0), (2,0), (3,0), (0,1), (1,1), (2,1), (3,1)]):
        df_subset = df.loc[(df["Sample"].isin([sample1, sample2])) & (df["True length"]==i)]

        ax[r,c].axvline(i, ymin=0, ymax=1, linestyle="--", color="black")

        sns.lineplot(df_subset, x="Homopolymer length", y="Rel. count", hue="Base", style="Sample", 
                    hue_order=["A", "C", "G", "U"], palette=["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"], linewidth=2, alpha=0.75, dashes={sample1: (2,1), sample2: (1, 0)},
                    ax=ax[r,c])
        ax[r,c].legend(title_fontsize=12, fontsize=11)
        if r==3:
            ax[r,c].set_xlabel("Observed homopolymer Length", fontsize=14)
        else:
            ax[r,c].set_xlabel(None)
        ax[r,c].set_ylabel("Frequency", fontsize=14)
        # ax.set_xlim(true_len-3,true_len+3)
        ax[r,c].set_xlim(0,15)
        ax[r,c].set_title(f"True homopolymer length: {i}", fontsize=16, ha="left", x=0)
    fig.tight_layout()
    return fig

def create_subplots(df: pd.DataFrame, sample1: str, sample2: str):

    fig, ax = plt.subplots(ncols=2, nrows=4, figsize=(15,20), sharey=True)

    for i, (r, c) in zip(range(3,11), [(0,0), (1,0), (2,0), (3,0), (0,1), (1,1), (2,1), (3,1)]):
        df_subset = df.loc[(df["Sample"].isin([sample1, sample2])) & (df["True length"]==i)]

        ax[r,c].axvline(i, ymin=0, ymax=1, linestyle="--", color="black")

        sns.lineplot(df_subset, x="Homopolymer length", y="Rel. count", hue="Sample", 
                    palette="tab10", linewidth=2, alpha=0.75, ax=ax[r,c])
        ax[r,c].legend(title_fontsize=12, fontsize=11)
        if r==3:
            ax[r,c].set_xlabel("Observed homopolymer Length", fontsize=14)
        else:
            ax[r,c].set_xlabel(None)
        ax[r,c].set_ylabel("Frequency", fontsize=14)
        # ax.set_xlim(true_len-3,true_len+3)
        ax[r,c].set_xlim(0,15)
        ax[r,c].set_title(f"True homopolymer length: {i}", fontsize=16, ha="left", x=0)
    fig.tight_layout()
    return fig

dirnames = ['A_RNA002_GRCh38_NCBIRefseq', 
            'A_RNA004_GRCh38_NCBIRefseq', 
            '1_RNA004_GRCh38_NCBIRefseq', 
            'B_RNA002_GRCh38_NCBIRefseq', 
            '2_RNA004_GRCh38_NCBIRefseq', 
            'Control_RNA002_GRCh38_NCBIRefseq',
            'B_RNA004_GRCh38_NCBIRefseq',
            '1_RNA002_GRCh38_NCBIRefseq']

df_subsets = []

for dirname in dirnames:
    d = np.load(f"/home/vincent/projects/charlotte_rna2_4/counterr_output/{dirname}/stats/hist_len_hp.npy")
    for i in range(4):
        for j in range(20):
            base = [i] * 20
            true_len = [j] * 20
            vals = list(d[i,j,:])
            df_subset = pd.DataFrame({"Base": base, "True length": true_len, "Count": vals})
            df_subset["Rel. count"] = df_subset["Count"] / df_subset["Count"].sum()
            df_subset["Sample"] = dirname
            df_subsets.append(df_subset)

data = pd.concat(df_subsets).reset_index().rename(columns={"index": "Homopolymer length"})
data["Base"] = data["Base"].replace({0:"A", 1:"C", 2:"G", 3:"U"})
data["Sample"] = data["Sample"].str.replace("_GRCh38_NCBIRefseq", "")

for a, b in [("A_RNA002", "A_RNA004"),
             ("B_RNA002", "B_RNA004"),
             ("1_RNA002", "1_RNA004")]:
    for i in range(3,11):
        fig = plot_true_len(i, data, a, b)
        fig.savefig(f"/home/vincent/projects/charlotte_rna2_4/plots/{a}_{b}_true_len_{i}.svg")

        fig = plot_true_len_by_base(i, data, a, b)
        fig.savefig(f"/home/vincent/projects/charlotte_rna2_4/plots/{a}_{b}_true_len_{i}_bybase.svg")

    fig = create_subplots_by_base(data, a, b)
    fig.savefig(f"/home/vincent/projects/charlotte_rna2_4/plots/{a}_{b}_true_len_combined_bybase.svg")
    fig = create_subplots(data, a, b)
    fig.savefig(f"/home/vincent/projects/charlotte_rna2_4/plots/{a}_{b}_true_len_combined.svg")