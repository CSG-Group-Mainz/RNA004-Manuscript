import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection

import argparse, os

FDR_THRESH = 0.001

def process(path_f: str, path_r: str, outdir: str) -> None:
    dml_f = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/differential_methylation_analysis/dss_output/forward_HEK293T_UHRR_dml.tsv", sep="\t")
    dml_f["strand"] = "+"
    dml_r = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/differential_methylation_analysis/dss_output/reverse_HEK293T_UHRR_dml.tsv", sep="\t")
    dml_r["strand"] = "-"

    dml = pd.concat([dml_f, dml_r])
    dml["log2FC"] = np.log2(dml["mu1"]/dml["mu2"])
    dml["fdr_comb"] = fdrcorrection(dml["pval"])[1]
    dml["neglog10FDR"] = -np.log10(dml["fdr_comb"])
    dml.loc[(dml["fdr"]>=FDR_THRESH) | ((dml["log2FC"]>=-1) & (dml["log2FC"]<=1)), "color"] = "#7f7f7f"
    dml.loc[(dml["fdr"]<FDR_THRESH) & (dml["log2FC"]<-1), "color"] = "#1f77b4"
    dml.loc[(dml["fdr"]<FDR_THRESH) & (dml["log2FC"]>1), "color"] = "#d62728"
    dml["dml_id"] = dml["chr"] + "_" + dml["pos"].astype(str)

    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(x=dml["log2FC"], y=dml["neglog10FDR"], c=dml["color"], alpha=0.75)
    ax.set_xlabel("log2(Fold change)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Differentially methylated loci (HEK293T vs. UHRR)")

    fig.savefig(os.path.join(outdir,"volcano_HEK_UHRR.svg"))
    fig.savefig(os.path.join(outdir,"volcano_HEK_UHRR.png"))

def process_annot(path_f: str, path_r: str, outdir: str) -> None:
    colnames = ["chr", "start", "end", "mean_methy_HEK293T", "mean_methy_UHRR", "diff_methy", "diff.se", "stat", "phi1", "phi2", "pval", "fdr", "postprob.overThreshold",
            "chr_annot", "start_annot", "end_annot", "strand", "feature", "gene_type", "gene_name", "gene_id", "count"]
    dml_f_annotated = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/differential_methylation_analysis/dss_output_annotated/forward_HEK293T_UHRR_dml.tsv", sep="\t", header=None, names=colnames).drop(["diff.se", "stat", "phi1", "phi2", "fdr", "postprob.overThreshold", "chr_annot", "start_annot", "end_annot", "count"], axis=1)
    dml_f_annotated = dml_f_annotated.loc[dml_f_annotated["feature"]=="gene"]

    dml_r_annotated = pd.read_csv("/home/vincent/projects/charlotte_rna2_4/differential_methylation_analysis/dss_output_annotated/reverse_HEK293T_UHRR_dml.tsv", sep="\t", header=None, names=colnames).drop(["diff.se", "stat", "phi1", "phi2", "fdr", "postprob.overThreshold", "chr_annot", "start_annot", "end_annot", "count"], axis=1)
    dml_r_annotated = dml_r_annotated.loc[dml_r_annotated["feature"]=="gene"]

    dml_annot = pd.concat([dml_f_annotated, dml_r_annotated])
    dml_annot["fdr"] = fdrcorrection(dml_annot["pval"])[1]
    dml_annot["neglog10FDR"] = -np.log10(dml_annot["fdr"])
    dml_annot["log2FC"] = np.log2(dml_annot["mean_methy_HEK293T"]/dml_annot["mean_methy_UHRR"])
    dml_annot["dml_id"] = dml_annot["chr"] + "_" + dml_annot["start"].astype(str) + "_" + dml_annot["end"].astype(str) + "_" + dml_annot["strand"]

    dml_annot_wo_duplicates = dml_annot.drop_duplicates(subset="dml_id")
    dml_annot_wo_duplicates.loc[(dml_annot_wo_duplicates["fdr"]>=FDR_THRESH) | ((dml_annot_wo_duplicates["log2FC"]>=-1) & (dml_annot_wo_duplicates["log2FC"]<=1)), "color"] = "#7f7f7f"
    dml_annot_wo_duplicates.loc[(dml_annot_wo_duplicates["fdr"]<FDR_THRESH) & (dml_annot_wo_duplicates["log2FC"]<-1), "color"] = "#1f77b4"
    dml_annot_wo_duplicates.loc[(dml_annot_wo_duplicates["fdr"]<FDR_THRESH) & (dml_annot_wo_duplicates["log2FC"]>1), "color"] = "#d62728"

    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(x=dml_annot_wo_duplicates["log2FC"], y=dml_annot_wo_duplicates["neglog10FDR"], c=dml_annot_wo_duplicates["color"], alpha=0.75)
    ax.set_xlabel("log2(Fold change)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Differentially methylated loci (HEK293T vs. UHRR)")

    dml_annot.to_csv(os.path.join(outdir, "dml_annotated.tsv"), sep="\t")
    fig.savefig(os.path.join(outdir, "volcano_HEK_UHRR_annotated_loci.svg"))
    fig.savefig(os.path.join(outdir, "volcano_HEK_UHRR_annotated_loci.png"))

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_f", type=str)
    parser.add_argument("--path_f_annot", type=str)
    parser.add_argument("--path_r", type=str)
    parser.add_argument("--path_r_annot", type=str)
    parser.add_argument("--outdir", type=str)

    args = parser.parse_args()
    process(args.path_f, args.path_r, args.outdir)
    process_annot(args.path_f_annot, args.path_r_annot, args.outdir)