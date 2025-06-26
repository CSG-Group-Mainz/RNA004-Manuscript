import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.offsetbox import AnchoredText
import functools as ft
import seaborn.objects as so


def modkitHIST(file:str, modification:str, color:str, outpath:str, title:str):

    #set plotting param
    matplotlib.rcParams['xtick.labelsize'] = 4

    #Read Table with probabilities 
    tsv = pd.read_csv(file+"probabilities.tsv", sep='\s+')
    tsv_mod = tsv.loc[tsv["code"] == modification]

    #Read Table with Modkit Thresholds
    thresh_tsv = pd.read_csv(file+"thresholds.tsv", sep='\s+')
    if modification == "a":
        thresh_tsv = thresh_tsv[(thresh_tsv["base"] == "A") & (thresh_tsv["percentile"] == 10)]

    if modification == "17802":
        thresh_tsv = thresh_tsv[(thresh_tsv["base"] == "T") & (thresh_tsv["percentile"] == 10)]

        

    #Merge Range Start and Range End buckets
    tsv_mod["bin"] = tsv_mod["range_start"].astype(str) + " - " + tsv_mod["range_end"].astype(str)

    #Find Bin where 0.98 is contained
    bin_098 = tsv_mod[(tsv_mod['range_start'] <= 0.98) & (tsv_mod['range_end'] >= 0.98)]
    bin_098 = bin_098["bin"].to_list()[0]

    #Find Bin where modkit 10th percentile value is contained
    threshold = thresh_tsv["threshold"].to_list()[0]
    bin_modkit = tsv_mod[(tsv_mod['range_start'] <= threshold) & (tsv_mod['range_end'] >= threshold)]
    bin_modkit = bin_modkit["bin"].to_list()[0]

    #Set Axis ticks for plotting
    Bin_labels = tsv_mod["bin"].to_list()
    Bin_labels_ticks = Bin_labels[0::16]
    Bin_labels_ticks.append(Bin_labels[-1])

    #PLOT as barplot
    plt.clf()
    plt.figure(figsize=(3,2))
    sns.set_context("paper")
    ax = sns.barplot(data=tsv_mod, x="bin", y="count", color=color)
    ax.axvline(bin_098, color='black')
    ax.axvline(bin_modkit, color='r')
    ax.set_title(title)
    ax.set_xticks(Bin_labels_ticks)
    plt.xticks(rotation=45)
    plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
    sns.despine()
    plt.savefig(outpath, format="pdf", bbox_inches="tight", dpi=600)


###### HEK

HEK_A_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_A_basecall.0.7.2.GRCh38/",
                       modification="a", color="#CD7D66", outpath="1_HEK_A_m6a.pdf", title="HEK293T_A m6A")

HEK_A_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_A_basecall.0.7.2.GRCh38/",
                       modification="17802", color="#CD7D66", outpath="2_HEK_A_pseu.pdf", title="HEK293T_A pseU")

HEK_B_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_B_basecall.0.7.2.GRCh38/",
                       modification="a", color="#CD7D66", outpath="3_HEK_B_m6a.pdf", title="HEK293T_B m6A")

HEK_B_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_B_basecall.0.7.2.GRCh38/",
                       modification="17802", color="#CD7D66", outpath="4_HEK_B_pseu.pdf", title="HEK293T_B pseU")


###### BLOOD

Blood_1_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_blood.0.7.2.GRCh38/",
                       modification="a", color="#0B3954", outpath="5_Blood_1_m6a.pdf", title="Blood_1 m6A")

Blood_1_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_blood.0.7.2.GRCh38/",
                       modification="17802", color="#0B3954", outpath="6_Blood_1_pseu.pdf", title="Blood 1 pseU")

Blood_2_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_S5_DRS_basecall.0.7.2.GRCh38/",
                       modification="a", color="#0B3954", outpath="7_Blood_2_m6a.pdf", title="Blood_2 m6A")

Blood_2_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_S5_DRS_basecall.0.7.2.GRCh38/",
                       modification="17802", color="#0B3954", outpath="8_Blood_2_pseu.pdf", title="Blood_2 pseU")

Blood_3_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_S6_DRS_basecall.0.7.2.GRCh38/",
                       modification="a", color="#0B3954", outpath="9_Blood_3_m6a.pdf", title="Blood_3 m6A")

Blood_3_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_S6_DRS_basecall.0.7.2.GRCh38/",
                       modification="17802", color="#0B3954", outpath="10_Blood_3_pseu.pdf", title="Blood_3 pseU")

Blood_IVT_m6A = modkitHIST(file="hist_dirs/hist_dirmerged_RNA004_S6_IVT/",
                       modification="a", color="#098C9A", outpath="11_Blood_IVT_m6A.pdf", title="Blood_IVT m6A")

Blood_IVT_pseu = modkitHIST(file="hist_dirs/hist_dirmerged_RNA004_S6_IVT/",
                       modification="17802", color="#098C9A", outpath="12_Blood_IVT_pseu.pdf", title="Blood_IVT pseU")


###### OLIGOS 

Oligo_1_m6A = modkitHIST(file="hist_dirs/new_sample_probs_oligo1/new_oligo1_",
                       modification="a", color="#8B5CF6", outpath="13_Oligo_1_m6A.pdf", title="Oligo_1 m6A")

Oligo_1_pseu = modkitHIST(file="hist_dirs/new_sample_probs_oligo1/new_oligo1_",
                       modification="17802", color="#8B5CF6", outpath="14_Oligo_1_pseu.pdf", title="Oligo_1 pseU")

Oligo_2_m6A = modkitHIST(file="hist_dirs/new_sample_probs_oligo2/new_oligo2_",
                       modification="a", color="#8B5CF6", outpath="15_Oligo_2_m6A.pdf", title="Oligo_2 m6A")

Oligo_2_pseu = modkitHIST(file="hist_dirs/new_sample_probs_oligo2/new_oligo2_",
                       modification="17802", color="#8B5CF6", outpath="16_Oligo_2_pseu.pdf", title="Oligo_2 pseU")


##### Vectors 

Vector_A_m6A = modkitHIST(file="hist_dirs/last_hist_dirA_vector_004.0.7.2.vector/",
                       modification="a", color="#8B5CF6", outpath="17_Vector_A_m6A.pdf", title="Vector_A m6A")

Vector_A_pseu = modkitHIST(file="hist_dirs/last_hist_dirA_vector_004.0.7.2.vector/",
                       modification="17802", color="#8B5CF6", outpath="18_Vector_A_pseu.pdf", title="Vector_A pseU")

Vector_B_m6A = modkitHIST(file="hist_dirs/last_hist_dirB_vector_004.0.7.2.vector/",
                       modification="a", color="#8B5CF6", outpath="19_Vector_B_m6A.pdf", title="Vector_B m6A")

Vector_B_pseu = modkitHIST(file="hist_dirs/last_hist_dirB_vector_004.0.7.2.vector/",
                       modification="17802", color="#8B5CF6", outpath="20_Vector_B_pseu.pdf", title="Vector_B pseU")

###### ONT Testdaten 

ONT_Ctrl_R1_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.CONTROL.R1.EPIref/",
                       modification="a", color="#2E0788", outpath="21_ONT_Ctrl_R1_m6A.pdf", title="ONT_Ctrl_R1_m6A")


ONT_Ctrl_R1_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.CONTROL.R1.EPIref/",
                       modification="17802", color="#2E0788", outpath="22_ONT_Ctrl_R1_pseu.pdf", title="ONT_Ctrl_R1_pseu")


ONT_Ctrl_R2_m6A = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.CONTROL.R2.EPIref/",
                       modification="a", color="#2E0788", outpath="23_ONT_Ctrl_R2_m6A.pdf", title="ONT_Ctrl_R2_m6A")


ONT_Ctrl_R2_pseu = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.CONTROL.R2.EPIref/",
                       modification="17802", color="#2E0788", outpath="24_ONT_Ctrl_R2_pseu.pdf", title="ONT_Ctrl_R2_pseu")





ONT_m6A_Ctrl_DRACH = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.M6A.CTRL_DRACH.M6Aref/",
                       modification="a", color="#2E0788", outpath="25_ONT_m6A_Ctrl_DRACH.pdf", title="ONT_m6A_Ctrl_DRACH")


ONT_m6A_DRACH = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.M6A.DRACH.M6Aref/",
                       modification="a", color="#2E0788", outpath="26_ONT_m6A_DRACH.pdf", title="ONT_m6A_DRACH")



ONT_m6A_Ctrl_DRACH = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.M6A.CTRL_DRACH.M6Aref/",
                       modification="a", color="#2E0788", outpath="27_ONT_m6A_Ctrl_DRACH.pdf", title="ONT_m6A_Ctrl_DRACH")


ONT_m6A_R1 = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.M6A.R1.EPIref/",
                       modification="a", color="#2E0788", outpath="28_ONT_m6A_R1.pdf", title="ONT_m6A_R1")


ONT_m6A_R2 = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.M6A.R2.EPIref/",
                       modification="a", color="#2E0788", outpath="29_ONT_m6A_R2.pdf", title="ONT_m6A_R2")



ONT_pseu_R1 = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.pseu.R1.EPIref/",
                       modification="17802", color="#2E0788", outpath="30_ONT_pseu_R1.pdf", title="ONT_pseu_R1")



ONT_pseu_R2 = modkitHIST(file="hist_dirs/hist_dirRNA004_ONT.pseu.R2.EPIref/",
                       modification="17802", color="#2E0788", outpath="31_ONT_pseu_R2.pdf", title="ONT_pseu_R2")
