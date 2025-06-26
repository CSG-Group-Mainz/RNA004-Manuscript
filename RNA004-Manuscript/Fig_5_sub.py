import upsetplot
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.offsetbox import AnchoredText
import functools as ft
from decimal import *
import statistics

#Function to prepare bedfile as df 
def BedToDf(bedfile:str):
    #read bedfile as pd.csv
    df = pd.read_csv(bedfile, sep='\t')
    print(df)
    #Merge chromosome, start position and strand info into one string
    df_merged = df.iloc[:,0].astype(str)+":"+df.iloc[:,1].astype(str)+":"+df.iloc[:,5].astype(str)
    #Get Ratios and coverages
    ratios = df.iloc[:,10].astype(float)
    coverages = df.iloc[:,4].astype(int)
    #Merge to DF
    dict_formatted = {"Position": list(df_merged), "Coverage": list(coverages), "Ratio": list(ratios)}
    formatted_df = pd.DataFrame(dict_formatted)
    return formatted_df

def GloriToDf(Glorifile:str):
    #Read Glori Sites
    df = pd.read_csv(Glorifile, delimiter='\t')
    #### REFORMAT #### All base positions must be shifted 1 to the left to compare with bedfiles
    df_1_base = list(df.iloc[:,1].astype(int))
    df_1_base = list(map(lambda x : x - 1, df_1_base))
    #Create concatenations
    pos = df.iloc[:,0].astype(str) + ":" + np.asarray(df_1_base).astype(str) + ":" + df.iloc[:,2].astype(str)
    pos = list(pos)

    coverages = df["AGcov"].tolist()
    ratios = df["Ratio"].tolist()
    ratios = list(map(lambda x : x * 100, ratios))

    df_dict = {"Position": pos, "Coverage": coverages, "Ratio": ratios}
    df_formatted = pd.DataFrame(df_dict)
    return df_formatted

def CorrPlot(A:str, B:str, A_name:str, B_name:str, ratio:int, cov:int, outpath:str, median:bool):
    #read and reformat bedfiles
    a = A
    a = a.loc[a["Ratio"] >= ratio]
    a = a.loc[a["Coverage"] >= cov]
    b = B
    b = b.loc[b["Ratio"] >= ratio]
    b = b.loc[b["Coverage"] >= cov]

    mergedata = a.merge(b, on='Position', how='inner', suffixes=('_A', '_B'))
    r = pearsonr(mergedata['Ratio_A'].tolist(), mergedata['Ratio_B'].tolist())
    r = r[0]
    r = round(r, 4)
    subtitle = "Correlation (r) = " + str(r) + ", n = " + str(mergedata.shape[0])
    
    if median:
      #find median
      median_a = round(statistics.median(mergedata["Ratio_A"].tolist()), 2)
      median_b = round(statistics.median(mergedata["Ratio_B"].tolist()), 2)
      A_name = A_name + ", median = " + str(median_a)
      B_name = B_name + ", median = " + str(median_b)


    h = sns.jointplot(x=mergedata["Ratio_A"].tolist(), y=mergedata["Ratio_B"].tolist(), kind="hex", height=3, ratio=5)
    h.set_axis_labels(A_name, B_name)
    anc = AnchoredText(subtitle, loc="upper left", frameon=False,prop=dict(fontsize=6.5, ))
    h.ax_joint.add_artist(anc)
    #h.ax_marg_x.axvline(x=0.010, linewidth=2, color="b")

    h.savefig(outpath)
    
    return None


def DeltaTable(A, B, A_name, B_name, ratio, coverage, outpath):

    #Filter data
    A = A.loc[A["Ratio"] >= ratio]
    A = A.loc[A["Coverage"] >= coverage]
    B = B.loc[B["Ratio"] >= ratio]
    B = B.loc[B["Coverage"] >= coverage]

    suffix_a = "_" + A_name
    suffix_b = "_" + B_name

    col_a = "Ratio" + suffix_a
    col_b = "Ratio" + suffix_b

    mergedata = A.merge(B, on='Position', how='inner', suffixes=(suffix_a, suffix_b))

    #Calculate Pearsons R
    r = pearsonr(mergedata[col_a].tolist(), mergedata[col_b].tolist())
    r = r[0]
    r = round(r, 4)

    #Calculate Delta of the ratio 
    mergedata['Delta'] = mergedata[col_a] - mergedata[col_b]
    pears = "Pearsons R: " + str(r) 
    mergedata[pears] = None

    mergedata.to_csv(outpath)


######### Upset Plot Main Figure 5 including IVT


# #Read bedfiles
# Blood_prev_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_blood.0.7.2.GRCh38_pseu.r1.mod.bed")
# Blood_S5_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_pseu.r1.mod.bed")
# Blood_S6_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S6_DRS_basecall.0.7.2.GRCh38_pseu.r1.mod.bed")
# Blood_S6_IVT_merged_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/merged_RNA004_S6_IVT_pseu.r1.mod.bed")

# ##### Filter Nanopore datasets for a minimum coverage of 10
# Blood_m6A_prev_cov10 = Blood_prev_m6A.loc[Blood_prev_m6A["Coverage"] >= 10]
# Blood_S5_m6A_cov10 = Blood_S5_m6A.loc[Blood_S5_m6A["Coverage"] >= 10]
# Blood_S6_m6A_cov10 = Blood_S6_m6A.loc[Blood_S6_m6A["Coverage"] >= 10]
# Blood_S6_IVT_merged_m6A_cov10 = Blood_S6_IVT_merged_m6A.loc[Blood_S6_IVT_merged_m6A["Coverage"] >= 10]

# ##### Save Coverage filtered datasets
# Blood_m6A_prev_cov10.to_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_1_pseu_cov10.csv")
# Blood_S5_m6A_cov10.to_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_2_pseu_cov10.csv")
# Blood_S6_m6A_cov10.to_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_3_pseu_cov10.csv")
# Blood_S6_IVT_merged_m6A_cov10.to_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_IVT_pseu_cov10.csv")



#Read coverage filtered datasets
Blood_m6A_prev_cov10 = pd.read_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_1_pseu_cov10.csv")
Blood_S5_m6A_cov10 = pd.read_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_2_pseu_cov10.csv")
Blood_S6_m6A_cov10 = pd.read_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_3_pseu_cov10.csv")
Blood_S6_IVT_merged_m6A_cov10 = pd.read_csv("/home/johannes/RNA004_Revision_1/Panel_5/Blood_IVT_pseu_cov10.csv")

#Get sites which are covered in all Blood samples
positions = [Blood_m6A_prev_cov10["Position"].to_list(), Blood_S5_m6A_cov10["Position"].tolist(), Blood_S6_m6A_cov10["Position"].tolist(),
             Blood_S6_IVT_merged_m6A_cov10["Position"].tolist()]

shared_positions = set.intersection(*map(set,positions))
print(len(shared_positions))

##Subset All dataset only to those joint covered sites
Blood_prev_m6A_subsetted = Blood_m6A_prev_cov10.loc[Blood_m6A_prev_cov10['Position'].isin(shared_positions)]
Blood_S5_m6A_subsetted = Blood_S5_m6A_cov10.loc[Blood_S5_m6A_cov10['Position'].isin(shared_positions)]
Blood_S6_m6A_subsetted = Blood_S6_m6A_cov10.loc[Blood_S6_m6A_cov10['Position'].isin(shared_positions)]
Blood_IVT_m6A_subsetted = Blood_S6_IVT_merged_m6A_cov10.loc[Blood_S6_IVT_merged_m6A_cov10['Position'].isin(shared_positions)]

#Sort by Position
Blood_prev_m6A_subsetted = Blood_prev_m6A_subsetted.iloc[Blood_prev_m6A_subsetted['Position'].map({v: k for k, v in enumerate(shared_positions)}).argsort()]
Blood_S6_m6A_subsetted = Blood_S6_m6A_subsetted.iloc[Blood_S6_m6A_subsetted['Position'].map({v: k for k, v in enumerate(shared_positions)}).argsort()]
Blood_S5_m6A_subsetted = Blood_S5_m6A_subsetted.iloc[Blood_S5_m6A_subsetted['Position'].map({v: k for k, v in enumerate(shared_positions)}).argsort()]
Blood_IVT_m6A_subsetted =  Blood_IVT_m6A_subsetted.iloc[Blood_IVT_m6A_subsetted['Position'].map({v: k for k, v in enumerate(shared_positions)}).argsort()]


#Only retain Positions with modification frequency over 10 Percent
Blood_prev_m6A_sig = Blood_prev_m6A_subsetted.loc[Blood_prev_m6A_subsetted["Ratio"] >= 10]
Blood_S6_m6A_sig = Blood_S6_m6A_subsetted.loc[Blood_S6_m6A_subsetted["Ratio"] >= 10]
Blood_S5_m6A_sig = Blood_S5_m6A_subsetted.loc[Blood_S5_m6A_subsetted["Ratio"] >= 10]
Blood_IVT_m6A_sig = Blood_IVT_m6A_subsetted.loc[Blood_IVT_m6A_subsetted["Ratio"] >= 10]

#####################################################################################################################

######### BOXPLOT ##########
#BOX PLOTS
blood_1_box = Blood_prev_m6A_sig
blood_2_box = Blood_S5_m6A_sig
blood_3_box = Blood_S6_m6A_sig
ivt_1_box = Blood_IVT_m6A_sig

plot_data = {
  'Blood 1': blood_1_box["Ratio"].tolist(),
  'Blood 2': blood_2_box["Ratio"].tolist(),
  'Blood 3': blood_3_box["Ratio"].tolist(),
  'Blood IVT': ivt_1_box["Ratio"].tolist(),
}

my_pal = {'Blood 1': '#1A87C7', 'Blood 2': '#1A87C7','Blood 3': '#1A87C7',
          'Glori 1': '#8B5CF6', 'Glori 2': '#8B5CF6',
          'Blood IVT': '#1A87C7'}

plt.clf()
plt.figure(figsize=(2,4))
ax = sns.violinplot(data=plot_data, palette=my_pal, linewidth=0.4)
#ax = sns.boxplot(data=plot_data, width=.3, palette=my_pal, boxprops={'zorder': 2}, ax=ax)

ax.set_box_aspect(1)
sns.despine()
ax.set(ylabel='pseU ratio')
ax.axhline(10, ls='--',  color = "#0B3954")

#add_median_labels(ax)
# category labels
plt.tick_params(axis='x', rotation=45)
plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
plt.savefig("Frequencies_significant.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.show()

#BOX PLOTS - frequency filter not for IVT
blood_1_box = Blood_prev_m6A_sig
blood_2_box = Blood_S5_m6A_sig
blood_3_box = Blood_S6_m6A_sig
ivt_1_box = Blood_IVT_m6A_sig
ivt_2_box = Blood_IVT_m6A_subsetted

plot_data = {
  'Blood 1': blood_1_box["Ratio"].tolist(),
  'Blood 2': blood_2_box["Ratio"].tolist(),
  'Blood 3': blood_3_box["Ratio"].tolist(),
  'Blood IVT filter': ivt_1_box["Ratio"].tolist(),
  'Blood IVT': ivt_2_box["Ratio"].tolist()
}

my_pal = {'Blood 1': '#1A87C7', 'Blood 2': '#1A87C7','Blood 3': '#1A87C7',
          'Glori 1': '#8B5CF6', 'Glori 2': '#8B5CF6',
          'Blood IVT': '#1A87C7', 'Blood IVT filter': '#1A87C7'}

plt.clf()
plt.figure(figsize=(2,4))
ax = sns.violinplot(data=plot_data, palette=my_pal, linewidth=0.4)
#ax = sns.boxplot(data=plot_data, width=.3, palette=my_pal, boxprops={'zorder': 2}, ax=ax)

ax.set_box_aspect(1)
sns.despine()
ax.set(ylabel='pseU ratio')
ax.axhline(10, ls='--',  color = "#0B3954")

#add_median_labels(ax)
# category labels
plt.tick_params(axis='x', rotation=45)
plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
plt.savefig("Frequencies_significant_IVTwithoutfreqfilter%_FINAL.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.show()

#####################################################################################################################

#### CORRELATION PLOTS ####

###### S5 vs S6

S5_vs_S6 = CorrPlot(A=BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_pseu.r1.mod.bed"),
              B= BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S6_DRS_basecall.0.7.2.GRCh38_pseu.r1.mod.bed"),
              A_name= "Blood 2", B_name="Blood 3", ratio=10, cov=10, outpath="FIG_5_S5_vs_S6_Corr_pseu.pdf")

S1_vs_S6 = CorrPlot(A=BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_blood.0.7.2.GRCh38_pseu.r1.mod.bed"),
              B= BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S6_DRS_basecall.0.7.2.GRCh38_pseu.r1.mod.bed"),
              A_name= "Blood 1", B_name="Blood 3", ratio=10, cov=10, outpath="FIG_5_S1_vs_S6_Corr_pseu.pdf")

S1_vs_S5 = CorrPlot(A=BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_blood.0.7.2.GRCh38_pseu.r1.mod.bed"),
              B= BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_pseu.r1.mod.bed"),
              A_name= "Blood 1", B_name="Blood 2", ratio=10, cov=10, outpath="FIG_5_S1_vs_S5_Corr_pseu.pdf")
