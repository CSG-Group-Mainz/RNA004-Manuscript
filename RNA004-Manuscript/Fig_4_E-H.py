import upsetplot
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.offsetbox import AnchoredText
import functools as ft

#Function to prepare bedfile as df 
def BedToDf(bedfile:str):
    #read bedfile as pd.csv
    df = pd.read_csv(bedfile, sep='\t', comment='t', header=None)
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

def CorrPlot(A:str, B:str, A_name:str, B_name:str, ratio:int, cov:int, outpath:str):
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
    
    h = sns.jointplot(x=mergedata["Ratio_A"].tolist(), y=mergedata["Ratio_B"].tolist(), kind="hex", height=3, ratio=5)
    h.set_axis_labels(A_name, B_name)
    anc = AnchoredText(subtitle, loc="upper left", frameon=False,prop=dict(fontsize=6.5, ))
    h.ax_joint.add_artist(anc)

    h.savefig(outpath)
    
    return None

######### Upset Plot Main Figure 4 D
Blood_m6A_prev = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_blood.0.7.2.GRCh38_m6A.r1.mod.bed")
Blood_S5_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed")
Blood_S6_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S6_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed")
Blood_Glori_1 = GloriToDf("/home/johannes/RNA004_Revision_1/new_data/AS-1466472-LR-78969_R1.totalm6A.FDR.csv")
Blood_Glori_2 = GloriToDf("/home/johannes/RNA004_Revision_1/new_data/AS-1466474-LR-78969_R1.totalm6A.FDR.csv")
Blood_S6_IVT_merged_m6A = BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/merged_RNA004_S6_IVT_m6A.r1.mod.bed")

##### Filter Nanopore datasets for a minimum coverage of 10
Blood_m6A_prev_cov10 = Blood_m6A_prev.loc[Blood_m6A_prev["Coverage"] >= 10]
Blood_S5_m6A_cov10 = Blood_S5_m6A.loc[Blood_S5_m6A["Coverage"] >= 10]
Blood_S6_m6A_cov10 = Blood_S6_m6A.loc[Blood_S6_m6A["Coverage"] >= 10]
Blood_S6_IVT_merged_m6A_cov10 = Blood_S6_IVT_merged_m6A.loc[Blood_S6_IVT_merged_m6A["Coverage"] >= 10]

#Get sites which are covered in all Blood samples
positions = [Blood_S5_m6A_cov10["Position"].tolist(), Blood_S6_m6A_cov10["Position"].tolist(),
             Blood_S6_IVT_merged_m6A_cov10["Position"].tolist()]


shared_positions = set.intersection(*map(set,positions))

##Subset All dataset only to those joint covered sites in blood samples, include Glori positions for this
Blood_S5_m6A_subsetted = Blood_S5_m6A_cov10.loc[Blood_S5_m6A_cov10['Position'].isin(shared_positions)]
Blood_S6_m6A_subsetted = Blood_S6_m6A_cov10.loc[Blood_S6_m6A_cov10['Position'].isin(shared_positions)]
Blood_Glori_1_subsetted = Blood_Glori_1.loc[Blood_Glori_1['Position'].isin(shared_positions)]
Blood_Glori_2_subsetted = Blood_Glori_2.loc[Blood_Glori_2['Position'].isin(shared_positions)]
Blood_IVT_m6A_subsetted = Blood_S6_IVT_merged_m6A_cov10.loc[Blood_S6_IVT_merged_m6A_cov10['Position'].isin(shared_positions)]

#Only retain Positions with modification frequency over 10 Percent 
Blood_S6_m6A_sig = Blood_S6_m6A_subsetted.loc[Blood_S6_m6A_subsetted["Ratio"] >= 10]
Blood_S5_m6A_sig = Blood_S5_m6A_subsetted.loc[Blood_S5_m6A_subsetted["Ratio"] >= 10]
Blood_Glori_1_sig = Blood_Glori_1_subsetted.loc[Blood_Glori_1_subsetted["Ratio"] >= 10]
Blood_Glori_2_sig = Blood_Glori_2_subsetted.loc[Blood_Glori_2_subsetted["Ratio"] >= 10]
Blood_IVT_m6A_sig = Blood_IVT_m6A_subsetted.loc[Blood_IVT_m6A_subsetted["Ratio"] >= 10]

#Merge Glori and Blood Replicates together for the upset plot
Blood_m6A_merged = set(Blood_S5_m6A_sig["Position"].tolist()).intersection(set(Blood_S6_m6A_sig["Position"].tolist()))
Glori_m6A_merged = set(Blood_Glori_1_sig["Position"].to_list()).intersection(set(Blood_Glori_2_sig["Position"].to_list()))


#Prepare for Plot
plotdata = upsetplot.from_contents({"Blood": Blood_m6A_merged,
                                    "Glori": Glori_m6A_merged,
                                    "Blood IVT": Blood_IVT_m6A_sig["Position"].to_list()})


plot_upset = upsetplot.UpSet(plotdata, show_counts=True)

plot_upset.style_subsets(present="Blood", facecolor="#0B3954", edgecolor="black")
plot_upset.style_subsets(present="Glori", facecolor='#8B5CF6', edgecolor="black")
plot_upset.style_subsets(present="Blood IVT", facecolor='#1A87C7', edgecolor="black")
plot_upset.style_subsets(present=["Blood", "Glori"], facecolor="#098C9A", edgecolor="black", hatch="")
plot_upset.style_subsets(present=["Blood", "Glori", "Blood IVT"], facecolor="grey", edgecolor="black", hatch="")
plot_upset.style_subsets(present=["Blood", "Blood IVT"], facecolor="grey", edgecolor="black", hatch="")
plot_upset.style_subsets(present=["Glori", "Blood IVT"], facecolor="grey", edgecolor="black", hatch="")

plot_upset.style_categories(
"Blood", bar_facecolor="#0B3954", bar_edgecolor="black"
)

plot_upset.style_categories(
"Glori", bar_facecolor='#8B5CF6', bar_edgecolor="black"
)

plot_upset.plot()
plt.suptitle("")
plt.savefig("Fig_4_D.pdf", format="pdf")
plt.show()

#####################################################################################################################

######### BOXPLOT ##########

#Get sites which are covered in ALL Blood samples including GLORI Sites 
positions = [Blood_m6A_prev_cov10["Position"].to_list(), Blood_S5_m6A_cov10["Position"].tolist(), Blood_S6_m6A_cov10["Position"].tolist(),
             Blood_S6_IVT_merged_m6A_cov10["Position"].tolist(), Blood_Glori_1["Position"].to_list(), Blood_Glori_2["Position"].to_list()]

shared_positions = set.intersection(*map(set,positions))
print(len(shared_positions))

##Subset All dataset only to those joint covered sites
Blood_prev_m6A_subsetted = Blood_m6A_prev_cov10.loc[Blood_m6A_prev_cov10['Position'].isin(shared_positions)]
Blood_S5_m6A_subsetted = Blood_S5_m6A_cov10.loc[Blood_S5_m6A_cov10['Position'].isin(shared_positions)]
Blood_S6_m6A_subsetted = Blood_S6_m6A_cov10.loc[Blood_S6_m6A_cov10['Position'].isin(shared_positions)]
Blood_Glori_1_subsetted = Blood_Glori_1.loc[Blood_Glori_1['Position'].isin(shared_positions)]
Blood_Glori_2_subsetted = Blood_Glori_2.loc[Blood_Glori_2['Position'].isin(shared_positions)]
Blood_IVT_m6A_subsetted = Blood_S6_IVT_merged_m6A_cov10.loc[Blood_S6_IVT_merged_m6A_cov10['Position'].isin(shared_positions)]

#Assign frequencies for plotting 
blood_1_box = Blood_prev_m6A_subsetted
blood_2_box = Blood_S5_m6A_subsetted
blood_3_box = Blood_S6_m6A_subsetted
ivt_1_box = Blood_IVT_m6A_subsetted
glori_1_box = Blood_Glori_1_subsetted
glori_2_box = Blood_Glori_2_subsetted

#dictionary for plotting
plot_data = {
  'Glori 1': glori_1_box["Ratio"].tolist(),
  'Glori 2': glori_2_box["Ratio"].tolist(),
  'Blood 1': blood_1_box["Ratio"].tolist(),
  'Blood 2': blood_2_box["Ratio"].tolist(),
  'Blood 3': blood_3_box["Ratio"].tolist(),
  'Blood IVT': ivt_1_box["Ratio"].tolist(),
}

#color palette for plotting
my_pal = {'Blood 1': '#1A87C7', 'Blood 2': '#1A87C7','Blood 3': '#1A87C7',
          'Glori 1': '#8B5CF6', 'Glori 2': '#8B5CF6',
          'Blood IVT': '#1A87C7'}

#Plot
plt.clf()
plt.figure(figsize=(2,4))
ax = sns.violinplot(data=plot_data, palette=my_pal, linewidth=0.4)
ax.set_box_aspect(1)
sns.despine()
ax.set(ylabel='methylation ratio')
ax.axhline(10, ls='--',  color = "#0B3954")
plt.tick_params(axis='x', rotation=45)
plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
plt.savefig("Fig_4_E.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.show()


#### CORRELATION PLOTS ####

GLORI_1_vs_S5 = CorrPlot(A=GloriToDf("/home/johannes/RNA004_Revision_1/new_data/AS-1466472-LR-78969_R1.totalm6A.FDR.csv"),
              B= BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed"),
              A_name= "Glori 1", B_name="Blood 2", ratio=10, cov=10, outpath="FIG_4_F.pdf")

S5_vs_S6 = CorrPlot(A=BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed"),
              B= BedToDf("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S6_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed"),
              A_name= "Blood 2", B_name="Blood 3", ratio=10, cov=10, outpath="FIG_4_G.pdf")


GLORI_CORR = CorrPlot(A=GloriToDf("/home/johannes/RNA004_Revision_1/new_data/AS-1466472-LR-78969_R1.totalm6A.FDR.csv"),
              B= GloriToDf("/home/johannes/RNA004_Revision_1/new_data/AS-1466474-LR-78969_R1.totalm6A.FDR.csv"),
              A_name= "Glori 1", B_name="Glori 2", ratio=10, cov=10, outpath="FIG_4_H.pdf")





























# #####################################################################################################################

# #Boxplot of frequencies of positions which are 1) exclusively found significant in GLORI and 2) exclusively sig. in Blood samples
# Blood_1_sig = Blood_prev_m6A_subsetted.loc[Blood_prev_m6A_subsetted["Ratio"] >= 10]
# Blood_2_sig = Blood_S5_m6A_subsetted.loc[Blood_S5_m6A_subsetted["Ratio"] >= 10]
# Blood_3_sig = Blood_S6_m6A_subsetted.loc[Blood_S6_m6A_subsetted["Ratio"] >= 10]

# #Get List of all significant Blood samples in all samples
# All_blood_significant = Blood_1_sig["Position"].to_list() + Blood_2_sig["Position"].to_list() + Blood_3_sig["Position"].to_list()
# All_Glori_significant = Blood_Glori_1_subsetted["Position"].to_list() + Blood_Glori_2_subsetted["Position"].to_list()

# #Get highly exclusive sites for glori
# Glori_1_excl = set(Blood_Glori_1_subsetted["Position"].to_list()) - set(All_blood_significant)
# Glori_2_excl = set(Blood_Glori_2_subsetted["Position"].to_list()) - set(All_blood_significant)

# #Get highly exclusive site for the respective Blood Samples
# print(len(All_Glori_significant))
# print(Blood_1_sig)
# Blood_1_excl = set(Blood_1_sig["Position"].to_list()) - set(All_Glori_significant)
# print(set(Blood_1_sig["Position"].to_list()).issubset(set(All_Glori_significant)))
# Blood_2_excl = set(Blood_2_sig["Position"].to_list()) - set(All_Glori_significant)
# Blood_3_excl = set(Blood_3_sig["Position"].to_list()) - set(All_Glori_significant)

# #Subset dataframes accordingly

# Glori_1_excl = Blood_Glori_1_subsetted.loc[Blood_Glori_1_subsetted['Position'].isin(Glori_1_excl)]
# Glori_2_excl = Blood_Glori_2_subsetted.loc[Blood_Glori_2_subsetted['Position'].isin(Glori_2_excl)]
# Blood_1_excl = Blood_prev_m6A_subsetted.loc[Blood_prev_m6A_subsetted['Position'].isin(Blood_1_excl)]
# Blood_2_excl = Blood_S5_m6A_subsetted.loc[Blood_S5_m6A_subsetted['Position'].isin(Blood_2_excl)]
# Blood_3_excl = Blood_S6_m6A_subsetted.loc[Blood_S6_m6A_subsetted['Position'].isin(Blood_3_excl)]



# plot_data = {
#   'Glori 1': Glori_1_excl["Ratio"].tolist(),
#   'Glori 2': Glori_2_excl["Ratio"].tolist(),
#   'Blood 1': Blood_1_excl["Ratio"].tolist(),
#   'Blood 2': Blood_2_excl["Ratio"].tolist(),
#   'Blood 3': Blood_3_excl["Ratio"].tolist()
# }

# my_pal = {'Blood 1': '#1A87C7', 'Blood 2': '#1A87C7','Blood 3': '#1A87C7',
#           'Glori 1': '#8B5CF6', 'Glori 2': '#8B5CF6',
#           'Blood IVT': '#1A87C7'}

# plt.clf()
# # almost verbatim from question
# plt.figure(figsize=(2,4))
# ax = sns.violinplot(data=plot_data, palette=my_pal, linewidth=0.4)
# #ax = sns.boxplot(data=plot_data, width=.3, palette=my_pal, boxprops={'zorder': 2}, ax=ax)

# ax.set_box_aspect(1)
# sns.despine()
# ax.set(ylabel='methylation ratio')
# ax.axhline(10, ls='--',  color = "#0B3954")


# plt.tick_params(axis='x', rotation=45)
# plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
# plt.savefig("Frequencies_EXCLUSIVE_m6A.pdf", format="pdf", bbox_inches="tight", dpi=600)
# plt.show()

# ##############################################################################################






###### S5 VS Glori 1


# ###### Delta Tables

# #HEK A and HEK B

# HEK_A_B_delta = DeltaTable(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"), B=BedToDf("data/RNA004_B_basecall.0.7.2.GRCh38_all.mod.chr.bed"), 
#                         A_name="HEK293T_A", B_name="HEK293T_B", ratio=0, coverage=0, outpath="HEK_A_B_Delta_m6a.csv")


# HEK_A_B_delta_filtered = DeltaTable(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"), B=BedToDf("data/RNA004_B_basecall.0.7.2.GRCh38_all.mod.chr.bed"), 
#                         A_name="HEK293T_A", B_name="HEK293T_B", ratio=10, coverage=10, outpath="HEK_A_B_Delta_m6a_filtered.csv")


# #### BLOOD VS HEK A


# HEK_A_Blood_S5_filtered = DeltaTable(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             B= BedToDf("data/RNA004_S5_DRS_basecall.0.7.2.GRCh38_all.m6A.mod.chr.bed"),
#             A_name= "HEK293T_A", B_name="Blood_S5", ratio=10, coverage=10, outpath="HEK_A_vs_Blood_S5_Delta_m6A_filtered.csv")


# HEK_A_Blood_S5_unfiltered = DeltaTable(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             B= BedToDf("data/RNA004_S5_DRS_basecall.0.7.2.GRCh38_all.m6A.mod.chr.bed"),
#             A_name= "HEK293T_A", B_name="Blood_S5", ratio=0, coverage=0, outpath="HEK_A_vs_Blood_S5_Delta_m6A_unfiltered.csv")

# ####### CORRELATION PLOTS 

# HEK_A_vs_B_m6A_filtered = CorrPlot(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             B= BedToDf("data/RNA004_B_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             A_name= "HEK293T_A", B_name="HEK293T_B", ratio=5, cov=20, outpath="HEK_A_vs_B_m6A_filtered.pdf")

# HEK_A_vs_B_m6A_unfiltered = CorrPlot(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             B= BedToDf("data/RNA004_B_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             A_name= "HEK293T_A", B_name="HEK293T_B", ratio=0, cov=0, outpath="HEK_A_vs_B_m6A_unfiltered.pdf")

# #### BLOOD VS HEK A

# HEK_A_Blood_S5_filtered = CorrPlot(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             B= BedToDf("data/RNA004_S5_DRS_basecall.0.7.2.GRCh38_all.m6A.mod.chr.bed"),
#             A_name= "HEK293T_A", B_name="Blood_S5", ratio=10, cov=10, outpath="HEK_A_vs_Blood_S5_m6A_filtered.pdf")


# HEK_A_Blood_S5_unfiltered = CorrPlot(A=BedToDf("data/RNA004_A_basecall.0.7.2.GRCh38_all.mod.chr.bed"),
#             B= BedToDf("data/RNA004_S5_DRS_basecall.0.7.2.GRCh38_all.m6A.mod.chr.bed"),
#             A_name= "HEK293T_A", B_name="Blood_S5", ratio=0, cov=0, outpath="HEK_A_vs_Blood_S5_m6A_unfiltered.pdf")

# ##UpSetPlot before and after subsampling, with filter, m6A PREVIOUS OLD BLOOD SAMPLE

# #Subsampled file is already filtered
# m6A_subsampled_old_blood = BedToDf("cov_sub_RNA004_blood.0.7.2.GRCh38_m6A.cov_sub_r1.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# m6A_subsampled_old_blood = m6A_subsampled_old_blood.loc[m6A_subsampled_old_blood["Coverage"] >= 10]
# m6A_subsampled_filter_old_blood = m6A_subsampled_old_blood.loc[m6A_subsampled_old_blood["Ratio"] >= 10]

# #Non-subsampled file needs to be filtered
# m6A_old_blood = BedToDf("data/RNA004_blood.0.7.2.GRCh38_all.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# m6A_old_blood = m6A_old_blood.loc[m6A_old_blood["Coverage"] >= 10]
# m6A_filter_old_blood = m6A_old_blood.loc[m6A_old_blood["Ratio"] >= 10]

# #UpsetPlot

# #Prepare for Plot
# plotdata = upsetplot.from_contents({"Blood_1_subsampled": m6A_subsampled_filter_old_blood["Position"].tolist(), "Blood_1": m6A_filter_old_blood["Position"].to_list()})
# print(plotdata)

# plot_upset = upsetplot.UpSet(plotdata, show_counts=True, show_percentages=True)

# plot_upset.style_subsets(present="Blood_1_subsampled", facecolor="#0B3954", edgecolor="black")
# plot_upset.style_subsets(present="Blood_1", facecolor='#1A87C7', edgecolor="black")
# plot_upset.style_subsets(present=["Blood_1_subsampled", "Blood_1"], facecolor="gray", edgecolor="black", hatch="")

# plot_upset.style_categories(
# "Blood_1_subsampled", bar_facecolor="#0B3954", bar_edgecolor="black"
# )

# plot_upset.style_categories(
# "Blood_1", bar_facecolor='#1A87C7', bar_edgecolor="black" 
# )

# plot_upset.plot()
# plt.suptitle("")
# plt.savefig("Blood_1_subsample_intersection.pdf", format="pdf")
# plt.show()




# ##UpSetPlot before and after subsampling, with filter, m6A NEW BLOOD SAMPLE

# #Subsampled file is already filtered
# m6A_subsampled_new_blood = BedToDf("cov_sub_RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.cov_sub_r1.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# m6A_subsampled_new_blood = m6A_subsampled_new_blood.loc[m6A_subsampled_new_blood["Coverage"] >= 10]
# m6A_subsampled_filter_new_blood = m6A_subsampled_new_blood.loc[m6A_subsampled_new_blood["Ratio"] >= 10]

# #Non-subsampled file needs to be filtered
# m6A_new_blood = BedToDf("data/RNA004_S5_DRS_basecall.0.7.2.GRCh38_all.m6A.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# m6A_new_blood = m6A_new_blood.loc[m6A_new_blood["Coverage"] >= 10]
# m6A_filter_new_blood = m6A_new_blood.loc[m6A_new_blood["Ratio"] >= 10]

# #UpsetPlot

# #Prepare for Plot
# plotdata = upsetplot.from_contents({"Blood_S5_subsampled": m6A_subsampled_filter_new_blood["Position"].tolist(), "Blood_S5": m6A_filter_new_blood["Position"].to_list()})
# print(plotdata)

# plot_upset = upsetplot.UpSet(plotdata, show_counts=True, show_percentages=True)

# plot_upset.style_subsets(present="Blood_S5_subsampled", facecolor="#0B3954", edgecolor="black")
# plot_upset.style_subsets(present="Blood_S5", facecolor='#1A87C7', edgecolor="black")
# plot_upset.style_subsets(present=["Blood_S5_subsampled", "Blood_S5"], facecolor="gray", edgecolor="black", hatch="")

# plot_upset.style_categories(
# "Blood_S5_subsampled", bar_facecolor="#0B3954", bar_edgecolor="black"
# )

# plot_upset.style_categories(
# "Blood_S5", bar_facecolor='#1A87C7', bar_edgecolor="black"
# )

# plot_upset.plot()
# plt.suptitle("")
# plt.savefig("Blood_S5_subsample_intersection.pdf", format="pdf")
# plt.show()


# #PSEUDOURIDINE


# ##UpSetPlot before and after subsampling, with filter, m6A PREVIOUS OLD BLOOD SAMPLE

# #Subsampled file is already filtered
# pseu_subsampled_old_blood = BedToDf("cov_sub_RNA004_blood.0.7.2.GRCh38_pseu.cov_sub_r1.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# pseu_subsampled_old_blood = pseu_subsampled_old_blood.loc[pseu_subsampled_old_blood["Coverage"] >= 10]
# pseu_subsampled_filter_old_blood = pseu_subsampled_old_blood.loc[pseu_subsampled_old_blood["Ratio"] >= 10]

# #Non-subsampled file needs to be filtered
# pseu_old_blood = BedToDf("data/RNA004_blood.0.7.2.GRCh38_all.pseu.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# pseu_old_blood = pseu_old_blood.loc[pseu_old_blood["Coverage"] >= 10]
# pseu_filter_old_blood = pseu_old_blood.loc[pseu_old_blood["Ratio"] >= 10]

# #UpsetPlot

# #Prepare for Plot
# plotdata = upsetplot.from_contents({"Blood_1_subsampled": pseu_subsampled_filter_old_blood["Position"].tolist(), "Blood_1": pseu_filter_old_blood["Position"].to_list()})
# print(plotdata)

# plot_upset = upsetplot.UpSet(plotdata, show_counts=True, show_percentages=True)

# plot_upset.style_subsets(present="Blood_1_subsampled", facecolor="#0B3954", edgecolor="black")
# plot_upset.style_subsets(present="Blood_1", facecolor='#1A87C7', edgecolor="black")
# plot_upset.style_subsets(present=["Blood_1_subsampled", "Blood_1"], facecolor="gray", edgecolor="black", hatch="")

# plot_upset.style_categories(
# "Blood_1_subsampled", bar_facecolor="#0B3954", bar_edgecolor="black"
# )

# plot_upset.style_categories(
# "Blood_1", bar_facecolor='#1A87C7', bar_edgecolor="black"
# )

# plot_upset.plot()
# plt.suptitle("")
# plt.savefig("Blood_1_subsample_intersection_PSEU.pdf", format="pdf")
# plt.show()


# ##UpSetPlot before and after subsampling, with filter, m6A NEW BLOOD SAMPLE

# #Subsampled file is already filtered
# pseu_subsampled_new_blood = BedToDf("cov_sub_RNA004_S5_DRS_basecall.0.7.2.GRCh38_pseu.cov_sub_r1.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# pseu_subsampled_new_blood = pseu_subsampled_new_blood.loc[pseu_subsampled_new_blood["Coverage"] >= 10]
# pseu_subsampled_filter_new_blood = pseu_subsampled_new_blood.loc[pseu_subsampled_new_blood["Ratio"] >= 10]

# #Non-subsampled file needs to be filtered
# pseu_new_blood = BedToDf("RNA004_S5_DRS_basecall.0.7.2.GRCh38_all.pseu.mod.chr.bed")

# #Filter by Coverage and Ratio >= 10
# pseu_new_blood = pseu_new_blood.loc[pseu_new_blood["Coverage"] >= 10]
# pseu_filter_new_blood = pseu_new_blood.loc[pseu_new_blood["Ratio"] >= 10]

# #UpsetPlot

# #Prepare for Plot
# plotdata = upsetplot.from_contents({"Blood_S5_subsampled": pseu_subsampled_filter_new_blood["Position"].tolist(), "Blood_S5": pseu_filter_new_blood["Position"].to_list()})
# print(plotdata)

# plot_upset = upsetplot.UpSet(plotdata, show_counts=True, show_percentages=True)

# plot_upset.style_subsets(present="Blood_S5_subsampled", facecolor="#0B3954", edgecolor="black")
# plot_upset.style_subsets(present="Blood_S5", facecolor='#1A87C7', edgecolor="black")
# plot_upset.style_subsets(present=["Blood_S5_subsampled", "Blood_S5"], facecolor="gray", edgecolor="black", hatch="")

# plot_upset.style_categories(
# "Blood_S5_subsampled", bar_facecolor="#0B3954", bar_edgecolor="black"
# )

# plot_upset.style_categories(
# "Blood_S5", bar_facecolor='#1A87C7', bar_edgecolor="black"
# )

# plot_upset.plot()
# plt.suptitle("")
# plt.savefig("Blood_pseu_S5_subsample_intersection.pdf", format="pdf")
# plt.show()
