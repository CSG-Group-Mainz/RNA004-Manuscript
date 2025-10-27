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


def BedToDf_2(bedfile:str):
    #read bedfile as pd.csv
    df = pd.read_csv(bedfile, sep='\t', comment='t', header=None)
    #Merge chromosome, start position and strand info into one string
    df_merged = df.iloc[:,0].astype(str)+":"+df.iloc[:,2].astype(str)
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


def GloriToDf_2(Glorifile:str):
    #Read Glori Sites
    df = pd.read_csv(Glorifile, delimiter='\t')
    #### REFORMAT #### All base positions must be shifted 1 to the left to compare with bedfiles
    df_1_base = list(df.iloc[:,1].astype(int))
    #Create concatenations
    pos = df.iloc[:,0].astype(str) + ":" + np.asarray(df_1_base).astype(str) 
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

#Only retain Positions with modification frequency over 10 Percent which are UNJOINT
Blood_1_sig_unjoint = Blood_m6A_prev_cov10.loc[Blood_m6A_prev_cov10["Ratio"] >= 10]
Blood_S6_m6A_sig_unjoint = Blood_S6_m6A_cov10.loc[Blood_S6_m6A_cov10["Ratio"] >= 10]
Blood_S5_m6A_sig_unjoint  = Blood_S5_m6A_cov10.loc[Blood_S5_m6A_cov10["Ratio"] >= 10]
Blood_Glori_1_sig_unjoint  = Blood_Glori_1.loc[Blood_Glori_1["Ratio"] >= 10]
Blood_Glori_2_sig_unjoint  = Blood_Glori_2.loc[Blood_Glori_2["Ratio"] >= 10]
Blood_IVT_m6A_sig_unjoint  = Blood_S6_IVT_merged_m6A_cov10.loc[Blood_S6_IVT_merged_m6A_cov10["Ratio"] >= 10]

#Merge Glori and Blood Replicates together for the upset plot
Blood_m6A_merged = set(Blood_S5_m6A_sig_unjoint["Position"].tolist()).intersection(set(Blood_S6_m6A_sig_unjoint["Position"].tolist())).intersection(set(Blood_1_sig_unjoint["Position"].tolist()))
Glori_m6A_merged = set(Blood_Glori_1_sig_unjoint["Position"].to_list()).intersection(set(Blood_Glori_2_sig_unjoint["Position"].to_list()))

#Prepare for Plot
plotdata = upsetplot.from_contents({"Blood dRNA": Blood_m6A_merged,
                                    "Blood GLORI": Glori_m6A_merged,
                                    "Blood IVT dRNA": Blood_IVT_m6A_sig_unjoint["Position"].to_list()})


plot_upset = upsetplot.UpSet(plotdata, show_counts=True)

plot_upset.style_subsets(present="Blood dRNA", facecolor="#0B3954", edgecolor="black")
plot_upset.style_subsets(present="Blood GLORI", facecolor='#8B5CF6', edgecolor="black")
plot_upset.style_subsets(present="Blood IVT dRNA", facecolor='#1A87C7', edgecolor="black")
plot_upset.style_subsets(present=["Blood dRNA", "Blood GLORI"], facecolor="#098C9A", edgecolor="black", hatch="")
plot_upset.style_subsets(present=["Blood dRNA", "Blood GLORI", "Blood IVT dRNA"], facecolor="grey", edgecolor="black", hatch="")
plot_upset.style_subsets(present=["Blood dRNA", "Blood IVT dRNA"], facecolor="grey", edgecolor="black", hatch="")
plot_upset.style_subsets(present=["Blood GLORI", "Blood IVT dRNA"], facecolor="grey", edgecolor="black", hatch="")

plot_upset.style_categories(
"Blood dRNA", bar_facecolor="#0B3954", bar_edgecolor="black"
)

plot_upset.style_categories(
"Blood GLORI", bar_facecolor='#8B5CF6', bar_edgecolor="black"
)

plot_upset.plot()
plt.suptitle("")
plt.savefig("Fig_4_D.pdf", format="pdf")
plt.show()

#####################################################################################################################

#color palette for plotting


my_pal = {'Blood 1': "#0B3954", 'Blood 2': "#177BB5",'Blood 3': "#38A6E5",
          'GLORI 1': '#8B5CF6', 'GLORI 2': '#8B5CF6',
          'Blood IVT': '#8B5CF6', 'Blood IVT filter': '#1A87C7'}



#Assign frequencies for plotting 
blood_1_box = Blood_m6A_prev_cov10.loc[Blood_m6A_prev_cov10["Ratio"] >= 10]
blood_2_box = Blood_S5_m6A_cov10.loc[Blood_S5_m6A_cov10["Ratio"] >= 10] 
blood_3_box = Blood_S6_m6A_cov10.loc[Blood_S6_m6A_cov10["Ratio"] >= 10]
ivt_1_box = Blood_S6_IVT_merged_m6A_cov10.loc[Blood_S6_IVT_merged_m6A_cov10["Ratio"] >= 10]
glori_1_box = Blood_Glori_1
glori_2_box = Blood_Glori_1

#dictionary for plotting
plot_data = {
  'GLORI 1': glori_1_box["Ratio"].tolist(),
  'GLORI 2': glori_2_box["Ratio"].tolist(),
  'Blood 1': blood_1_box["Ratio"].tolist(),
  'Blood 2': blood_2_box["Ratio"].tolist(),
  'Blood 3': blood_3_box["Ratio"].tolist(),
  'Blood IVT': ivt_1_box["Ratio"].tolist(),
}

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
plt.savefig("Fig_4_E_MAIN.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.show()


######### READ BEDFILES FOR FREQUENCY PLOTS
Blood_m6A_prev = BedToDf_2("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_blood.0.7.2.GRCh38_m6A.r1.mod.bed")
Blood_S5_m6A = BedToDf_2("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed")
Blood_S6_m6A = BedToDf_2("/home/johannes/RNA004_Revision_1/new_data/chr/RNA004_S6_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed")
Blood_Glori_1 = GloriToDf_2("/home/johannes/RNA004_Revision_1/new_data/AS-1466472-LR-78969_R1.totalm6A.FDR.csv")
Blood_Glori_2 = GloriToDf_2("/home/johannes/RNA004_Revision_1/new_data/AS-1466474-LR-78969_R1.totalm6A.FDR.csv")
Blood_S6_IVT_merged_m6A = BedToDf_2("/home/johannes/RNA004_Revision_1/new_data/chr/merged_RNA004_S6_IVT_m6A.r1.mod.bed")

##### Filter Nanopore datasets for a minimum coverage of 10
Blood_m6A_prev_cov10 = Blood_m6A_prev.loc[Blood_m6A_prev["Coverage"] >= 10]
Blood_S5_m6A_cov10 = Blood_S5_m6A.loc[Blood_S5_m6A["Coverage"] >= 10]
Blood_S6_m6A_cov10 = Blood_S6_m6A.loc[Blood_S6_m6A["Coverage"] >= 10]
Blood_S6_IVT_merged_m6A_cov10 = Blood_S6_IVT_merged_m6A.loc[Blood_S6_IVT_merged_m6A["Coverage"] >= 10]



######### BOXPLOT ##########
Blood_1_path = '/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/DORADO_0_7_2_folder/RNA002_blood.0.7.2.GRCh38.bam'
Blood_2_path = '/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S5_analysis_5_12_24/RNA004_S5_DRS_basecall.0.7.2.GRCh38.bam'
Blood_3_path = '/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S6_IVT_analysis/RNA004_S6_DRS_basecall.0.7.2.GRCh38.bam'
Blood_IVT_1_path = '/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S6_IVT_analysis/RNA004_S6_IVT_basecall.0.7.2.GRCh38.bam'
Blood_IVT_2_path = '/raid/chhewel_analysis/RNA004_Manuscript_reanalysis/S6_IVT_analysis/RNA004_S6_IVT_LI_basecall.0.7.2.GRCh38.bam'


Glori_1_m6a_sites = pd.read_table("/home/johannes/RNA004_Revision_1/BED_DEPTH/PAPER_GLORI_blood_1.SAMPLES_DEPTH.txt")
Glori_2_m6a_sites = pd.read_table("/home/johannes/RNA004_Revision_1/BED_DEPTH/PAPER_GLORI_blood_2.SAMPLES_DEPTH.txt")
Blood_1_m6a_sites = pd.read_table("/home/johannes/RNA004_Revision_1/BED_DEPTH/RNA004_blood.0.7.2.GRCh38_m6A.r1.mod.SAMPLES_DEPTH.txt")
Blood_2_m6a_sites = pd.read_table("/home/johannes/RNA004_Revision_1/BED_DEPTH/RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.SAMPLES_DEPTH.txt")
Blood_3_m6a_sites = pd.read_table("/home/johannes/RNA004_Revision_1/BED_DEPTH/RNA004_S6_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.SAMPLES_DEPTH.txt")
Blood_IVT_m6a_sites = pd.read_table("/home/johannes/RNA004_Revision_1/BED_DEPTH/merged_RNA004_S6_IVT_m6A.r1.mod.SAMPLES_DEPTH.txt")

# Subset for joint covered sites 
#GLORI 1 SUBSETTIING
Glori_1_m6a_sites_joint = Glori_1_m6a_sites[(Glori_1_m6a_sites[Blood_1_path] >= 10) & (Glori_1_m6a_sites[Blood_2_path] >= 10) & (Glori_1_m6a_sites[Blood_3_path] >= 10) & 
                                            (Glori_1_m6a_sites[Blood_IVT_1_path] + Glori_1_m6a_sites[Blood_IVT_2_path] >= 10)]

Glori_1_m6a_sites_joint_concat = Glori_1_m6a_sites_joint["#CHROM"].astype(str) + ":" + Glori_1_m6a_sites_joint["POS"].astype(str)

#GLORI 2 SUBSETTING
Glori_2_m6a_sites_joint = Glori_2_m6a_sites[(Glori_2_m6a_sites[Blood_1_path] >= 10) & (Glori_2_m6a_sites[Blood_2_path] >= 10) & (Glori_2_m6a_sites[Blood_3_path] >= 10) &
                                            (Glori_2_m6a_sites[Blood_IVT_1_path] + Glori_2_m6a_sites[Blood_IVT_2_path] >= 10)]

Glori_2_m6a_sites_joint_concat = Glori_2_m6a_sites_joint["#CHROM"].astype(str) + ":" + Glori_2_m6a_sites_joint["POS"].astype(str)

#BLOOD 1 SUBSETTING 

Blood_1_m6a_sites_joint = Blood_1_m6a_sites[(Blood_1_m6a_sites[Blood_2_path] >= 10) & (Blood_1_m6a_sites[Blood_3_path] >= 10) &
                                            (Blood_1_m6a_sites[Blood_IVT_1_path] + Blood_1_m6a_sites[Blood_IVT_2_path] >= 10)]

Blood_1_m6a_sites_joint_concat = Blood_1_m6a_sites_joint["#CHROM"].astype(str) + ":" + Blood_1_m6a_sites_joint["POS"].astype(str)


#BLOOD 2 SUBSETTING 

Blood_2_m6a_sites_joint = Blood_2_m6a_sites[(Blood_2_m6a_sites[Blood_1_path] >= 10) & (Blood_2_m6a_sites[Blood_3_path] >= 10) &
                                            (Blood_2_m6a_sites[Blood_IVT_1_path] + Blood_2_m6a_sites[Blood_IVT_2_path] >= 10)]

Blood_2_m6a_sites_joint_concat = Blood_2_m6a_sites_joint["#CHROM"].astype(str) + ":" + Blood_2_m6a_sites_joint["POS"].astype(str)

#BLOOD 3 SUBSETTING 

Blood_3_m6a_sites_joint = Blood_3_m6a_sites[(Blood_3_m6a_sites[Blood_1_path] >= 10) & (Blood_3_m6a_sites[Blood_2_path] >= 10) &
                                            (Blood_3_m6a_sites[Blood_IVT_1_path] + Blood_3_m6a_sites[Blood_IVT_2_path] >= 10)]

Blood_3_m6a_sites_joint_concat = Blood_3_m6a_sites_joint["#CHROM"].astype(str) + ":" + Blood_3_m6a_sites_joint["POS"].astype(str)


#BLOOD IVT SUBSETTING 

Blood_IVT_m6a_sites_joint = Blood_IVT_m6a_sites[(Blood_IVT_m6a_sites[Blood_1_path] >= 10) & (Blood_IVT_m6a_sites[Blood_2_path] >= 10) &
                                                (Blood_IVT_m6a_sites[Blood_3_path] >= 10)]

Blood_IVT_m6a_sites_joint_concat =Blood_IVT_m6a_sites_joint["#CHROM"].astype(str) + ":" + Blood_IVT_m6a_sites_joint["POS"].astype(str)


##########
# Violinplot of joint covered m6A sites (Supplementary Figure)

################################################################################################
Blood_1_violin_sig = Blood_m6A_prev_cov10[(Blood_m6A_prev_cov10["Position"].isin(Blood_1_m6a_sites_joint_concat)) & (Blood_m6A_prev_cov10["Ratio"] >= 10)]
Blood_2_violin_sig = Blood_S5_m6A_cov10[(Blood_S5_m6A_cov10["Position"].isin(Blood_2_m6a_sites_joint_concat)) & (Blood_S5_m6A_cov10["Ratio"] >= 10)]
Blood_3_violin_sig = Blood_S6_m6A_cov10[(Blood_S6_m6A_cov10["Position"].isin(Blood_3_m6a_sites_joint_concat)) & (Blood_S6_m6A_cov10["Ratio"] >= 10)]
Blood_ivt_violin_sig = Blood_S6_IVT_merged_m6A_cov10[(Blood_S6_IVT_merged_m6A_cov10["Position"].isin(Blood_IVT_m6a_sites_joint_concat)) & (Blood_S6_IVT_merged_m6A_cov10["Ratio"] >= 10)]
Glori_1_violin_sig = Blood_Glori_1[(Blood_Glori_1["Position"].isin(Glori_1_m6a_sites_joint_concat)) & (Blood_Glori_1["Ratio"] >= 10)]
Glori_2_violin_sig = Blood_Glori_2[(Blood_Glori_2["Position"].isin(Glori_2_m6a_sites_joint_concat)) & (Blood_Glori_2["Ratio"] >= 10)]

#Assign frequencies for plotting 
blood_1_box = Blood_1_violin_sig 
blood_2_box = Blood_2_violin_sig  
blood_3_box = Blood_3_violin_sig
ivt_1_box = Blood_ivt_violin_sig
glori_1_box = Glori_1_violin_sig
glori_2_box = Glori_2_violin_sig

#dictionary for plotting
plot_data = {
  'GLORI 1': glori_1_box["Ratio"].tolist(),
  'GLORI 2': glori_2_box["Ratio"].tolist(),
  'Blood 1': blood_1_box["Ratio"].tolist(),
  'Blood 2': blood_2_box["Ratio"].tolist(),
  'Blood 3': blood_3_box["Ratio"].tolist(),
  'Blood IVT': ivt_1_box["Ratio"].tolist(),
}


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
plt.savefig("Fig_SUPP_4B_violin.pdf", format="pdf", bbox_inches="tight", dpi=600)
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