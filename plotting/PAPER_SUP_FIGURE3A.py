import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from decimal import *


####### Function to iterate over varying filterthresholds of bedfiles and comparing with validated modifications sites 
#Input:
#bedfile: String of Path to bedfile for threshold iterations (Bedfile obtained after using modkit)
#validation: List of validated modification sites 
#modification: String of modification code to be analyzed in the bedfile 
#Output: Name of saved SVG

def filter_iteration(bedfile: str, validation: list, modification: str, output: str, Mode):

    #Clear potential previous plots
    plt.clf()

    #Read bedfile and reformat for further use
    bed =  pd.read_csv(bedfile, delimiter='\t')
    bed.columns = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"]
    bed = bed.loc[bed["4"] == modification]    

    unfiltered_concat = bed.iloc[:,0].astype(str) + ":" + bed.iloc[:,1].astype(str) + ":" + bed.iloc[:,5].astype(str)

    #####INIT PLOT PARAMS
    plotdata = {}
    ratios = np.arange(0, 36, 1)
    coverages = np.arange(0, 55, 5)

    #Calculate Validation Proportions with different threshold settings
    for coveragefilter in coverages:
        print(coveragefilter)
        current= []
        for ratiofilter in ratios:
            #Filter by coverage and ratio according to current iteration
            currentbed = bed.loc[(bed["5"] >= coveragefilter) & (bed["11"] >= ratiofilter)]
            #Obtain current position concats to compare with validated positions
            concat = currentbed.iloc[:,0].astype(str) + ":" + currentbed.iloc[:,1].astype(str) + ":" + currentbed.iloc[:,5].astype(str)

            #Calculate number of true positives and false positives by checking intersections with the validation positions
            true_pos = list(set(concat) & set(validation))
            false_pos = list(set(concat) - set(validation))

            #Get the length of TP and FP
            tp = len(true_pos)
            fp = len(false_pos)


            Recall = Decimal(tp/len(Glori_consensus_pos))
            Precision = Decimal(tp/(tp + fp))

            F1 = 2 * ((Precision*Recall)/(Precision+Recall))
            #Here I still need to figure out which metric to use ... Absolute numbers, ratios, proper metrics .... ask charlotte
            #current.append((tp/len(Glori_consensus_pos)))
            if Mode == "Recall":
                current.append((tp/len(Glori_consensus_pos)))

            if Mode == "FPR":
                current.append((len(concat)/len(unfiltered_concat)))

            if Mode == "Precision":
                current.append(tp/len(concat))

            if Mode == "F1":
                current.append(F1)

        plotdata[str(coveragefilter)] = current

    #Convert data to long Format
    data_complete = pd.DataFrame()
    for key in plotdata:
        data_current = pd.DataFrame()
        data_current[Mode] = plotdata[key]
        data_current['Ratio'] = ratios
        data_current['Coverage'] = [str(key)] * len(plotdata[key]) 
        data_complete = pd.concat([data_complete, data_current])

    pdf_title = output.split(".")[0]
    pdf_title = pdf_title + ".csv"

    data_complete.to_csv(pdf_title)

    #Plot 
    ticks = list(np.arange(0, 36, 1))
    plot = sns.lineplot(data=data_complete, x = "Ratio", y = Mode, hue="Coverage", palette="colorblind")
    plot.tick_params(axis='x', rotation=45)
    plot.set_xticks(ticks)
    sns.move_legend(plot, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig(output, bbox_inches='tight')



###################################################################################################################################

##### BLOOD M6A #####

#Read Glori Sites of both blood samples
Glori_1 = pd.read_csv("new_data/AS-1466474-LR-78969_R1.totalm6A.FDR.csv", delimiter='\t')
Glori_2 = pd.read_csv("new_data/AS-1466472-LR-78969_R1.totalm6A.FDR.csv", delimiter='\t')

#### REFORMAT #### All base positions must be shifted 1 to the left to compare with bedfiles
Glori_1_base = list(Glori_1.iloc[:,1].astype(int))
Glori_1_base = list(map(lambda x : x - 1, Glori_1_base))
Glori_2_base = list(Glori_2.iloc[:,1].astype(int))
Glori_2_base = list(map(lambda x : x - 1, Glori_2_base))

#Create concatenations
Glori_1_pos = Glori_1.iloc[:,0].astype(str) + ":" + np.asarray(Glori_1_base).astype(str) + ":" + Glori_1.iloc[:,2].astype(str)
Glori_2_pos = Glori_2.iloc[:,0].astype(str) + ":" + np.asarray(Glori_2_base).astySpe(str)  + ":" + Glori_2.iloc[:,2].astype(str)

Glori_1_pos = list(Glori_1_pos)
Glori_2_pos = list(Glori_2_pos)

#Find Consensus positions between both bedfiles
Glori_consensus_pos = list(set(Glori_1_pos) & set(Glori_2_pos))

#Blood Sample vs GLORI Sample as ground truth 
filter_iteration(bedfile="new_data/RNA004_S5_DRS_basecall.0.7.2.GRCh38_m6A.r1.mod.bed", validation=Glori_consensus_pos,
                  modification="a", output="SUPP_FIG_3_A.svg", Mode="F1")