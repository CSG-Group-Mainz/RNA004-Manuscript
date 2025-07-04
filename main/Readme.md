This script does process the raw pod5 data for all of the main samples in this manuscript. The output are aligned bam files (GRCh38), basic QC statistics (NanoComp) and gene counts (FeatureCounts). 

##### Some notes for the basecaller options
Dorado 0.7.2 in the version for this paper had to be run with a number of custom and legacy options, below you will find basic examples for RNA004.
First, the models had to be downloaded and then specified as filepaths, because dorado can not extract the chemistry version for all samples reliably (note: this is not an issue for newer runs, just runs that were produced before the chemistry left the beta testing stage). Example command to download a canonical model:
```
dorado download --model rna004_130bps_sup@v5.0.0
```  
###### Example command for RNA004 chemistry
```
dorado basecaller ./rna004_130bps_sup@v5.0.0 ${pod5_input_folder} --estimate-poly-a -r --emit-moves --modified-bases m6A pseU -b 320 --device "cuda:0,cuda:1,cuda:2,cuda:3" > ${unaligned_output_bam}
```  
Option **'-b 320'**:  
This dorado version had a known bug, where the automatic batch size detection of the basecaller was not working, causing runs to crash when specifying a sup model and two modbase models [link](https://github.com/nanoporetech/dorado/issues/866). Thus, batchsize had to be adjusted manually, for our graphic card (A100) 320 was working.  
Option **'--device "cuda:0,cuda:1,cuda:2,cuda:3"'**:  
In newer dorado versions the synthax has slightly changed and running like this will cause the run to crash. In addition, this option can only be selected, if 4 graphic cards of the same type are available.

