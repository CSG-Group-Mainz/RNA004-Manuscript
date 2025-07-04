# Custom vector pseudouridine calling

The script in this section was written to analyse the data in the section ***"RNA004 accurately reads the pseudouridylation stochiometry in a targeted reporter system"*** in the main paper. It will take raw *pod5 data as input and emit aligned bam files + modified bed files + modification probability histograms. Of additional note in this section is that the standard dorado pseudouridine calling was edited to call C bases instead of U bases, in order to evaluate how the pseudouridine modbase model does perform on miscalled C bases from the canonical model.


In order to prepare the custom modbase calling, the normal models have to be downloaded first:  
```
dorado download --model rna004_130bps_sup@v5.0.0
dorado download --model rna004_130bps_sup@v5.0.0_pseU@v1
```  
Then, one does navigate to the model folder:  
```
cd rna004_130bps_sup@v5.0.0_pseU@v1/
```
And edit the config file in your favorite text editor, here we chose vim:  
```
vim config.toml
```
Now the relevant line is line 19 starting with motif=T:  
![image](https://github.com/user-attachments/assets/69eee054-79b3-49dc-a2ac-b6fc9173070c)  
We change line 19 to C and - for good measure - line 11 mod_bases to another modcode that is not associated with pseu. So there is less chance of mixing up our results when evaluating the mod bed files later on. That means, the resulting file will now look like this:
![image](https://github.com/user-attachments/assets/4095586e-4921-423b-afc9-a0f13a1c7525)  
Now, after the file has been closed and saved, we are able to run dorado with the custom modbase model, by specifying the exact filepaths. Please keep in mind that, depending on the other models run, the dorado version and your graphic card model the batch size needs to be adjusted. For us (A100 + dorado 0.7.2) that was --b 320.
```
dorado basecaller ./rna004_130bps_sup@v5.0.0 {pod5_input_data}   -r --emit-moves --modified-bases-models ./rna004_130bps_sup@v5.0.0_pseU@v1 > {output.unaligned.bam}
```
After this step the bam files can be aligned and processed with modkit pileup as normal.


