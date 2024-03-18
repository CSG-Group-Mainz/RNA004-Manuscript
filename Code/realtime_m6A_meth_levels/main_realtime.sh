#!/usr/bin/bash

i=10

#input paths where the POD5 files need to be located
IVT_PATH=/mnt/ssd_share_01/new_folder/RESEARCH/RNA004_runs/RNA004_peripheral_blood/2024_02_28_401_24_IVT/20240228_2051_2F_PAS25093_1dc40d55/pod5_pass
DIR_PATH=/mnt/ssd_share_01/new_folder/RESEARCH/RNA004_runs/RNA004_peripheral_blood/2024_02_28_401_24/20240228_2051_2G_PAS24854_4617597d/pod5_pass
#output path where analysis results should be stored
OUT_IVT_PATH=/raid/fehof_analysis/parallel_test/realtime_output/IVT
OUT_DIR_PATH=/raid/fehof_analysis/parallel_test/realtime_output/Direct
#location of the script to perform analyses
SCRIPT=/raid/fehof_analysis/parallel_test/Direct_RNA004/parallel_analysis_rt.sh

#carry out while loop until 100 files are processed
while [ $i -le 200 ]; do

	#first check for Direct run output and see if already i files are there
	cd $DIR_PATH
	if [ $(ls -1 | wc -l) -ge $i ]; then
		echo "$i Direct POD5 files are available, check for IVT POD5 files"
		#next go to IVT path and verify that there are also already i files available
		cd $IVT_PATH
		if [ $(ls -1 | wc -l) -ge $i ]; then
			echo "$i IVT POD5 files are available, start copying files"
			#copy the first i IVT pod5 files to a new subdirectory
			files=$( ls $IVT_PATH/*.pod5 | head -n $i )
			cd $OUT_IVT_PATH
        		mkdir $i.basecall
        		cp $files $OUT_IVT_PATH/$i.basecall
			echo "$i IVT files copied"

			#copy the first i Direct pod5 files to a new subdirectory
			cd $DIR_PATH
                	files=$( ls $DIR_PATH/*.pod5 | head -n $i )
                	cd $OUT_DIR_PATH
                	mkdir $i.basecall
                	cp $files $OUT_DIR_PATH/$i.basecall
			echo "$i Direct files copied"

			#analyse the data with parallel_analysis.sh script
			bash $SCRIPT $OUT_IVT_PATH/$i.basecall $OUT_DIR_PATH/$i.basecall $i

			#fuse the two metagene plots for Direct and IVT together
			/raid/fehof_analysis/software/magick convert +append \
			$OUT_IVT_PATH/$i.basecall/$i.genebody_methylation.RNA004_GRCh38_aligned.png \
			$OUT_DIR_PATH/$i.basecall/$i.genebody_methylation.RNA004_GRCh38_aligned.png \
			$OUT_DIR_PATH/$i.metagene_m6A_combined.png

			#create gif out of the merged metagene plots
			/raid/fehof_analysis/software/magick convert -delay 50 -loop 0 $(ls -v $OUT_DIR_PATH/*.png) $OUT_DIR_PATH/metagene_m6A.gif 
			echo "created corrected gif for $i files"
		#increment i by 10
	           	((i+=10))
		#sleep 20 seconds if not enough IVT files are there
		else
			echo "sleeping 20s - not enough IVT POD5 files available (target: $i files)"
	               	sleep 20
		fi
	#sleep 20 seconds if not enough Direct files are there
	else
		echo "sleeping 20s - not enough Direct POD5 files available (target: $i files)"
		sleep 20
	fi
done
echo "done"
