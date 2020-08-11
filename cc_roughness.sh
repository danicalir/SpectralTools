#!/bin/bash

for slope in 14 20 24 39
do
	for r in 0.1 0.2 0.3 0.4 0.5 #1 2.5 5
	do
		startmessage="[`date`] STARTING $r m roughness kernel on $slope degree point cloud..."
		echo $startmessage
		echo -e "\n\n$startmessage \n" >> $slope.roughnesslog.txt
		C:/Program\ \Files/CloudCompare/CloudCompare -silent -auto_save off -log_file roughtemp.txt -o "$slope-SegmentRoughness.bin" -rough $r -no_timestamp -save_clouds -clear
		endmessage="[`date`] $r m roughness kernel COMPLETE!"
		cat roughtemp.txt >> $slope.roughnesslog.txt
		echo -e "\n$endmessage" >> $slope.roughnesslog.txt
		echo $endmessage
	done
done

exit 0
