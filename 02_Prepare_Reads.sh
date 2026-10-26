#Must run using "bash -login" for conda to work
echo "Started at" $(date)

CURR_PATH="/run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/"

#Iterate over both sequenced markers
for MARKER in {"16S","ITS"}
do

	#Remove primer sequences with cutadapt
	#First four bases of primers were already removed as part of the barcode
	forward_16S="CAGCMGCCGCGGTAA" #515F
	reverse_16S="TACHVGGGTWTCTAAT" #806R
	forward_ITS="RTCATCGARTCTTTG" #gITS7
	reverse_ITS="CCGCTTATTGATATGC" #ITS4
	if [[ $(ls ${CURR_PATH}${MARKER}_Cleaned/*_R1_Cut.fastq | wc -l) -eq $(cat ${CURR_PATH}${MARKER}_Cleaned/*_R1_names.txt | wc -l) ]]; then
			echo "Primers of ${MARKER} already cut, skipping cutting."
		else
		echo "Cutting primers from ${MARKER} sequences."
		conda activate cutadapt-5.2
		for FILE in ${CURR_PATH}${MARKER}_Cleaned/*_R1.fastq
		do
			FILE_REV=$( echo $FILE | sed s/"_R1"/"_R2"/ )
			FILE_OUT=$( echo $FILE | sed s/".fastq"/"_Cut.fastq"/ )
			FILE_REV_OUT=$( echo $FILE_REV | sed s/".fastq"/"_Cut.fastq"/ )
			
			if [[ ${MARKER} == "16S" ]]; then
			cutadapt -b $forward_16S -B $reverse_16S -o $FILE_OUT -p $FILE_REV_OUT $FILE $FILE_REV --minimum-length 200 --quiet 
			else
			cutadapt -b $forward_ITS -B $reverse_ITS -o $FILE_OUT -p $FILE_REV_OUT $FILE $FILE_REV --minimum-length 200 --quiet 
			fi
		done
		conda deactivate
		echo "Finished cutting ${MARKER} at" $(date)
	fi

	#Use Trimmomatic to remove low quality reads, sliding window of 6 bases and quality 20 SLIDINGWINDOW:6:20
	#Remove sequences shorter than 100 bp MINLEN:100
	if [[ $(ls ${CURR_PATH}${MARKER}_Cleaned/*_Cut.fastq | wc -l) -eq $(ls ${CURR_PATH}${MARKER}_Cleaned/*_Cut_Trimmed.fastq | wc -l) ]]; then
			echo "Expected number of trimmed ${MARKER} sequences, skipping trimming."
		else
			echo "Trimming ${MARKER} sequences..."
			for FILE in ${CURR_PATH}${MARKER}_Cleaned/*_R1.fastq
			do
				SAMPLE=$(echo $FILE | sed s/"_R1.fastq"/""/ )	
				java -jar /run/media/SharedStorage/MBU_Postdoc/04_Software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 7 -summary ${SAMPLE}_trimsummary.txt ${SAMPLE}_R1_Cut.fastq ${SAMPLE}_R2_Cut.fastq ${SAMPLE}_R1_Cut_Trimmed.fastq ${SAMPLE}_R1_Cut_Trimmed_unpaired.fastq ${SAMPLE}_R2_Cut_Trimmed.fastq ${SAMPLE}_R2_Cut_Trimmed_unpaired.fastq SLIDINGWINDOW:6:20 MINLEN:100
			done
			echo "Finished trimming ${MARKER} at" $(date)
	fi


	if [[ ${MARKER} == "ITS" ]]; then 
	
		#Join the reads
		if [[ $(ls ${CURR_PATH}${MARKER}_Cleaned/*_join.fastq | wc -l) -eq $(ls ${CURR_PATH}${MARKER}_Cleaned/*R1.fastq | wc -l) ]]; then
				echo "Expected number of joined ${MARKER} sequences, skipping joining."
			else
				echo "Joining ${MARKER} sequences..."
				for FILE in ${CURR_PATH}${MARKER}_Cleaned/*R1.fastq
				do
					SAMPLE=$(echo $FILE | sed s/"_R1.fastq"/""/)
					
					/run/media/SharedStorage/MBU_Postdoc/04_Software/ea-utils-1.04.807/clipper/fastq-join -v " " -p 15 -m 40 ${SAMPLE}_R1_Cut_Trimmed.fastq ${SAMPLE}_R2_Cut_Trimmed.fastq -o ${SAMPLE}_%.fastq
				done
				echo "Joining ITS reads finished at" $(date)
		fi

		#Extract ITS region
		if [[ $(ls ${CURR_PATH}${MARKER}_Cleaned/*_join.fastq | wc -l) -eq $(ls ${CURR_PATH}${MARKER}_Cleaned/*_Joined_ITSxpress.fastq | wc -l) ]]; then
				echo "Expected number of ${MARKER} extracted sequences, skipping ITS extraction."
			else
				echo "Extracting ${MARKER} sequences..."
				conda activate ITSxpress_2.1.4
				for FILE in ${CURR_PATH}${MARKER}_Cleaned/*_join.fastq
				do
					SAMPLE=$(echo $FILE | sed s/"_join.fastq"/""/ )
					
					#Extract the ITS region
					itsxpress --threads 7 --region ITS2 --taxa Fungi --single_end --fastq $FILE --log ${SAMPLE}_ITSxpress_logfile.txt --outfile ${SAMPLE}_Joined_ITSxpress.fastq
				done
				conda deactivate
			echo "ITS extracted at" $(date)
		fi
	fi

	#Create a summary of the process
	echo "Creating ${MARKER} Summary"

	if [[ ${MARKER} == "ITS" ]]; then 
			echo "Sample,Forward,Forward_Trimmed,Merged,ITS-Extracted" > ${CURR_PATH}${MARKER}_Cleaned/${MARKER}_Processing_Summary.csv
			for FILE in ${CURR_PATH}${MARKER}_Cleaned/*_Joined_ITSxpress.fastq
			do
				SAMPLE=$(echo $FILE | sed s/"_Joined_ITSxpress.fastq"/""/ | sed s/"\/run\/media\/SharedStorage\/MBU_Postdoc\/02_Projects\/13_Narsarsuaq\/02_Cleaned_Data\/${MARKER}_Cleaned\/"/""/ )
				INPUTFORWARD=$(( $(cat ${CURR_PATH}${MARKER}_Cleaned/${SAMPLE}_R1.fastq | wc -l) / 4 ))
				TRIMFORWARD=$(( $(cat ${CURR_PATH}${MARKER}_Cleaned/${SAMPLE}_R1_Cut_Trimmed.fastq | wc -l) / 4 ))
				MERGED=$(( $(cat ${CURR_PATH}${MARKER}_Cleaned/${SAMPLE}_join.fastq | wc -l) / 4 ))
				EXTRACTED=$(( $(cat ${CURR_PATH}${MARKER}_Cleaned/${SAMPLE}_Joined_ITSxpress.fastq | wc -l) / 4 ))
				echo "${SAMPLE},$INPUTFORWARD,$TRIMFORWARD,$MERGED,$EXTRACTED" >> ${CURR_PATH}${MARKER}_Cleaned/${MARKER}_Processing_Summary.csv
			done
		else
			echo "Sample,Forward,Forward_Trimmed" > ${CURR_PATH}${MARKER}_Cleaned/${MARKER}_Processing_Summary.csv
			for FILE in ${CURR_PATH}${MARKER}_Cleaned/*_R1_Cut_Trimmed.fastq
			do
				SAMPLE=$(echo $FILE | sed s/"_R1_Cut_Trimmed.fastq"/""/ | sed s/"\/run\/media\/SharedStorage\/MBU_Postdoc\/02_Projects\/13_Narsarsuaq\/02_Cleaned_Data\/${MARKER}_Cleaned\/"/""/ )
				INPUTFORWARD=$(( $(cat ${CURR_PATH}${MARKER}_Cleaned/${SAMPLE}_R1.fastq | wc -l) / 4 ))
				TRIMFORWARD=$(( $(cat ${CURR_PATH}${MARKER}_Cleaned/${SAMPLE}_R1_Cut_Trimmed.fastq | wc -l) / 4 ))
				echo "${SAMPLE},$INPUTFORWARD,$TRIMFORWARD" >> ${CURR_PATH}${MARKER}_Cleaned/${MARKER}_Processing_Summary.csv
			done
	fi
done

echo "Finished at" $(date)
