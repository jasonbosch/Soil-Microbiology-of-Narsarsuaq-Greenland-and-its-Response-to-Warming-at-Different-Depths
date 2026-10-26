MARKERS=("16S" "ITS")

echo "Started at" $(date)

for MARKER in ${MARKERS[@]}
do
	
	#Prepare directory
	if test -d "/run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned"; then
		echo "Directory for cleaned $MARKER reads exists, skipping directory creation."
	else
		echo "Creating directory for cleaned $MARKER reads..."
		mkdir /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned
		echo "Done at" $(date)
	fi
	
	#Set libraries
	if [[ ${MARKER} == "16S" ]]; then
		echo "Using ${MARKER} libraries."
		LIBRARIES=("LibD28" "LibD35")
	else
		echo "Using ${MARKER} libraries."
		LIBRARIES=("LibD28" "LibD37")
	fi
	
	#Set reads
	READS=("R1" "R2")
	
	for LIBRARY in ${LIBRARIES[@]}
	do
	
		for READ in ${READS[@]}
		do
			#Identify samples
			if test -f "/run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned/${LIBRARY}_${READ}_names.txt"; then
				echo "Samples for $MARKER $LIBRARY ${READ}already identified, skipping sample identification."
			else
				echo "Identifying $MARKER $LIBRARY ${READ}samples..."
				grep '^@M' /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/01_Raw_Data/original_fastq_${MARKER}/${LIBRARY}_CAM_${READ}.fastq /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/01_Raw_Data/original_fastq_${MARKER}/${LIBRARY}_ALG_${READ}.fastq | awk -F '|' '{print $2}' | sort | uniq | awk -F ' ' '{print $1}' > /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned/${LIBRARY}_${READ}_names.txt
			fi
			
			#Split fastq file
			if [[ $(ls /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned/${LIBRARY}*${READ}.fastq | wc -l) -eq $(cat /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned/${LIBRARY}_${READ}_names.txt | wc -l) ]]; then
				echo "$MARKER $LIBRARY ${READ} already split, skipping split."
			else
				echo "Splitting $MARKER $LIBRARY ${READ}..."
				for i in $(cat /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned/${LIBRARY}_${READ}_names.txt); do grep --no-group-separator -A3 $i /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/01_Raw_Data/original_fastq_${MARKER}/${LIBRARY}_*_${READ}.fastq > /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/${MARKER}_Cleaned/${LIBRARY}_${i}_${READ}.fastq ;done
				echo "Done at" $(date)
			fi
		done
	done
done
