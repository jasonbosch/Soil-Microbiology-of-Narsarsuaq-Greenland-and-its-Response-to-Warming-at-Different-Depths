##Clustering ITS Sequences into 97% OTUs with VSEARCH
#Largely based on https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline

#Set the path
CURR_PATH="/run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/ITS_Cleaned"

#Dereplicate the sequences
if [[ $(ls ${CURR_PATH} | grep -c "_Joined_ITSxpress.fastq") -eq $(ls ${CURR_PATH} | grep -c ".derep.fasta") ]]; then
	echo "Expected number of files, skipping dereplication."
	else
	echo "Dereplicating the sequences."
	for FILE in ${CURR_PATH}/*_Joined_ITSxpress.fastq
	do
		SAMPLE=$(echo $FILE | sed s/".\+\/"/""/ | sed s/"_Joined_ITSxpress.fastq"/""/)
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--derep_prefix $FILE \
		--output $(echo $FILE | sed s/".fastq"/""/).derep.fasta \
		--relabel $SAMPLE. \
		--sizeout \
		--fasta_width 0
	done
fi

#Combine all dereplicated files
echo "Combining dereplicated files"
cat $CURR_PATH/*.derep.fasta > $CURR_PATH/ITS_Combined.fasta

#Dereplicate across samples
if test -f ${CURR_PATH}/"ITS_Combined_Dereplicated_Sequences.fasta"; then
	echo "Combined seqeunces already dereplicated, skipping combined dereplication."
	else
	echo "Dereplicating the combined seqeunces."
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--derep_prefix $CURR_PATH/ITS_Combined.fasta \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--uc $CURR_PATH/ITS_Combined_Dereplicated_Table.tsv \
		--output $CURR_PATH/ITS_Combined_Dereplicated_Sequences.fasta
fi

#Cluster at 97% similarity
if test -f ${CURR_PATH}/"ITS_Clusters_0.97.fasta"; then
	echo "Samples already clustered, skipping clustering."
	else
	echo "Clustering the sequences."
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--cluster_size $CURR_PATH/ITS_Combined_Dereplicated_Sequences.fasta \
		--threads 6 \
		--id 0.97 \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--centroids $CURR_PATH/ITS_Clusters_0.97.fasta
fi

#Remove singletons
if test -f ${CURR_PATH}/"ITS_Clusters_0.97_NoSingletons.fasta"; then
	echo "Singletons already removed, skipping step."
	else
	echo "Removing singletons."
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--sortbysize $CURR_PATH/ITS_Clusters_0.97.fasta \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--minsize 2 \
		--output $CURR_PATH/ITS_Clusters_0.97_NoSingletons.fasta
fi

#Identify non-chimeras with de novo chimaera detection
if test -f ${CURR_PATH}/"ITS_Clusters_0.97_NoSingletons_NoChimeras.fasta"; then
	echo "Chimeras already removed, skipping step."
	else
	echo "Removing chimeras."
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--uchime_denovo $CURR_PATH/ITS_Clusters_0.97_NoSingletons.fasta \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--qmask none \
		--nonchimeras  $CURR_PATH/ITS_Clusters_0.97_NoSingletons_NoChimeras.fasta
fi



#Relabel OTUs
if test -f ${CURR_PATH}/"ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs.fasta"; then
	echo "OTUs already relabeled, skipping step."
	else
	echo "Relabelling OTUs."
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--fastx_filter $CURR_PATH/ITS_Clusters_0.97_NoSingletons_NoChimeras.fasta \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--relabel OTU_ \
		--fastaout $CURR_PATH/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs.fasta
fi

#Map sequences to OTUs
if test -f ${CURR_PATH}/"ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs_Table.txt"; then
	echo "Sequences already mapped to OTUs, skipping step."
	else
	echo "Mapping sequences to OTUs..."
		/run/media/SharedStorage/MBU_Postdoc/04_Software/vsearch-2.30.4-linux-x86_64/bin/vsearch \
		--usearch_global $CURR_PATH/ITS_Combined.fasta \
		--threads 6 \
		--db $CURR_PATH/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs.fasta \
		--id 0.97 \
		--strand plus \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--qmask none \
		--dbmask none \
		--otutabout $CURR_PATH/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs_Table.txt
fi

#Sort the OTU table
sort -k1.5n $CURR_PATH/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs_Table.txt > /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs_Table_Sorted.txt

#Remove the comment
sed -i s/"#"/""/ /run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs_Table_Sorted.txt

