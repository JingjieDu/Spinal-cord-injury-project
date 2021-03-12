module load blast/2.4.0+
cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/Prophage_blast/
blastn -query /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/Prophage_blast/all_maxbin.fa -db /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/Viral_contigs_5KB -out Prophage_blast_all_maxbin.txt -num_threads 12 -evalue 0.001 -perc_identity 100 -outfmt "6 qseqid sseqid qlen slen pident length evalue bitscore stitle"
