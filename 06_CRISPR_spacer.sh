#get the CRISPR spacers of microbial genomes
cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/CRISPR_spacer/
/fs/project/PAS1117/bioinformatic_tools/minced-master/minced -spacers -minNR 2 all_maxbin.fa microbial_output.crisprs

#use viral_contigs blast against MAGs spacers
cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/CRISPR_spacer/
module load blast
blastn -num_threads 7 -db Viral_contigs_5KB -query microbial_output_spacers.fa -qcov_hsp_perc 95 -word_size 5 -reward 1 -penalty -3 -outfmt 6 -out Viral_contigs_vs_microbial_spacers.bln
