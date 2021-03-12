module use /fs/project/PAS1117/modulefiles
module load tRNAscan-SE/1.23
cd  /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/tRNA_scan/
#look the seceond structures of viral_contigs
tRNAscan-SE -G -f /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/tRNA_scan/Viral_contigs_5KB_tRNA.fa /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/Viral_contigs_5KB.fasta
#look for the second structures of host sequences
tRNAscan-SE -B -f /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/tRNA_scan/all_maxbin_second_tructure_tRNA.fa /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/tRNA_scan/all_maxbin.fa
#remove the second structures and only keep the sequences
grep -e 'maxbin' -e 'Seq: ' all_maxbin_second_structure_tRNA.fa | sed 's/maxbin/>maxbin/g' | sed 's/Seq: //g' > all_maxbin_second_structure_tRNA_clean.fa
