cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/BAM_files/

module load bowtie2

bowtie2-build --threads 12 /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/00_Viral_contigs/Viral_contigs_5KB_95-80.fna SCI_viral_populations

