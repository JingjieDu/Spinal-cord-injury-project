module load singularity/current

cd /fs/project/PAS1117/00_SCI_VconTACT_Jingjie/
/fs/project/PAS1117/bioinformatic_tools/Prodigal2.6.1/prodigal -i Viral_contigs_10KB_95-80.fna -p meta -a Viral_contigs_10KB_95-80_proteins_new.faa

singularity run /fs/project/PAS1117/modules/singularity/vContact-Gene2Contig-1.0.0.img -p SCI_viral_contigs_All_GVD_protein.faa -o SCI_viral_contigs_All_GVD_protein.csv -s Prodigal-FAA

singularity run /fs/project/PAS1117/modules/singularity/vConTACT2-0.9.9.simg -r SCI_viral_contigs_All_GVD_protein.faa --rel-mode BLASTP -p SCI_viral_contigs_All_GVD_protein.csv --db ProkaryoticViralRefSeq88-Merged --pcs-mode MCL --vcs-mode ClusterONE -t 28 -f
