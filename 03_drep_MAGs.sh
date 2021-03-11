cd /fs/project/PAS1117/00_SCI_drep_MAGs/
#source /users/PAS1117/osu9433/ProjectDir_Cronin/Bioinformatic_tools/miniconda3/bin/activate unitem
#checkm lineage_wf -t 20 -x fasta /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/maxbin/ checkm_out
#checkm qa --tab_table -t 20 checkm_out/lineage.ms checkm_out/ > checkm_out/qa_report.txt
if cd checkm_out/ ; then
python /users/PAS1117/osu9433/ProjectDir_Cronin/Bioinformatic_tools/miniconda3/envs/unitem/bin/get_genome_info_fasta.py
cd ..
fi
source /users/PAS1117/osu9433/ProjectDir_Cronin/Bioinformatic_tools/miniconda3/bin/activate dRep
dRep dereplicate dRep_97 -g /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/maxbin/*.fasta -comp 60 -con 10 -p 20 --genomeInfo checkm_out/genomeinfo.csv -sa 0.97
dRep dereplicate dRep_95 -g /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/maxbin/*.fasta -comp 60 -con 10 -p 20 --genomeInfo checkm_out/genomeinfo.csv -sa 0.95
