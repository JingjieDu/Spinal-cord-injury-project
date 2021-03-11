module load python 
source /users/PAS1117/osu8392/local/src/anaconda3/bin/activate /users/PAS1117/osu8392/local/src/anaconda3
 
cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/ 
python /fs/project/PAS1117/viral_hhm/tools/DeepVirFinder/dvf.py -i /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/01.SCI_ion_all_coassembly.fasta -o /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/ -l 300 -c 12
