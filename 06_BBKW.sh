module use /fs/project/PAS1117/modulefiles
module load WrightonPipeline

cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/03.BBKW_annotations
Wrighton-Sullivan_Pipeline.py -i Viral_contigs_5KB.fasta -t nucl -m 1500 -o /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/00_Viral_contigs --tool pfam
