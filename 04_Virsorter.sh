module load bioperl/1.6.1
module use /fs/project/PAS1117/modulefiles
module load blast/2.2.26
module use /users/PAS1117/osu9156/local/share/modulefiles/

conda activate virsorter
export PERL5LIB=/home/jdu7/perl5/lib/perl5/
cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/
wrapper_phage_contigs_sorter_iPlant.pl SCI_ion_all_coassembly 01.SCI_ion_all_coassembly.fasta 2 /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/ 12
