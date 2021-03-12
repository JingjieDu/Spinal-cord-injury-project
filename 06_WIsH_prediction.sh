cd /users/PAS1117/osu10035/Bins_lt_60_perc/
module use /fs/project/PAS1117/modulefiles
module load WIsH
WIsH -c predict -g /users/PAS1117/osu10035/Viral_contigs/ -m /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/WIsH_model/ -r /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/WIsH_results_control/ -b -n /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Host_predictions/WIsH_results_control/nullParameters.tsv -t 28 -b
