source /users/PAS1117/osu9433/ProjectDir_Cronin/Bioinformatic_tools/miniconda3/bin/activate py27

/users/PAS1117/osu9433/ProjectDir_Cronin/Bioinformatic_tools/coverm/coverm2 genome --single /fs/project/PAS1117/SCI/01_Clean_reads/*clean.fq.gz -f /fs/project/PAS1117/00_SCI_drep_MAGs/dereplicated_genomes_60_percent/*.fasta --bam-file-cache-directory /fs/scratch/Sullivan_Lab/Jingjie_SCI/SCI_Ion_torrent_read_mapping_all/ --min-read-percent-identity .95 --min-read-aligned-percent .75 -t 20 > /fs/project/PAS1117/00_SCI_drep_MAGs/read_mapping_MAGs_new.txt
