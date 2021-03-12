cd /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/Protein_clustering/

# create database
/users/PAS1117/osu8392/local/src/mmseqs2/bin/mmseqs createdb /fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/03.SCI_ion_all_coassembly.faa targetDB

# index the database
/users/PAS1117/osu8392/local/src/mmseqs2/bin/mmseqs createindex targetDB
# cluster
/users/PAS1117/osu8392/local/src/mmseqs2/bin/mmseqs cluster targetDB DB_clu /fs/scratch/Sullivan_Lab/SCI_tmp --min-seq-id 0.3 --cov-mode 1 -c 0.7 -e 0.00001
# create tsv
/users/PAS1117/osu8392/local/src/mmseqs2/bin/mmseqs createtsv targetDB targetDB DB_clu DB_clu2.tsv


cut -f 1 DB_clu2.tsv > 00.first_column.tsv
awk '!seen[$0]++' 00.first_column.tsv > 01.unique.tsv
nl -n rz 01.unique.tsv > 02.unique_with_numbers.tsv
sed -e 's/^/cluster/' 02.unique_with_numbers.tsv > 03.unique_with_clusters.tsv
awk '{ print $2 " " $1}' 03.unique_with_clusters.tsv > 04.unique_with_clusters_rearranged.tsv
