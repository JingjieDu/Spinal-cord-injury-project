module load python/2.7

# export 3rd party softwares to the environment

export PATH=$PATH:/users/PAS1117/osu8392/local/src
export PATH=$PATH:/users/PAS1117/osu8392/local/src/hmmer-3.2.1/src
export PATH=$PATH:/users/PAS1117/osu8392/local/src/pplacer-Linux-v1.1.alpha19
export PATH=$PATH:/users/PAS1117/osu8392/local/src/FastANI
export PATH=$PATH:/fs/project/PAS1117/biobin/Prodigal-2.6.2

# export the required perl and python libraries to the environment

export PERL5LIB=/users/PAS1117/osu8392/perl5/lib/perl5/
export PYTHONPATH=$PYTHONPATH:/users/PAS1117/osu8392/.local/lib/python2.7/site-packages/

# Start the process!
cd /fs/project/PAS1117/00_SCI_drep_MAGs/
/users/PAS1117/osu8392/local/src/GTDBTk/bin/gtdbtk classify_wf --cpus 40 --genome_dir /fs/project/PAS1117/00_SCI_drep_MAGs/dereplicated_genomes_60_percent --extension fasta --out_dir Bins_taxonomy
