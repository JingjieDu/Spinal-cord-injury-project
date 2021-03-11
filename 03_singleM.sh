
source /users/PAS1117/osu8392/local/src/anaconda3/bin/activate /users/PAS1117/osu8392/local/src/anaconda3/envs/cocacola_env/
export PATH=$PATH:/users/PAS1117/osu8392/local/src/bin/
export PATH=$PATH:/users/PAS1117/osu8392/local/src/fxtract/
export PATH=$PATH:/users/PAS1117/osu8392/local/src/diamond_v0.9.22
export PATH=$PATH:/users/PAS1117/osu8392/local/src/orfm-0.7.1_Linux_x86_64
export PATH=$PATH:/users/PAS1117/osu8392/local/src/Krona/bin
export PATH=$PATH:/users/PAS1117/osu8392/local/src/pplacer-Linux-v1.1.alpha19
export PATH=$PATH:/users/PAS1117/osu8392/local/src/vsearch/bin
export PATH=$PATH:/users/PAS1117/osu8392/local/src/smafa_0.1


# export a required perl library to the environment

export PERL5LIB=/users/PAS1117/osu8392/perl5/lib/perl5/
export PYTHONPATH=$PYTHONPATH:/users/PAS1117/osu8392/.local/lib/python2.7/site-packages/

# Start the process!
cd /fs/project/PAS1117/SCI/
perl /users/PAS1117/osu8392/local/src/squeezeM/scripts/00.run_singleM.pl SCI_ion_all_coassembly
