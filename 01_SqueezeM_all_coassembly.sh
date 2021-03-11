# export a required perl library to the environment

export PERL5LIB=/users/PAS1117/osu8392/perl5/lib/perl5/
export PYTHONPATH=$PYTHONPATH:/users/PAS1117/osu8392/.local/lib/python2.7/site-packages/

# Start the process!
cd /fs/project/PAS1117/SCI/
perl /users/PAS1117/osu8392/local/src/squeezeM/scripts/squeezeM.pl -m coassembly -p SCI_ion_all_coassembly -s samples.Ion.Torrent.all -f raw.Ion.Torrent.all/ -t 48 -a spades --spades_options "-m 1500 --iontorrent -k 21,33,55,77,99,127"
