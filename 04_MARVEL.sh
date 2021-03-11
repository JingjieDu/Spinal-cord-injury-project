module load python/3.6
export PERL5LIB=/users/PAS1117/osu8392/perl5/lib/perl5/
export PYTHONPATH=$PYTHONPATH:/users/PAS1117/osu8392/.local/lib/python3.6/site-packages
cd /users/PAS1117/osu8392/local/src/MARVEL/
python marvel_bins.py -i /fs/scratch/Sullivan_Lab/SCI/input_for_MARVEL -t 28
