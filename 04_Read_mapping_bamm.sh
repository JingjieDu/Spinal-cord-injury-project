module load samtools
module use /users/PAS1117/osu9156/local/share/modulefiles/
module load python
module load pysam
module load bamm
module load bwa

cd /fs/project/PAS1117/SCI/01_Clean_reads/

bamm make -d /fs/project/PAS1117/SCI/03_renamed_scaffolds/Lam_scaffolds.fasta -s Lam_16_IonXpress_007_SINGLETONS_clean.fq.gz Lam_20_IonXpress_008_SINGLETONS_clean.fq.gz Lam_32_IonXpress_009_SINGLETONS_clean.fq.gz Lam_34_IonXpress_005_SINGLETONS_clean.fq.gz Lam_4_IonXpress_006_SINGLETONS_clean.fq.gz -t 12 -m 190 -f -o /fs/project/PAS1117/SCI/03_renamed_scaffolds/
bamm make -d /fs/project/PAS1117/SCI/03_renamed_scaffolds/T10_scaffolds.fasta -s T10_12_IonXpress_003_SINGLETONS_clean.fq.gz T10_14_IonXpress_004_SINGLETONS_clean.fq.gz T10_56_IonXpress_005_SINGLETONS_clean.fq.gz T10_58_IonXpress_006_SINGLETONS_clean.fq.gz T10_60_IonXpress_007_SINGLETONS_clean.fq.gz -t 12 -m 190 -f -o /fs/project/PAS1117/SCI/03_renamed_scaffolds/
bamm make -d /fs/project/PAS1117/SCI/03_renamed_scaffolds/T4_scaffolds.fasta -s T4_10_IonXpress_007_SINGLETONS_clean.fq.gz T4_36_IonXpress_008_SINGLETONS_clean.fq.gz T4_38_IonXpress_009_SINGLETONS_clean.fq.gz T4_52_IonXpress_010_SINGLETONS_clean.fq.gz T4_6_IonXpress_006_SINGLETONS_clean.fq.gz -t 12 -m 190 -f -o /fs/project/PAS1117/SCI/03_renamed_scaffolds/
