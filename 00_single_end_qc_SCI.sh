cd /fs/project/PAS1117/SCI/00_Raw_reads/
module load bowtie2
module load samtools
module load bamtools

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=Lam_16_IonXpress_007.fastq.gz ftl=10 ftr=249 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=Lam_16_IonXpress_007_SINGLETONS.fq.gz outm=Lam_16_IonXpress_007_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=Lam_16_IonXpress_007_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=Lam_16_IonXpress_007_SINGLETONS.fq.gz outm=Lam_16_IonXpress_007_SINGLETONS_mouse.fq.gz outu=Lam_16_IonXpress_007_SINGLETONS_clean.fq.gz nodisk 
rm Lam_16_IonXpress_007_SINGLETONS.fq.gz
module load fastqc
fastqc Lam_16_IonXpress_007_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=Lam_20_IonXpress_008.fastq.gz ftl=10 ftr=249 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=Lam_20_IonXpress_008_SINGLETONS.fq.gz outm=Lam_20_IonXpress_008_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=Lam_20_IonXpress_008_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=Lam_20_IonXpress_008_SINGLETONS.fq.gz outm=Lam_20_IonXpress_008_SINGLETONS_mouse.fq.gz outu=Lam_20_IonXpress_008_SINGLETONS_clean.fq.gz nodisk 
rm Lam_20_IonXpress_008_SINGLETONS.fq.gz
module load fastqc
fastqc Lam_20_IonXpress_008_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=Lam_32_IonXpress_009.fastq.gz ftl=10 ftr=249 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=Lam_32_IonXpress_009_SINGLETONS.fq.gz outm=Lam_32_IonXpress_009_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=Lam_32_IonXpress_009_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=Lam_32_IonXpress_009_SINGLETONS.fq.gz outm=Lam_32_IonXpress_009_SINGLETONS_mouse.fq.gz outu=Lam_32_IonXpress_009_SINGLETONS_clean.fq.gz nodisk 
rm Lam_32_IonXpress_009_SINGLETONS.fq.gz
module load fastqc
fastqc Lam_32_IonXpress_009_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=Lam_34_IonXpress_005.fastq.gz ftl=10 ftr=239 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=Lam_34_IonXpress_005_SINGLETONS.fq.gz outm=Lam_34_IonXpress_005_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=Lam_34_IonXpress_005_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=Lam_34_IonXpress_005_SINGLETONS.fq.gz outm=Lam_34_IonXpress_005_SINGLETONS_mouse.fq.gz outu=Lam_34_IonXpress_005_SINGLETONS_clean.fq.gz nodisk 
rm Lam_34_IonXpress_005_SINGLETONS.fq.gz
module load fastqc
fastqc Lam_34_IonXpress_005_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T10_14_IonXpress_004.fastq.gz ftl=10 ftr=244 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T10_14_IonXpress_004_SINGLETONS.fq.gz outm=T10_14_IonXpress_004_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T10_14_IonXpress_004_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T10_14_IonXpress_004_SINGLETONS.fq.gz outm=T10_14_IonXpress_004_SINGLETONS_mouse.fq.gz outu=T10_14_IonXpress_004_SINGLETONS_clean.fq.gz nodisk 
rm T10_14_IonXpress_004_SINGLETONS.fq.gz
module load fastqc
fastqc T10_14_IonXpress_004_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T10_56_IonXpress_005.fastq.gz ftl=10 ftr=249 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T10_56_IonXpress_005_SINGLETONS.fq.gz outm=T10_56_IonXpress_005_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T10_56_IonXpress_005_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T10_56_IonXpress_005_SINGLETONS.fq.gz outm=T10_56_IonXpress_005_SINGLETONS_mouse.fq.gz outu=T10_56_IonXpress_005_SINGLETONS_clean.fq.gz nodisk 
rm T10_56_IonXpress_005_SINGLETONS.fq.gz
module load fastqc
fastqc T10_56_IonXpress_005_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T10_58_IonXpress_006.fastq.gz ftl=10 ftr=244 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T10_58_IonXpress_006_SINGLETONS.fq.gz outm=T10_58_IonXpress_006_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T10_58_IonXpress_006_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T10_58_IonXpress_006_SINGLETONS.fq.gz outm=T10_58_IonXpress_006_SINGLETONS_mouse.fq.gz outu=T10_58_IonXpress_006_SINGLETONS_clean.fq.gz nodisk 
rm T10_58_IonXpress_006_SINGLETONS.fq.gz
module load fastqc
fastqc T10_58_IonXpress_006_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T10_60_IonXpress_007.fastq.gz ftl=10 ftr=244 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T10_60_IonXpress_007_SINGLETONS.fq.gz outm=T10_60_IonXpress_007_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T10_60_IonXpress_007_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T10_60_IonXpress_007_SINGLETONS.fq.gz outm=T10_60_IonXpress_007_SINGLETONS_mouse.fq.gz outu=T10_60_IonXpress_007_SINGLETONS_clean.fq.gz nodisk 
rm T10_60_IonXpress_007_SINGLETONS.fq.gz
module load fastqc
fastqc T10_60_IonXpress_007_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T4_10_IonXpress_007.fastq.gz ftl=10 ftr=249 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T4_10_IonXpress_007_SINGLETONS.fq.gz outm=T4_10_IonXpress_007_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T4_10_IonXpress_007_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T4_10_IonXpress_007_SINGLETONS.fq.gz outm=T4_10_IonXpress_007_SINGLETONS_mouse.fq.gz outu=T4_10_IonXpress_007_SINGLETONS_clean.fq.gz nodisk 
rm T4_10_IonXpress_007_SINGLETONS.fq.gz
module load fastqc
fastqc T4_10_IonXpress_007_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T4_36_IonXpress_008.fastq.gz ftl=10 ftr=244 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T4_36_IonXpress_008_SINGLETONS.fq.gz outm=T4_36_IonXpress_008_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T4_36_IonXpress_008_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T4_36_IonXpress_008_SINGLETONS.fq.gz outm=T4_36_IonXpress_008_SINGLETONS_mouse.fq.gz outu=T4_36_IonXpress_008_SINGLETONS_clean.fq.gz nodisk 
rm T4_36_IonXpress_008_SINGLETONS.fq.gz
module load fastqc
fastqc T4_36_IonXpress_008_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T4_38_IonXpress_009.fastq.gz ftl=10 ftr=224 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T4_38_IonXpress_009_SINGLETONS.fq.gz outm=T4_38_IonXpress_009_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T4_38_IonXpress_009_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T4_38_IonXpress_009_SINGLETONS.fq.gz outm=T4_38_IonXpress_009_SINGLETONS_mouse.fq.gz outu=T4_38_IonXpress_009_SINGLETONS_clean.fq.gz nodisk 
rm T4_38_IonXpress_009_SINGLETONS.fq.gz
module load fastqc
fastqc T4_38_IonXpress_009_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T4_52_IonXpress_010.fastq.gz ftl=10 ftr=229 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T4_52_IonXpress_010_SINGLETONS.fq.gz outm=T4_52_IonXpress_010_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T4_52_IonXpress_010_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T4_52_IonXpress_010_SINGLETONS.fq.gz outm=T4_52_IonXpress_010_SINGLETONS_mouse.fq.gz outu=T4_52_IonXpress_010_SINGLETONS_clean.fq.gz nodisk 
rm T4_52_IonXpress_010_SINGLETONS.fq.gz
module load fastqc
fastqc T4_52_IonXpress_010_SINGLETONS_clean.fq.gz

/fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=T4_6_IonXpress_006.fastq.gz ftl=10 ftr=249 out=stdout.fq| /fs/project/PAS1117/bioinformatic_tools/bbmap/bbduk.sh -Xmx1G threads=1 overwrite=t in=stdin.fq fastawrap=10000 out=T4_6_IonXpress_006_SINGLETONS.fq.gz outm=T4_6_IonXpress_006_qc_failed.fq.gz minlength=30 qtrim=rl maq=10 maxns=0 overwrite=t stats=T4_6_IonXpress_006_qualTrimming.stats statscolumns=5 trimq=10 
/fs/project/PAS1117/bioinformatic_tools/bbmap/bbmap.sh -Xmx48G threads=12 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 overwrite=t ref=/fs/project/PAS1117/SCI/Mouse_genome/GCF_000001635.26_GRCm38.p6_genomic.fna in=T4_6_IonXpress_006_SINGLETONS.fq.gz outm=T4_6_IonXpress_006_SINGLETONS_mouse.fq.gz outu=T4_6_IonXpress_006_SINGLETONS_clean.fq.gz nodisk 
rm T4_6_IonXpress_006_SINGLETONS.fq.gz
module load fastqc
fastqc T4_6_IonXpress_006_SINGLETONS_clean.fq.gz

mv *_SINGLETONS_clean.fq.gz /fs/project/PAS1117/SCI/01_Clean_reads
mv *.html /fs/project/PAS1117/SCI/00_Raw_reads/QC_results
mv *.zip /fs/project/PAS1117/SCI/00_Raw_reads/QC_results

