setwd("/fs/project/PAS1117/SCI/01_Clean_reads/")
wd <- getwd()
outputd <- "/fs/project/PAS1117/SCI/02_Assemblies"
pattern <- "(.*)_IonXpress"
###################################################
library(gsubfn)
read_files <- list.files()

for (i in 1:length(read_files)){
   sample <- strapplyc(read_files[i], pattern, simplify = TRUE)
  assembly_Script <- paste("#PBS -N 2_assembly_", sample, "\n",
    "#PBS -m bea
    #PBS -j oe
    #PBS -A PAA0034
    #PBS -l nodes=1:ppn=28
    #PBS -l mem=124GB
    #PBS -l walltime=100:0:0 
    cd ", wd, "\n",
    "module load bowtie2
    module load samtools
    module load bamtools 
    
    mkdir ", sample, "_assembly 
    /fs/project/PAS1117/bioinformatic_tools/SPAdes-3.11.1-Linux/bin/spades.py -t 28 -m 124 --iontorrent -s ", read_files[i] , " -o ", sample, "_assembly
    cd ", sample, "_assembly
    mv ", "scaffolds.fasta ", outputd, "/" , sample, "_scaffolds.fasta", "\n",
    "mv ", "contigs.fasta ", outputd, "/" , sample, "_contigs.fasta"
    ,sep="")
  write(assembly_Script, file = paste("02_Assemble_",sample,".sh", sep = ""))
}
