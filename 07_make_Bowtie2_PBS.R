setwd("/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/data/raw_fastq/")
wd <- getwd()
outputd <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/BAM_files"
pattern <- "(.*).fq.gz"
###################################################
library(gsubfn)
#read_files <- list.files(pattern = "\\21_F.fq.gz$")
read_files <- Sys.glob("*21_F.fq.gz")

for (i in 1:length(read_files)){
   sample <- strapplyc(read_files[i], pattern, simplify = TRUE)
  BBKW_Script <- paste("#PBS -N 7_Bowtie2_", sample, "\n",
    "#PBS -m bea
    #PBS -j oe
    #PBS -A PAA0034
    #PBS -l nodes=1:ppn=28
    #PBS -l walltime=168:0:0 

    module load samtools/1.3.1
    cd ", wd, "\n",
    "/fs/project/PAS1117/bioinformatic_tools/bowtie2-2.3.3.1-linux-x86_64/bowtie2 --no-unal --non-deterministic -p 28 -X 2000 -x " , outputd, "/SCI_viral_populations -U ", read_files[i] , " | samtools view -@ 28 -Sb - > ", sample, ".bam" , sep="")
  write(BBKW_Script, file = paste("07_Bowtie2_",sample,".sh", sep = ""))
}
