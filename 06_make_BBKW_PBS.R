setwd("/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/03.BBKW_annotations/")
wd <- getwd()
outputd <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results/03.BBKW_annotations"
pattern <- "01.SCI_ion_all_coassembly_(.*).fasta"
###################################################
library(gsubfn)
read_files <- list.files(pattern = "\\.fasta$")

for (i in 1:length(read_files)){
   sample <- strapplyc(read_files[i], pattern, simplify = TRUE)
  BBKW_Script <- paste("#PBS -N 6_BBKW_", sample, "\n",
    "#PBS -m bea
    #PBS -j oe
    #PBS -A PAA0034
    #PBS -l nodes=1:ppn=12
    #PBS -l walltime=100:0:0 
    module use /fs/project/PAS1117/modulefiles
    module load WrightonPipeline
    cd ", wd, "\n",
    "mkdir ", outputd, "/", sample, "\n", 
    "Wrighton-Sullivan_Pipeline.py -i " , outputd, "/", read_files[i] , " -t nucl -m 500 -o ", outputd, "/", sample, " --tool pfam", sep="")
  write(BBKW_Script, file = paste("06_BBKW_",sample,".sh", sep = ""))
}
