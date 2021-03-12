# This script combines all the results from MARVEL, VirSorter, DeepVirFinder and CAT in one big table.
# The table will be written to viral_hhm/05_combined_results/combined_results.csv
### Libraries ##################################################################
library(data.table)
library(dplyr)
library(stringr)
library(seqinr)


### Directories ################################################################
proj_dir <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/"
setwd(proj_dir)

CAT_dir <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/CAT_results/"
DeepVirFinder_dir <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/DVF_results/"
MARVEl_outfiles_dir <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/MARVEL_Results/out_files/"
VirSorter_dir <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/SCI_ion_all_coassembly/"


### Grabbing results from DeepVirFinder ########################################
dfv_files <- as.vector(list.files(DeepVirFinder_dir, full.names = T))

dfv_DT <- data.table(name = character(),
                     score = numeric(),
                     pvalue = numeric())

for(i in 1:length(dfv_files)){
  local_DT <- fread(input = dfv_files[i])
  local_DT <- select(local_DT, name, score, pvalue)
  dfv_DT <- rbind(dfv_DT, local_DT)
  rm(local_DT)
}

contigs_DT <- dfv_DT
rm(dfv_files, dfv_DT)



### Grabbing results from MARVEL ###############################################
MARVEL_files <- list.files(MARVEl_outfiles_dir, full.names = T)
# MARV_map <- fread(input = "~/viral_hhm/marv_names_mapped.csv")        # here the mapped names are stored. 
outfiles_curated <- "/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/viral_results/MARVEL_Results/curated_out_files/"

# Curating the output
for(i in 1:length(MARVEL_files)){
  outfile <- as.data.table(readLines(con = MARVEL_files[i]))
  outfile$stars <- lapply(outfile$V1, str_detect, pattern = "\\*\\*\\*")
  outfile <- outfile[outfile$stars == TRUE]
  outfile <- select(outfile, V1)
  
  fwrite(outfile, file = paste0(outfiles_curated, "/", i, ".csv"), sep = " ")
}
rm(outfile, i, MARVEL_files)

# Reading the curated output
curated_files <- list.files(outfiles_curated, full.names = T)
file_list <- lapply(curated_files, fread)
marv_DT <- rbindlist(file_list, fill = T)
rm(file_list)

# some data tweaking, matching names with the mapped names
marv_DT <- select(marv_DT, V2, V4)
# MARV_map$marv_ID <- str_replace(string = MARV_map$marv_ID,
#                                pattern = ".fa",
#                                replacement = "")

# changing the names
# ChangeNames <- function(marv_name){
#  original_ID <- MARV_map[marv_ID == marv_name]$original_ID
#  return(original_ID)
# }
# this takes a few minutes
# marv_DT$V2 <- lapply(marv_DT$V2, ChangeNames)
marv_DT$V2 <- lapply(marv_DT$V2, unlist)

# convert all to data.table
marv_DT <- as.data.table(marv_DT)
contigs_DT <- as.data.table(contigs_DT)

# add new column
contigs_DT$marv_percentage <- NA

#change the name of columns in marv_DT
names(marv_DT) <- c("name", "marv_percentage")
marv_DT$name <- unlist(marv_DT$name)
marv_DT <- select(marv_DT, name, marv_percentage)

# join the two DTs
contigs_DT <- merge(x = contigs_DT, y = marv_DT, by = "name", all.x = T)

# select only relevant cols
contigs_DT <- select(contigs_DT, name, score, pvalue, marv_percentage.y)

# rename the cols
names(contigs_DT) <- c("ID", "dfv_score", "dfv_pvalue", "marv_percentage")

# save contigs_DT and remove the other DTs
fwrite(x = contigs_DT, file = paste0(proj_dir, "/contigs_DT.csv"), sep = "\t")
rm(marv_DT, MARV_map, curated_files, marv_name, MARVEl_outfiles_dir, outfiles_curated)



### Grabbing results from CAT ##################################################
CAT_subdirs <- list.dirs(CAT_dir, recursive = F)
CAT_subdirs_short <- list.dirs(CAT_dir, recursive = F, full.names = F)
CAT_files_contigs_class <- paste0(CAT_subdirs, "/", CAT_subdirs_short, "_contigs_classification.txt")

# reading contigs_classification.txt-files
CAT_ls <- lapply(CAT_files_contigs_class, fread, fill = TRUE)
CAT_DT <- rbindlist(CAT_ls)
CAT_DT <- select(CAT_DT, "#Contig_ID", No_genes, No_annotated_genes,
                 Contribution_of_annotated_ORFs, "taxid_superkingdom:contribution")
names(CAT_DT) <- c("ID", "Genes", "Annotated_Genes", "Contribution_of_annotated_ORFs", "Kingdom")

# calculation of number_of_genes per contig
# CAT_faa_files <- paste0(CAT_subdirs, "/", CAT_subdirs_short, "_proteins.faa")
CAT_faa_files <- paste0("/fs/project/PAS1117/SCI/SCI_ion_all_coassembly/results", "/", "03.SCI_ion_all_coassembly.faa")

# APPROACH 5
fasta_files <- lapply(CAT_faa_files, read.fasta)
fasta_names <- lapply(fasta_files, getName)
fasta_names <- lapply(fasta_names, as.data.table)
rm(fasta_files)

random_ls <- list()

for(i in 1:length(fasta_names)){
  # select only viable contigs (reducing the number of searches)
  local_genes <- fasta_names[[i]]
  local_genes_strip <- sub('[_][^_]+$', '', local_genes$V1)
  local_contigs <- filter(CAT_DT, str_detect(CAT_DT$ID, pattern = "NODE")) %>%
    select(ID)
  # local_contigs <- filter(CAT_DT, str_detect(CAT_DT$ID, pattern = CAT_subdirs_short[i])) %>%
  #  select(ID)
  # create local function
#  find_genes <- function(ID){
#    no_of_genes <- filter(local_genes, str_detect(string = local_genes$V1,
#                                                  pattern = ID)) %>%
#      nrow()
#    return(no_of_genes)
#  }
  gene_frequency <- as.data.table(as.data.frame(table(local_genes_strip)))
  names(gene_frequency) <- c("ID", "no_of_genes")
  #local_contigs$no_of_genes <- lapply(local_contigs$ID, find_genes)
  local_contigs <- merge(x = local_contigs, gene_frequency, by = "ID", all.x = T)
  random_ls[[i]] <- local_contigs
}

random_DT <- rbindlist(random_ls)

# add info to contigs_DT
copy_DT <- contigs_DT
contigs_DT <- merge(x = contigs_DT, random_DT, by = "ID", all.x = T)
setnames(contigs_DT, old = "no_of_genes", new = "CAT_genes")
contigs_DT$CAT_genes <- as.numeric(as.character(contigs_DT$CAT_genes)) # converts NULL into NA
fwrite(x = contigs_DT, file = paste0(proj_dir, "/contigs_DT.csv"), sep = "\t")

# Add CAT_Taxonomy to contigs_DT
tax <- select(CAT_DT, ID, Kingdom)
contigs_DT <- merge(contigs_DT, tax, by = "ID", all.x = T)
setnames(contigs_DT, old = "Kingdom", new = "CAT_Kingdom")
rm(tax)

# Add Contribution of annotated ORFs
orfs <- select(CAT_DT, ID, "Contribution_of_annotated_ORFs")
contigs_DT <- merge(contigs_DT, orfs, by = "ID", all.x = T)
setnames(contigs_DT, old = "Contribution_of_annotated_ORFs", new = "CAT_annot_ORFs")
setnames(contigs_DT, old = "CAT_annot_ORFs", new = "CAT_contributing_genes")


# Add "Contribution_percentage"
contigs_DT$CAT_Contribution_Percentage <- str_extract(string = contigs_DT$CAT_Kingdom,
                                                      pattern = "[\\d\\.]{1,}$")
contigs_DT$CAT_Contribution_Percentage <- as.numeric(contigs_DT$CAT_Contribution_Percentage)

# Calculate and add "No_Annotated_genes"
contigs_DT$CAT_annotated_genes <- contigs_DT$CAT_contributing_genes/contigs_DT$CAT_Contribution_Percentage*100
fwrite(x = contigs_DT, file = paste0(proj_dir, "/contigs_DT.csv"), sep = "\t")

# Some redoing of the percentage, since the current perce is based on a bitscore
# the calc of the annot_genes is therefore wrong
contigs_DT$CAT_my_percentage <- contigs_DT$CAT_contributing_genes/contigs_DT$CAT_genes*100
contigs_DT$CAT_annotated_genes <- NA


### Grabbing results from VirSorter ############################################
VS_dirs <- list.dirs(VirSorter_dir, recursive = F, full.names = F)
VS_files <- paste0(VirSorter_dir, VS_dirs, "/", VS_dirs,
                   "_global-phage-signal_with-prophage-coordinates.csv")

# VS_cat1_files <- paste0(VirSorter_dir, "/", VS_dirs, "/", VS_dirs,
#                        "_cat-1.fasta")
# VS_cat2_files <- paste0(VirSorter_dir, "/", VS_dirs, "/", VS_dirs,
#                        "_cat-2.fasta")
# VS_cat3_files <- paste0(VirSorter_dir, "/", VS_dirs, "/", VS_dirs,
#                        "_cat-3.fasta")
# VS_cat4_files <- paste0(VirSorter_dir, "/", VS_dirs, "/", VS_dirs,
#                        "_prophages_cat-4.fasta")
# VS_cat5_files <- paste0(VirSorter_dir, "/", VS_dirs, "/", VS_dirs,
#                        "_prophages_cat-5.fasta")
# VS_cat6_files <- paste0(VirSorter_dir, "/", VS_dirs, "/", VS_dirs,
#                        "_prophages_cat-6.fasta")

VS_cat1_files <- paste0(VirSorter_dir, VS_dirs,"/", VS_dirs,
                        "_cat-1.fasta")
VS_cat2_files <- paste0(VirSorter_dir, VS_dirs,"/", VS_dirs,
                        "_cat-2.fasta")
VS_cat3_files <- paste0(VirSorter_dir, VS_dirs,"/", VS_dirs,
                        "_cat-3.fasta")
VS_cat4_files <- paste0(VirSorter_dir, VS_dirs,"/", VS_dirs,
                        "_prophages_cat-4.fasta")
VS_cat5_files <- paste0(VirSorter_dir, VS_dirs,"/", VS_dirs,
                        "_prophages_cat-5.fasta")
VS_cat6_files <- paste0(VirSorter_dir, VS_dirs,"/", VS_dirs,
                        "_prophages_cat-6.fasta")

VS1_ls <- list()
VS2_ls <- list()
VS3_ls <- list()
VS4_ls <- list()
VS5_ls <- list()
VS6_ls <- list()

for(i in 1:length(VS_dirs)){
  
  if(file.size(VS_cat1_files[i]) > 0){
    VS1_ls[[i]] <- fread(VS_cat1_files[i], header = F) %>% 
      filter(str_detect(string = V1, pattern = ">"))    
  }else{VS1_ls[[i]] <- NA}
  
  if(file.size(VS_cat2_files[i]) > 0){
  VS2_ls[[i]] <- fread(VS_cat2_files[i], header = F) %>% 
    filter(str_detect(string = V1, pattern = ">"))
  }else{VS2_ls[[i]] <- NA}
  
  if(file.size(VS_cat3_files[i]) > 0){
  VS3_ls[[i]] <- fread(VS_cat3_files[i], header = F) %>% 
    filter(str_detect(string = V1, pattern = ">"))
  }else{VS3_ls[[i]] <- NA}
  
  if(file.size(VS_cat4_files[i]) > 0){
  VS4_ls[[i]] <- fread(VS_cat4_files[i], header = F) %>% 
    filter(str_detect(string = V1, pattern = ">"))
  }else{VS4_ls[[i]] <- NA}
  
  if(file.size(VS_cat5_files[i]) > 0){
  VS5_ls[[i]] <- fread(VS_cat5_files[i], header = F) %>% 
    filter(str_detect(string = V1, pattern = ">"))
  }else{VS5_ls[[i]] <- NA}
  
  if(file.size(VS_cat6_files[i]) > 0){
  VS6_ls[[i]] <- fread(VS_cat6_files[i], header = F) %>% 
    filter(str_detect(string = V1, pattern = ">"))
  }else{VS6_ls[[i]] <- NA}
  
  print(paste0("Worked for ", i))
}

na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

VS1_ls <- na.omit.list(VS1_ls)
VS1_DT <- rbindlist(VS1_ls, fill = T)

VS2_ls <- na.omit.list(VS2_ls)
VS2_DT <- rbindlist(VS2_ls, fill = T)

VS3_ls <- na.omit.list(VS3_ls)
VS3_DT <- rbindlist(VS3_ls, fill = T)

VS4_ls <- na.omit.list(VS4_ls)
VS4_DT <- rbindlist(VS4_ls, fill = T)

VS5_ls <- na.omit.list(VS5_ls)
VS5_DT <- rbindlist(VS5_ls, fill = T)

VS6_ls <- na.omit.list(VS6_ls)
VS6_DT <- rbindlist(VS6_ls, fill = T)

VS_DT <- rbind(VS1_DT, VS2_DT, VS3_DT, VS4_DT, VS5_DT, VS6_DT)
VS_DT$Category <- as.numeric(str_extract(VS_DT$V1, pattern = "\\d$"))

names(VS_DT) <- c("ID", "VS_Category")
#VS_DT$ID <- str_replace(VS_DT$ID, ">", "")
VS_DT$ID <- str_replace(VS_DT$ID, ">SCI_ion_all_coassembly_", "")
VS_DT$ID <- str_replace(VS_DT$ID, "-cat_", "")
VS_DT$ID <- str_replace(VS_DT$ID, "\\d$", "")
VS_DT$ID <- str_replace(VS_DT$ID, "_gene.{1,}$", "")
VS_DT$ID <- str_replace(VS_DT$ID, "-circular", "")
VS_DT$ID <- str_replace(VS_DT$ID, "-\\d{1,100}to\\d{1,100}-prophage$", "")



rm(VS1_DT, VS1_ls, VS2_DT, VS2_ls, VS3_DT, VS3_ls,VS5_DT, VS5_ls, VS4_ls, VS4_DT, VS6_DT, VS6_ls)
rm(VS_cat1_files, VS_cat2_files, VS_cat3_files, VS_cat4_files, VS_cat5_files, VS_cat6_files)

# Change ID in Category from "." to "_"
contigs_DT$ID <- str_replace_all(string = contigs_DT$ID,
                                 pattern = "\\.",
                                 replacement = "_")
# Add Category
contigs_DT <- merge(x = contigs_DT, y = VS_DT, by = "ID", all.x = T)
contigs_DT <- select(contigs_DT, "ID", "dfv_score", "dfv_pvalue", "marv_percentage",
                     "CAT_genes", "CAT_Kingdom", "CAT_Contribution_Percentage",
                     "CAT_contributing_genes", "CAT_my_percentage", "VS_Category")


fwrite(x = contigs_DT, file = paste0(proj_dir, "/contigs_DT.csv"), sep = ",")

