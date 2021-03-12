module use /fs/project/PAS1117/modulefiles
module load singularityImages
cd /fs/project/PAS1117/00_SCI_VconTACT_Jingjie/cluster_genomes/
singularity run /fs/project/PAS1117/modules/singularity/ClusterGenomes-1.1.0.img -f SCI_GVD_refseq_all_10kb.faa -c 80 -i 95 -o outputs_10kb 
