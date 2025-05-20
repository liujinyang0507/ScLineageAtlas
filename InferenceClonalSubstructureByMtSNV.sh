#!/bin/bash

input_bam=$1
barcode=$2
sample=$3
out_dir=$4
patient=$5

#conda install -c bioconda cellsnp-lite
#pip install -U mquad -i http://pypi.douban.com/simple/ --trusted-host pypi.douban.com
#pip install -U vireoSNP -i http://pypi.douban.com/simple/ --trusted-host pypi.douban.com

### SMART-SEQ2 ####
#for sequencing data with individual .bam files for each cell + no UMItags/barcodes

#run cellsnp-lite mode2a on bam list
#change --chrom= to whatever reference genome you aligned to - in #this case we use hg19
#ls *.bam > bam.lst
#cellsnp-lite -S bam.lst -i bam.lst -o cellsnp --UMItag None #--genotype --gzip --chrom=chrM -p 10

### 10X/UMI-BASED ###
#for 10x and UMI-based sequencing data where there is only 1 big .bam + barcodes run cellsnp-lite on the .bam directly
#change --chrom= to whatever reference genome you aligned to, in most cases 10x data #are aligned to GrCh38 so the chr name is MT

echo "barcode:$(realpath $barcode)"
echo "******$sample analysis start******"


if [ ! -d  ${out_dir} ];then
    mkdir $out_dir
fi

if [ ! -d ${out_dir}/out_cellsnp ];then
    cellsnp-lite -s $input_bam -b $barcode -O ${out_dir}/out_cellsnp --chrom=M --UMItag Auto  --genotype --gzip -p 20
    if [ $? -eq 1 ];then
        echo "The Cellsnp-lite tool run error, please check if the $barcode is empty, exit 1."
        exit 2
    fi
fi

n_vars_cellsnp=$(less ${out_dir}/out_cellsnp/cellSNP.tag.DP.mtx|awk 'NR==3{print $3}')

if [ "a$n_vars_cellsnp" = "a" ];then
    echo "The Cellsnp-lite tool was unable to detect any variants, exit 1."
    exit 1
fi

if [ $n_vars_cellsnp -eq 0 ];then
    echo "The Cellsnp-lite tool was unable to detect any variants, exit 1."
    exit 1
fi

if [ ! -d ${out_dir}/out_mquad ];then
    mquad -c ${out_dir}/out_cellsnp -o ${out_dir}/out_mquad -p 20 --minDP 1
    if [ $? -eq 1 ];then
        echo "The Mquad tool was unable to identify the knee or elbow point, exit 2."
        exit 2
    fi
fi

if [ ! -f ${out_dir}/out_mquad/passed_dp.mtx ];then
    echo "The Mquad tool had run, but it was not unable to identify the knee or elbow point, exit 2."
    exit 2
fi

n_vars_mquad=$(less ${out_dir}/out_mquad/passed_dp.mtx|awk 'NR==3{print $3}')

if [ "a$n_vars_mquad" = "a" ];then
    echo "The Mquad tool was unable to detect any variants, exit 2."
    exit 2
fi

if [ $n_vars_mquad -eq 0 ];then
    echo "The Mquad tool was unable to detect any variants with clone information, exit 2."
    exit 2
fi
variants=$(less ${out_dir}/out_mquad/passed_variant_names.txt|awk '{print NR}'|tail -n1)
#Identifitation the optimal number of clones by setting the ELBO growth rate cutoff value to 0.04
echo "$sample start vireoSNP analysis"
if [ ! -f ${out_dir}/out_vireoSNP/${sample}_cell_assignment_probability_to_each_donor.png ];then
    python3.8 /home/data/CellLineage/code/VireoSNPNumDonorCutoff.py -ad ${out_dir}/out_mquad/passed_ad.mtx -dp ${out_dir}/out_mquad/passed_dp.mtx\
    -gr 0.08 -mc $variants -mi 30 -ni 30 -o ${out_dir}/out_vireoSNP -r $sample
    
    n_clone=$(awk 'NR>1{print $0}' ${out_dir}/out_vireoSNP/${sample}_the_optimal_number_of_clones.csv)
    python3.8 /home/data/CellLineage/code/IdentificationAndVisualizationOfClonesUsingVireoSNP.py -ad ${out_dir}/out_mquad/passed_ad.mtx -dp ${out_dir}/out_mquad/passed_dp.mtx\
    -var ${out_dir}/out_mquad/passed_variant_names.txt -b $barcode -n $n_clone -mi 30 -ni 300 -o ${out_dir}/out_vireoSNP -r $sample
    
    if [ $? -eq 1 ];then
        echo "IdentificationAndVisualizationOfClonesUsingVireoSNP.py run error, exit 3."
        exit 3
    fi

   # echo $sample >>  downstream_analysis_sample_list.txt
    #echo $patient >> success_sample.txt
    
fi
    
mkdir ${out_dir}/out_diff_analysis
