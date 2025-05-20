#!/bin/bash

obs_files=$1
sample_info_file=$2
hdf5_file=$3
umap_distribution_file=$4
analysis_folder=$5
cd $analysis_folder

if [ $# -eq 5 ];then
    echo "Check if necessary input files exist."
    for file in $obs_files $sample_info_file $hdf5_file $umap_distribution_file
    do
        if [ -f $file ];then
            echo "$file:$(realpath $file)" 
        else
            echo "$file file does not exists, program exit!"
            exit
        fi
    done

    #eval "$(conda shell.bash hook)"
    #conda activate CellLineage2 && echo "env activates successfully"


    echo Program Start Running $(date +%F%n%T)

    awk -F ',' -v OFS="," '{if(NR==1){print $0}else{if($1 !~ /-1$/){$1=$1"-1"}; print}}' $obs_files > new_obs_file_tmp.csv
    less new_obs_file_tmp.csv|awk -F "," '{print $1}'|awk 'BEGIN{print "barcodes"}NR>1{print $0}'> barcodes.csv
    paste -d "," new_obs_file_tmp.csv barcodes.csv > new_obs_file.csv

    less $umap_distribution_file|awk -F "," '{print $1","$2","$3","$4}' > umap.csv

    tree -if|grep bam|awk -F "->" '{print $1"\t"$2}'|awk '{if ($NF~/bam$/){print $0}}' > Bam_path.txt
    less Bam_path.txt|awk -F "/" '{for (i=1;i<=NF;i++){ if ($i~/SRR/)print $i}}'|sort|uniq > Bam_sample.txt
    awk -F "," 'FNR==NR{a[$1];next} ($1 in a)' Bam_sample.txt $sample_info_file > Bam_sample_merge.txt
    less Bam_sample_merge.txt |awk -F "," '!a[$2]++{print}' > merge_sample_info_file.txt
    sed -i "1i run_accession,sample_title" merge_sample_info_file.txt

    sample_info_list=$(less merge_sample_info_file.txt|awk -F "," 'NR>1{print $1}'|sort|uniq)
    for i in $sample_info_list;do less Bam_path.txt|grep $i;done > merge_bam_path.txt
    rm -rf new_obs_file_tmp.csv barcodes.csv Bam_path.txt Bam_sample.txt Bam_sample_merge.txt

    if [ -d out_diff_analysis ];then
        rm -rf out_diff_analysis/*
    fi

    less merge_bam_path.txt|while read link_bam_path true_bam_path
    do
        sample=$(echo $link_bam_path| awk -F "/" '{for (i=1;i<=NF;i++){ if ($i~/SRR/)print $i}}')
        patient=$(less merge_sample_info_file.txt|grep $sample|awk -F "," '{print $2}')
        dir_sample=$(dirname $(echo $link_bam_path))
        out_dir=${dir_sample}/mtSNV
        echo "out_dir:$(realpath $out_dir)"
        python /home/data/CellLineage/code/Merge_obs_barcode.py new_obs_file.csv ${dir_sample}/barcodes.tsv $patient $sample ${dir_sample}/${sample}_merge_barcodes.csv
        
        #rm -rf $out_dir/out_vireoSNP $out_dir/out_diff_analysis $out_dir/out_scite ${out_dir}/out_cellsnp ${out_dir}/out_mquad
        rm -rf $out_dir
        
        sh /home/data/CellLineage/code/InferenceClonalSubstructureByMtSNV.sh $true_bam_path ${dir_sample}/${sample}_merge_barcodes.csv $sample\
        $out_dir >${dir_sample}/${sample}_run_out.txt  2>${dir_sample}/${sample}_run_error.txt 

        mv $out_dir/out_vireoSNP/${sample}_cells_clones_infomation.csv ./out_diff_analysis
        
        less ./out_diff_analysis/${sample}_cells_clones_infomation.csv|awk -F "," -v p=$patient 'BEGIN{print "barcodes,clone_id"}NR>1{print p"_"$1","$2}' >\
        ./out_diff_analysis/${sample}_cells_clones_updata_barcodes_infomation.csv
    done

    less downstream_analysis_sample_list.txt|sort|uniq > uniq_downstream_analysis_sample_list.txt
    awk -F "," 'FNR==NR{a[$1];next} ($1 in a)' uniq_downstream_analysis_sample_list.txt $sample_info_file > downstream_analysis_sample_information.txt
    sed -i "1i run_accession,sample_title" downstream_analysis_sample_information.txt
    rm -rf downstream_analysis_sample_list.txt
    
    transfer_h5ad=$(echo $hdf5_file|awk -F "." '{print $1}')_transfer.h5ad
    if [ ! -f $transfer_h5ad ];then
        python3.8 /home/data/CellLineage/code/Transfer_CountH5ad.py $hdf5_file
    fi
    Rscript /home/data/CellLineage/code/nebula_clonetype.R $transfer_h5ad new_obs_file.csv  downstream_analysis_sample_information.txt ./out_diff_analysis ./out_diff_analysis 
    rm -rf $(echo $transfer_h5ad|cut -d "." -f1).h5seurat
    
    less downstream_analysis_sample_information.txt|awk -F "," 'NR>1{print $1"\t"$2}'|while read i j
    do
        Rscript /home/data/CellLineage/code/GO_KEGG_Aanlysis.R ./out_diff_analysis/${i}_DEG_among_clonetypes.csv ./out_diff_analysis $i
        python3.8 /home/data/CellLineage/code/Overview_Clone.py umap.csv out_diff_analysis/${i}_cells_clones_updata_barcodes_infomation.csv ./out_diff_analysis
        mv ./out_diff_analysis/${i}* ./${i}/mtSNV/out_diff_analysis/
    done
    echo Program End Running $(date +%F%n%T)
else
    echo "This script requires 5 parameters to run."
fi
