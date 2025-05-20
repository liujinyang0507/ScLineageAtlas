#!/bin/bash

project=$1
sample_info_file=$2

echo $(find ./ -name "*cell_assignment_probability_to_each_donor.png" -type f|awk -F "/" '{for (i=1;i<=NF;i++){ if ($i~/SRR/)print $i}}')|tr " " "\n"|grep -v png>success_sample_tmp.txt
if [ "a$(cat success_sample_tmp.txt)" = "a" ];then
	#bam_path=$(tree -if |grep possorted_genome_bam.bam)
	#if [ "a$bam_path" = "a" ];then
#		echo "$project:Input bam does not exist. Please be prepared before starting the analysis task."
#		exit
#	fi

	for file in GSE*_obs.csv $sample_info_file *count.h5ad distribution.csv
	do
		if [ ! -f $file ];then
			echo "$project:Input $file does not exist. Please be prepared before starting the analysis task." 
			exit
		fi
	done
fi

awk -F "," 'FNR==NR{a[$1];next} ($1 in a){print $2}' success_sample_tmp.txt $sample_info_file > success_sample.txt
less $sample_info_file |awk -F "," '{print $2 }'|sort|uniq > Total_Sample.txt
awk -F "," 'FNR==NR{a[$1];next} !($1 in a)' success_sample.txt Total_Sample.txt > Missing_Samples_tmp.txt
awk -F "," 'FNR==NR{a[$1];next} ($2 in a)' Missing_Samples_tmp.txt $sample_info_file > Missing_Samples.txt
sed -i "1i run_accession,sample_title" Missing_Samples.txt

if [ -f project_check_result.log ];then
	rm -rf project_check_result.log
fi

if [ -s Missing_Samples.txt ];then
	echo "The sample analysis is incomplete, the missing sample information is as follows." >> project_check_result.log
	cat Missing_Samples.txt >> project_check_result.log
	echo 'You also view the samples details are missing from the "Missing_Samples.txt" file' >> project_check_result.log
else
	echo "The sample analysis of this project is complete." >> project_check_result.log
fi
echo "\n" >> project_check_result.log

for sample in $(cat success_sample_tmp.txt)
do	
	echo "**************************************************************************" >> project_check_result.log
	echo "$sample: Start of result verification." >> project_check_result.log
	echo "1. check barcode " >> project_check_result.log

	if [ ! -f ./$sample/${sample}_merge_barcodes.csv ];then
		echo " ./$sample/${sample}_merge_barcodes.csv \n\
		input barcode file does not exist, barcode check failed." >> project_check_result.log
	fi

	
	echo "2. checking directory out_vireoSNP, listing missing files." >> project_check_result.log
	for file in allele_frequency.csv allele_frequency_heatmap.png cell_assignment_probability_to_each_donor.csv cell_assignment_probability_to_each_donor.png\
		ELBO.csv finding_the_optimal_number_of_clones_using_the_ELBO_curve.png genotype_matrix.csv mean_allelic_ratio.csv mean_allelic_ratio.png\
		the_evidence_lower_bound_ELBO_values_and_the_number_of_clone_information.csv the_optimal_number_of_clones.csv
	do
		if [ ! -f ${sample}/mtSNV/out_vireoSNP/${sample}_${file} ];then
			echo "${sample}/mtSNV/out_vireoSNP/${sample}_${file} " >> project_check_result.log
		fi
	done

	echo "3. checking directory out_diff_analysis, listing missing files." >> project_check_result.log
	for file in cells_clones_infomation.csv cells_clones_updata_barcodes_infomation.csv DEG_among_clonetypes.csv\
	fraction_of_each_clone_of_each_cell_type.csv fraction_of_each_clone_of_each_cell_type.png GO.png GO_total.csv\
	GO_visualization.csv KEGG.png KEGG_total.csv KEGG_visualization.csv multiple_volcano.pdf multiple_volcano.png\
	umap_celltype_clone_information.csv umap_clones.html 
	do
		if [ ! -f ${sample}/mtSNV/out_diff_analysis/${sample}_${file} ];then
			echo "${sample}/mtSNV/out_diff_analysis/${sample}_${file}" >> project_check_result.log
		fi

	done
	
	#echo "4. checking directory out_scite, listing missing files." >> project_check_result.log
	#file="scite_ml0.newick"
	#if [ ! -f ${sample}/mtSNV/out_scite/${sample}_${file} ];then
#		echo "${sample}/mtSNV/out_scite/${sample}_${file}" >> project_check_result.log
#	fi
	echo "**************************************************************************" >> project_check_result.log
	echo "\n" >> project_check_result.log
done
