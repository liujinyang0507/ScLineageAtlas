#!/bin/bash
rm -rf run_sample_tmp.csv
for sample in $(ls|grep SRR)
do
	python /home/data/CellLineage/code/Get_Sample_Info_By_Barcodes_Obs.py $sample GSE1*_obs.csv ./$sample/barcodes.tsv ./${sample}/${sample}_sample.csv
	cat ./${sample}/${sample}_sample.csv>>run_sample_tmp.csv
done
