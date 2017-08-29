#!/usr/bin/env bash

dir=/home/labs/amit/diklag/BaseCalls/
cell_barcodes=/home/labs/amit/diklag/BaseCalls/5_CBC PCR1-F_well_IDs_V1.csv
B_primer=AGCGACCTCGGGTGGGAAC
A_primer=GTACACGGCAGGGTCAGG


OUTPUT_FILE=${dir}/global_table.csv
temp_file=${dir}/temp.csv
echo I1,I2,R1,plate > ${OUTPUT_FILE} # Copy header if it is the first file
for plate in "$(ls ${dir}*fastq.gz | cut -d'_' -f1,2,3 | sort | uniq | paste -s -d' ')";
do
    echo ${plate}
    #I1="$(gunzip -c ${plate}_I1_001.fastq.gz | awk 'NR%4==2')"
    #I2="$(gunzip -c ${plate}_I2_001.fastq.gz | awk 'NR%4==2')"
    #R1="$(gunzip -c ${plate}_R1_001.fastq.gz | awk 'NR%4==2')"
    #paste -d"," ${I1} ${I2} ${R1} > ${temp_file} # Append from the 2nd line each file
    #awk -v a="$plate" 'BEGIN { FS="," ; OFS="," } FNR > 0 { $4=a; print }' ${temp_file} >> ${OUTPUT_FILE}

done