#!/usr/bin/env bash
# running the commannd: submit_plates.sh will start the whole pipeline
# this script creates new job for running gelSeq.py code for each plate
#### please change the platesrc file before running this script!!! ###
source platesrc.cfg
for dir in $plates_dirs;
do
    for fastq1 in ${dir}*R1_001.fastq.gz;
    do
        fastq2=${fastq1%%R1_*.fastq.gz}R2_001.fastq.gz
        base_name="$(basename $fastq1)"
        base_name=${base_name%%_S*}
        base_name=${base_name%%B}
        #base_name="$(echo ${base_name} | sed 's/[^0-9]//g')"
        plate_name="$(grep ^${base_name} ${table}  | cut -d',' -f2)"
        echo ${plate_name}
        if [ "${plate_name}" != "" ]
        then
            out_path=${output_dir}/${plate_name}.log
            error_path=${output_dir}/${plate_name}.err
            bsub -q new-short -R rusage[mem=1000] -o ${out_path} -e ${error_path} python3.5 ${code_dir}/gelSeq.py plate ${fastq1} ${fastq2} ${plate_name} ${output_dir} --loci=B
        fi

    done
done
