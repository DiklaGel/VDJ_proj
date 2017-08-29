#!/usr/bin/env bash
source platesrc.cfg

i=0
OUTPUT_FILE=${output_dir}/result.csv
for plate_dir in ${output_dir}/AB*;
do
    if [ -a ${plate_dir}/final_table.csv ]
    then
        ORIG_FILE=${plate_dir}/final_table.csv
        NEW_FILE=${plate_dir}/full_table.csv
        amp_batch="$(basename $plate_dir)"
        patient_name="$(grep ${amp_batch} ${table}  | cut -d',' -f4 |  sed 's/\s$//g' | sed 's/\s/_/g')"
        echo $amp_batch
        awk -v a="$amp_batch" -v p="$patient_name" 'BEGIN { FS="," ; OFS="," } FNR == 1 { $22="Amp Batch"; $23="Patient"; print } FNR > 1 { $22=a; $23=p; print }' ${ORIG_FILE} > ${NEW_FILE}
        if [[ $i -eq 0 ]] ; then
            head -1 ${NEW_FILE} > ${OUTPUT_FILE} # Copy header if it is the first file
        fi
        tail -n +2 ${NEW_FILE} >> ${OUTPUT_FILE} # Append from the 2nd line each file
        i=$(( $i + 1 ))                        # Increase the counter
    fi
done