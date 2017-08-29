#!/usr/bin/env bash
#python3.5 ../gelSeq.py plate /home/labs/amit/weiner/Work/HANJAY/MiSeq/20170504_MiSeq/fastq/104_S11_L001_R1_001.fastq.gz /home/labs/amit/weiner/Work/HANJAY/MiSeq/20170504_MiSeq/fastq/104_S11_L001_R2_001.fastq.gz 104_S11_L001 /home/labs/amit/diklag/PycharmProjects/VDJ_proj/plates --loci=B
#plate /home/labs/amit/weiner/Work/HANJAY/MiSeq/20170504_MiSeq/fastq/108_S12_L001_R1_001.fastq.gz /home/labs/amit/weiner/Work/HANJAY/MiSeq/20170504_MiSeq/fastq/108_S12_L002_R1_001.fastq.gz 108_S12_L001 /home/labs/amit/diklag/PycharmProjects/VDJ_proj/plates --loci=B

#plate /home/labs/amit/weiner/Work/HANJAY/MiSeq/20170504_MiSeq/fastq/85_S5_L001_R1_001.fastq.gz /home/labs/amit/weiner/Work/HANJAY/MiSeq/20170504_MiSeq/fastq/85_S5_L002_R1_001.fastq.gz AB2516 /home/labs/amit/diklag/PycharmProjects/VDJ_proj/plates --loci=B


#cell -s Hsap --loci=B /home/labs/amit/diklag/PycharmProjects/VDJ_proj/plates/AB2512/A3.fasta A3 /home/labs/amit/diklag/PycharmProjects/VDJ_proj/plates/AB2512
# todo - please change all the paths to your directories !!!!!!
for cell in /home/labs/amit/diklag/PycharmProjects/VDJ_proj/plates/109_S13_L001/*
do
    cell_name="$(basename $cell)"
    cell_name=${cell_name%%.fasta}
    echo $cell_name
    if [ -f ${cell}/reads/${cell_name}.fasta ];
    then
        fasta=${cell}/reads/${cell_name}.fasta
        out_path=${cell}/out.out
        error_path=${cell}/err.err
        if [ -f ${out_path} ];
        then
            rm ${out_path}
        fi
        if [ -f ${error_path} ];
        then
            rm ${error_path}
        fi
        echo $out_path
        echo $error_path
        bsub -q new-short -o ${out_path} -e ${error_path} "python3.5 /home/labs/amit/diklag/PycharmProjects/VDJ_proj/gelSeq.py cell -s Hsap --loci=B ${fasta} ${cell_name} ${cell}"
        fi
done