!/bin/bash

base=/path/to/some_folder

# Trimmomatic

for IN in [/path/to/raw_data_folder*_R1.fastq.gz];
    do name1=${IN/R1/R2};
    base=$(basename ${IN} _R1.fastq.gz);
    trimmomatic PE -phred33 ${IN} ${name1} ${base}_f_paired.fastq.gz ${base}_f_unpaired.fastq.gz ${base}_r_paired.fastq.gz ${base}_r_unpaired.fastq.gz ILLUMINACLIP:$TRIMH/adapters/TruSeq3-PE-2.fa:2:30:10 TRAILING:30 MINLEN:30;
done

# Quality control of the trimmed data

find $trimmomatic_data_folder/* -type f -print0 | xargs -0 -P 0 -I @ bash -c 'fastqc -o ./QC_trimmed_all_even @'

# bowtie2 transcriptome alignment

# Make the Index

BT2_HOME=/path/to/bowtie2_folder

$BT2_HOME/bowtie2-build /path/to/reference bowtie_transc_index

# Create a loop over all files

basedata=/path/to/trimmomatic_data_folder

for file in ${basedata}*_{1,2,3,4,5,6,7,8,9,10,11,12}_[1234]*f_paired.fastq.gz;
    do name=${file/f/r}; name2=$(basename -s _f_paired.fastq.gz ${file});
    $BT2_HOME/bowtie2 -x ${base}group_example/bowtie_transc_index -1 ${file} -2 ${name} |
    /path/to/samtools-1.10/samtools view -b - > ${name2}.bam ;
done

# Create a for loop to count how many reads align to each transcript

bamfiles=/path/to/bowtie_data_folder

for file in ${bamfiles}*_{1,2,3,4,5,6,10,11,12}_[13]*.bam;
    do name=$(basename ${file});
    /path/to/samtools-1.10/samtools view -q 35 -f 2 ${bamfiles}${name} | cut -f 3 | sort | uniq -c >
    /path/to/readcounts_from_bowtie2_bamfiles/${name/.bam/_readcount_q35_f2.txt};
done
