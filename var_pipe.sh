#! /usr/bin/bash

### run like ./var_pipe.sh srr.ids &

if [ $# -lt 2 ]; then
    echo "Need srr.ids and refGene_exons.bed files"
    exit 1
fi

#export FREEBAYES=$TOOLS/freebayes-v1.3.0-1
#export TRIMMOMATIC=/export/home/bwiley4/data/unit10/Trimmomatic-0.39/trimmomatic-0.39.jar

# get reads from SRA
parallel -a $1 fastq-dump --split-files --origfmt --gzip;

readarray -t mySRRs < $1

for id in "${mySRRs[@]}"; do
	echo "Processing ${id} **********************************************************************************";
	java -jar $TRIMMOMATIC PE ${id}_1.fastq.gz ${id}_2.fastq.gz \
		${id}_1.trim.fastq.gz S${id}_1un.trim.fastq.gz \
		${id}_2.trim.fastq.gz ${id}_2un.trim.fastq.gz \
		SLIDINGWINDOW:4:20 \
		MINLEN:20;


	bowtie2 -x hg19 -1 ${id}_1.trim.fastq.gz -2 ${id}_2.trim.fastq.gz --threads 32 --fast -S ${id}_paired.hg19.trim.sam;
	samtools sort ${id}_paired.hg19.trim.sam -o ${id}_paired.hg19.trim.sorted.bam -O bam;
	samtools index ${id}_paired.hg19.trim.sorted.bam;

	$FREEBAYES -f hg19_chr17.fa ${id}_paired.hg19.trim.sorted.bam > ${id}_paired.hg19.trim.vcf;

	vcffilter -f "DP > 10 & AF = 0.5" ${id}_paired.hg19.trim.vcf > ${id}_paired.hg19.trim.filtered.vcf;
	vcfannotate -b $2 \
		 -k REFGENE \
		 ${id}_paired.hg19.trim.filtered.vcf > ${id}_paired.hg19.trim.filtered.annot.vcf;

	echo ${id};
	count=$(bedtools intersect -a ${id}_paired.hg19.trim.filtered.vcf -b $2 | wc -l);	
	echo -e "${id}\t${count}" >> results.txt

done;
