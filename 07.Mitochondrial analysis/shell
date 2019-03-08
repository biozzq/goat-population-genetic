#alignment the clean data to the modified mitochondrial geneome (Given that mitochondrial genomes are circular, we added 300 of the first base pairs to the end of the reference to assure equal coverage of the sequences across the mtDNA).
export sample=$1
export prefix=$2
export fastqgz1=$3
export fastqgz2=$4
bwa aln -t 4 -f $prefix\_1.sai $modified_reference.fa $fastqgz1
bwa aln -t 4 -f $prefix\_2.sai $modified_reference.fa $fastqgz2
bwa sampe -r '@RG\tID:'$prefix'\tLB:'$prefix'\tPL:ILLUMINA\tSM:'$sample'' $modified_reference.fa $prefix\_1.sai $prefix\_2.sai $fastqgz1 $fastqgz2 | samtools view -F4 -b -S | samtools sort -@ 4 -o $prefix.sort.bam
java -Xmx50g -jar picard.jar MarkDuplicates INPUT=$prefix.sort.bam OUTPUT=$prefix.sort.dedup.bam METRICS_FILE=$prefix\_dedup REMOVE_DUPLICATES=true CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES=2000

#assembly using MIA
#!/bin/sh
if [ $# -ne 2 ];then
 echo "error.. need args"
 echo "command:$0 <bam> <outprefix>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export bam=$1
export out=$2
samtools view -q 20 -b -f 2 -F 1024 $bam | java -jar picard.jar SamToFastq INPUT=/dev/stdin FASTQ=$out.fq INCLUDE_NON_PRIMARY_ALIGNMENTS=false INTERLEAVE=true VALIDATION_STRINGENCY=LENIENT
mia -H 1 -F -i -c -r $reference.fa -f $out.fq -m $out.result
fi

#generate the final results
ma -M $out.result.iteration -f 5 -I $sample >$sample.fa

#multiple sequence alignment
mafft --thread 10 --globalpair --maxiterate 16 --inputorder $merge.fa > $merge.aln
