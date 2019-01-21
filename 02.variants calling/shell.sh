#read alignment, sort, and remove PCR duplication
bwa mem -t 4 -M -R '@RG\tID:$group\tLB:$group\tPL:ILLUMINA\tSM:$sample' $reference.fa $fastq1 $fastq2 | samtools view -b -S -o $sample.bam
java -jar picard.jar SortSam I=$sample.bam O=$sample.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
java -jar picard.jar MarkDuplicates I=$sample.sort.bam O=$sample.dup.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=$sample.dup.txt VALIDATION_STRINGENCY=LENIENT
java -jar picard.jar SetNmMdAndUqTags I=$sample.dup.bam O=$sample.dup.sort.fix.bam CREATE_INDEX=true R=$reference.fa VALIDATION_STRINGENCY=LENIENT

#variants calling using ANGSD
angsd -bam $bam.list -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -minQ 20 -minMapQ 30 -C 50 -ref $reference.fa -r $chr -out $chr -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth $max

angsd -bam $bam.list -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 -remove_bads 1 -minQ 20 -minMapQ 30 -C 50 -ref $reference.fa -r $chr -out $chr -doMaf 1 -doMajorMinor 1 -GL 1 -setMinDepth $mindepth -setMaxDepth $maxdepth -doCounts 1 -dosnpstat 1 -SNP_pval 1

#variants calling using GATK
java -jar GenomeAnalysisTK.jar -R $reference.fa -T HaplotypeCaller -L $chr -ERC GVCF -I $bam -o $chr.g.vcf.gz
java -jar GenomeAnalysisTK.jar -R $reference.fa -T CombineGVCFs --variant $chr.list -o $chr.gvcf.gz
java -jar GenomeAnalysisTK.jar -R $reference.fa -T GenotypeGVCFs --variant $chr.gvcf.gz --includeNonVariantSites -o $chr.vcf.gz
