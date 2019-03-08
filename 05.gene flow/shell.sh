#generate the genotype posterior probability. Only positions sequenced at least once in each individual were considered
#!/bin/sh
export list=$1
export chr=$2
export mindepth=$3  #600
export maxdepth=$4  #2500
angsd -bam $list -only_proper_pairs 0 -uniqueOnly 1 -remove_bads 1 -minQ 20 -minMapQ 25 -ref $reference.fa -r $chr -out $chr -doMaf 1 -minInd $num -skipTriallelic 1 -doMajorMinor 1 -GL 1 -setMinDepth $mindepth -setMaxDepth $maxdepth -doCounts 1 -P 4 -doGlf 2 -SNP_pval 1e-6 -setMinDepthInd 1

#randomly sampled genotypes according to their posterior probabilities.
python beagle2treemix.py --beagle $beagle.gz --outfile $treemix.gz

#onverts beagle2treemix data to the format used by ADMIXTOOLS
perl 04.treemix2geno.pl <beagle2treemix output> <outputprefix.geno>

###########inferring the introgressed segments#####
#sliding window (non-overlapping) D-statistic
export size=$1
export bam=$2
export chr=$3
export out=$4
export maxdepth=$5
export blocksize=$6
angsd -doAbbababa2 1 -sizeFile $size -remove_bads 1 -only_proper_pairs 1 -rmTrans 0 -bam $bam -doCounts 1 -out $out -anc $outgroup.fa -useLast 0 -minQ 20 -minMapQ 25 -p 4 -maxDepth $maxdepth -uniqueOnly 1 -r $chr -blockSize $blocksize

#calculate the pairwise IBS distances
#!/bin/sh
export bam=$1
export region=$2
export out=$3
angsd -bam $bam -minMapQ 20 -minQ 20 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -doIBS 1 -doCounts 1 -makeMatrix 1 -minMaf 0.05 -P 4 -rf $region -out $out

#inferring the introgressed alleles. The 24 modern bezoars are used as outgroup. 
export chr=$1
export start=$2
export end=$3
export out=$4
java -jar sprime.jar gt=$phased.vcf.gz outgroup=$outgroup_ID map=$chr.map out=$out maxfreq=0.01 mu=4.32e-9 minscore=150000 chrom=$chr:$start-$end excludesamples=$exclude_ID
