##sliding window Fst analysis
#calculate per pop saf for each populatoin
export bamlist=$1
export out=$2
export mindepth=$3
export maxdepth=$4
export minind=$5
angsd -bam $bamlist -only_proper_pairs 1 -baq 1 -uniqueOnly 1 -remove_bads 1 -minQ 20 -minMapQ 30 -C 50 -ref $reference.fa -P 8 -r $chr -out $out -doMajorMinor 1 -GL 1 -doMaf 1 -doCounts 1 -doSaf 1 -skipTriallelic 1 -setMinDepth $mindepth -setMaxDepth $maxdepth -minInd $minind -anc $outgroup.fa

#calculate the 2dsfs prior
export idx1=$1
export idx2=$2
export out=$3
realSFS -P 8 $idx1 $idx2 >$out #To reduce memory usage, one can limit the number of sites being loaded into memory using the -nSites argument.

#performs a sliding-window analysis
if [ $# -ne 6 ]; then
 echo "error.. need args"
 echo "command:$0 <saf1.index> <saf2.index> <sfs> <winsize> <stepsize> <out>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export saf1_idx=$1
export saf2_idx=$2
export sfs=$3
export win=$4
export step=$5
export out=$6
realSFS fst index -P 8 $saf1_idx $saf2_idx -sfs $sfs -fstout $out -whichFst 1
realSFS fst stats2 $out.fst.idx -win $win -step $step -whichFST 1 -type 2 >$out.$win.$step.fst
fi

##sliding window Thetas and Tajima tests
#estimate the site allele frequency likelihood by chromosome
export bamlist=$1
export out=$2
export mindepth=$3
export maxdepth=$4
export minind=$5
angsd -bam $bamlist -only_proper_pairs 1 -baq 1 -uniqueOnly 1 -remove_bads 1 -minQ 20 -minMapQ 30 -C 50 -ref $reference.fa -P 8 -r $chr -out $out -doMajorMinor 1 -GL 1 -doMaf 1 -doCounts 1 -doSaf 1 -skipTriallelic 1 -setMinDepth $mindepth -setMaxDepth $maxdepth -minInd $minind -anc $outgroup.fa

#Obtain the maximum likelihood estimate of the SFS using the realSFS
export pop=$1
realSFS cat 1.saf.idx 2.saf.idx 3.saf.idx 4.saf.idx 5.saf.idx 6.saf.idx 7.saf.idx 8.saf.idx 9.saf.idx 10.saf.idx 11.saf.idx 12.saf.idx 13.saf.idx 14.saf.idx 15.saf.idx 16.saf.idx 17.saf.idx 18.saf.idx 19.saf.idx 20.saf.idx 21.saf.idx 22.saf.idx 23.saf.idx 24.saf.idx 25.saf.idx 26.saf.idx 27.saf.idx 28.saf.idx 29.saf.idx -outnames $pop
realSFS $pop.saf.idx -P 8 >$pop.sfs #To reduce memory usage, one can limit the number of sites being loaded into memory using the -nSites argument.


#Calculate the statistics(-doThetas) using the SFS as prior information (-pest) for each site
export bamlist=$1
export chr=$2
export sfs=$3
angsd -bam $bamlist -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -baq 1 -minQ 20 -minMapQ 30 -C 50 -ref $reference.fa -P 8 -r $chr -sites $chr.site -out $chr -GL 1 -doCounts 1 -doSaf 1 -anc $outgroup.fa -doThetas 1 -pest $pop.sfs

#Estimate Tajimas D and other statistics using sliding window
if [ $# -ne 4 ]; then
 echo "error.. need args"
 echo "command:$0 <thetas.idx> <windowsize> <stepsize> <output.prefix>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export thetas=$1
export win=$2
export step=$3
export out=$4
thetaStat do_stat $thetas -win $win -step $step -type 2 -outnames $out
fi

##sliding window XPEHH analysis
#Estimate the raw XP-EHH value for each site 
export chr=$1
vcftools --gzvcf $chr.phase.vcf.gz --keep wild --out wild.$chr --recode
vcftools --gzvcf $chr.phase.vcf.gz --keep domestic --out dom.$chr --recode
perl -e 'while(<>){next if /^#/; chomp; @tmp=split/\t/,$_,4; print "$tmp[0]\t$tmp[0]:$tmp[1]\t$tmp[1]\t$tmp[1]\n"}' wild.$chr.recode.vcf >$chr.map
selscan --xpehh --vcf wild.$chr.recode.vcf --vcf-ref dom.$chr.recode.vcf --map $chr.map --cutoff 0.01 --out $chr --threads 4
#norm
norm --xpehh --files *xpehh.out --qbins 10 --bp-win --winsize 20000
#Estimate XP-EHH using sliding window
perl XPEHH.window.pl <chr fai> <normalized XPEHH> <window size> <step size> <output>
