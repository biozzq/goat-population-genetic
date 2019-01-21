##########################NJ tree analysis###################################
# generate the genomic distance proportion matrix using all SNPs
plink --bfile goat --chr-set 29 --distance-matrix --out distance
# format the .mdist file for MEGA to read it
perl plink.distance.matrix.to.mega.pl <mdist.id> <mdist> <individual number> <out prefix>

##########################admixture analysis#################################
#!/bin/sh
# generate random number: shuf -i 1-10000 -n 1
export K=$1
export seed=$2
admixture -s $seed --cv $prune.bed $K -j4 | tee log${K}.out

##########################PCA analysis#######################################

#Lines(part of lines) starting with symbol "#" are comments. smartpca accepts these lines as comments and ignores when reading the script.
#Command to run PCA in linux command line: smartpca -p Goat_184.par
genotypename:	 $prune.bed
snpname:		 $prune.bim
indivname:		 BEZ_EUR_AFR_SWA_SA_EA.ind
evecoutname:     BEZ_EUR_AFR_SWA_SA_EA.evec  #output file of eigenvectors.
evaloutname:     BEZ_EUR_AFR_SWA_SA_EA.eval  #output file of all eigenvalues
#optional commands:
numoutevec:		 10  #number of eigenvectors to output.  Default is 10.
numoutlieriter:  0 #maximum number of outlier removal iterations.  Default is 5.  To turn off outlier removal, set this parameter to 0.
numoutlierevec:  10 #number of principal components along which to remove outliers during each outlier removal iteration.  Default is 10.
nsnpldregress:   0 #If set to a positive integer, then LD correction is turned on,  and input to PCA will be the residual of a regression involving that many  previous SNPs, according to physical location.  See Patterson et al. 2006.   Default is 0 (no LD correction).  If desiring LD correction,  recommend value is: 2.
outliersigmathresh: 6.0 #number of standard deviations which an individual must exceed, along one of the top (numoutlierevec) principal components, in order for that individual to be removed as an outlier.  Default is 6.0.
outlieroutname: Goat_184_PCA.log #output logfile of outlier individuals removed. If not specified,  smartpca will print this information to stdout, which is the default.
usenorm: YES #Whether to normalize each SNP by a quantity related to allele freq.  Default is YES.  (When analyzing microsatellite data, should be set to NO.
phylipoutname: BEZ_EUR_AFR_SWA_SA_EA.fst.matrix  #output file containing an fst matrix which can be used as input to programs in the PHYLIP package, such as the "fitch" program for constructing phylogenetic trees.

#########################treemix analysis#####################################
#!/bin/sh
if [ $# -ne 3 ]; then
 echo "error.. need args"
 echo "command:$0 <input> <migration events> <LD>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export input=$1
export migration=$2
export LD=$3 # recommend using a value of n that far exceeds the known extent of LD in the organism.
treemix -i $input -m $migration -k $LD -root Outgroup -o migration_$migration
fi

#########################LD analysis###########################################
#!/bin/sh
if [ $# -ne 2 ]; then
 echo "error.. need args"
 echo "command:$0 <popname> <chr>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export pop=$1
export chr=$2
java -jar Haploview.jar -n -minMAF 0.05 -hwcutoff 0.001 -pedfile $pop.chr-$chr.ped -info $pop.chr-$chr.info -dprime -out $pop.$chr -memory 4000
gzip $pop.$chr.LD
fi
