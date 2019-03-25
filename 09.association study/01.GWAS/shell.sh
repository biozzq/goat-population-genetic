#### the first three PCA values (eigenvectors) derived from whole-genome SNPs represented the fixed effect in the mixed model to correct for stratification
#/bin/bash
if [ $# -ne 1 ]; then
 echo "error.. need args"
 echo "command:$0 <PCA file[eigenvec]>"
 exit 1
fi
IN=$1
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' OFS="\t" ${IN} > covariate


#### The kinship derived from whole-genome SNPs of 16 animals was set as a random effect to control for family effects.
#/bin/bash
if [ $# -ne 3 ]; then
 echo "error.. need args"
 echo "command:$0 <VCF file> <Output name> <autosome number>"
 exit 1
fi
VCF=$1
OUT=$2
CHR=$3
plink --vcf ${VCF} --output-missing-genotype 0 --geno 0.1 --maf 0.01 --recode 12 transpose --out ${OUT} --allow-extra-chr --chr-set ${CHR}
######### Add SNP IDs
cut -d " " -f5- ${OUT}.tped > ${OUT}.tped.tmp
cut -d " " -f1-4 ${OUT}.tped > ${OUT}.tped.headcol
awk '{print $1," ",$1,":",$4," ",$3," ",$4}' OFS="" ${OUT}.tped.headcol > ${OUT}.tped.headcol.tmp
paste -d " " ${OUT}.tped.headcol.tmp ${OUT}.tped.tmp > ${OUT}.tped
mv ${OUT}.tped.headcol.tmp ${OUT}.tped.headcol
rm ${OUT}.tped.tmp
#############
emmax-kin -v -h -s -d 10 ${OUT}


#### A mixed linear model program, EMMAX, was used for the association analysis.
#!/bin/sh
if [ $# -ne 4 ]; then
  echo "error.. need args"
  echo "command:$0 <Input> <Trait file> <covariate file> <output>"
  exit 1
fi
IN=$1
Trait=$2
covariate=$3
OUT=$4
#emmax -v -d 10 -t ${IN} -p ${Trait} -k ${IN}.hIBS.kinf -c ${covariate} -o ${OUT}
emmax -v -d 10 -t ${IN} -p ${Trait} -k ${IN}.hIBS.kinf -o ${OUT}
paste ${IN}.tped.headcol ${OUT}.ps > ${OUT}.result.txt
awk '{$3=null;$5=null;$6=null;print}' ${OUT}.result.txt > ${OUT}.assoc.txt
sed -i '1i\CHR SNPID BP P' ${OUT}.assoc.txt
python3 ~/script/GWAS/Emmax/EmmaxFilter.py -i ${OUT}.assoc.txt -o ${OUT}.assoc.plot.txt
Rscript ~/script/GWAS/manqqplot.r ${OUT}.assoc.plot.txt ${OUT}
