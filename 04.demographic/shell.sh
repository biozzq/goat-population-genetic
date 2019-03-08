#######generater the regions of reference genome that should be considered to be uniquely alignable
#!/bin/sh
export SV_DIR=/svtoolkit
export LD_LIBRARY_PATH=/svtoolkit/bwa:${LD_LIBRARY_PATH}
classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
mx="-Xmx40g"
chroms=`cat chr | cut -f 1`
for chr in ${chroms}
do
java -cp ${classpath} ${mx} org.broadinstitute.sv.apps.ComputeGenomeMask -R ASM.fa  -O `pwd`/${chr}.fasta -readLength 36 -sequence ${chr}
python2 softmasked_to_bed.py ${chr}.fasta >${chr}.alignable.region
done

#########MSMC analysis#############################
#generate mask files containing the sites with depth between half and twice of mean depth  with custom scripts from msmc-tools (https://github.com/stschiff/msmc-tools)
#!/bin/sh
bam=$1
depth=$2
out=$3
mkdir -p $out
cd $out
for i in $bam
do
	for j in {1..29}
	do
	samtools mpileup -q 20 -Q 20 -C 50 -g -r $j -u -f $reference.fa $bam | bcftools call -c -V indels | python3.5 /msmc-tools/bamCaller.py $depth $out\_$j.mask.bed.gz | gzip -c > $out\_$j.vcf.gz
	done
done

#generate input files used in effective population size analysis for each chromosome
if [ $# -ne 3 ]; then
 echo "error.. need args"
 echo "command:$0 <ind1> <ind2> <pop>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export ind1=$1
export ind2=$2
export pop=$3
mkdir -p $pop
for i in {1..29}
do
	python3.5 /msmc-tools/generate_multihetsep.py --chr $i $ind1.$i.vcf.gz $ind2.$i.vcf.gz --mask $i.alignable.region --mask $ind1\_$i.mask.bed.gz --mask $ind2\_$i.mask.bed.gz >$ind1.$ind2.$i.multihetsep
done
fi

#effective population size analysis
#!/bin/sh
if [ $# -ne 2 ]; then
 echo "error.. need args"
 echo "command:$0 <out.prefix> <mul.prex>"
 exit 1
else
  for arg in "$@"
  do
   echo $arg
  done
  echo -n "msmc2 -i 50 -t 10 -p 1*2+25*1+1*2+1*3 -o $1 -I 0,1,2,3" >>$1.sh
  for i in {1..29}
   do
    echo -n " $2.$i.multihetsep" >>$1.sh
   done
fi

#generate input files used in split history analysis for each chromosome
#!/bin/sh
if [ $# -ne 6 ]; then
 echo "error.. need args"
 echo "command:$0 <ind1> <ind2> <ind3> <ind4> <pop1> <pop2>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export ind1=$1
export ind2=$2
export ind3=$3
export ind4=$4
export pop1=$5
export pop2=$6
for i in {1..29}
do
python3.5 /msmc-tools/generate_multihetsep.py --chr $i $ind1.$i.vcf.gz $ind2.$i.vcf.gz $ind3.$i.vcf.gz $ind4.$i.vcf.gz --mask $i.alignable.region --mask $ind1\_$i.mask.bed.gz --mask $ind2\_$i.mask.bed.gz --mask $ind3\_$i.mask.bed.gz --mask $ind4\_$i.mask.bed.gz >$ind1.$ind2.$ind3.$ind4.$i.multihetsep
done
fi
#split history analysis
#!/bin/sh
if [ $# -ne 3 ]; then
 echo "error.. need args"
 echo "command:$0 <pop1> <pop2> <mulsep.prefix>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
echo -n "msmc2 -i 30 -t 10 -p 1*2+20*1+1*2+1*3 -o $1 -I 0,1,2,3" >>msmccommand1.sh
for i in {1..29}
 do
  echo -n " $3.$i.multihetsep" >>msmccommand1.sh
 done
echo -n "msmc2 -i 30 -t 10 -p 1*2+20*1+1*2+1*3 -o $2 -I 4,5,6,7" >>msmccommand2.sh
for i in {1..29}
 do
  echo -n " $3.$i.multihetsep" >>msmccommand2.sh
 done
echo -n "msmc2 -i 30 -t 10 -P 0,0,0,0,1,1,1,1 -s -p 1*2+20*1+1*2+1*3 -o $1_$2" >>msmccommand3.sh
for i in {1..29}
 do
  echo -n " $3.$i.multihetsep" >>msmccommand3.sh
 done
fi
#generate a combined msmc output file using the combineCrossCoal.py script from the msmc-tools repository
python3 /msmc-tools/combineCrossCoal.py $pop1\_$pop2.final.txt $pop1.final.txt $pop2.final.txt >$results


#########SMC++ analysis############################
#converts VCF data to the format used by SMC++
#!/bin/sh
if [ $# -ne 3 ]; then
 echo "error.. need args"
 echo "command:$0 <poplabel> <distinguished> <popname|sep by comma>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done

export poplabel=$1
export distinguished=$2
export popname=$3
for chr in {1..29}
do
	smc++ vcf2smc -d $distinguished $distinguished --mask $chr.mask.bed.gz $chr.vcf.gz $distinguished.$chr.smc.gz $chr $poplabel:$popname
done
fi
#bootstrapping by breaking up the genome into 5-Mb segments and then randomly sampling with replacement
#!/bin/sh
export prefix=$1
bootstrap_smcpp.py --nr_chromosomes 29 $prefix $prefix.1.smc.gz $prefix.2.smc.gz $prefix.3.smc.gz $prefix.4.smc.gz $prefix.5.smc.gz $prefix.6.smc.gz $prefix.7.smc.gz $prefix.8.smc.gz $prefix.9.smc.gz $prefix.10.smc.gz $prefix.11.smc.gz $prefix.12.smc.gz $prefix.13.smc.gz $prefix.14.smc.gz $prefix.15.smc.gz $prefix.16.smc.gz $prefix.17.smc.gz $prefix.18.smc.gz $prefix.19.smc.gz $prefix.20.smc.gz $prefix.21.smc.gz $prefix.22.smc.gz $prefix.23.smc.gz $prefix.24.smc.gz $prefix.25.smc.gz $prefix.26.smc.gz $prefix.27.smc.gz $prefix.28.smc.gz $prefix.29.smc.gz
#estimate population size history
#!/bin/sh
source ~/.bashrc
if [ $# -ne 1 ]; then
 echo "error.. need args"
 echo "command:$0 <output>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export out=$1
mkdir -p $out
/stor9000/apps/users/NWSUAF/2012010954/Software/smcpp_v1.11.1dev/bin/smc++ estimate --cores 4 -o $out --knots 32 --ftol 1e-4 --regularization-penalty 8 --em-iterations 20 4.32e-9 *.gz
fi

#generate the input format used in inferring divergence times jointly with population size histories. The current version of our implementation assumes a 'clean split' model, in which no gene flow occurs after the populations split
#!/bin/sh
source ~/.bashrc
if [ $# -ne 6 ]; then
 echo "error.. need args"
 echo "command:$0 <pop1label> <pop1name|sep by comma> <distinguished_pop1> <pop2label> <pop2name|sep by comma> <distinguished_pop2>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export pop1label=$1
export pop1name=$2
export distinguished_pop1=$3
export pop2label=$4
export pop2name=$5
export distinguished_pop2=$6
for chr in {1..29}
do
smc++ vcf2smc --cores 4 -d $distinguished_pop1 $distinguished_pop1 --mask $chr.mask.bed.gz $chr.clean.vcf.gz $pop1label.$pop2label.$distinguished_pop1.$chr.smc.gz $chr $pop1label:$pop1name $pop2label:$pop2name
smc++ vcf2smc --cores 4 -d $distinguished_pop2 $distinguished_pop2 --mask $chr.mask.bed.gz $chr.clean.vcf.gz $pop2label.$pop1label.$distinguished_pop2.$chr.smc.gz $chr $pop2label:$pop2name $pop1label:$pop1name
done
fi

#bootstrapping by breaking up the genome into 5-Mb segments and then randomly sampling with replacement
#!/bin/sh
export prefix=$1
bootstrap_smcpp.py --nr_chromosomes 29 $prefix $prefix.1.smc.gz $prefix.2.smc.gz $prefix.3.smc.gz $prefix.4.smc.gz $prefix.5.smc.gz $prefix.6.smc.gz $prefix.7.smc.gz $prefix.8.smc.gz $prefix.9.smc.gz $prefix.10.smc.gz $prefix.11.smc.gz $prefix.12.smc.gz $prefix.13.smc.gz $prefix.14.smc.gz $prefix.15.smc.gz $prefix.16.smc.gz $prefix.17.smc.gz $prefix.18.smc.gz $prefix.19.smc.gz $prefix.20.smc.gz $prefix.21.smc.gz $prefix.22.smc.gz $prefix.23.smc.gz $prefix.24.smc.gz $prefix.25.smc.gz $prefix.26.smc.gz $prefix.27.smc.gz $prefix.28.smc.gz $prefix.29.smc.gz

#estimate divergence times
#!/bin/sh
if [ $# -ne 3 ]; then
 echo "error.. need args"
 echo "command:$0 <pop1model> <pop2model> <outdir>"
 exit 1
else
for arg in "$@"
do
 echo $arg
done
export pop1model=$1
export pop2model=$2
export out=$3
mkdir -p $out
smc++ split --cores 8 -o $out --em-iterations 20 $pop1model $pop2model *.smc.gz
fi
