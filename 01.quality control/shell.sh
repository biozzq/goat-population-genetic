#######trimmomatic-0.38.jar
java -jar trimmomatic-0.38.jar PE -threads 2 -trimlog $leftname.log $directory/$leftname$ARGV[1] $directory/$leftname$ARGV[2] $leftname\_1.fq.gz $leftname\_1.single.fastq.gz $leftname\_2.fq.gz $leftname\_2.single.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33

#######fastp
#!/bin/sh
export fastq1=$1
export fastq2=$2
export prefix=$3
~/bin/fastp -i $fastq1 -I $fastq2 -o $prefix\_1.fq.gz -O $prefix\_2.fq.gz --cut_front --cut_front_window_size=1 --cut_front_mean_quality=3 --cut_tail --cut_tail_window_size=1 --cut_tail_mean_quality=3 --cut_right --cut_right_window_size=4 --cut_right_mean_quality=15 --length_required 40 --json $prefix.json --html $prefix.html --dont_overwrite
#~/bin/fastp -i $fastq -o $prefix.fq.gz --cut_front --cut_front_window_size=1 --cut_front_mean_quality=3 --cut_tail --cut_tail_window_size=1 --cut_tail_mean_quality=3 --cut_right --cut_right_window_size=4 --cut_right_mean_quality=15 --dont_overwrite --length_required 40
