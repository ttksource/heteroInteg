REF=
REPIDX=
FILELIST=
#FILELIST=
for i in $(cat $FILELIST)
do
	filename=$(basename "$i" .fastq.gz)
	echo "$filename"
	echo "$i"
	echo "#!/bin/bash
#SBATCH --job-name=hifi.map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20g
date
winnowmap -t 8 -H --MD -W $REPIDX -ax map-pb $REF $i | samtools view -bh -@ 8 - > $filename.bam
samtools sort -@ 8 $filename.bam > $filename.sort.bam
date " > job.$filename.sh
	sbatch job.$filename.sh
done
