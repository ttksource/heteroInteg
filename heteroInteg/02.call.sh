p=
mkdir -p $p
# deepvariant
mkdir -p 
echo '#!/bin/bash
#SBATCH --job-name=dp$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=35g
INPUT=
OUTPUT=/dpvar/$b
REF=/Reference/
singularity run -B '${INPUT}':'/input','${OUTPUT}':'/output','${REF}':'/ref' /deepvariant_1.8.0.sif /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/ref/CHM13v2m.fasta --reads=/input/$b.hifi.sort.bam --output_vcf=/output/$b.deepvariant.vcf.gz --sample_name $b --intermediate_results_dir /output/intermediate_results_dir --vcf_stats_report --num_shards=24' > $p/run.$b_deepvariant.sh
sbatch $p/run.$b_deepvariant.sh
# cuteSV hifi
mkdir -p /$b/hifi/cuteSV
echo '#!/bin/bash
#SBATCH --job-name=cuH$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=3g
REF=/Reference/CHM13v2m.fasta
BAM=$b.hifi.sort.bam
OUTPUT=/$b/hifi/cuteSV
cuteSV --genotype -l 50 -t 16 -s 5 -S $b.hifi --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 $BAM $REF $b.hifi.cuteSV.vcf ${OUTPUT}' > $p/run.$b.hifi_cuteSV.sh
sbatch $p/run.$b.hifi_cuteSV.sh
# cuteSV ont
mkdir -p /$b/ont/cuteSV
echo '#!/bin/bash
#SBATCH --job-name=cuO$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=15g
REF=/Reference/CHM13v2m.fasta
BAM=$b.ont.sort.bam
OUTPUT=/$b/ont/cuteSV
cuteSV --genotype -l 50 -t 16 -s 5 -S $b.ont --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 $BAM $REF $b.ont.cuteSV.vcf ${OUTPUT}' > $p/run.$b.ont_cuteSV.sh
sbatch $p/run.$b.ont_cuteSV.sh
# svim hifi
mkdir -p /$b/hifi/svim
echo '#!/bin/bash
#SBATCH --job-name=svmH$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30g
REF=/Reference/CHM13v2m.fasta
BAM=$b.hifi.sort.bam
OUTPUT=/$b/hifi/svim
conda init
conda activate svim_env
svim alignment ${OUTPUT} ${BAM} ${REF}' > $p/run.$b.hifi_svim.sh
sbatch $p/run.$b.hifi_svim.sh
# svim ont
mkdir -p /$b/ont/svim
echo '#!/bin/bash
#SBATCH --job-name=svmO$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30g
REF=/Reference/CHM13v2m.fasta
BAM=$b.ont.sort.bam
OUTPUT=/$b/ont/svim
conda init
conda activate svim_env
svim alignment ${OUTPUT} ${BAM} ${REF}' > $p/run.$b.ont_svim.sh
sbatch $p/run.$b.ont_svim.sh
# sniffle hifi
mkdir -p /$b/hifi/sniffle
echo '#!/bin/bash
#SBATCH --job-name=snfH$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30g
REF=/Reference/CHM13v2m.fasta
BAM=$b.hifi.sort.bam
OUTPUT=/$b/hifi/sniffle
sniffles -i ${BAM} -v ${OUTPUT}/$b.hifi.sniffle.vcf.gz --reference ${REF} -t 24' > $p/run.$b.hifi_sniffle.sh
sbatch $p/run.$b.hifi_sniffle.sh
# sniffle ont
mkdir -p /$b/ont/sniffle
echo '#!/bin/bash
#SBATCH --job-name=snfO$b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30g
REF=/Reference/CHM13v2m.fasta
BAM=$b.ont.sort.bam
OUTPUT=/$b/ont/sniffle
sniffles -i ${BAM} -v ${OUTPUT}/$b.ont.sniffle.vcf.gz --reference ${REF} -t 24' > $p/run.$b.ont_sniffle.sh
sbatch $p/run.$b.ont_sniffle.sh
