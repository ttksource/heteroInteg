sample=$1
for i in 
$sample.deepvariant.*.filter.vcf.gz;do
	name=$(basename $i .vcf.gz)
	echo $i
	bcftools annotate -x ID,QUAL,FILTER,^INFO/SVTYPE,^INFO/SVLEN,^FMT/GT $i |bcftools reheader -h header.2.txt |bgzip > $name.std.vcf.gz
	bgzip -dc $name.std.vcf.gz > $name.std.vcf
	rm $name.std.vcf.gz.tbi
	tabix -p vcf $name.std.vcf.gz
done
