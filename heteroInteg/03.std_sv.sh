rm *.tbi
for i in *.vcf.gz;do
	name=$(basename $i .vcf.gz)
	echo $i
	if [[ $i =~ cute ]];then
		# bcftools annotate -c INFO/SUPPORT:=INFO/RE $i| bcftools annotate -x ID,QUAL,FILTER,^INFO/END,^INFO/SVTYPE,^INFO/SVLEN,^INFO/SUPPORT,^FMT/GT |bcftools reheader -h header.txt |bgzip > $name.std.vcf.gz
		bcftools annotate -c INFO/SUPPORT:=INFO/RE $i| bcftools annotate -x ID,QUAL,FILTER,^INFO/END,^INFO/SVTYPE,^INFO/SVLEN,^FMT/GT |bcftools reheader -h header.2.txt |bgzip > $name.std.vcf.gz
	else
		# bcftools annotate -x ID,QUAL,FILTER,^INFO/END,^INFO/SVTYPE,^INFO/SVLEN,^INFO/SUPPORT,^FMT/GT $i |bcftools reheader -h header.txt |bgzip > $name.std.vcf.gz
		bcftools annotate -x ID,QUAL,FILTER,^INFO/END,^INFO/SVTYPE,^INFO/SVLEN,^FMT/GT $i |bcftools reheader -h header.2.txt |bgzip > $name.std.vcf.gz
	fi
done
for f in ./*.vcf.gz;do
	name=$(basename $f .vcf.gz)
	bgzip -dc $f > $name.vcf
	bcftools index -t $f
done
rm *pav*
