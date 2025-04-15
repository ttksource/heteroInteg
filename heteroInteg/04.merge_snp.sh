sample=$1
ref=
p=
cen=
telo=

echo "------ $sample ------" >> records.cnt
for svtype in SNP INDEL;do
	name=$sample.$svtype
	f2=$p/$sample.deepvariant.$svtype.filter.std.vcf.gz
	f1=$p/$sample.pav.$svtype.filter.std.vcf.gz
	zcat $f1 | awk -v name="pav" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> records.cnt
	zcat $f2 | awk -v name="deepvariant" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> records.cnt
	mkdir $name.bcfIsec_temp
	# ------ norm ------
	bcftools norm -f $ref $f1 --check-ref s -O z > $name.bcfIsec_temp/f1.norm.vcf.gz
	bcftools norm -f $ref $f2 --check-ref s -O z > $name.bcfIsec_temp/f2.norm.vcf.gz
	tabix -p vcf $name.bcfIsec_temp/f1.norm.vcf.gz
	tabix -p vcf $name.bcfIsec_temp/f2.norm.vcf.gz
	# ------ bcftools intersection ------
	bcftools isec -n=2 -c all -p $name.bcfIsec_temp $f1 $f2
	bcftools merge --force-samples -m all -i SVTYPE:join,SVLEN:join -O z -o $name.bcfIsec_temp/merge.draft.vcf.gz  $name.bcfIsec_temp/f1.norm.vcf.gz $name.bcfIsec_temp/f2.norm.vcf.gz
	tabix -p vcf $name.bcfIsec_temp/merge.draft.vcf.gz
	bcftools view -T $name.bcfIsec_temp/sites.txt $name.bcfIsec_temp/merge.draft.vcf.gz -O v -o $name.bcfIsec_temp.vcf
	echo "bcftools intersection done"
	cnt=$name.bcfIsec_temp.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $cnt >> records.cnt
	# ------ select same GT ------
	file=$name.bcfIsec_temp.vcf
	awk '$1 ~ /^#/{print "True"}' $file > ifGTsame_temp.bool
	awk '$1 !~ /^#/{print $10,$11}' $file | awk '{sub(/\|/, ""); sub(/\//,""); print}' > GT_temp.list
	echo "with open(\"ifGTsame_temp.bool\",\"a\") as wri:
    with open(\"GT_temp.list\" ,\"r\") as read:
        for line in read:
            l = line.strip().split()
            a, b = l
            if set(list(a))==set(list(b)):
                wri.write(\"True\n\")
            else:
                wri.write(\"False\n\")" > ifGTsame_temp.py
	python ifGTsame_temp.py
	paste ifGTsame_temp.bool $file|awk '$1 ~ /True/{print}'| cut -f2-11 |bcftools sort -O v > $name.f_temp.vcf
	#tabix -p vcf $name.vcf.gz
	#zcat $name.vcf.gz | awk -v name="$name.sameGT" '$0 !~ /^#/{cnt+=1} END{print name,cnt}'>> records.cnt
	
	awk '$1 ~ /^#/{print;next}' $name.f_temp.vcf > $name.cen_temp.vcf
	bedtools subtract -a $name.f_temp.vcf -b $cen >> $name.cen_temp.vcf
	awk '$1 ~ /^#/{print;next}' $name.f_temp.vcf > $name.t_temp.vcf
	bedtools subtract -a $name.cen_temp.vcf -b $telo >> $name.t_temp.vcf
	awk '$1 ~ /^#/{print;next} $10 ~ /1\|0/ || $10 ~ /0\|1/ {print}' $name.t_temp.vcf > $name.vcf
	bgzip $name.vcf
	tabix -p vcf $name.vcf.gz

	rm -r *bcfIsec_temp/
	rm *_temp*
done

