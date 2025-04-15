sample=$1
n_isec=2
sv_distance=1000
remark="none"
ref=
p=
cen=
telo=

echo "------ $sample, isec=$n_isec, sv_distance=$sv_distance, remark=$remark ------" >> records.cnt
for seqtech in ont hifi;do
	name=$sample.$seqtech
	echo -e "$name"
	f1=$p/$name.sniffle.filter.std.vcf.gz
	f2=$p/$name.svim.filter.std.vcf.gz
	f3=$p/$name.cuteSV.filter.std.vcf.gz
	mkdir $name.bcfIsec
	# ------ norm ------
	bcftools norm -f $ref $f1 --check-ref s -O z > $name.bcfIsec/f1.norm.vcf.gz
	bcftools norm -f $ref $f2 --check-ref s -O z > $name.bcfIsec/f2.norm.vcf.gz
	bcftools norm -f $ref $f3 --check-ref s -O z > $name.bcfIsec/f3.norm.vcf.gz
	bcftools index -t $name.bcfIsec/f1.norm.vcf.gz
	bcftools index -t $name.bcfIsec/f2.norm.vcf.gz
	bcftools index -t $name.bcfIsec/f3.norm.vcf.gz
	# ------ bcftools intersection ------
	bcftools isec -n+$n_isec -c all -p $name.bcfIsec $f1 $f2 $f3
	bcftools merge --force-samples -m all -i 'SVTYPE:join,SVLEN:join' -O z -o $name.bcfIsec/merge.draft.vcf.gz  $name.bcfIsec/f1.norm.vcf.gz $name.bcfIsec/f2.norm.vcf.gz $name.bcfIsec/f3.norm.vcf.gz
	bcftools index -t $name.bcfIsec/merge.draft.vcf.gz
	bcftools view -T $name.bcfIsec/sites.txt $name.bcfIsec/merge.draft.vcf.gz -O v -o $name.bcfIsec_temp.vcf
	echo "bcftools intersection done"
	cnt=$name.bcfIsec_temp.vcf
	awk -v name="$name.bcfIsec" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $cnt >> records.cnt
	rm -r $name.bcfIsec
	# ------ SURVIVOR intersection ------
	echo -e "$p/$name.sniffle.filter.std.vcf\n$p/$name.svim.filter.std.vcf\n$p/$name.cuteSV.filter.std.vcf" > surList_temp
	echo -e "$p/$name.sniffle.filter.std.vcf\n$p/$name.svim.filter.std.vcf\n$p/$name.cuteSV.filter.std.vcf\n"
	SURVIVOR merge surList_temp $sv_distance $n_isec 1 -1 -1 50 $name.surIsec_temp.vcf
	cnt=$name.surIsec_temp.vcf
	awk -v name="$name.surIsec" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $cnt >> records.cnt
	echo "SURVIVOR intersection done"
	# ------ bcftools/SURVIVOR union ------
	echo -e "$name.bcfIsec_temp.vcf\n$name.surIsec_temp.vcf" > surList_temp
	echo -e "$name.bcfIsec_temp.vcf\n$name.surIsec_temp.vcf\n"
	SURVIVOR merge surList_temp $sv_distance 1 1 -1 -1 50 $name.union_temp.vcf
	cnt=$name.union_temp.vcf
	awk -v name="$name.union" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $cnt >> records.cnt
	echo "bcftools/SURVIVOR union done"
done
# ------ hifi/ont union ------
echo -e "$sample.hifi.union_temp.vcf\n$sample.ont.union_temp.vcf" > surList_temp
echo -e "$sample.hifi.union_temp.vcf\n$sample.ont.union_temp.vcf"
SURVIVOR merge surList_temp $sv_distance 1 1 -1 -1 50 $sample.union_temp.vcf
cnt=$sample.union_temp.vcf
awk -v name="$sample.union" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $cnt >> records.cnt
# ------ intersect with pav ------
pavPath=/share/home/zhanglab/user/qiulingxin/projects/HV/02.varcall/pav/02.grep_filter/$sample.pav.SV.filter.vcf.gz
bgzip -dc $pavPath > $sample.pav_temp.vcf
awk -v name="pav" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $sample.pav_temp.vcf >> records.cnt
echo -e "$sample.pav_temp.vcf\n$sample.union_temp.vcf" > surList_temp
echo -e "$sample.pav_temp.vcf\n$sample.union_temp.vcf\n"
SURVIVOR merge surList_temp $sv_distance 2 1 -1 -1 50 $sample.merge_orig_temp.vcf
cnt=$sample.merge_orig_temp.vcf
awk -v name="$sample.merge" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $cnt >> records.cnt
bcftools annotate -x ^INFO/END,^INFO/SVTYPE,^INFO/SVLEN,^FMT/GT -o $sample.merge.rawGT_temp.vcf $sample.merge_orig_temp.vcf
# ------ select same GT ------
file=$sample.merge.rawGT_temp.vcf
awk '$1 ~ /^#/{print "True"}' $file > ifGTsame_temp.bool
awk '$1 !~ /^#/{print $10,$11}' $file | awk '{sub(/\|/, ""); sub(/\//,""); print}' > GT_temp.list
echo "
with open(\"ifGTsame_temp.bool\",\"a\") as wri:
    with open(\"GT_temp.list\" ,\"r\") as read:
        for line in read:
            l = line.strip().split()
            a, b = l
            if set(list(a))==set(list(b)):
                wri.write(\"True\n\")
            else:
                wri.write(\"False\n\")" > ifGTsame_temp.py
python ifGTsame_temp.py
paste ifGTsame_temp.bool $file | awk '$1 ~ /True/{print}'| cut -f2-11 |bcftools sort -O v > $sample.SV_temp.vcf
#tabix -p vcf $sample.SV_temp.vcf
#bgzip -dc $sample.SV.vcf.gz > ${sample}_withXY.SV.vcf
#awk '$1 !~ /^#/ && ($1 ~/chrX/ || $1~/chrY/){print}' $sample.merge.rawGT_temp.vcf |cut -f1-10  >> ${sample}_withXY.SV.vcf
#bcftools sort ${sample}_withXY.SV.vcf -O z > ${sample}_withXY.SV.vcf.gz
#rm ${sample}_withXY.SV.vcf
#tabix -p vcf ${sample}_withXY.SV.vcf.gz
#zcat ${sample}_withXY.SV.vcf.gz| awk -v name="${sample}_withXY.sameGT" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> records.cnt
awk -v name="$sample.sameGT" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' $sample.SV_temp.vcf >> records.cnt

awk '$1 ~ /^#/{print;next}'  $sample.SV_temp.vcf > $sample.cen_temp.vcf
bedtools subtract -a $sample.SV_temp.vcf -b $cen >> $sample.cen_temp.vcf
awk '$1 ~ /^#/{print;next}'  $sample.SV_temp.vcf  > $sample.SV.vcf
bedtools subtract -a $sample.cen_temp.vcf -b $telo >> $sample.SV.vcf
bgzip $sample.SV.vcf
tabix -p vcf $sample.SV.vcf.gz

rm *_temp*
