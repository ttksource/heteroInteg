sample=$1
ref=
input=
outdir=
rm $sample*
# ------ grep ------
# grep snp, indel_sv
gatk SelectVariants \
-R $ref \
-V $input/$sample.deepvariant.vcf.gz \
--select-type-to-include SNP \
-O $outdir/$sample.deepvariant.SNP.vcf.gz
gatk SelectVariants \
-R $ref \
-V $input/$sample.deepvariant.vcf.gz \
--select-type-to-include INDEL \
-O $outdir/$sample.deepvariant.INDEL_SV_temp.vcf.gz
# grep indel
zcat  $outdir/$sample.deepvariant.INDEL_SV_temp.vcf.gz | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {svlen=length($5)-length($4);if(svlen>-50 && svlen<50){print}}' | bgzip > $outdir/$sample.deepvariant.INDEL.vcf.gz
tabix -p vcf $outdir/$sample.deepvariant.INDEL.vcf.gz
rm *_temp*
# ------ filter ------
echo " ------ $sample ------ " >> reads.cnt
# filter SNP
depth=8
support=3
file=$outdir/$sample.deepvariant.SNP.vcf.gz
name=$(basename $file .vcf.gz)
gatk VariantFiltration \
        -V $file \
        -O $outdir/$name.label.vcf.gz \
        -G-filter "AD[1] < $support" -G-filter-name "lowSUP" \
        -G-filter "DP < $depth" -G-filter-name "DP" \
        -G-filter "isHomRef == 1" -G-filter-name "Ref"
zcat $outdir/$name.label.vcf.gz |awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 !~ /NC_012920.1/ && $7 ~ /PASS/ && $9 !~ /FT/ {split($10, a, ":"); if (a[1] !~ /\./) print}' | bgzip >  $outdir/$name.filter.vcf.gz
tabix -p vcf $outdir/$name.filter.vcf.gz
zcat $i |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> reads.cnt
zcat $outdir/$name.filter.vcf.gz |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name".filter",cnt}' >> reads.cnt
# filter INDEL
file=$outdir/$sample.deepvariant.INDEL.vcf.gz
name=$(basename $file .vcf.gz)
gatk VariantFiltration \
        -V $file \
        -O $outdir/$name.label.vcf.gz \
        -G-filter "AD[1] < $support" -G-filter-name "lowSUP" \
        -G-filter "DP < $depth" -G-filter-name "DP" \
        -G-filter "isHomRef == 1" -G-filter-name "Ref"
zcat $outdir/$name.label.vcf.gz |awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 !~ /NC_012920.1/ && $7 ~ /PASS/ && $9 !~ /FT/{split($10, a, ":"); if (a[1] !~ /\./) print}' | bgzip >  $outdir/$name.filter.vcf.gz
tabix -p vcf $outdir/$name.filter.vcf.gz
zcat $i |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> reads.cnt
zcat $outdir/$name.filter.vcf.gz |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name".filter",cnt}' >> reads.cnt
