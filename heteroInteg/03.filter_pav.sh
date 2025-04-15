sample=$1
ref=
input=
outdir=
rm $sample*
# ------ grep ------
# grep snp, indel_sv
gatk SelectVariants \
-R $ref \
-V $input/$sample.pav.vcf.gz \
--select-type-to-include SNP \
-O $outdir/$sample.pav.SNP.vcf.gz
gatk SelectVariants \
-R $ref \
-V $input/$sample.pav.vcf.gz \
--select-type-to-include INDEL \
-O $outdir/$sample.pav.INDEL_SV_temp.vcf.gz
# grep indel
zcat  $outdir/$sample.pav.INDEL_SV_temp.vcf.gz | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {svlen=length($5)-length($4);if(svlen>-50 && svlen<50){print}}' | bgzip > $outdir/$sample.pav.INDEL.vcf.gz
tabix -p vcf $outdir/$sample.pav.INDEL.vcf.gz
# grep sv
zcat $outdir/$sample.pav.INDEL_SV_temp.vcf.gz | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $8 ~ /SVLEN/{print}' | bgzip > $outdir/$sample.pav.withSVLEN_temp.vcf.gz
tabix -p vcf $outdir/$sample.pav.withSVLEN_temp.vcf.gz
gatk VariantFiltration \
-V $outdir/$sample.pav.withSVLEN_temp.vcf.gz \
-O $outdir/$sample.pav.SVlabel_temp.vcf.gz \
-filter "SVLEN < 50 && SVLEN > -50" --filter-name "lowLEN"
zcat $outdir/$sample.pav.SVlabel_temp.vcf.gz |awk 'BEGIN{OFS="\t"} /^#/ {print; next} $7 !~ /lowLEN/{print}' | bgzip > $outdir/$sample.pav.SV.vcf.gz
tabix -p vcf $outdir/$sample.pav.SV.vcf.gz
rm *_temp*
# ------ filter ------
echo " ------ $sample ------ " >> reads.cnt
# filter SNP
file=$outdir/$sample.pav.SNP.vcf.gz
name=$(basename $file .vcf.gz)
zcat $file | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 ~ /chrX/ || $1 ~ /chrY/{print;next} $1 !~ /NC_012920.1/ && $10 !~ /\./ && $10 !~ /0\|0/{print}'| bgzip > $outdir/$name.filter.vcf.gz
tabix -p vcf $outdir/$name.filter.vcf.gz
zcat $file |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> reads.cnt
zcat $outdir/$name.filter.vcf.gz |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name".filter",cnt}' >> reads.cnt
# filter INDEL
file=$outdir/$sample.pav.INDEL.vcf.gz
name=$(basename $file .vcf.gz)
zcat $file | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 ~ /chrX/ || $1 ~ /chrY/{print;next} $1 !~ /NC_012920.1/ && $10 !~ /\./ && $10 !~ /0\|0/{print}'| bgzip > $outdir/$name.filter.vcf.gz
tabix -p vcf $outdir/$name.filter.vcf.gz
zcat $file |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> reads.cnt
zcat $outdir/$name.filter.vcf.gz |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name".filter",cnt}' >> reads.cnt
# filter SV
file=$outdir/$sample.pav.SV.vcf.gz
name=$(basename $file .vcf.gz)
zcat $file | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 ~ /chrX/ || $1 ~ /chrY/{print;next} $1 !~ /NC_012920.1/ && $10 !~ /\./ && $7 ~ /PASS/ && $10 !~ /0\|0/{print}'| bgzip > $outdir/$name.filter.vcf.gz
tabix -p vcf $outdir/$name.filter.vcf.gz
zcat $file |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name,cnt}' >> reads.cnt
zcat $outdir/$name.filter.vcf.gz |awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print name".filter",cnt}' >> reads.cnt
