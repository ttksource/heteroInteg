outdir=
indir=

vlo=0.3
vhi=0.7
support=3
depth=8

for i in $indir/*.vcf;do
    name=$(basename $i .vcf)
	echo "------ $name ------" >> reads.cnt
	# 0: original
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $i >> reads.cnt
    # 1: variant type: BND
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $8 !~ /BND/' $i > $outdir/$name.1_type_temp.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $outdir/$name.1_type_temp.vcf >> reads.cnt
	# 2 : autosome: mitochondria, XY
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 !~ /NC_012920.1/ && $1 !~ /chrX/ && $1 !~ /chrY/{print}' $outdir/$name.1_type_temp.vcf > $outdir/$name.2_autosome_temp.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $outdir/$name.2_autosome_temp.vcf >> reads.cnt
	# 3: length: no SVLEN, SVLEN
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $8 ~ /SVLEN/' $outdir/$name.2_autosome_temp.vcf > $outdir/$name.3_withlen_temp.vcf
	gatk VariantFiltration \
	-V $outdir/$name.3_withlen_temp.vcf \
	-O $outdir/$name.3_withlen_label_temp.vcf \
	-filter "SVLEN < 50 && SVLEN > -50" --filter-name "lowLEN"
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 !~ /lowLEN/{print}' $outdir/$name.3_withlen_label_temp.vcf >  $outdir/$name.3_res_temp.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $outdir/$name.3_res_temp.vcf >> reads.cnt
	# 4: depth: DP, DV
	# cuteSV
    if [[ $name =~ cute  ]]; then
        # flag: SVLEN, DP, DV, HomRef, VAF
        gatk VariantFiltration \
        -V $outdir/$name.3_res_temp.vcf \
        -O $outdir/$name.4_label_temp.vcf \
        -filter "RE < $support" --filter-name "lowSUP" \
        -G-filter "DV + DR < $depth" -G-filter-name "DP"
    fi
        # svim
    if [[ $name =~ svim  ]]; then
        # flag: SVLEN, DP, DV, HomRef, VAF
        gatk VariantFiltration \
        -V $outdir/$name.3_res_temp.vcf \
        -O $outdir/$name.4_label_temp.vcf \
        -filter "SUPPORT < $support" --filter-name "lowSUP" \
        -G-filter "DP < $depth" -G-filter-name "DP"
    fi
        # sniffle
    if [[ $name =~ sniffle  ]]; then
        # flag: SVLEN, DP, DV, HomRef, VAF
        gatk VariantFiltration \
        -V $outdir/$name.3_res_temp.vcf \
        -O $outdir/$name.4_label_temp.vcf \
        -filter "SUPPORT < $support" --filter-name "lowSUP" \
        -G-filter "DV + DR < $depth" -G-filter-name "DP"
    fi
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 !~ /lowSUP/ && $9 !~ /FT/{print}' $outdir/$name.4_label_temp.vcf >  $outdir/$name.4_res_temp.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $outdir/$name.4_res_temp.vcf >> reads.cnt
	# 5: genotype: no GT, HomRef, VAF
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {split($10, a, ":"); if (a[1] !~ /\./) print}' $outdir/$name.4_res_temp.vcf > $outdir/$name.5_gt_temp.vcf
	# cuteSV
    if [[ $name =~ cute  ]]; then
        # flag: SVLEN, DP, DV, HomRef, VAF
        gatk VariantFiltration \
        -V $outdir/$name.5_gt_temp.vcf \
		-O $outdir/$name.5_label_temp.vcf \
        -G-filter "isHomRef == 1" -G-filter-name "Ref" \
        -G-filter "1.0*DV/(DV+DR)>$vlo && 1.0*DV/(DV+DR)<$vhi && isHet == 0" -G-filter-name "Gf" \
        -G-filter "1.0*DV/(DV+DR)>$vhi && isHomVar == 0" -G-filter-name "Gf1"
    fi
    # svim
    if [[ $name =~ svim  ]]; then
		# flag: SVLEN, DP, DV, HomRef, VAF
		gatk VariantFiltration \
        -V $outdir/$name.5_gt_temp.vcf \
        -O $outdir/$name.5_label_temp.vcf \
        -G-filter "isHomRef == 1" -G-filter-name "Ref" \
        -G-filter "1.0*AD[1]/(AD[1]+AD[0])>$vlo && 1.0*AD[1]/(AD[1]+AD[0])<$vhi && isHet == 0" -G-filter-name "Gf" \
		-G-filter "1.0*AD[1]/(AD[1]+AD[0])>$vhi && isHomVar == 0" -G-filter-name "Gf1"
    fi
    # sniffle
    if [[ $name =~ sniffle  ]]; then
        # flag: SVLEN, DP, DV, HomRef, VAF
        gatk VariantFiltration \
        -V $outdir/$name.5_gt_temp.vcf \
        -O $outdir/$name.5_label_temp.vcf \
        -G-filter "isHomRef == 1" -G-filter-name "Ref" \
        -G-filter "1.0*DV/(DV+DR)>$vlo && 1.0*DV/(DV+DR)<$vhi && isHet == 0" -G-filter-name "Gf" \
        -G-filter "1.0*DV/(DV+DR)>$vhi && isHomVar == 0" -G-filter-name "Gf1"
    fi
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $9 !~ /FT/{print}' $outdir/$name.5_label_temp.vcf >  $outdir/$name.5_res_temp.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $outdir/$name.5_res_temp.vcf >> reads.cnt
	# 6: flag
	awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 ~ /PASS/{print}' $outdir/$name.5_res_temp.vcf >  $outdir/$name.filter.vcf
	awk -v name="$name" '$0 !~ /^#/{cnt+=1} END{print cnt}' $outdir/$name.filter.vcf >> reads.cnt
	rm *_temp*
done
