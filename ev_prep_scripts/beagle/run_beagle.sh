#!/bin/bash

_JAVA_OPTIONS="-Xmx8g"

vcf_file="$1"

while read chrom
do
java -Xmx8g -jar beagle.r1399.jar gl=$vcf_file chrom="$chrom" nthreads=4 burnin-its=10 phase-its=10 out=beagle/"$chrom"_beagle
vcftools --gzvcf beagle/"$chrom"_beagle.vcf.gz --hap-r2 --out hap_ld_stats --hap-r2-positions hap_ld_positions
done < beagle_chrom.txt