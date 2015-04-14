#!/bin/bash

_JAVA_OPTIONS="-Xmx8g"

while read chrom
do
java -Xmx8g -jar beagle.r1399.jar gt=marines_pac.vcf.recode.vcf chrom="$chrom" nthreads=4 out=beagle/"$chrom"_beagle
done < beagle_chrom.txt