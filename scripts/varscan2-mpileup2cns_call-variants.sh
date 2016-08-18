#! /usr/bin/env bash
# script varscan2-mpileup2cns_call-variants.sh
 
# required
# samtools version: 1.x (htslib)
# varscan2 version 2.4.1
# vcftools 0.1.14
# vcfutils.pl (samtools)
# vcf2index custom function
 
infolder=BDGP5.78_bwa-mapping
infile=SRR833244_dm5-pe.bam
outfolder=BDGP5.78_bwa_varscan2_variants
mkdir -p ${outfolder}
outfile=BDGP5.78_mpileup2cns.vcf

ref=$REFERENCE/dm5/Drosophila_melanogaster.BDGP5.dna.toplevel.fa.gz

# call SNV and InDels variants from samtools mpileup
maxmem="24G"

samtools mpileup -f ${ref} ${infolder}/${infile} | \
	java -Xmx${maxmem} -jar $VARSCAN/VarScan.v2.4.1.jar mpileup2cns \
	--variants --output-vcf 1 --p-value 0.01 \
	>${outfolder}/${outfile} \
	2>${outfolder}/varscan2_mpileup2cns.err

## compress and index result
bgzip -c ${outfolder}/${outfile} \
	> ${outfolder}/${outfile}.gz &&
	tabix -p vcf ${outfolder}/${outfile}.gz

echo
