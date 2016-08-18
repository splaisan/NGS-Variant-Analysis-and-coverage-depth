#! /usr/bin/env bash
 
# script: vcf-compare-2.sh
 
# required:
# vcf-compare (vcftools)
# two vcf files from the same reference genome
#
# SP@NC, 2016-03-30

outfolder=${3:-"vcf_compare"}
mkdir -p ${outfolder}

ts=$(date +%s)

header=$(cat <<-END_HEREDOC
## venn comparison results for:
# 1:	${1}
# 2:	${2}
END_HEREDOC
)

echo "${header}" > ${outfolder}/vcf-compare-${ts}.txt

vcf-compare ${1} ${2} | \
	grep ^VN | cut -f 2- | \
	sed -e "s:${1}:1:g" |
	sed -e "s:${2}:2:g" |
	column -t \
	>> ${outfolder}/vcf-compare-${ts}.txt

cat ${outfolder}/vcf-compare-${ts}.txt
