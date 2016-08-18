#! /usr/bin/env bash
 
# script: bedtools_genomecoverage.sh
 
# required:
# bedtools Version: v2.25.0
# IGVtools v2.3.4
# a genome text table with chromosome names and their lengths
#
# Summary: Compute the coverage of a feature file across a genome.
#
# SP@NC, 2016-03-29

# we use here the genome size file generated from the Picard .dict file
genome=reference/Drosophila_melanogaster.BDGP5.dna.toplevel.chrom.sizes

# sample or full data
infolder=BDGP5.78_bwa-mapping

outfolder="BDGP5.78_bedtools_genomeCvg"
mkdir -p ${outfolder}

# input and output
inbam=${infolder}/*.bam
name=${infolder%%.bam}

echo
echo "# processing ${inbam}"
 
# compute a coverage summary used to plot the coverage distribution with R
bedtools genomecov -ibam ${inbam} \
	-g ${genome} > ${outfolder}/${name}-histo.txt

# compute full genome coverage for IGV track
# -bga: Report depth in BedGraph format like with -bg. 
# However with this option, regions with zero coverage are also reported. 
# This allows one to quickly extract all regions of a genome with 0 coverage 
# by applying: “grep -w 0$” to the output (see below).
# the output is fairly large and can be reduced with '''igvtools toTDF'''
bedtools genomecov -bga \
	-ibam ${inbam} \
	-g ${genome} > ${outfolder}/${name}.bg
 
# convert .bg to .igv and compress with igvtools on the fly
(gawk 'BEGIN{FS="\t"; OFS="\t"}{print $1, $2, $3, "genomecov", $4}' \
	${outfolder}/${name}.bg > ${outfolder}/${name}-cvg.igv) && \
	(igvtools toTDF -z 4 -f min,max,mean ${outfolder}/${name}-cvg.igv \
	${outfolder}/${name}-cvg.tdf ${genome})
 
# extract 0-coverage regions from the bg file to a second file with '1' as score
grep -w 0$ ${outfolder}/${name}.bg | \
	( gawk 'BEGIN{FS="\t"; OFS="\t"}{print $1, $2, $3, "no-coverage", 1}' \
	> ${outfolder}/nocvg-${name}-cvg.igv ) && \
	( igvtools toTDF -z 4 -f min,max,mean ${outfolder}/nocvg-${name}-cvg.igv \
	${outfolder}/nocvg-${name}-cvg.tdf ${genome})
	
