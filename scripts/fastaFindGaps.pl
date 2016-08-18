#!/usr/bin/perl -w
use strict;
use File::Basename;
use Bio::SeqIO;
use Getopt::Std;
use File::Tee qw(tee);

# title: fastaFindGaps.pl
# Search N-regions in multifasta (genomes)
# and produce a BED file with found GAP locations
#
# adapted from http://stackoverflow.com/questions/10319696
#
# Stephane Plaisance (VIB-NC+BITS) 2015/04/02; v1.01
# small edits 2015/04/21; v1.02
# support for archives 2015/06/30 v1.03
#
# handle complex fasta headers including description
# visit our Git: https://github.com/BITS-VIB

# disable buffering to get output during long process (loop)
$|=1;

getopts('i:o:l:h');
our($opt_i, $opt_o, $opt_l, $opt_h);

my $usage="## Usage: fastaFindGaps.pl <-i fasta-file> 
# Additional optional parameters are:
# <-o BED output (optional, deduced from input file)>
# <-l minsize in bps (default to 100bps)>
# <-h to display this help>";

####################
# declare variables
####################

my $fastain = $opt_i || die $usage."\n";
my $minlen = $opt_l || 100;
defined($opt_h) && die $usage."\n";

# counters
our $total = 0;
our $present = 0;
our $absent = 0;
our $totallen = 0;
our $presentlen = 0;
our $absentlen = 0;
our $nlength = 0;

# open stream from BED file
my $outpath = dirname($fastain);
my $basename = basename($fastain);
(my $outbase = $basename) =~ s/\.[^.]+$//;

my $outfile;
if ( defined($opt_o) ) {
	$outfile = $outpath."/".$opt_o;
} else {
	$outfile = $outpath."/".$outbase."-".$minlen."_Gaps.bed";
}

# include size limit and max intensity in file names
open OUT, "> $outfile" || die $!;

# keep log copy of STDOUT (comment out if you do not have 'File::Tee' installed
tee STDOUT, '>', $outfile."_log.txt" or die $!;

# create parser for multiple fasta files
my $parser = OpenArchiveFile($fastain);

# look for $minlen N's in a row
my $motif="[N]{".$minlen.",}";
my $totcnt = 0;

############################################
# loop over records and return hits to BED #
############################################

while( my $seq_obj = $parser->next_seq() ) {
	my $counter=0;

	# load id, and description into strings and merge into header
	my $seqid = $seq_obj->id;
	my $seqdesc = defined $seq_obj->desc ? $seq_obj->desc : "";
	my $seqheader = join(" ", $seqid, $seqdesc);
	$seqheader =~ s/\s+$//;
	# get sequence length
	my $len = $seq_obj->length;

	# count
	$total ++;
	$totallen += $len;
	$present += 1;
	$presentlen += $len;
	print STDOUT "## Searching sequence $seqid for $motif\n";
	my $sequence = $seq_obj->seq();

	# scan for motif and report hits
	while ( $sequence =~ m/$motif/gi ) {
		$counter++;
		my $match_start = $-[0]+1; # BED is zero-based !
		my $match_end = $+[0];
		my $match_seq = $&;
		$nlength += length($&);
		# print in BED5 format when present in cmap
		print OUT join("\t", $seqid, $match_start,
			$match_end, "N-region", length($&), "+")."\n";
		}

	# end for this sequence
	print STDOUT "# found $counter matches for $seqid\n";
	$totcnt += $counter;
}

# close filehandle
close OUT;

# report absent maps and absent length
# reformat lengths with thousand separator
my $percentpresent = sprintf '%.1f%%', 100*$presentlen/$totallen;
my $percentn = sprintf '%.1f%%', 100*$nlength/$presentlen;
$totallen =~ s/\d{1,3}(?=(\d{3})+(?!\d))/$&,/g;
$presentlen =~ s/\d{1,3}(?=(\d{3})+(?!\d))/$&,/g;
$absentlen =~ s/\d{1,3}(?=(\d{3})+(?!\d))/$&,/g;
$nlength =~ s/\d{1,3}(?=(\d{3})+(?!\d))/$&,/g;

print STDOUT "\n############################ summary #####################################\n";
print STDOUT "# $present fasta entries ($presentlen bps)\n";
print STDOUT "# reported a total of $totcnt N-regions of $minlen bps or more\n";
print STDOUT "# for a total N-length of $nlength bps ($percentn of represented sequences)\n";

exit 0;

#### Subs ####
sub OpenArchiveFile {
    my $infile = shift;
    my $FH;
    if ($infile =~ /.fa$|.fasta$/) {
    $FH = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');
    }
    elsif ($infile =~ /.bz2$/) {
    $FH = Bio::SeqIO -> new(-file => "bgzip -c $infile |", -format => 'Fasta');
    }
    elsif ($infile =~ /.gz$/) {
    $FH = Bio::SeqIO -> new(-file => "gzip -cd $infile |", -format => 'Fasta');
    }
    elsif ($infile =~ /.zip$/) {
    $FH = Bio::SeqIO -> new(-file => "unzip -p $infile |", -format => 'Fasta');
    } else {
	die ("$!: do not recognise file type $infile");
	# if this happens add, the file type with correct opening proc
    }
    return $FH;
}

