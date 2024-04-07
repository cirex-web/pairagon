#!/usr/bin/perl

use strict;
use Getopt::Long;
use Alignment;

my ($infile, $intype, $outtype);
my ($chr, $program, $species, $gid, $tid);
my ($targetfn, $queryfn);

GetOptions("i:s"=>\$infile,
	   "x:s"=>\$intype,
	   "o:s"=>\$outtype,
	   "c:s"=>\$chr,
	   "p:s"=>\$program,
	   "s:s"=>\$species,
	   "gid:s"=>\$gid,
	   "tid:s"=>\$tid,
	   "t:s"=>\$targetfn,
	   "q:s"=>\$queryfn);

if(!$outtype || (!$infile && !$intype) || ($infile && $intype)){
    print stderr
	"Usage: $0 {-i input_file | -x input_type} \n".
	"          {-o output_type} [other args] \n\n".
	"   -i     The input file (type inferred from filename)\n".
	"   -x     The input type (contents from stdin)\n".
	"   -o     The output type can be vulgar, gtf, gff, psl, seed, or alignment\n\n".
	"Certain output types can use special arguments:\n".
	" -o gtf  [-c chromosome] [-p program]\n".
	"         [-gid gene_id] [-tid transcript_id]\n".
	" -o gff  {-t target_file} {-q query_file}\n".
	"         [-c chromosome] [-s species]\n".
	" -o psl  {-t target_file} {-q query_file}\n".
	" -o seed {-t target_file} {-q query_file}\n".
	" -o alignment {-t target_file} {-q query_file}\n\n".
	"   -c     Chromosome identifier (default: target id)\n".
	"   -p     Program identifier (default: \"Alignment\")\n".
	"   -g     Gene id (default: query id)\n".
	"   -t     Transcript id (default: gene id)\n".
	"   -s     Species id (default: \"??\")\n";
    exit;
}

my $param = "$infile$intype";
my $aln = Alignment::read($param);
if(!$chr){ $chr=$aln->{tid}; }

if(scalar @{$aln->{sectionsref}} != 0){
    $_ = $outtype;
    if(/vulgar/i){
	print "vulgar: " . $aln->vulgar();
    }elsif(/gtf/i){
	if(!$program) { $program="Alignment"; }
	if(!$gid){ $gid = $aln->{qid}; }
	if(!$tid){ $tid = $gid; }

	print $aln->gtf({chr=>$chr, program=>$program,
			 gid=>$gid, tid=>$tid}); 
    }elsif(/psl/i || /gff/i || /alignment/i || /seed/i){
	if(!$queryfn || !$targetfn){
	    die "Need to specify -q and -t options";
	}
	my ($qhead, $query) = getContent("$queryfn");
	my ($thead,$target)=getContent("$targetfn");

	if(/psl/i){
	    print $aln->psl($query, $target);
	    print "\n";
	}elsif(/gff/i){
	    if(!$species){ $species = "??"; }
	    print $aln->gff($query, $target, {chr=>$chr, sp=>$species});
	}elsif(/alignment/i){
	    print $aln->print_alignment($query, $target);
	}elsif(/seed/i){
	    print $aln->getSeed($query, $target, $qhead);
	}
    }
}

sub getContent
{
    my ($filename) = @_;
    open(CONTENT, $filename) or die("Cannot open $filename!");
    my $head = <CONTENT>;
    chomp($head);
    my $body = "";
    while((my $line=<CONTENT>)){
	chomp($line);
	$body.=$line;
    }
    return ($head, $body);
}
