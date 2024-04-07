#! /usr/bin/perl -w

use strict;
use Getopt::Long;

my ($targetfn, $queryfn, $outputdir);
my ($prefix, $cross, $seed);
my ($vulgar, $gtf, $gff, $psl, $alignment, $help);
my ($chr, $species, $gid, $tid);

GetOptions("t:s"=>\$targetfn,
	   "q:s"=>\$queryfn,
	   "o:s"=>\$outputdir,
	   "p:s"=>\$prefix,
	   "a:s"=>\$seed,
	   "cross!"=>\$cross,
	   "vulgar!"=>\$vulgar,
	   "gtf!"=>\$gtf,
	   "gff!"=>\$gff,
	   "psl!"=>\$psl,
	   "alignment!"=>\$alignment,
	   "help!"=>\$help,

	   #output specific parameters:
	   "c:s"=>\$chr,
	   "s:s"=>\$species,
	   "gid:s"=>\$gid,
	   "tid:s"=>\$tid,
    );

if(!$targetfn || !$queryfn || !$outputdir || $help){
    print 
	"Usage: $0 {-t target_file} {-q query_file} \n".
	"          {-o output_dir} [-p prefix] \n".
	"          [-a (GMap | seed_alignment)]\n".
	"          [-cross] [-vulgar] [-gtf] [-gff]\n".
	"          [-psl] [-alignment]\n".
	"          [output params]\n\n".
	"    -t         Target/Genome fasta filename\n".
	"    -q         Query/cDNA fasta filename\n".
	"    -o         Output directory\n".
	"    -p         Prefix for output files. Default is query filename.\n".
	"    -a GMap    Use GMap to generate seed alignment\n".
	"    -a aln     Use the alignment with filename aln for a seed.\n".
	"    -cross     Use cross species parameter\n".
	"    -vulgar    Outputs vulgar format\n".
	"    -gtf       Outputs GTF \n".
	"    -gff       Outputs GFF \n".
	"    -psl       Outputs PSL \n".
	"    -alignment Outputs base by base alignment\n\n".
	"    (for output params, see ./bin/alignmentConvert.pl)\n";
    
    
    exit;
}

if(!$prefix){
    $queryfn =~ /\/?([^\/]*)\Z/;
    $prefix = $1;
}

my $seedfile = "$outputdir/$prefix.seed";
my $pairfile = "$outputdir/$prefix.pair";
my $seedflag = "";

if($seed){
    my $seedCmd;
    $_ = $seed;
    if(/GMap/i){
	$seedCmd = "gmap -g $targetfn $queryfn --format=3 | ./bin/alignmentConvert.pl -x GMap -o seed -t $targetfn -q $queryfn > $seedfile";
    }else{
	$seedCmd = "./bin/alignmentConvert.pl -i $seed -o seed -t $targetfn -q $queryfn > $seedfile";
    }
    print "Creating seed...\n";
    print "$seedCmd\n";
    `$seedCmd`;
    $seedflag = "--seed=$seedfile";
}

my $paramfile = "parameterfs/";
if($cross){
    $paramfile .= "pairagonx.zhmm";
}else{
    $paramfile .= "pairagon.zhmm";
}

my $pCmd = "./bin/pairagon $paramfile $queryfn $targetfn $seedflag --splice_mode=cdna -o -i > $pairfile";
print "Running Pairagon\n";
`$pCmd`;

my $extraparams="";
if($chr){ $extraparams.="-c $chr "; }
if($species){ $extraparams.="-s $species "; }
if($gid){ $extraparams.="-gid $gid "; }
if($tid){ $extraparams.="-tid $tid "; }


if($vulgar){
    print "Outputting vulgar\n";
    `./bin/alignmentConvert.pl -i $pairfile -o vulgar > $outputdir/$prefix.vulgar`;
}

if($gtf){
    print "Outputting GTF\n";
    my $program;
    if($cross){
	$program = "PairagonX";
    }else{
	$program = "Pairagon";
    }
    `./bin/alignmentConvert.pl -i $pairfile -o gtf -p $program $extraparams > $outputdir/$prefix.gtf`;
}
if($gff){
    print "Outputting GFF\n";
    `./bin/alignmentConvert.pl -i $pairfile -o gff -q $queryfn -t $targetfn $extraparams > $outputdir/$prefix.gff`;
}
if($psl){
    print "Outputting PSL\n";
    `./bin/alignmentConvert.pl -i $pairfile -o psl -q $queryfn -t $targetfn > $outputdir/$prefix.psl`;
}
if($alignment){
    print "Outputting Alignment\n";
    `./bin/alignmentConvert.pl -i $pairfile -q $queryfn -t $targetfn -o alignment > $outputdir/$prefix.alignment`;
}

print "Done!\n";
