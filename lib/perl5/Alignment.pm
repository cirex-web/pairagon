package Alignment;

use strict;

sub new
{
    my ($self, 
	$qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref)=@_;

    if(!defined($qstart)||$qstart==-1){
	$qstart=0;
	$qend=0;
	$tstart=0;
	$tend=0;
	$qstrand = "+";
	$tstrand = "+";
    }
    if(!defined($score)){
	$score = 0;
    }
    $self = {
	qid => $qid,
	qstart=>$qstart,
	qend=>$qend,
	qstrand=>$qstrand,
	tid => $tid,
	tstart=>$tstart,
	tend=>$tend,
	tstrand=>$tstrand,
	score=>$score,
	sectionsref=>$sectionsref,
    };
    bless $self, 'Alignment';
    return $self;
}

#
# Usage:
# Alignment::read("filename") 
#     read from file, infer type from extension
# Alignment::read("type")
#     read from stdin, use given type
# Alignment::read("filename", "type") 
#     read from file, use given type
#
sub read
{
    my ($param1, $param2) = @_;
    my ($file, $type);
    if(defined($param1)){
	if(defined($param2)){
	    $file = $param1;
	    $type = $param2;
	}elsif(-f $param1){
	    $file = $param1;
	    $type = getType($file);
	}else{
	    $type = $param1;
	}
    }else{
	die "No parameters specified";
    }

    my $fh;
    if(defined($file)){
	open(FH, $file) or die "Cannot open $file!";
	$fh = \*FH;
    }else{
	$fh = \*STDIN;
    }
    
    my $self;
    
    if($type eq "Est2Genome"){
	$self = parseEstgen($fh);
    }elsif($type eq "Pairagon" || $type eq "SPairagon" || $type eq "Paragorn"){
	$self = parsePairagon($fh);
    }elsif($type eq "GMap"){
	$self = parseGFF($fh);
    }elsif($type eq "Exalin"){
	$self = parseExalin($fh);
    }elsif($type eq "Spidey"){
	$self = parseSpidey($fh);
    }elsif($type eq "Blat" || $type eq "Palma"){
	$self = parseBlat($fh);
    }elsif($type eq "Sim4"){
	$self = parseSim4($fh);
    }elsif($type eq "Exonerate" || $type eq "vulgar"){
	$self = parseVulgar($fh);
    }elsif($type eq "Splign"){
	$self = parseSplign($fh);
    }elsif($type eq "Spaln"){
	$self = parseSpaln($fh);
    }elsif($type eq "Blast"){
	$self = parseBlast($fh);
    }elsif($type eq "Xat"){
	$self = parseXat($fh);
    }elsif($type eq "Geneseqer"){
	$self = parseGeneseqer($fh);
    }else{
	die "'$type' not recognized!";
    }

    close($fh);
    
    return $self;	
}

sub getType
{
    $_ = shift(@_);
    if(/pair\Z/) { return "Pairagon" }
    elsif(/estgen\Z/) { return "Est2Genome" }
    elsif(/gff\Z/) { return "GMap" }
    elsif(/exalin\Z/) { return "Exalin" }
    elsif(/sim4\Z/) { return "Sim4" }
    elsif(/spidey\Z/) { return "Spidey" }
    elsif(/psl\Z/) { return "Blat" }
    elsif(/splign\Z/) { return "Splign" }
    elsif(/spaln\Z/) { return "Spaln" }
    elsif(/vulgar\Z/) { return "vulgar" }
    elsif(/lab\Z/) { return "Blast" }
    elsif(/xat\Z/) { return "Xat" }
    elsif(/geneseqer\Z/) { return "Geneseqer" }
    else{
	die "Cannot recognize type '$_'";
    }
}

sub parseEstgen
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;
    my $parsingAlignment = 0;
    my ($s1, $s2, $a);
   
    my $qlen = 0;

    while(<$fh>){
	if(/\# cDNA.*, (\d+)bp/){
	    $qlen = $1;
	}elsif(/between .* est and .* genome/){
	    ($qstrand, $tstrand) = read_note($_);
	}elsif(/Span/){
	    my @x = split(/\s+/);
	    $score = $x[1];
	    
	    $tstart = $x[3]-1;
	    $tend = $x[4];
	    $tid = $x[5];
	
	    if($tstrand eq "+" && $qstrand eq "+"){
		$qstart = $x[6]-1;
		$qend = $x[7];
	    }elsif($tstrand eq "-" && $qstrand eq "-"){
		$qstart = $x[6]-1;
		$qend = $x[7];
	    }else{
		$qstart = $qlen - $x[7];
		$qend = $qlen - $x[6]+1;
	    }
	    $qid = $x[8];
	    
	    
	}elsif($parsingAlignment){
	    if(/\A(\s*$tid\s*\d+\s)(\S+)\s*\d*/){
		$s1 .= $2;
		my $tab = length($1);
		my $len = length($2);
		$a .= substr(<$fh>, $tab, $len);
		$s2.= substr(<$fh>, $tab, $len);
	    }
	}else{
	    if(/\A.* vs .*:/){
		$parsingAlignment = 1;
	    }
	}
    }
    close(FILE);
    
    my $site = "";
    my $donor;
    
    if($qstrand eq "-" && $tstrand ne "-"){
	$s1 = reverse($s1);
	$s2 = reverse($s2);
	$a = reverse($a);
    }

    my $count = 0;
    my $label = ""; 
    my $ilength = "";

    for(my $i=0;$i<length($s1);$i++){
	my $c1 = substr($s1, $i, 1);
	my $c2 = substr($s2, $i, 1);
	my $ca = substr($a, $i, 1);
	
	my $newlabel = "";

	if($ca eq "|"){
	    $newlabel = "M";
	}elsif($c2 eq "-"){
	    $newlabel = "G2";
	}elsif($c1 eq "-"){
	    $newlabel = "G1";
	}elsif($c1 ne "." && $c2 ne "."){
	    $newlabel = "M";
	}else{
	    $_ = $ca;
	    
	    if(/\d/){
		$ilength.=$ca;
	    }
	    $newlabel = "I";
	}

	if($label eq $newlabel && $newlabel eq "I"){

	}elsif($label eq $newlabel || $label eq ""){
	    $count++;
	}else{
	    my @triple;
	    if($label eq "G1"){
		add_section(\@sections, "G", $count, 0);
		
	    }elsif($label eq "G2"){
		add_section(\@sections, "G", 0, $count);
	    }elsif($label eq "I"){
		if($qstrand eq "-" && $tstrand ne "-"){
		    $ilength = reverse $ilength;
		}
		add_intron(\@sections, $ilength, $tstrand);
		$ilength = "";
	    }else{
		add_section(\@sections, $label, $count, $count);
	    }
	    	    
	    $count = 1;
	}
	$label = $newlabel;
    }
    
    if($label eq "G1"){
	add_section(\@sections, "G", $count, 0);
    }elsif($label eq "G2"){
	add_section(\@sections, "G", 0, $count);
    }elsif($label eq "I"){
	add_section(\@sections, $label, 0, $count);
    }else{
	add_section(\@sections, $label, $count, $count);
    }

    if($tstrand eq "-" || $qstrand eq "-"){
	@sections = reverse(@sections);
    }
    
    $sectionsref = \@sections;

    return new Alignment($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
}

sub parsePairagon
{
    my $fh = shift(@_);
    
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;
    $qstart = -1;
    my $ibound = -1;
    my $qlen = 0;
    while(<$fh>){
	
	if(/\A[^\#]\S+\s+([-\d]*)/){
	    $score+=$1;
	}
	
	if(/\# Genomic Sequence: >(\w+)\s*/){
	    $tid = $1;
	}elsif(/\# Note Best alignment is between (.*) est and (.*) genome/){
	    ($qstrand, $tstrand) = read_note($_);
	}elsif(/cDNA    Sequence: >(\S+), (\d+)bp/){
	    $qid = $1;
	    $qlen = $2;
	}elsif(/Match/){
	    
	    my @fields = split(/[\t ]+/);
	    my $length = $fields[4]-$fields[3]+1;
	    add_section(\@sections, "M", $length, $length);
	    
	    if($tstrand eq "+"){
		if($qstart < 0){
		    $tstart = $fields[3]-1;
		    $qstart = $fields[6]-1;
		}
		$tend = $fields[4];
		$qend = $fields[7];
	    }else{
		if($qstart < 0){
		    $tstart = $fields[3]-1;
		    $qend = $qlen - $fields[6]+1;
		}
		$tend = $fields[4];
		$qstart = $qlen - $fields[7];
	    }
	    
	
	}elsif(/\AGenomic/){
	    my @fields = split(/[\t ]+/);
	    my $length = $fields[4]-$fields[3]+1;
	    add_section(\@sections, "G", 0, $length);
	    $tend = $fields[4];	    
	}elsif(/\ACDna/){
	    my @fields = split(/[\t ]+/);
	    my $length = $fields[4]-$fields[3]+1;
	    add_section(\@sections, "G", $length, 0);
	    $qend = $fields[4];
	}elsif(/Donor/){
	    my @fields = split(/[\t ]+/);
	    if($tstrand eq "+"){
		$ibound = $fields[3];
	    }else{
		my $iend = $fields[4];
		my $length = $iend - $ibound + 1;
		add_intron(\@sections, $length, $tstrand);		
	    }
	}elsif(/\AAcc/){
	    my @fields = split(/[\t ]+/);

	    if($tstrand eq "+"){
		my $iend = $fields[4];
		my $length = $iend - $ibound + 1;
		add_intron(\@sections, $length, $tstrand);	
	    }else{
		$ibound = $fields[3];
	    }
	}
	
    }
    if($tstrand eq "-"){
	@sections = reverse(@sections);
    }
    
    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub parseGFF
{
    my $fh = shift(@_);
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    $qstart = $qend = $tstart = $tend = -1;

    my @lines;
    while(<$fh>){
	if(!/\A\#/){
	    push(@lines, $_);
	}
    }
    
    $lines[0]=~ /\S+\t(\S+).*\t(\+|\-)\t.*Name=([^;]+);/;
    my $tfile = $1;
    $qstrand = "+";
    $tstrand = $2;
    $qid = $3;
    $tid = get_header_text($tfile);

    my $n = scalar @lines;
    
    for(my $i=0;$i<$n;$i++){
	my $line;			
	if($tstrand eq "+"){
	    $line = $lines[$i];
	}else{
	    $line = $lines[$n-$i-1];
	}
	
	$line=~ /cDNA_match\t(\d+)\t(\d+)\t(\d+).*\s(\d+) (\d+);Gap=(.*)\Z/;
	my ($s1, $f1, $segmentscore, $s2, $f2) = ($1, $2, $3, $4, $5);
	
	my @codes = split(/\s+/, $6);

	if($tstart < 0){
	    $tstart = $s1 -1;
	    if($tstrand eq "+"){
		$qstart = $s2 - 1;
	    }else{
		$qend = $f2;
	    }
	}else{
	    my $tgap = $s1-$tend-1;
	    add_intron(\@sections, $tgap, $tstrand);
	    
	    my $qgap;
	    if($tstrand eq "+"){
		$qgap = $s2 - 1 - $qend;
	    }else{
		$qgap = $qstart - $f2;
	    }
	   
	    if($qgap > 0){
		add_section(\@sections, "G", $qgap, 0);
	    }
	}
	
	$tend = $f1;
	
	if($tstrand eq "+"){
	    $qend = $f2;
	}else{
	    $qstart = $s2-1;
	    @codes = reverse @codes;
	}

	$score += $segmentscore;

	for my $code (@codes){
	    $_ = $code;
	    if(/M(\d+)/){
		add_section(\@sections, "M", $1, $1);
	    }elsif(/I(\d+)/){
		add_section(\@sections, "G", $1, 0);
	    }elsif(/D(\d+)/){
		add_section(\@sections, "G", 0, $1);
	    }
	}	
    }

    if($tstrand eq "-"){
	@sections = reverse @sections;
    }
    $sectionsref = \@sections;
    
    my $self =  new Alignment($qid, $qstart, $qend, $qstrand,
				   $tid, $tstart, $tend, $tstrand,
				   $score, $sectionsref);
    
    
    return $self;
}

sub parseExalin
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;
    my ($s1, $s2);
    my @intronlengths;
    my $intronindex = 0;
    $qstart = $qend = -1;
    $tstart = $tend = -1;

    while(<$fh>){
	
	if(/T-SEQUENCE:\s+(.*)/){
	    $qid = $1;
	}elsif(/G-SEQUENCE:\s+(\S*)/){
	    $tid = $1;
	}elsif(/SCORE: ([\d\.]*) nats/){
	    $score = $1;
	}elsif(/SPLICE-DIR: (.)/){
	    $qstrand = $1;
	}elsif(/T-STRAND: (.)/){
	    $tstrand = $1;
	}elsif(/EXON: (\d+) - (\d+) \((\d+) - (\d+)\)/){
	    
	    if($tstrand eq "+"){
		if($qstart < 0){
		    $qstart = $1-1;
		    $tstart = $3-1;
		}else{
		    my $length = $3-1-$tend;
		    my $diff = $1-1-$qend;
		    if($diff != 0){
			$length = $length-$diff-1;
		    }
		    push(@intronlengths, $length);
		}
		$qend = $2;
		$tend = $4;
	    }else{
		
		if($qend < 0){
		    $qend = $1;
		    $tstart = $3-1;
		}else{
		    my $length = $3-1-$tend;
		    
		    my $diff = $qstart-$1;
		    if($diff != 0){
			$length = $length-$diff-1;
		    }
		    push(@intronlengths, $length);
		}
		$qstart = $2-1;
		$tend = $4;
	    }
	}elsif(/T:\s+(\S+)/){
	    $s1 .= $1;
	}elsif(/G:\s+(\S+)/){
	    $s2 .= $1;
	}
    }
    my $label;
    my $count = 0;
    for(my $i=0;$i<length($s1);$i++){
	my $c1 = substr($s1, $i, 1);
	my $c2 = substr($s2, $i, 1);
		
	my $newlabel = "";

	if($c2 eq "-"){
	    $newlabel = "G2";
	}elsif($c1 eq "-"){
	    $newlabel = "G1";
	}elsif($c1 ne "." && $c2 ne "."){
	    $newlabel = "M";
	}else{
	    $newlabel = "I";
	}

	if($label eq $newlabel && $newlabel eq "I"){

	}elsif($label eq $newlabel || $label eq ""){
	    $count++;
	}else{
	    if($label eq "G2"){
		add_section(\@sections, "G", $count, 0);
	    }elsif($label eq "G1"){
		add_section(\@sections, "G", 0, $count);
	    }elsif($label eq "I"){
		add_intron(\@sections, $intronlengths[$intronindex++], $tstrand);
	    }else{
		add_section(\@sections, $label, $count, $count);
	    }
	    	    
	    $count = 1;
	}
	$label = $newlabel;
    }
    
    if($label eq "G1"){
	add_section(\@sections, "G", $count, 0);
    }elsif($label eq "G2"){
	add_section(\@sections, "G", 0, $count);
    }elsif($label eq "I"){
	add_section(\@sections, $label, 0, $count);
    }else{
	add_section(\@sections, $label, $count, $count);
    }

    if($tstrand eq "-"){
	@sections = reverse(@sections);
    }

    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub parseSim4
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;
    my ($tseq, $qseq, $aseq) = ("", "", "");
    my @intronlengths;
    my $intronindex = 0;
    my $tlength;
    $qstart = $qend = -1;
    $tstart = $tend = -1;
    $qstrand = $tstrand = "+";
    $score = 0;
    
    my $label = "";
    my $count = 0;
    
    while(<$fh>){
	if(/\A\>(\S+)\s/){
	    $qid = $1;
	    <$fh> =~ /\A\>(\S+)\s/;
	    $tid = $1;
	}elsif(/(\d+)-(\d+)\s+\((\d+)-(\d+)\)/){
	    
	    if($tstrand eq "+"){
		if($qstart < 0){
		    $qstart = $1-1;
		    $tstart = $3-1;
		}else{
		    my $tlength = $3-1-$tend;
		    my $qlength = $1-1-$qend;
		    if($qlength > 0){
			my @x = ($qlength, $tlength);
			push(@intronlengths, \@x);
		    }else{
			push(@intronlengths, $tlength);
		    }
		}
		$qend = $2;
		$tend = $4;
	    }else{
		if($qstart < 0){
		    $qstart = $1-1;
		    $tend = $tlength-($3-1);
		    
		}else{
		    my $ts2 = $tlength-($3-1);
		    my $tlength = $tstart-$ts2;
		    my $qlength = $1-1-$qend;
		    if($qlength > 0){
			my @x = ($qlength, $tlength);
			push(@intronlengths, \@x);
		    }else{
			push(@intronlengths, $tlength);
		    }
		}
		$qend = $2;
		$tstart = $tlength-$4;
	    }
	}elsif(/\(complement\)/){
	    $tstrand = "-";
	}elsif(/\A(\s*\d+ )[\s\.:]*\Z/){
	    my $tab = length($1);
	    my $starti = $1 + 0;

	    if($starti == 0 && 
	       (scalar @sections > 0 || $count > 0)){
		if($label eq "GT"){
		    add_section(\@sections, "G", $count, 0);
		}elsif($label eq "GQ"){
		    add_section(\@sections, "G", 0, $count);
		}else{
		    add_section(\@sections, $label, $count, $count);
		}
		
		$count = 0;
		my $arrref = $intronlengths[$intronindex++];
		my @arr = @$arrref;
		add_section(\@sections, "G", $arr[0], 0);
		add_intron(\@sections, $arr[1], $qstrand);
	    }
	    
	    $qseq = substr(<$fh>, $tab);
	    $aseq = substr(<$fh>, $tab);
	    $tseq = substr(<$fh>, $tab);
	    chomp($qseq);
	    chomp($aseq);
	    chomp($tseq);

	    for(my $i=0;$i<length($tseq);$i++){
		my $tc = substr($tseq, $i, 1);
		my $qc = substr($qseq, $i, 1);
		my $ac = substr($aseq, $i, 1);
		
		my $newlabel = "";
		
		if($ac eq "-" && $qc eq " "){
		    $newlabel = "GQ";
		}elsif($ac eq "-" && $tc eq " "){
		    $newlabel = "GT";
		}elsif($tc ne " " && $qc ne " "){
		    $newlabel = "M";
		}else{
		    $newlabel = "I";
		}
		
		if($label eq $newlabel && $newlabel eq "I"){
		    
		}elsif($label eq $newlabel || $label eq ""){
		    $count++;
		}else{
		    if($label eq "GT"){
			add_section(\@sections, "G", $count, 0);
		    }elsif($label eq "GQ"){
			add_section(\@sections, "G", 0, $count);
		    }elsif($label eq "I"){
			add_intron(\@sections, $intronlengths[$intronindex++], $qstrand);
		    }else{
			add_section(\@sections, $label, $count, $count);
		    }
	    	    
		    $count = 1;
		}
		$label = $newlabel;
	    }
	}elsif(/seq2 = .*, (\d+) bp/){
	    $tlength = $1;
	}
    }
        
    if($label eq "GT"){
	add_section(\@sections, "G", $count, 0);
    }elsif($label eq "GQ"){
	add_section(\@sections, "G", 0, $count);
    }elsif($label eq "I"){
	add_section(\@sections, $label, 0, $count);
    }else{
	add_section(\@sections, $label, $count, $count);
    }

    if($qstrand eq "-"){
	@sections = reverse @sections;
    }

    $sectionsref = \@sections;
    
    my $self =  new Alignment($qid, $qstart, $qend, $qstrand,
				   $tid, $tstart, $tend, $tstrand,
				   $score, $sectionsref);
    if($qstrand eq "-"){
	#$self->reverse();
    }
    return $self;
}
sub parseBlat
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    my @fields;

    while(<$fh>){
	chomp($_);
	@fields = split(/\t/, $_);
	if(scalar @fields > 16 && !/match/){
	    
	    last;
	}
    }

    
    $qid = $fields[9];
    $tid = $fields[13];
    $qstart = $fields[11];
    $qend = $fields[12];
    my $tsize = $fields[14];
    $tstart=$fields[15];
    $tend = $fields[16];
    $tstrand = "+";
    $qstrand = $fields[8];
    
    $score = 0;
    my @sizes = split(/,/, $fields[18]);
    my @qstarts = split(/,/, $fields[19]);
    my @tstarts = split(/,/, $fields[20]);
    
    my $qlast = -1;
    my $tlast = -1;
    
    # for palma formatting
    if($tstarts[0] > $tstarts[-1]){
	for(my $i=0;$i<scalar @tstarts;$i++){
	    my $ts = $tstarts[$i];
	    $ts = $tsize - $ts - $sizes[$i];
	    $tstarts[$i] = $ts;
	}
	
	my $temp = $tstart;
	$tstart = $tsize - $tend;
	$tend = $tsize - $temp;
	
    }


    for(my $i=0;$i<scalar @sizes;$i++){
	my $size = $sizes[$i];
	my $qs = $qstarts[$i];
	my $ts = $tstarts[$i];
		
	if($tlast > 0){
	    my $tgap = $ts-$tlast;
	    my $qgap = $qs-$qlast;
	    
	    if($tgap > 0){
		add_intron(\@sections, $tgap, $qstrand);
	    }

	    if($qgap > 0){
		add_section(\@sections, "G", $qgap, 0);
	    }
	}

	add_section(\@sections, "M", $size, $size);
	$tlast = $ts + $size;
	$qlast = $qs + $size;
    }

    $sectionsref = \@sections;
    my $self =  new Alignment($qid, $qstart, $qend, $qstrand,
				   $tid, $tstart, $tend, $tstrand,
				   $score, $sectionsref);
    if($qstrand eq "-"){
	$self->reverse();
    }
    return $self;

    
}
sub parseSpidey
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;
    
    my ($tseq, $qseq) = ("", "");
    my @intronlengths;
    my $intronindex = 0;
    $tid = "";
    $qstart = $qend = -1;
    $tstart = $tend = -1;
    $qstrand = $tstrand = "+";
    $score = 0;
    
    #<$fh> =~ /Genomic: .*\|(\S+)\s/;
    #$tid = $1;
    #<$fh> =~ /mRNA: .*\|(\S+)\s/;
    #$qid = $1;

    my $state = 0;
    while(<$fh>){
	if(/Genomic: .*\|(\S+)\s/){
	    $tid = $1;
	}elsif(/mRNA: .*\|(\S+)\s/){
	    $qid = $1;
	}elsif(/Exon \d+: (\d+)-(\d+) \(gen\)\s+(\d+)-(\d+) \(mRNA\)\Z/){
	    
	    if($qstart < 0){
		if($1 > $2){
		    $tstrand = "-";
		}

		if($3 > $4){
		    $qstrand = "-";
		}
	    }

	    my ($qgap, $tgap);
	    
	    if($tstrand eq "+" && $qstrand eq "+"){
		if($qstart < 0){
		    $tstart = $1 - 1;
		    $qstart = $3 - 1;
		}else{
		    $tgap = $1-1-$tend;
		    $qgap = $3-1-$qend;
		}
		$tend = $2;
		$qend = $4;
	    }elsif($tstrand eq "+" && $qstrand eq "-"){
		if($qstart < 0){
		    $tstart = $1;
		    $qend = $3;
		}else{
		    $tgap = $1-$tend;		    
		    $qgap = $qstart-$3;
		}
		$tend = $2+1;
		$qstart = $4 - 1;
	    }elsif($tstrand eq "-" && $qstrand eq "+"){
		if($qstart < 0){
		    $tend = $1;
		    $qstart = $3 - 1;
		}else{
		    $tgap = $tstart-$1;
		    $qgap = $3-$qend-1;
		}
		$tstart = $2-1;
		$qend = $4;
	    }elsif($tstrand eq "-" && $qstrand eq "-"){
		if($qstart < 0){
		    $tend = $1;
		    $qend = $3;
		}else{
		    $tgap = $tstart-$1;
		    $qgap = $qstart-$3;
		}
		$tstart = $2-1;
		$qstart = $4-1;
	    }
	    
	    if($tgap < 0){
		$tseq = substr($tseq, 0, length($tseq)+$tgap);
		$qseq = substr($qseq, 0, length($qseq)+$tgap);
		$qgap -= $tgap;
	    }
	    
	    if($qgap < 0){
		$tseq = substr($tseq, 0, length($tseq)+$qgap);
		$qseq = substr($qseq, 0, length($qseq)+$qgap);
		$tgap -= $qgap;
	    }

	    if($tgap > 0){
		$tseq = substr($tseq, 0, length($qseq));
		$tseq .= "I";
		$qseq .= "I";
		push(@intronlengths, $tgap);
	    }

	    if($qgap > 0){
		$tseq = substr($tseq, 0, length($qseq));
		$tseq .= "-"x$qgap;
		$qseq .= "N"x$qgap;
	    }

	    $state = 0;
	}elsif(/Genomic: /){
	    last;
	}elsif($state==0 &&
	       /(\A\s*[ACTG\-]+)\Z/){
	    $tseq .= $1;
	    $state = 1;
	}elsif($state==1){
	    $state = 2;
	}elsif($state==2 && /(\A\s*[ACTG\-]+)\Z/){
	    $qseq .= $1;
	    $state = 3;
	}elsif($state==3 && /\A\s*\Z/){
	    $state = 0;
	}

    }

    my $label;
    my $count = 0;
    for(my $i=0;$i<length($tseq);$i++){
	my $tc = substr($tseq, $i, 1);
	my $qc = substr($qseq, $i, 1);

	my $newlabel = "";
	if($qc eq "-"){
	    $newlabel = "GQ";
	}elsif($tc eq "-"){
	    $newlabel = "GT";
	}elsif($tc eq "I"){
	    $newlabel = "I";
	}elsif($tc ne " " && $qc ne "" && $qc ne " "){
	    $newlabel = "M";
	}else{
	    $newlabel = "X";
	}

	if($label eq $newlabel && $newlabel eq "X"){
	    $label = $newlabel;
	}elsif($label eq $newlabel){
	    $count++;
	}else{
	    if($label eq "GT"){
		add_section(\@sections, "G", $count, 0);
	    }elsif($label eq "GQ"){
		add_section(\@sections, "G", 0, $count);
	    }elsif($label eq "M"){
		add_section(\@sections, $label, $count, $count);
	    }

	    if($newlabel eq "I"){
		add_intron(\@sections, $intronlengths[$intronindex++], $qstrand);
	    }
	    	    
	    $count = 1;
	}
	$label = $newlabel;
    }
    
    if($label eq "GT"){
	add_section(\@sections, "G", $count, 0);
    }elsif($label eq "GQ"){
	add_section(\@sections, "G", 0, $count);
    }elsif($label eq "I"){
	add_section(\@sections, $label, 0, $count);
    }elsif($label eq "M"){
	add_section(\@sections, $label, $count, $count);
    }

    $sectionsref = \@sections;
   
    my $self =  new Alignment($qid, $qstart, $qend, $qstrand,
				   $tid, $tstart, $tend, $tstrand,
				   $score, $sectionsref);
    if($qstrand eq "-"){
	$self->reverse();
    }
    return $self;
}

sub parseVulgar
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    my $line;
    my $only;
    while(<$fh>){
	if(/\Avulgar: /){
	    $line = $_;
	    last;
	}
	$only = $_;
    }

    if(!$line){
	$line = $only;
    }

    my @parts = split(/\s+/, $line);
    if($parts[0] eq "vulgar:"){
	shift(@parts);
    }
    $qid = shift(@parts);
    $qstart = shift(@parts);
    $qend = shift(@parts);
    $qstrand = shift(@parts);
    $tid = shift(@parts);
    $tstart = shift(@parts);
    $tend = shift(@parts);
    $tstrand = shift(@parts);
    $score = shift(@parts);

    if($qstrand eq "-"){
	my $temp = $qstart;
	$qstart = $qend;
	$qend = $temp;
    }
    if($tstrand eq "-"){
	my $temp = $tstart;
	$tstart = $tend;
	$tend = $temp;
    }
    

    while(@parts){
	my $label = shift(@parts);
	my $n1 = shift(@parts);
	my $n2 = shift(@parts);
	my @triple = ($label, $n1, $n2);
	push(@sections, \@triple);
    }

    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub parseSplign
{
    my $fh = shift(@_);
        
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    $qstart = -1;
    my $saved = "";
    while(<$fh>){
	if(/\A\#/){

	}else{
	    my @parts = split(/\t/);
	    if($parts[9] eq "<L-Gap>"){
		next;
	    }
	    if($qstart < 0){
		$qid = $parts[1];
		$tid = $parts[2];
		$qstart = $parts[5]-1;
		$qstrand = "+";
		$tstart = $parts[7]-1;
		$tend = $parts[8];
		if($tstart < $tend){
		    $tstrand = "+";
		}else{
		    $tstrand = "-";
		}
		$score = 0; # fix?
		$saved = $parts[0];
		
	    }elsif($parts[0] ne $saved){
		last;
	    }elsif($parts[9] eq "<M-Gap>" ||
		   $parts[9] eq "<R-Gap>"){
		add_section(\@sections, "G", $parts[4], 0);
		next;
	    }elsif($tstrand eq "+"){
		my $length = $parts[7]-1-$tend;
		add_intron(\@sections, $length, $qstrand);
	    }else{
		my $length = $tend - 1 - $parts[7];
		add_intron(\@sections, $length, $qstrand);
	    }
	    
	    $qend = $parts[6];
	    $tend = $parts[8];

	    my $count = 0;
	    my @x = split(/([MRID]\d*)/, $parts[10]);
	    for my $code (@x){
		$_ = $code;
		if(/M(\d*)/ || /R(\d*)/){
		    if($1){
			$count += $1;
		    }else{
			$count++;
		    }
		}elsif(/I(\d*)/){
		    if($count>0){
			add_section(\@sections, "M", $count, $count);
			$count = 0;
		    }
		    if($1){
			$count = $1;
		    }else{
			$count = 1;
		    }
		    add_section(\@sections, "G", 0, $count);
		    $count = 0;
		    
		}elsif(/D(\d*)/){
		    if($count>0){
			add_section(\@sections, "M", $count, $count);
			$count = 0;
		    }
		    if($1){
			$count = $1;
		    }else{
			$count = 1;
		    }
		    add_section(\@sections, "G", $count, 0);
		    $count = 0;
		}
	    }
	    
	    if($count>0){
		add_section(\@sections, "M", $count, $count);
	    }
	}
    }
    if($tstrand eq "-"){
	my $temp = $tstart;
	$tstart = $tend -1;
	$tend = $temp + 1;
    }
    
    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub parseSpaln
{
    my $fh = shift(@_);
    
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    my $line = <$fh>;
    $line = <$fh>;
    $line =~ /([\<\>])(\S*).*([\<\>])(\S*)/;
    $qid = $4;
    $tid = $2;
    $qstrand = get_sign($3);
    $tstrand = get_sign($1);
    $tstart = $qstart = -1;
    
    my $tseq = "";
    my ($tlen, $qlen) = (0,0);

    my $label = "M";
    my $count = 0;
    while(<$fh>){

	if(/Score = ([\d\.]+)/){
	    $score = $1;
	}elsif(/(\d+) ([^\|]*)\| $tid\Z/){
	    if($tstart < 0){
		$tstart = $1 - 1;
	    }
	    
	    $tseq = $2;
	}elsif(/(\d+) (.*)\| $qid\Z/){
	    if($qstart < 0){
		$qstart = $1 - 1;
	    }
	    my $qseq = $2;
	    for(my $i=0;$i<length($tseq);$i++){
		my $tc = substr($tseq, $i, 1);
		my $qc = substr($qseq, $i, 1);
		my $newlabel ;
		if($tc eq "-" || ($tc eq " " && $qc ne " ")){
		    $newlabel = "G1";
		    $qlen++;
		}elsif($qc eq "-"){
		    $newlabel = "G2";
		    $tlen++;
		}elsif($tc eq " "){
		    $newlabel = "X";
		}elsif($qc eq " "){
		    $newlabel = "I";
		    $tlen++;
		}else{
		    $newlabel = "M";
		    $tlen++;
		    $qlen++;
		}

		if($label eq $newlabel && $newlabel eq "I"){
		    $count++;
		}elsif($label eq $newlabel || $label eq ""){
		    $count++;
		}elsif($label ne "X"){
		    if($label eq "G1"){
			add_section(\@sections, "G", $count, 0);
		    }elsif($label eq "G2"){
			add_section(\@sections, "G", 0, $count);
		    }elsif($label eq "I"){
			add_intron(\@sections, $count, $qstrand);
		    }else{
			add_section(\@sections, $label, $count, $count);
		    }
	    	    
		    $count = 1;
		}
		$label = $newlabel;
	    }
	    
	}elsif(/;; skip (\d+) nt\'s/){
	    $count += $1;
	    $tlen += $1;
	}
    }

    if($count > 1 && $label eq "M"){
	add_section(\@sections, $label, $count, $count);
    }

    if($qstrand eq "+"){
	$qend = $qstart + $qlen;
    }else{
	$qend = $qstart+1;
	$qstart = $qend - $qlen;
    }

    if($tstrand eq "+"){
	$tend = $tstart + $tlen;
    }else{
	$tend = $tstart + 1;
	$tstart = $tend - $tlen;
    }

    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub parseBlast
{
    my $fh = shift(@_);
    
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    my $aflag = 0;
    my $qlast = -1;
    my $tlast = -1;
    while(<$fh>){
	if(/h \{/){
	    $_ = <$fh>;
	    
	    if(/>(\S+)\s/){
		$tid = $1;
	    }
	    if(/reverse/){
		$tstrand = "-";
	    }else{
		$tstrand = "+";
	    }
	    
	    $_ = <$fh>;
	    if(/>([^\s\"]+)/){
		$qid = $1;
	    }

	    if(/reverse/){
		$qstrand = "-";
	    }else{
		$qstrand = "+";
	    }
	}elsif(/a \{/ && $aflag == 0 ){
	    $aflag = 1;
	}elsif($aflag == 1){
	    if(/s (\d+)/){
		$score = $1;
	    }elsif(/b (\d+) (\d+)/){
		$tstart = $1-1;
		$qstart = $2-1;
		
	    }elsif(/e (\d+) (\d+)/){
		$tend = $1;
		$qend = $2;
	    }elsif(/l (\d+) (\d+) (\d+) (\d+)/){
		my ($ts, $qs, $te, $qe) = ($1, $2, $3, $4);
		my $size = $qe-$qs;
		
		if($tlast > 0){

		    if(1==1){
			my $qlen = $qs-$qlast;
			my $tlen = $ts-$tlast;
			if($qlen>0){
			    add_section(\@sections, "G", $qlen, 0);
			}

			if($tlen>0){
			    add_section(\@sections, "G", 0, $tlen);
			}
		    }else{
		    if($tlast == $ts){
			my $length = $qs-$qlast;
			add_section(\@sections, "G", $length, 0);
		    }elsif($ts-$tlast > 30){
			my $length = $ts-$tlast;
			add_intron(\@sections, $length, $qstrand);
		    }else{
			my $length = $ts - $tlast;
			add_section(\@sections, "G", 0, $length);
		    }}
		}
		
		add_section(\@sections, "M", $size, $size);
		$tlast = $ts + $size;
		$qlast = $qs + $size;
		
	    }elsif(/\}/){
		$aflag = 2;
	    }
	    
	}
    }

    $sectionsref = \@sections;
    my $self =  new Alignment($qid, $qstart, $qend, $qstrand,
				   $tid, $tstart, $tend, $tstrand,
				   $score, $sectionsref);
    if($qstrand eq "-"){
	$self->reverse();
    }
    return $self;
}

sub parseXat
{
    my $fh = shift(@_);
    
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;
    
    $qstart = -1;
    $qid = -1;
    
    while(<$fh>){
	if(/seq1 = .* \((.*)\)/){
	    if($qid < 0){
		$qid = $1;
	    }else{
		last;
	    }
	}elsif(/seq2 = .* \((.*)\).*score=(.*)\Z/){
	    $tid = $1;
	    $score = $2;
	}elsif(/\[(.)\] (\d+) - (\d+)\s+(\d+) - (\d+).*, (.*)\Z/){
	    if($qstart < 0){
		if($4 < $5){
		    $tstrand = "+";
		}else{
		    $tstrand = "-";
		}
		$qstrand = "+";
	    }

	    if($tstrand eq "+"){
		if($qstart < 0){
		    $qstart = $2-1;
		    $tstart = $4-1;
		}else{
		    my $qgap = $2-1-$qend;
		    my $tgap = $4-1-$tend;
		    if($tgap>0){
			add_intron(\@sections, $tgap, $qstrand);
		    }
		    
		    if($qgap>0){
			add_section(\@sections, "G", $qgap, 0);
		    }
		}
		$qend = $3;
		$tend = $5;
	    }else{
		
		if($qstart < 0){
		    $qstart = $2-1;
		    $tend = $4;
		}else{
		    my $qgap = $2-1-$qend;
		    my $tgap = $tstart-$4;
		    
		    if($tgap>0){
			add_intron(\@sections, $tgap, $qstrand);
		    }
		    
		    if($qgap>0){
			add_section(\@sections, "G", $qgap, 0);
		    }
		}
		$qend = $3;
		$tstart = $5 - 1;
	    }

	    my @parts = split(/\s+/, $6);
	    for my $part (@parts){
		$part =~ /(.)-(\d+)/;
		if($1 eq "S"){
		    add_section(\@sections, "M", $2, $2);
		}elsif($1 eq "D"){
		    add_section(\@sections, "G", 0, $2);
		}elsif($1 eq "I"){
		    add_section(\@sections, "G", $2, 0);
		}
	    }
	    
	    
	}
    }

    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub parseGeneseqer
{
    my $fh = shift(@_);
    
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my ($tsx, $qsx, $tex, $qex);
    my $flag = 0;
    
    my @sections;

    while(<$fh>){
	if(/Exon\s+(\d+)\s+(\d+)\s+(\d+) .* cDNA\s+(\d+)\s+(\d+)/){
	    if($1==1){
		$tsx = $2-1;
		$qsx = $4-1;
	    }
	    $tex = $3;
	    $qex = $5;
	}elsif(/\MATCH\s+(\S*)(\+|\-)\s*(\S*)(\+|\-)\s+([\d\.]+)\s/){
	    if(!defined($score) || $5 > $score){
		$tid = "chr$1";
		$tstrand = $2;
		$qid = "chr$3";
		$qstrand = $4;
		$score = $5;
		$qstart = $qsx;
		$qend = $qex;

		if($tstrand eq "+"){
		    $tstart = $tsx;
		    $tend = $tex;
		}else{
		    $tstart = $tex-1;
		    $tend = $tsx+1;
		}
		@sections = ();
		$flag = 1;
	    }else{
		$flag = 0;
	    }
	}elsif(/\AAlignment/ && $flag){
	    my $label;
	    my $count = 0;

	    while(1==1){
		<$fh>; #blank line;
		my $tseq = <$fh>;
		my $aseq = <$fh>;
		my $qseq = <$fh>;
		<$fh>; #blank line;

		$_ = $tseq;
		if(/PGS/){
		    last;
		}

		my $n = length($aseq);
		for(my $i=0;$i<$n;$i++){
		    my $tc = substr($tseq, $i,1);
		    my $qc = substr($qseq, $i,1);
		    my $newlabel;

		    if($tc eq " "){
			$newlabel = $label;
		    }elsif($tc eq "-"){
			$newlabel = "G1";
		    }elsif($qc eq "-"){
			$newlabel = "G2";
		    }elsif($qc eq "."){
			$newlabel = "I";
		    }else{
			$newlabel = "M";
		    }

		    if($tc eq " "){

		    }elsif($label eq $newlabel && $newlabel eq "I"){
			$count++;
		    }elsif($label eq $newlabel || $label eq ""){
			$count++;
		    }elsif($label ne "X"){
			if($label eq "G1"){
			    add_section(\@sections, "G", $count, 0);
			}elsif($label eq "G2"){
			    add_section(\@sections, "G", 0, $count);
			}elsif($label eq "I"){
			    add_intron(\@sections, $count, $qstrand);
			}else{
			    add_section(\@sections, $label, $count, $count);
			}
			
			$count = 1;
		    }
		    $label = $newlabel;
		}
		
	    }
	    if($count > 1 && $label eq "M"){
		add_section(\@sections, $label, $count, $count);
	    }
	}
    }

    $sectionsref = \@sections;
    my $self =  new Alignment($qid, $qstart, $qend, $qstrand,
				   $tid, $tstart, $tend, $tstrand,
				   $score, $sectionsref);
    if($qstrand eq "-"){
	$self->reverse();
    }
    return $self;
}

sub parseSomething
{
    my $fh = shift(@_);
    
    my ($qid, $qstart, $qend, $qstrand,
	$tid, $tstart, $tend, $tstrand,
	$score, $sectionsref);
    my @sections;

    $sectionsref = \@sections;
    
    return new Alignment($qid, $qstart, $qend, $qstrand,
			      $tid, $tstart, $tend, $tstrand,
			      $score, $sectionsref);
}

sub fix
{
    my ($self) = @_;
    my @sections = @{$self->{sectionsref}};
    my @newsections;
    my $lastmatch = 0;
    for my $sectionref (@sections){
	my @section = @$sectionref;
	if($section[0] eq "M"){
	    $lastmatch+=$section[1];
	}else{
	    if($lastmatch > 0){
		my @newsection = ("M", $lastmatch, $lastmatch);
		push(@newsections, \@newsection);
		$lastmatch = 0;
	    }
	    
	    push(@newsections, $sectionref);
	}
    }
    if($lastmatch > 0){
	my @newsection = ("M", $lastmatch, $lastmatch);
	push(@newsections, \@newsection);
    }
    $self->{sectionsref} = \@newsections;
}

sub add_section
{
    my ($sectionsref, $label, $n1, $n2) = @_;
    my @triple = ($label, $n1, $n2);
    push(@$sectionsref, \@triple);
}

sub add_intron
{
    my ($sectionsref, $length, $strand) = @_;

    if($length >= 30){
	my ($l1, $l2);
	if($strand eq "+"){
	    ($l1, $l2) = ("5", "3");
	}else{
	    ($l1, $l2) = ("3", "5");
	}
	
	add_section($sectionsref, $l1, 0, 2);
	add_section($sectionsref, "I", 0, $length-4);
	add_section($sectionsref, $l2, 0, 2);
    }else{
	add_section($sectionsref, "G", 0, $length);
    }
}

sub read_note
{
    shift =~ /Best alignment is between (.*) est and (.*) genome, .* imply (\S*) /;
    my $es = get_sign($1);
    my $gs = get_sign($2);
    if($3 eq "REVERSED"){
	$es = flip($es);
	$gs = flip($gs);
    }
    return ($es, $gs);
}

sub get_header_text
{
    my $file = shift(@_);
    open(FILE, $file) or return "";
    my $line = <FILE>;
    $line =~ />(\S+)/;
    close(FILE);
    if(defined($1)){
	return $1;
    }else{
	die "Error parsing target header";
    }
}

sub get_sign
{
    $_ = shift;
    if(/forward/ || /\>/){
	return "+";
    }else{
	return "-";
    }
}

sub flip
{
    $_ = shift;
    if(/\+/){
	return "-";
    }else{
	return "+";
    }
}

sub reverse
{
    my ($self) = @_;
    my @sections = @{$self->{sectionsref}};
    @sections = reverse @sections;
    $self->{sectionsref} = \@sections;
    $self->{qstrand} = flip($self->{qstrand});
    $self->{tstrand} = flip($self->{tstrand});
}

sub vulgar
{
    my ($self) = @_;
    my $base = "$self->{qid} ";
    $base .= vulgar_head($self->{qstart}, $self->{qend}, $self->{qstrand});
    $base .= "$self->{tid} ";
    $base .= vulgar_head($self->{tstart}, $self->{tend}, $self->{tstrand});
    $base .= "$self->{score}\t";
    for my $trip (@{$self->{sectionsref}}){
	my @triple = @$trip;
	
	$base .= join(" ", @triple);
	$base .= " ";
    }
    return "$base\n";
}

sub vulgar_head
{
    my ($start, $end, $strand) = @_;
    if($strand eq "+"){
	return "$start $end $strand ";
    }else{
	return "$end $start $strand ";
    }
}

sub offset
{
    my ($self, $offset) = @_;
    $self->{tstart}+=$offset;
    $self->{tend}+=$offset;
}

sub segments
{
    my ($self) = @_;
    my @segments;
    my ($s1, $e1, $s2, $e2, $m1, $m2);
    if($self->{qstrand} eq "+"){
	$s1 = $e1 = $self->{qstart};
	$m1 = 1;
    }else{
	$s1 = $e1 = $self->{qend};
	$m1 = -1;
    }
    
    if($self->{tstrand} eq "+"){
	$s2 = $e2 = $self->{tstart};
	$m2 = 1;
    }else{
	$s2 = $e2 = $self->{tend};
	$m2 = -1;
    }

    my $flag = 1;

    for my $trip (@{$self->{sectionsref}}){
	my @triple = @$trip;
	my $label = $triple[0];
	
	if($label eq "5" || $label eq "3" || $label eq "I"){
	    if($flag){
		push(@segments, make_segment($s1, $e1, $s2, $e2, $m1, $m2));
		$flag = 0;
		
	    }
	    $e1 += $triple[1]*$m1;
	    $e2 += $triple[2]*$m2;
	}else{
	    if(!$flag){
		$flag = 1;
		$s1 = $e1;
		$s2 = $e2;
	    }
	    $e1 += $triple[1]*$m1;
	    $e2 += $triple[2]*$m2;
	    
	}
	
    }
        
    push(@segments, make_segment($s1, $e1, $s2, $e2, $m1, $m2));
    
    return @segments;
    
}

sub highscoringpairs
{
    my ($self) = @_;
    my @segments;
    my ($s1, $e1, $s2, $e2, $m1, $m2);
    if($self->{qstrand} eq "+"){
	$s1 = $e1 = $self->{qstart};
	$m1 = 1;
    }else{
	$s1 = $e1 = $self->{qend};
	$m1 = -1;
    }
    
    if($self->{tstrand} eq "+"){
	$s2 = $e2 = $self->{tstart};
	$m2 = 1;
    }else{
	$s2 = $e2 = $self->{tend};
	$m2 = -1;
    }

    my $flag = 1;

    for my $trip (@{$self->{sectionsref}}){
	my @triple = @$trip;
	my $label = $triple[0];
	
	if($label eq "M"){
	    $s1 = $e1;
	    $s2 = $e2;
	    $e1 += $triple[1]*$m1;
	    $e2 += $triple[2]*$m2;
	    push(@segments, make_segment($s1, $e1, $s2, $e2, $m1, $m2));
	    
	}else{
	    $e1 += $triple[1]*$m1;
	    $e2 += $triple[2]*$m2;
	}
	
    }
        
    return @segments;
    
}

sub make_segment
{
    my ($s1, $e1, $s2, $e2, $m1, $m2) = @_;
    my ($a, $b, $c, $d);
    if($m1==1){
	$a = $s1;
	$b = $e1 - 1;
    }else{
	$a = $s1 - 1;
	$b = $e1;
    }

    if($m2==1){
	$c = $s2;
	$d = $e2 - 1;
    }else{
	$c = $s2 - 1;
	$d = $e2;
    }
    my @arr = ($a,$b,$c,$d);
    return \@arr;
}

sub gtf
{
    my ($self, $mapref)=@_;
    my %map = %$mapref;
    my $s = "";
    
    my $strand;
    if($self->{tstrand} eq "-" ^ 
       $self->{qstrand} eq "-"){
	$strand = "-";
    }else{
	$strand = "+";
    }

    my @segments = $self->segments();
    if($self->{tstrand} eq "-"){
	@segments = reverse @segments;
    }
    
    for my $segmentref (@segments){
	my @segment = @$segmentref;
	
	$s.="$map{chr}\t$map{program}\tCDS\t";
	
	if($segment[2] <= $segment[3]){
	    $s.="$segment[2]\t";
	    $s.="$segment[3]\t";
	}else{
	    $s.="$segment[3]\t";
	    $s.="$segment[2]\t";
	}

	$s.=".\t$strand\t.\t";
	$s.="gene_id \"$map{gid}\"; transcript_id \"$map{tid}\";";
	$s.="\n";
    }
    return $s;
}

sub psl
{
    my ($self, $query, $target) = @_;
    my @fields;
    my @sections = @{$self->{sectionsref}};
    if(scalar @sections == 0){
	return "";
    }
    my $qsize = length($query);

    my ($tseq, $aseq, $qseq) = $self->alignment($query, $target);
    my ($matches, $mismatches, $ncount, 
	$qNumInserts, $qBaseInserts,
	$tNumInserts, $tBaseInserts) = (0,0,0,0,0,0,0);
    
    for(my $i=0;$i<length($tseq);$i++){
	my $tc = substr($tseq,$i,1);
	my $qc = substr($qseq,$i,1);

	if($qc eq "." || $qc eq "-" || $tc eq "-"){
	    
	}elsif($tc eq $qc){
	    $matches++;
	}else{
	    $mismatches++;
	}
    }

    my (@sizes, @qstarts, @tstarts);
    
    my ($tcount, $qcount) = (0,0);
    my ($ti, $qi) = ($self->{tstart},
		     $self->{qstart});

    my @sections = @{$self->{sectionsref}};
    if($self->{tstrand} eq "-"){
	@sections = reverse @sections;
	$qi += $qsize - $self->{qend} ;
    }

    for my $sectionref (@sections){
	my @section = @$sectionref;
	if($section[0] eq "M"){
	    if($qcount > 0){
		$qNumInserts++;
		$qBaseInserts+=$qcount;
		$qcount=0;
	    }

	    if($tcount>0){
		$tNumInserts++;
		$tBaseInserts+=$tcount;
		$tcount=0;
	    }

	    my $nx;
	    if($self->{tstrand} eq "+"){
		$nx = $qi;
	    }else{
		$nx = $qi - $self->{qstart};
	    }

	    push(@sizes, $section[1]);
	    push(@tstarts, $ti);
	    push(@qstarts, $nx);
	}else{
	    $qcount += $section[1];
	    $tcount += $section[2];
	}

	$qi += $section[1];
	$ti += $section[2];
    }


    $fields[0] = $matches;
    $fields[1] = $mismatches;
    $fields[2] = 0;
    $fields[3] = $ncount;
    $fields[4] = $qNumInserts;
    $fields[5] = $qBaseInserts;
    $fields[6] = $tNumInserts;
    $fields[7] = $tBaseInserts;
    $fields[8] = $self->{tstrand};
    $fields[9] = $self->{qid};
    $fields[10] = $qsize;
    $fields[11] = $self->{qstart};
    $fields[12] = $self->{qend};
    $fields[13] = $self->{tid};
    $fields[14] = length($target);
    $fields[15] = $self->{tstart};
    $fields[16] = $self->{tend};
    $fields[17] = scalar @sizes;
    $fields[18] = join(",", @sizes);
    $fields[19] = join(",", @qstarts);
    $fields[20] = join(",", @tstarts);

    for(my $i=18;$i<=20;$i++){
	if(length($fields[$i])>0){
	    $fields[$i].=",";
	}
    }

    return join("\t", @fields);
}

sub exon_identities
{
    my ($self, $query, $target) = @_;
    my @ids;
    my ($match, $mismatch, $ins, $dels)= (0,0,0,0);
    
    my ($tseq, $aseq, $qseq) = $self->alignment($query, $target);
    
    my $flag = 0;

    for(my $i=0;$i<length($tseq);$i++){
	my $tc = substr($tseq,$i,1);
	my $qc = substr($qseq,$i,1);

	if($qc eq "."){
	    if($flag==0){
		push(@ids, calculate_identity($match, $mismatch, $ins, $dels));
		$flag = 1;
		($match, $mismatch, $ins, $dels) = (0,0,0,0);
	    }
	}else{
	    $flag = 0;
	}

	if($qc eq "."){

	}elsif($qc eq "-"){
	    $ins++;
	}elsif($tc eq "-"){
	    $dels++;
	}elsif($tc eq $qc){
	    $match++;
	}else{
	    $mismatch++;
	}
    }
	
    push(@ids, calculate_identity($match, $mismatch, $ins, $dels));
    return @ids;
}

sub identity
{
    my ($self, $query, $target, $stats) = @_;
    my ($match, $mismatch, $ins, $dels) = (0,0,0,0);
    my ($tseq, $aseq, $qseq) = $self->alignment($query, $target);
    
    for(my $i=0;$i<length($tseq);$i++){
	my $tc = substr($tseq,$i,1);
	my $qc = substr($qseq,$i,1);

	if($qc eq "."){
	    #intron
	}elsif($qc eq "-"){
	    $ins++;
	}elsif($tc eq "-"){
	    $dels++;
	}elsif($tc eq $qc){
	    $match++;
	}else{
	    $mismatch++;
	}
    }

    my $id = calculate_identity($match, $mismatch, $ins, $dels);
    if($stats){
	return ($match, $mismatch, $ins, $dels, $id);
    }else{
	return $id;
    }
}

sub calculate_identity
{
    my ($match, $mismatch, $ins, $dels) = @_;
    
    if($match>0){
	my $id = $match*100/
	    ($match+$mismatch+$ins+$dels);
	return sprintf("%.1f", $id);
    }else{
	return "0.0";
    }
}

sub compare_alignment
{
    my ($self, $other)=@_;
    my @segments1 = $self->segments();
    my @segments2 = $other->segments();
    
    my %exonmap;
    my ($correctexons, $correctbases,
	$totalbases1, $totalbases2) = (0,0,0,0);
    
    for my $segmentref (@segments1){
	my @segment = @$segmentref;
	my $exon = $segment[2] . " " . $segment[3];
	print "$exon\n";
	$exonmap{$exon} = 1;
	$totalbases1 += abs($segment[3]-$segment[2])+1;
    }
    print "\n";
    for my $segmentref (@segments2){
	my @segment = @$segmentref;
	my ($start1, $end1) = ($segment[2], $segment[3]);

	my $exon = "$start1 $end1";
	print "$exon";
	if($exonmap{$exon} == 1){
	    $correctexons++;
	    print "x";
	}
	print "\n";
	
	if($start1>$end1){
	    $end1 = $start1;
	    $start1 = $segment[3];	    
	}
	
	$totalbases2 += $end1-$start1+1;

	
	for my $segmentref2 (@segments1){
	    my @segment2 = @$segmentref2;
	    my ($start2, $end2);
	    if($segment2[2]<@segment2[3]){
		($start2, $end2) = ($segment2[2], $segment2[3]);
	    }else{
		($start2, $end2) = ($segment2[3], $segment2[2]);
	    }
	    
	    if($end2 > $start1 && $start2 < $end1){
		my ($start, $end);
		if($start1 > $start2){
		    $start = $start1;
		}else{
		    $start = $start2;
		}

		if($end1 < $end2){
		    $end = $end1;
		}else{
		    $end = $end2;
		}

		$correctbases+= $end - $start+1;
	    }
	}
    }

    #print "$correctbases\n$totalbases1 $totalbases2\n";
    print scalar @segments1 . " " . scalar @segments2 . "\n";
    print "$correctexons\n";

    my ($exsn, $exsp, $nusn, $nusp);
    if(scalar @segments1 > 0){
	$exsn = 100*$correctexons/scalar @segments1;
    }else{
	$exsn = 0;
    }

    if(scalar @segments2 > 0){
	$exsp = 100*$correctexons/scalar @segments2;
    }else{
	$exsp = 0;
    }

    if($totalbases1 > 0){
	$nusn = 100*$correctbases/$totalbases1;
    }else{
	$nusn = 0;
    }

    if($totalbases2 > 0){
	$nusp = 100*$correctbases/$totalbases2;
    }else{
	$nusp = 0;
    }

    return ($exsn, $exsp, $nusn, $nusp);
    
}

sub btab
{
    my ($self, $mapref)=@_;
    my %map = %$mapref;
    my $s = "";
    my @segments = $self->segments();
    my @identities = @{$map{identities}};
    
    my @args;
    $args[0] = $map{chr};
    $args[3] = $map{program};
    $args[5] = $map{id};

    $args[10] = $map{identity};
    $args[14] = "1";
    $args[1] = $args[2] = $args[4] = $args[11] = $args[12]="";

    my $lineno = $map{lineno};
    my $id_index = 0;

    if($self->{tstrand} eq "-"){
	@segments = reverse @segments;
	@identities = reverse @identities;
    }

    for my $segmentref (@segments){
	my @segment = @$segmentref;

	$args[6] = $segment[2]+1;
	$args[7] = $segment[3]+1;
	$args[8] = $segment[0]+1;
	$args[9] = $segment[1]+1;
	$args[10] = $identities[$id_index++];
	$args[13] = $lineno++;
	$s.= join("\t", @args);
	$s.= "\n";
    }
    return $s;
}

sub gff
{
    my ($self, $query, $target, $mapref) = @_;
    my %map = %$mapref;
    my $gff = "##gff-version   3\n";
    
    my ($chr, $sp);
    if($map{"chr"}){
	$chr = $map{"chr"};
    }else{
	$self->{tid} =~ /chr([A-Za-z\d]+)[^A-Za-z\d]?/;
	$chr = $1;
    }

    if($map{"sp"}){
	$sp = $map{"sp"};
    }else{
	$sp = ".";
    }
    
    my @fields;
    $fields[0] = $chr;
    $fields[1] = $sp;
    $fields[2] = "cDNA_match";
    $fields[3] = $self->{tstart}+1;
    $fields[4] = $self->{tend};
    $fields[5] = 0;
    $fields[6] = $self->{tstrand};
    $fields[7] = ".";
    
    $gff .= join("\t", @fields);

    my @info;
    $info[0] = "ID=$self->{qid}";
    $info[1] = "Target=$self->{qid} " . 
	($self->{qstart}+1)." ".
	($self->{qend});
    $gff .= "\t";
    $gff .= join(";", @info);

    $gff .= "\n";
    my @ids = $self->exon_identities($query, $target);

    my @sections = @{$self->{sectionsref}};
    if(scalar @sections == 0){
	return "";
    }

    my $baseField = "Parent=$self->{qid};Target=$self->{qid} ";
    my $ti = $self->{tstart}+1;
    my $qi = $self->{qstart}+1;
    my $qstart = $qi;
    my $tstart = -1;
    my $qm = 1;
    my $idi = 0;

    if($self->{tstrand} eq "-"){
	@sections = reverse @sections;
	$qi = $self->{qend} ;
	$qstart = $qi;
	$qm = -1;
    }

    $fields[2] = "match_part";

    my @bits;
    
    for my $sectionref (@sections){
	my @section = @$sectionref;
	if($section[0] eq "M"){
	    my $size = $section[1];
	    if($tstart < 0){
		$tstart = $ti;
	    }
	    push(@bits, "M$size");
	    $qi += $size * $qm;
	    $ti += $size;
	}elsif($section[0] eq "G"){
	    if($section[1] > 0){
		push(@bits, "I$section[1]");
	    }else{
		push(@bits, "D$section[2]");
	    }
	    $qi += $section[1] * $qm;
	    $ti += $section[2];
	}elsif($tstart >= 0 && 
	       ($section[0] eq "5" ||
	        $section[0] eq "3")){
	    my $qrange;
	    if($self->{tstrand} eq "+"){
		$qrange = "$qstart " . ($qi-1) . ";Gap=";
		$qrange .= join(" ", @bits);
	    }else{
		$qrange = ($qi+1) . " $qstart;Gap=";
		@bits = reverse @bits;
		$qrange .= join(" ", @bits);
	    }
	    $fields[3] = $tstart;
	    $fields[4] = $ti-1;
	    $fields[5] = $ids[$idi];
	    $idi++;
	    $fields[8] = "$baseField$qrange";
	    $gff .= join("\t", @fields);
	    $gff .= "\n";
	    @bits = ();
	    $tstart = -1;
	    $ti += $section[2];
	    $qstart = $qi;
	}else{
	    $qi += $section[1];
	    $ti += $section[2];
	}
    }
    my $qrange;
    if($self->{tstrand} ne "-"){
	$qrange = "$qstart " . ($qi-1) . ";Gap=";
	$qrange .= join(" ", @bits);
    }else{
	$qrange = ($qi+1) . " $qstart;Gap=";
	@bits = reverse @bits;
	$qrange .= join(" ", @bits);
    }
    $fields[3] = $tstart;
    $fields[4] = $ti-1;
    $fields[5] = $ids[$idi];
    $idi++;
    $fields[8] = "$baseField$qrange";
    $gff .= join("\t", @fields);
    $gff .= "\n";
    
    return $gff;
}

sub getValidation
{
     my ($self, $query, $target) = @_;
     my ($introns, $validsplices, $splicemismatches, 
	 $incontiguities, $pctAligned) =
	     $self->getValidationStats($query, $target);
     my $s = "";

     my ($pctIntrons, $pctSplices);
     
     if($introns == 0){
	 $pctIntrons = 100;
	 $pctSplices = 0;
     }else{
	 $pctIntrons = $validsplices*100/$introns;
	 $pctSplices = $splicemismatches*50/$introns;
     }

     if($incontiguities>0 && $pctAligned>0){
	 $s.= "Incontiguous alignment ";
     }
     if($pctAligned < 90){
	 $s.= "less than 90% ";
     }
     if($pctIntrons<100){
	 $s.= "Splice site validations failed ";
     }
     if($pctSplices>0){
	 $s.= "found mismatch at splice boundary ";
     }

     if(length($s)==0){
	 return "Validated";
     }else{
	 return $s;
     }
     
}

sub getSeed
{
    my ($self, $query, $target, $header) = @_;

    my $seed = "";
    my @hsps = $self->highscoringpairs();

    $seed .= $header;
    $seed .= "\n";
    $seed .= "genomic_boundary_start=";
    $seed .= "1";
    $seed .= " genomic_boundary_end=";
    $seed .= length($target);
    $seed .= " strand=$self->{tstrand}\n";
    $seed .= "count=".scalar @hsps."\n";
    
    
    my $strand = $self->{tstrand};
    if($strand eq "-"){
	@hsps = reverse @hsps;
    }
    
    for my $hsp (@hsps){
	my @arr = @$hsp;
	$arr[2]++;
	$arr[3]++;
	
	if($strand eq "+"){
	    $arr[0]++;
	    $arr[1]++;
	    $seed .= "($arr[2], $arr[0]) ($arr[3], $arr[1])\n";
	}else{
	    $arr[0] = length($query) - $arr[0];
	    $arr[1] = length($query) - $arr[1];
	    $seed .= "($arr[3], $arr[1]) ($arr[2], $arr[0])\n";
	}
    }

    return $seed;

}

sub getValidationStats
{
    my ($self, $query, $target) = @_;
        
    my ($qstrand, $tstrand) = 
	($self->{qstrand}, $self->{tstrand});
    
    #if($qstrand eq "-"){
#	$query = reverseCompliment($query);
#    }
#    if($tstrand eq "-"){
#	$target = reverseCompliment($target);
#    }

    # valid splice sites
    my @splices = $self->getSpliceSites($query, $target);
    my $introns = scalar @splices;
    my $validsplices=0;
    for my $splice (@splices){
	if($splice eq "gt-ag" ||
	   $splice eq "gc-ag" ||
	   $splice eq "at-ac"){
	    $validsplices++;
	}
    }

    # incontiguities
    my $ascope = $self->{qend}-$self->{qstart};
    my $qlength = length($query);
    my $incontiguities = $qlength-$ascope;
    
    # percent aligned and splice site mismatches
    my $splicemismatches = 0;
    my $aligned = 0;
    my ($tseq, $aseq, $qseq) = $self->alignment($query, $target);
    
    my $state = 0;
    my $count = 0;
    for(my $i=0;$i<length($tseq);$i++){
	my $tc = substr($tseq, $i, 1);
	my $qc = substr($qseq, $i, 1);
	
	$_ = $tc;
	
	if($state==0){
	    if($tc eq $qc){
		$count++;
		$aligned++;
	    }elsif($tc eq "-"){
		$state = 1;
		$count = 1;
	    }elsif($qc eq "."){
		$state = 2;
		if($count==0){
		    $splicemismatches++;
		}
	    }elsif($tc eq "-" || 
		   $qc eq "-"){
		$count = 0;
	    }else{
		$aligned++;
		$count = 0;
		
	    }
	}elsif($state==1){
	    if($tc eq "-"){
		$count++;
	    }elsif($qc eq "."){
		$splicemismatches++;
		$state = 2;
	    }else{
		$aligned+=$count+1;
		$count = 1;
		$state = 0;
	    }
	}elsif($state==2){
	    if($qc ne "."){
		if($tc ne $qc){
		    $splicemismatches++;
		    $count = 0;
		}else{
		    $count = 1;
		}
		$state = 0;
	    }
	}
    }

    my $pctAligned = $aligned*100/$qlength;
    
    return ($introns, $validsplices, $splicemismatches, 
	    $incontiguities, $pctAligned);
}

sub getSpliceSites
{
    my ($self, $query, $target) = @_;
    my ($tseq, $aseq, $qseq) = $self->alignment($query, $target);
    my @splices = split(/[ACTG\.\-]+/, $tseq);
    
    my @pairs;
    my $last = "";
    for my $splice (@splices){
	if($last eq ""){
	    $last = $splice;
	}else{
	    push(@pairs, "$last-$splice");
	    $last = "";
	}
    }
    
    return @pairs;
}

sub reverseCompliment
{
    my ($body) = @_;
    $body =~ tr/ACTGactg/TGACTGAC/;
    $body = reverse $body;
    return $body;
}

sub alignment
{
    my ($self, $query, $target) = @_;

    my ($qstrand, $tstrand) = 
	($self->{qstrand}, $self->{tstrand});

    if($qstrand eq "-"){
	$query = reverseCompliment($query);
    }
    if($tstrand eq "-"){
	$target = reverseCompliment($target);
    }
       
    my $tseq = "";
    my $aseq = "";
    my $qseq = "";
    
    my ($query_index, $target_index);
    my ($mq, $mt);

    if($qstrand eq "+"){
	$query_index = $self->{qstart};
	$mq = 1;
    }else{
	$query_index = $self->{qend};
	$mq = -1;
    }
    if($tstrand eq "+"){
	$target_index = $self->{tstart};
	$mt = 1;
    }else{
	$target_index = $self->{tend};
	$mt = -1;
    }
    
    for my $trip (@{$self->{sectionsref}}){
	my @triple = @$trip;
	my $label = $triple[0];
	

	if($label eq "5" || $label eq "3"){
	    my $splice = ssubstr($target, $target_index, 2, $tstrand);
	    $splice =~ tr/ACTG/actg/;
	    $tseq .= $splice;
	    $qseq .= "..";
	    $aseq .= "  ";
	    
	}elsif($label eq "I"){
	    my $ilen = 20;

	    $tseq .= "." x $ilen;
	    $qseq .= "." x $ilen;
	    my $ns = " $triple[2] ";
	    if(length($ns) % 2 == 1){
		$ns .= ">";
	    }
	    my $start = ($ilen - length($ns))/2;
	    	    
	    $aseq .= ">" x $start;
	    $aseq .= $ns;
	    $aseq .= ">" x $start;
	    
	}elsif($label eq "M"){
	    
	    for(my $i=0;$i<$triple[1];$i++){
		my $qbase = ssubstr($query, $query_index+$i*$mq, 1, $qstrand);
		my $tbase = ssubstr($target, $target_index+$i*$mt, 1, $tstrand);
		
		$qseq .= $qbase;
		$tseq .= $tbase;
		
		if($qbase eq $tbase){
		    $aseq .= "|";
		}else{
		    $aseq .= " ";
		}
		
	    }
	}elsif($label eq "G"){
	    if($triple[2] == 0){
		$qseq .= ssubstr($query, $query_index, $triple[1], $qstrand);
		$tseq .= "-" x $triple[1];
		$aseq .= " " x $triple[1];
	    }else{
		$tseq .= ssubstr($target, $target_index, $triple[2], $tstrand);
		$qseq .= "-" x $triple[2];
		$aseq .= " " x $triple[2];
	    }	    
	}
	
	$query_index += $triple[1]*$mq;
	$target_index += $triple[2]*$mt;	
    }

    return ($tseq, $aseq, $qseq);
}

sub print_alignment
{
    my ($self, $query, $target) = @_;
    my ($tseq, $aseq, $qseq) = $self->alignment($query, $target);
    
    my $s = "";
    my $limit = 60;
    for(my $i=0;$i<length($qseq);$i+=$limit){
	$s .= substr($tseq, $i, $limit);
	$s .= "\n";
	$s .= substr($aseq, $i, $limit);
	$s .= "\n";
	$s .= substr($qseq, $i, $limit);
	$s .= "\n\n";
    }
    return $s;
}

sub ssubstr
{
    my ($seq, $start, $length, $strand) = @_;
    
    if(!$strand || $strand eq "+"){
	return substr($seq, $start, $length);
    }else{
	my $s = substr($seq, length($seq) - $start, $length);
	#$s =~ tr/ACTG/TGAC/;
	#$s = reverse $s;
	return $s;
    }
}


1;
