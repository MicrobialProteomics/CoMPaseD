#!/usr/bin/perl
#
#############################################################################
# Program       : promast                                                   #
# Author        : Luis Mendoza                                              #
# Date          : 07.07.18                                                  #
# SVN Info      : $Id$
#                                                                           #
# Find all instances of input peptide sequence(s) by index lookup           #
# Copyright (C) 2018 Luis Mendoza                                           #
#                                                                           #
# This library is free software; you can redistribute it and/or             #
# modify it under the terms of the GNU Lesser General Public                #
# License as published by the Free Software Foundation; either              #
# version 2.1 of the License, or (at your option) any later version.        #
#                                                                           #
# This library is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         #
# General Public License for more details.                                  #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public          #
# License along with this library; if not, write to the Free Software       #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
#                                                                           #
# Institute for Systems Biology                                             #
# 401 Terry Avenue North                                                    #
# Seattle, WA  98109  USA                                                   #
# lmendoza@isb                                                              #
#                                                                           #
#############################################################################
#
use strict;
use warnings;
no warnings ('uninitialized', 'substr');
use threads;
use threads::shared;

use Getopt::Std;
use File::Basename;
use Benchmark qw(:all);
my @btimes;
push @btimes, Benchmark->new; #0 :: exec start time

my $VersionInfo = "stand-alone 1.0";
eval "require tpplib_perl";
unless ($@) {  # ToDo: fix this
#    use tpplib_perl; # exported TPP lib function points
#    $VersionInfo = tpplib_perl::getTPPVersionInfo();
}

$|++; # autoflush

my $USAGE=<<"EOU";
-------------------------------------------------------------------------------------
Program:    $0 ($VersionInfo)
Purpose:    Maps peptide sequences to all proteins using indexed segments.
            Specify either a single <peptide> (use X for wildcard), a list of peptides
             (one per line) contained in <file>, or a <pepxml_file>.

Usage:      $0 [options]  <fasta_file> <peptide>|<file>|<pepxml_file>
Requires:   File <fasta_file>.pep.idx  (unless using -i option)
Options:
            -U        omit UNMAPPED sequences from report
            -P        omit position of mapped sequence; only report protein
            -I        collapse isoforms to a single reported entry [pepx output option]
            -c        provide sequence context of mapped peptide (unmod protein sequence, and flanking AAs)
            -u        convert input sequences to uppercase
            -o <fmt>  output format, one of: tsv, text, pepx, html, json  [default:tsv]
            -f <num>  fuzzy sequence mapping, with <num> unknown aminoacids  [max:3]
                       note: when setting this to 3, only consecutive AAs are considered
            -m <tol>  only consider "fuzzy" peptides of mass within +/- <tol> of original peptide
            -t <num>  number of threads to use for (potentially) faster processing  [default:1]
            -i        input file is index file, not source fasta
			-n <str>  output filename
For Developers:
            -z        print performance metrics
            -D        print debug information
-------------------------------------------------------------------------------------
EOU
#            -V        do not consider PEFF variants [default:use them] -- future option!

my %options;
# getopts('VUPIcuo:f:m:t:izD', \%options);
getopts('VUPIcuo:f:m:t:in:zD', \%options);

my $DEBUG      = $options{'D'} || '';
my $index      = shift || die "Error: Please provide an input fasta file.\n".$USAGE;
$index        .= $options{'i'} ? '' : ".pep.idx";
if (! -e $index) { die  "Error: could not find index file $index!\n".$USAGE; }

my $pepsrc     = shift || die "I need a peptide sequence or pepXML file!\n".$USAGE;

my %keys;
my %prots;
my %protseqs;
my %segidx;
my %index;
my %peps : shared;
my %modpeps;
my %aamass;
my %aamods;
my $pepkeysize = 5; # default in case not in index file
my $ItoL       = 0; # default in case not in index file
my $usepeff    = $options{'V'} ? 0:1;  #ToDo: implement this
$options{'o'}||= 'tsv';

my $aa_list    = 'ACDEFGHIKLMNPQRSTVWY';

####
my $filename      = $options{'n'};
####
&readIdxHdr();   # populate %prots and capture index settings
$aa_list =~ s/I// if $ItoL;
push @btimes, Benchmark->new; #1 :: setup and read index

&collectPeps($pepsrc);
push @btimes, Benchmark->new; #2 :: read peptide list / parse pepXML

&getPepKeys();
push @btimes, Benchmark->new; #3 :: split peptides into segments

&extractKeys();
push @btimes, Benchmark->new; #4 :: read segments index entries from file

$options{'t'} ? &dispatchPeps($options{'t'}) : &processPeps();
push @btimes, Benchmark->new; #5 :: actual computation

&printReport();
push @btimes, Benchmark->new; #6 :: sorting / output

&perfReport() if $options{'z'};
exit(0);


####################################################################################
sub perfReport {
####################################################################################
  my $told = shift(@btimes);
  my $step = 0;
  for (@btimes) {
    my $td = timediff($_, $told);
    print STDERR ++$step,": ",timestr($td),"\n";
    $told = $_;
  }
}


####################################################################################
sub collectPeps {
####################################################################################
  my $p = shift;

  my $parse_xml = 0;  # 0 seems faster

  if (! -e $p) {  # assume it is a peptide sequence
    print STDERR "Processing single peptide sequence $p";
    $p = uc $p if $options{'u'};
    $p =~ s/I/L/g if $ItoL;

    if ($p =~ /X/) {
      print STDERR "...with wildcards";
      print STDERR " (ignoring fuzzy matching)" if $options{'f'};
      &expand_wildcards('',$p);
    }
    elsif ($options{'f'}) {
      print STDERR "...with fuzzy mapping of ".$options{'f'}." sites";
      print STDERR "(consecutive)" if $options{'f'} == 3;

      if ($options{'f'} !~ /[123]/ || $options{'f'} > 3) {
	print STDERR "Error: There is a maximum of 3 sites when specifying fuzzy matching.  Please re-run with a different number.\n";
	exit 2;
      }

      if ($options{'m'}) {
	&initMasses();
	$options{'pepmass'} = &calcMass($p);
	print STDERR "...and mass range within: ".$options{'pepmass'}."+/-".$options{'m'};
      }

      my @q;
      if ($options{'f'} == 3) {
	for my $i (0..length($p)-3) {
	  my $q = $p;
	  substr($q, $i,3) = 'XXX';
	  push @q, $q;
	}
      }
      else {
	for my $i (0..length($p)-1) {
	  my $q = $p;
	  substr($q, $i,1) = 'X';

	  if ($options{'f'} == 2) {
	    for my $j ($i+1..length($p)-1) {
	      my $r = $q;
	      substr($r, $j,1) = 'X';
	      push @q, $r;
	    }
	  }
	  else {
	    push @q, $q;
	  }
	}
      }


      for (@q) {
	&expand_wildcards('',$_);
      }
    }
    else {
      $peps{$p} = 1;
    }
    print STDERR "\n";

  }
  elsif ($p !~ /pep.xml$/) {  # assume it is a peptide list
    open(PEPS, $p) || die "\n\nError: cannot open $p for reading peptides. Exiting";
    while(<PEPS>) {
      s/[\r\n]//g; # cheap chop
      $_ = uc if $options{'u'};
      s/I/L/g if $ItoL;
      $peps{$_} = 1;
    }
    close (PEPS);
    print STDERR "Read ". scalar(keys %peps) ." peptide sequences from $p\n";
  }
  elsif ($parse_xml) {
    use XML::Parser;
    my $parser = XML::Parser->new(
      Handlers => { Start => \&proc_pepxml_element },
      ErrorContext => 2 );
    eval { $parser->parsefile( $p ); };
    die "ERROR_PARSING_XML:$@" if($@);

    print STDERR "Read ". scalar(keys %peps) ." peptide sequences from $p\n";
  }
  else {
    open(PXML, $p) || die "\n\nError: cannot open $p for reading peptides. Exiting";
    while(<PXML>) {
      if (index($_, '<search_hit') != -1) {
	/peptide="(\w+)"/;
	my $peptide =  $1;
	/hit_rank="(\d+)"/;
	my $rank =  $1;

	$peptide =~ s/I/L/g if $ItoL;
	$peps{$peptide} = 1 if ($rank == 1);
      }
    }
    close (PXML);
    print STDERR "Read ". scalar(keys %peps) ." peptide sequences from $p\n";
  }


  if (scalar(keys %peps) < 100 && $options{'t'}) {
    print STDERR "(too few spectra; setting numthreads to 1)\n";
    $options{'t'} = 1;
  }
}

####################################################################################
sub proc_pepxml_element {
####################################################################################
  my ($handler, $element, %atts) = @_;

  if ($element eq 'search_hit') {
    my $rank    = $atts{'hit_rank'} || 1;
    my $peptide = $atts{'peptide'};

    $peptide =~ s/I/L/g if $ItoL;
    $peps{$peptide} = 1 if ($rank == 1);
  }

}


####################################################################################
sub expand_wildcards {    # recurse me...recurse me, my friend...
####################################################################################
  my $str = shift;
  my $pep = shift;

  my $aas = substr $pep, length($str), 1;
  if ($aas eq 'X') { $aas = $aa_list };

  for my $aa (split("", $aas)) {
    if (1+length($str) == length($pep)) {

      if (!$options{'m'}) {
	$peps{"$str$aa"} = 1;
      }
      elsif (&isMassWithinRange("$str$aa")) {
	$peps{"$str$aa"} = 1;
      }

    }
    else {
      &expand_wildcards("$str$aa", $pep);
    }
  }
}



####################################################################################
sub isMassWithinRange {
####################################################################################
  my $pepseq = shift;

#  if ($pepseq eq 'CEFGM') {
#    print STDERR "\n================= $pepseq\n";
#  }

  my $base_mass = &calcMass($pepseq);

  return 0 if $base_mass > $options{'pepmass'}+$options{'m'};  #assumes all possible mods are > 0
  return 1 if $base_mass > $options{'pepmass'}-$options{'m'};

  return 1 if ($peps{$pepseq});

  ## block below... return 0;

  my $delta = $base_mass - $options{'pepmass'};
  my @pot_mods;
  for my $aa (split //, "nc$pepseq") {
    push @pot_mods, $aamods{$aa} if $aamods{$aa}; # add!  && $mods{$aa} < ($delta+$options{'m'};
  }
  for (0..$#pot_mods) {
    $pot_mods[$_] .= ",0";
  }

#  if ($pepseq eq 'CEFGM') {
#    print STDERR "================= delta :: $delta\n";
#    print STDERR "================= mods  :: @pot_mods\n";
#    my $retc = &assemble_mods(0, 0, 'CEFGM', $delta, \@pot_mods);
#    print STDERR "================= valid :: $retc\n";
#    return $retc;
#  }


  my $found = &assemble_mods(0, 0, '', $delta, \@pot_mods);
  if ($found) {
    $modpeps{$pepseq} = "$pepseq$found";
    return 1;
  }
  else {
    return 0;
  }
}


####################################################################################
sub assemble_mods {    # recurse me...recurse me, my friend...
####################################################################################
  my $aapos = shift;
  my $extramass = shift;
  my $modstr = shift;
  my $deltamass = shift;
  my $modaas = shift;

  for my $maa (split(",", @$modaas[$aapos])) {
    my $newextra = $extramass + $maa;

    next if $newextra + $deltamass > $options{'m'};

    my $newmodstr = $maa ? "$modstr(+$maa)" : $modstr;

    if ($newextra + $deltamass > -($options{'m'})) {
      print STDERR "pepmod:$newmodstr\n" if $DEBUG;
      return $newmodstr;
    }
    elsif ($aapos+1 == scalar @$modaas) {
      return 0;
    }

    my $fm = &assemble_mods($aapos+1, $newextra, $newmodstr, $deltamass, $modaas);
    return $fm if $fm;

    # next
  }

  return 0;
}


####################################################################################
sub extractProtList {
####################################################################################
  my %findprots;
  for my $pep (keys %peps) {
    if ($peps{$pep}) {

      for (split /:/, $peps{$pep}) {
	my ($prot, $pos) = split /,/;
	my $pseq = (split /\s/,$prots{$prot})[0];

	print "Found prot entry: $prot :: [$pseq] $prots{$prot}\n" if $DEBUG;

	$findprots{$pseq} = 1;
      }
    }
  }

  my $fasta = dirname($index).'/'.$options{'fasta'};

  unless (-e $fasta) {
    print STDERR "Source sequence file ($fasta) not found. Cannot provide peptide mapping context.\n";
    return;
  }
  open(FASTA, $fasta) || do {
    print STDERR "Could not open file for input ($fasta). Cannot provide peptide mapping context.\n";
    return;
  };

  print STDERR "Extracting original protein sequences from fasta file...\n";
  my $readseq = 0;
  my $protacc = '';
  my $protseq = '-';
  while(<FASTA>) {
    if (/^>(\S*)/) {
      if ($readseq) {
	$protseqs{$protacc} = $protseq.'-';
	delete $findprots{$protacc};
	$readseq = 0;
      }

      $protacc = $1;
      $protseq = '-';

      for (keys %findprots) {
	if ($_ eq $protacc) {
	  $readseq = 1;
	  last;
	}
      }
    }
    elsif ($readseq) {
      s/[\r\n]//g; # cheap chop
      s/I/L/g if $ItoL;
      $protseq .= $_;
    }


    last if keys %findprots < 1;
  }
  # capture last sequence in file, if still looking
  $protseqs{$protacc} = "${protseq}-" if $readseq;


  close FASTA;
}


####################################################################################
sub printReport {
####################################################################################

# open(Outfile, ">", "outfile.tsv");
open(Outfile, ">", $filename);

  &extractProtList() if $options{'c'};

  my $report_type = $options{'o'};   # tsv text pepx html json

  if ($report_type eq 'tsv') {
    print Outfile "peptide\tprotein";
    print Outfile "\tlocation" unless $options{'P'};
    print Outfile "\tprevAA\tin_fasta\tnextAA" if $options{'c'};
    print Outfile "\n";
  }
  elsif ($report_type eq 'html') {
    print "<table>\n<tr><th>peptide</th><th>protein</th>";
    print "<th>location</th>" unless $options{'P'};
    print "<th>prevAA</th><th>in_fasta</th><th>nextAA</th>" if $options{'c'};
    print "</tr>\n";
  }
  elsif ($report_type eq 'json') {
    print "{\n\t\"mappings\" : [";
  }

  my $found = 0;
  my $jsoncomma = '';
  for my $pep (sort keys %peps) {
    if ($peps{$pep}) {

      if ($report_type eq 'text') {
	print "$pep found at following proteins";
	print ",locations" unless $options{'P'};
	print ":";
      }

      my $comma = '';
      my $prevprot = '';

      for (sort {$prots{(split /,/,$a)[0]} cmp $prots{(split /,/,$b)[0]} or (split /,/,$a)[1]<=>(split /,/,$b)[1]} split /:/, $peps{$pep}) {
	my ($prot, $pos) = split /,/, $_;

	## for (sort {(split /\./,$a)[1]<=>(split /\./,$b)[1] or (split /\./,$a)[0]<=>(split /\./,$b)[0]} split /:/, $peps{$pep}) {   # use for reverse decimal indices
	## my ($pos, $prot) = split /\./, $_;   # use for reverse decimal indices

	# $pos++;  # since processPeps_good records zero-based position

	$found++;

	if ($report_type eq 'tsv') {
	  my $mod = $modpeps{$pep} ? $modpeps{$pep} : $pep;

	  if ($options{'P'}) {
	    if ($prots{$prot} ne $prevprot) {
	      $prevprot = $prots{$prot};
	      print Outfile "$mod\t$prots{$prot}\n";
	    }
	  }
	  else {
	    print Outfile "$mod\t$prots{$prot}\t$pos";

	    if ($options{'c'}) {
	      my $seq = $protseqs{(split /\s/,$prots{$prot})[0]};
	      my $prev = substr $seq, ($pos-1), 1;
	      my $orig = substr $seq, $pos, length $pep;
	      my $foll = substr $seq, ($pos+length $pep), 1;
	      print Outfile "\t$prev\t$orig\t$foll";
	    }
	    print Outfile "\n";
	  }
	}

	elsif ($report_type eq 'html') {
	  my $mod = $modpeps{$pep} ? $modpeps{$pep} : $pep;

	  if ($options{'P'}) {
	    if ($prots{$prot} ne $prevprot) {
	      $prevprot = $prots{$prot};
	      print "<tr><td>$mod</td><td>$prots{$prot}</td></tr>\n";
	    }
	  }
	  else {
	    print "<tr><td>$mod</td><td>$prots{$prot}</td><td>$pos</td>";

	    if ($options{'c'}) {
	      my $seq = $protseqs{(split /\s/,$prots{$prot})[0]};
	      my $prev = substr $seq, ($pos-1), 1;
	      my $orig = substr $seq, $pos, length $pep;
	      my $foll = substr $seq, ($pos+length $pep), 1;
	      print "<td>$prev</td><td>$orig</td><td>$foll</td>";
	    }
	    print "</tr>\n";
	  }
	}

	elsif ($report_type eq 'json') {
	  my $mod = $modpeps{$pep} ? $modpeps{$pep} : $pep;

	  if ($options{'P'}) {
	    if ($prots{$prot} ne $prevprot) {
	      $prevprot = $prots{$prot};
	      print "$jsoncomma\n\t{\n";
	      print "\t\t\"protein\" : \"$prots{$prot}\",\n";
	      print "\t\t\"peptide\" : \"$mod\"\n";
	      print "\t},\n";
	    }
	  }
	  else {
	    print "$jsoncomma\n\t{\n";
	    print "\t\t\"protein\" : \"$prots{$prot}\",\n";
	    print "\t\t\"peptide\" : \"$mod\",\n";
	    print "\t\t\"location\" : $pos,\n";

	    if ($options{'c'}) {
	      my $seq = $protseqs{(split /\s/,$prots{$prot})[0]};
	      my $prev = substr $seq, ($pos-1), 1;
	      my $orig = substr $seq, $pos, length $pep;
	      my $foll = substr $seq, ($pos+length $pep), 1;

	      print "\t\t\"prevAA\" : \"$prev\",\n";
	      print "\t\t\"orig_seq\" : \"$orig\",\n";
	      print "\t\t\"follAA\" : \"$foll\"\n";
	    }
	    print "\t}";
	    $jsoncomma = ',';
	  }
	}

	elsif ($report_type eq 'text') {
	  if ($options{'P'}) {
	    if ($prots{$prot} ne $prevprot) {
	      $prevprot = $prots{$prot};
	      print " $prots{$prot}";
	    }
	  }
	  else {
	    print " $prots{$prot},$pos";
	  }
	}

	elsif ($report_type eq 'pepx') {  # matches -n option in pepx
	  $prots{$prot} =~ s/-.*// if $options{'I'};  # ignore isoforms

	  if ($prots{$prot} ne $prevprot) {
	    $prevprot = $prots{$prot};
	    $prots{$prot} =~ s/.*\|(.*)\|.*/$1/;
	    print "$comma$prots{$prot}";
	    $comma = ',';
	  }
	}

      }
      print " $pep\n" if $report_type eq 'pepx';
      print "\n" if $report_type eq 'text';
    }
    elsif (!$options{'U'}) {
      print "$pep NOT found!\n" if $report_type eq 'text';
      print "NO_MATCH $pep\n" if $report_type eq 'pepx';
      if ($report_type eq 'tsv') {
	print "$pep\tUNMAPPED\t0";
	print "\t-\t".lc($pep)."\t-" if $options{'c'};
	print "\n";
      }
      elsif ($report_type eq 'html') {
	print "<tr><td>$pep</td><td>UNMAPPED</td><td>0</td>";
	print "<td>-</td><td>".lc($pep)."</td><td>-</td>" if $options{'c'};
	print "</tr>\n";
      }
      elsif ($report_type eq 'json') {
	print "$jsoncomma\n\t{\n";
	print "\t\t\"protein\" : \"UNMAPPED\",\n";
	print "\t\t\"peptide\" : \"$pep\",\n";
	print "\t\t\"location\" : 0,\n";

	if ($options{'c'}) {
	  print "\t\t\"prevAA\" : \"-\",\n";
	  print "\t\t\"orig_seq\" : \"". lc($pep) . "\",\n";
	  print "\t\t\"follAA\" : \"-\"\n";
	}
	print "\t}";
	$jsoncomma = ',';
      }
    }
  }
  print "</table>\n" if $report_type eq 'html';
  print "\n\t],\n\t\"message\" : \"$found mappings found\",\n\t\"status\" : \"OK\"\n}\n" if $report_type eq 'json';
}


####################################################################################
sub processPep_faster {   # NOTE: final peptide position within protein is ONE-based
####################################################################################
  my $pep = shift;
  print STDERR "="x80 . "\n$pep (". int(length($pep)/$pepkeysize+0.9999) .")\n" if $DEBUG;

  my %all_locs;
  my %cur_locs;
  keys(%all_locs) = 128;
  keys(%cur_locs) = 128;

  my $prot;
  my $pos;

  # initialize these:
  my $peppos = length($pep)-$pepkeysize;
  my $offset = ($peppos<=$pepkeysize ? $peppos : $pepkeysize);

  # start at end of peptide...
  my $key = substr $pep, $peppos, $pepkeysize;
  for (split /:/, $index{$key}) {
    ($prot, $pos) = split /,/;
    $all_locs{"$prot,".($pos-$offset)} = 1; # unless $offset > $pos;  # slower...
    ##$all_locs{--$_} = 1;   # use for reverse decimal indices
  }
  print STDERR "\t$key :: ". scalar keys(%all_locs) ."\n" if $DEBUG;
#    for (keys %all_locs) { print STDERR "\t\t$_\n"; } if $DEBUG;


  # ...and now the rest
  while ($peppos > 0) {
    $peppos -= ($peppos<=$pepkeysize) ? $peppos : $pepkeysize;

#	if ($peppos>0) { $offset = ($peppos<=$pepkeysize ? $peppos : $pepkeysize); }
#	else {$offset = 0; }

    $offset = ($peppos>0) * ($peppos<=$pepkeysize ? $peppos : $pepkeysize);


    $key = substr $pep, $peppos, $pepkeysize;

    %cur_locs = ();
    for (split /:/, $index{$key}) {
      if ($all_locs{$_}) {
	($prot, $pos) = split /,/;
	$cur_locs{"$prot,".($pos-$offset)} = 1; # unless $offset > $pos;;  # slower...
	## $cur_locs{--$_} = 1;   # use for reverse decimal indices
      }
    }

    %all_locs = %cur_locs;

    print STDERR "\t$key :: ". scalar keys(%all_locs) ."\n" if $DEBUG;
  }

#    for (keys %all_locs) { print "\t\t$_\n"; }
#    $peps{$pep} = $mapped;

  $peps{$pep} = join(':', keys %all_locs);
}

####################################################################################
sub processPep_good_UNUSED {   # NOTE: final peptide position within protein is ZERO-based
####################################################################################
  my $pep = shift;
  print STDERR "="x80 . "\n$pep (". (length($pep)-$pepkeysize+1) .")\n" if $DEBUG;

  my %all_locs;
  my %cur_locs;
  keys(%all_locs) = 128;
  keys(%cur_locs) = 128;

  my $prot;
  my $pos;

  # start at end of peptide...
  my $key = substr $pep, -$pepkeysize;
  for (split /:/, $index{$key}) {
    ($prot, $pos) = split /,/;
    $all_locs{"$prot,".($pos-1)} = 1;
    ##$all_locs{--$_} = 1;   # use for reverse decimal indices
  }
  print STDERR "\t$key :: ". scalar keys(%all_locs) ."\n" if $DEBUG;

  # ...and now the rest
  for my $i (1..length($pep)-$pepkeysize) {
    $key = substr $pep, -($pepkeysize+$i), $pepkeysize;

    %cur_locs = ();
    for (split /:/, $index{$key}) {
      if ($all_locs{$_}) {
	($prot, $pos) = split /,/;
	$cur_locs{"$prot,".($pos-1)} = 1;
	## $cur_locs{--$_} = 1;   # use for reverse decimal indices
      }
    }

    %all_locs = %cur_locs;
    print STDERR "\t$key :: ". scalar keys(%all_locs) ."\n" if $DEBUG;
  }

#    for (keys %all_locs) {
#	$mapped .= "$_:";
#    }
#    $peps{$pep} = $mapped;

  $peps{$pep} = join(':', keys %all_locs);
}

####################################################################################
sub processPeps {
####################################################################################
  my $thn = shift || 0;
  print STDERR "Starting thread $thn...\n" if $thn;

  my $processed = 0;
  for my $p (keys %peps) {
    if ($peps{$p} eq '1') {
      $peps{$p} = ''; # also blocks other threads from grabbing this entry
      &processPep_faster($p);
      $processed++;
    }
  }
  print STDERR '' . ($thn?"Thread $thn p":'P') ."rocessed $processed peptides\n";
  return;
}

####################################################################################
sub dispatchPeps {
####################################################################################
  my $maxthreads = shift || 1;

  if ($maxthreads < 2) {
    &processPeps();
  }
  else {
    for my $i (1..$maxthreads) {
      threads->new(\&processPeps, $i);
    }
    # Catch when these threads end
    foreach my $onthr (threads->list()) {
      $onthr->join();
    }
  }
}


####################################################################################
sub getALLPepKeys_UNUSED {   # no longer used (left for debugging purposes)
####################################################################################
  keys(%keys) = 10 * (scalar keys %peps);

  for my $pep (keys %peps) {
    if (length($pep) < $pepkeysize) {
      print STDERR "Warning: peptide $pep is too short; ignoring...\n";
      next;
    }
    for (my $i = 0; $i <= length($pep) - $pepkeysize; $i++) {
      my $key = substr $pep, $i, $pepkeysize;
      $keys{$key} = 1;
    }
  }
  print STDERR "Retrieved ". scalar(keys %keys) ." segments from peptide list\n";
}

####################################################################################
sub getPepKeys {
####################################################################################
  keys(%keys) = 5 * (scalar keys %peps);

  for my $pep (keys %peps) {
    my $p = length($pep);

    if ($p < $pepkeysize) {
      print STDERR "Warning: peptide $pep is too short; ignoring...\n";
      delete $peps{$pep};
      next;
    }
    if ($pep =~ /[JUZBOX]/ ||
	$pep !~ /^[A-Z]+$/ ) {
      print STDERR "Warning: peptide $pep contains non-AA characters; ignoring...\n";
      delete $peps{$pep};
      next;
    }

    while ($p>0) {
      $p -= $pepkeysize;
      if ($p<0)  { $p = 0; }  # always get the start of the peptide
      $keys{substr $pep, $p, $pepkeysize} = 1;
    }
  }
  print STDERR "Retrieved ". scalar(keys %keys) ." segments from peptide list\n";
}


####################################################################################
sub readIdxHdr {
####################################################################################
  open(INDX, $index) || die "\n\nError: cannot open $index for reading header. Exiting";

  print STDERR "Reading index file header";
  my $inalias = 0;

  while(<INDX>) {
    last if /BeginIndex/;

    $inalias++ if /BeginProteins/;

    if (/KeyLength=(\d*)/) {
      $pepkeysize = $1;
      print STDERR "...setting keysize to $pepkeysize";
    }
    if (/AASub.*I->L/) {
      $ItoL = 1;
      print STDERR "...using I->L substitution";
    }
    if (/OriginalFile=(.*)/) {
      $pepkeysize = $1;
      $options{'fasta'} = $1;
    }

    if (/::/) {
      s/[\r\n]//g; # cheap chop
      my ($L, $R) = split /::/;

      if ($inalias) {
	$prots{$L} = $R;
      }
      else {
	$segidx{$L} = $R;
	$segidx{'LandmarkLen'} ||= length $L;
      }
    }

  }
  close (INDX);
  print STDERR "\n";
}


####################################################################################
sub extractKeys_noSegIdx_UNUSED {  # original
####################################################################################
  if (!%keys) {
    print STDERR "Error: no keys were generated; cannot look up peptides.\n";
    exit 1;
  }

  open(INDX, $index) || die "\n\nError: cannot open $index. Exiting";
  binmode INDX;
  my $last_match_byte_offset = 0;

  print STDERR "Reading segments from index file...\n";

  my @k = sort keys %keys;
  keys(%index) = scalar @k;  # presize hash for performance

  my $k = shift(@k);
  print STDERR "Looking for key $k..." if $DEBUG;

  while (1) {
    while(<INDX>) {
      if (s/^$k\://) {
	s/[\r\n]//g; # cheap chop

	$index{$k} = $_;
#	$index{$k} = [ split /:/, $_ ];

	print STDERR "found: $_\n" if $DEBUG;
	$last_match_byte_offset = tell INDX;

	$k = shift(@k) || last;
	print STDERR "Looking for key $k..." if $DEBUG;
      }
    }

    if (!$index{$k}) {
      $index{$k} = '';
      print STDERR "$k NOT found! ";# if $DEBUG;
      $k = shift(@k) || last;

      print STDERR "Going back to $last_match_byte_offset\n";# if $DEBUG;
      seek INDX, $last_match_byte_offset, 0;

      print STDERR "Looking for key $k..." if $DEBUG;
    }
    else {
      last;
    }
  }

  close (INDX);
}



####################################################################################
sub extractKeys {
####################################################################################
  if (!%keys) {
    print STDERR "Error: no keys were generated; cannot look up peptides.\n";
    exit 1;
  }
  my $sgsize = $segidx{'LandmarkLen'} || 2;


  open(INDX, $index) || die "\n\nError: cannot open $index. Exiting";
  binmode INDX;

  print STDERR "Reading segments from index file...\n";

  my @k = sort keys %keys;
  keys(%index) = scalar @k;  # presize hash for performance

  my $k = shift(@k);
  my $f2 = substr $k, 0, $sgsize;
  seek (INDX, $segidx{$f2}, 0) if $segidx{$f2};
  print STDERR "-->$f2 [$segidx{$f2}]\n" if $DEBUG;
  print STDERR "Looking for key $k..." if $DEBUG;

  while (1) {
    while(<INDX>) {
      if (s/^$k\://) {
	s/[\r\n]//g; # cheap chop

	$index{$k} = $_;
#	$index{$k} = [ split /:/, $_ ];

	print STDERR "found: $_\n" if $DEBUG;

	$k = shift(@k) || last;

	if ($f2 ne substr $k, 0, $sgsize) {
	  $f2 = substr $k, 0, $sgsize;
	  seek (INDX, $segidx{$f2}, 0) if $segidx{$f2};
	  print STDERR "-->$f2 [$segidx{$f2}]\n" if $DEBUG;
	}
	print STDERR "Looking for key $k..." if $DEBUG;
      }
    }

    if (!$index{$k}) {
      $index{$k} = '';
      print STDERR "-- segment $k NOT found! \n";# if $DEBUG;
      $k = shift(@k) || last;

      $f2 = substr $k, 0, $sgsize;
      seek (INDX, $segidx{$f2}, 0) if $segidx{$f2};
      print STDERR "-->$f2 [$segidx{$f2}]\n" if $DEBUG;
      print STDERR "Looking for key $k..." if $DEBUG;
    }
    else {
      last;
    }
  }

  close (INDX);
}


####################################################################################
sub calcMass {
####################################################################################
  my $seq = shift || '';

  my $mass = $aamass{'c'} + $aamass{'n'} + $aamass{'p'};
  for (split("", $seq)) {
    $mass += $aamass{$_} ? $aamass{$_} : 0;
  }

  return $mass;
}

####################################################################################
sub initMasses {
####################################################################################
  # AA mono masses (from comet.params)
  $aamass{'G'} = 57.02146;
  $aamass{'A'} = 71.03711;
  $aamass{'S'} = 87.03203;
  $aamass{'P'} = 97.05276;
  $aamass{'V'} = 99.06841;
  $aamass{'T'} = 101.04768;
  $aamass{'C'} = 103.00918;
  $aamass{'L'} = 113.08406;
  $aamass{'I'} = 113.08406;
  $aamass{'N'} = 114.04293;
  $aamass{'D'} = 115.02694;
  $aamass{'Q'} = 128.05858;
  $aamass{'K'} = 128.09496;
  $aamass{'E'} = 129.04259;
  $aamass{'M'} = 131.04048;
  $aamass{'O'} = 132.08988;
  $aamass{'H'} = 137.05891;
  $aamass{'F'} = 147.06841;
  $aamass{'U'} = 150.95363;
  $aamass{'R'} = 156.10111;
  $aamass{'Y'} = 163.06333;
  $aamass{'W'} = 186.07931;
  $aamass{'n'} = 1.007825;
  $aamass{'c'} = 17.00274;
  $aamass{'p'} = 1.00727649;

  while (<DATA>) {
    s/\s//g;
    my @mm = split /,/;
    if ($aamods{$mm[0]}) {
      $aamods{$mm[0]} .= ",$mm[1]";
    }
    else {
      $aamods{$mm[0]} = $mm[1];
    }
  }

  if ($DEBUG) {
    print "\n\n";
    for (sort keys %aamods) {
      print "Found mod on $_ = $aamods{$_}\n";
    }
  }

}


# default mods
# aminoacid (or term), mass diff, [optional] description
__DATA__
K,42.0106, Acetyl
n,42.0106, Acetyl
C,57.0215, Carbamidomethyl
C,58.0055, Carboxymethyl
Q,0.9840, Deamidation
N,0.9840, Deamidation
c,14.0157, Methylation
M,15.9949, Oxidation
S,79.9663, Phosphorylation
T,79.9663, Phosphorylation
Y,79.9663, Phosphorylation
