#!/usr/bin/perl
#
# Description: parses Groenen table with accumulated recombination rate.  It is used to calculate recombination rates over chromosomes, with the output of delta recombination rate over centered intervals.
# 	       Specific format for the recombination table: 
# 	       		Descriptional rows that are labeled with "#" will be skipped (i.e. not parsed).
# 	       		Necessary columns for the table:
#					1. chromosome name;
# 	       			2. accumulated genetic distance (cM);
# 	       			3. physical position of marker. #??: Note that since it is a bed file, the position is represented by two consecutive integers.
# 	       		**IMPORTANT: table is presorted genetic position.  Any rows with discordant physical distances are skipped/deleted.
#
# Usage: $0 <tab-delimited.txt with comment-line styled txt file>,<column-1>,<column-2>,<column-3>
#		(column numbers are 1-based)
#
# YJ 12/10/2014, modified from Groenen-tbl_parser_v2.pl
# YJ 3/27/2015, modified to add an option to print in BedGraph format

use strict;

my ( $n_seq, $o_seq, $n_cm, $o_cm, $n_pos, $o_pos, );
my ( $delta_cm, $pos );
my ( $counter, $discarded, $used, $total, $n_line, $marked, $n_chr );


my @par = split ",", $ARGV[0];

if (  @par != 4  and  @par != 5 ){
	die "Usage: $0 <tab-delimited.txt>,<column-1>,<column-2>,<column-3>,<(optional)m>\n",
	    "column-1. chromosome name; \n",
 	    "column-2. accumulated genetic distance (cM); \n",
 	    "column-3. physical position of marker. \n",
	    "m: print in marker format (print uniform values instead of recomb. rates\n";
}

open (INFILE, "<", $par[0] ) or die ("Cannot open input file:#!\n");

my $markers_used = "used_markers.list";
my $markers_discarded = "discarded_markers.list";

open (USED, ">", $markers_used ) or die ("Cannot create output file: #!\n");
open (DISCARDED, ">", $markers_discarded ) or die ("Cannot create output file: #!\n");


while ( my $line = <INFILE> ){
	chomp $line;
	$n_line ++;

	if ( $line =~ /^#/ ){
		print DISCARDED "$line\n"; 
		$marked ++;
		next; 
	}
	
	$counter ++;
	my @field = split "\t", $line;
#	print STDOUT "$counter\t$field[0]\n";
	if ( $counter == 1 ){
		$o_seq = $field[$par[1]-1];
		$o_cm = $field[$par[2]-1];
		$o_pos = $field[$par[3]-1];
		next;
	}

### if position does not keep increasing, skip the row without saving any information.
	if ($field[$par[1]-1] eq $o_seq && $field[$par[3]-1] < $o_pos ){
#		print STDERR "Warning: physical positions are not properly sorted!  This row is skipped!\n";
		print DISCARDED "$line\n";
		$discarded ++;
		$total ++;
		next;
	}	
	$n_seq = $field[$par[1]-1]; 
	$n_cm = $field[$par[2]-1];
	$n_pos = $field[$par[3]-1];

	if ( $n_seq eq $o_seq ){
		if ( $n_cm < $o_cm ){
			die ( "Program aborted: genetic distances are not properly sorted at line $n_line.\n$line\n" ); # change to "die"?
		}
		if ( $par[4] eq "m" ){
			my $temp_pos = $o_pos + 1;
			print STDOUT "$n_seq\t$o_pos\t$temp_pos\t100\n"; # 100 is just an arbituary number.
		}else{
			$delta_cm = $n_cm - $o_cm;
			$pos = ( $n_pos + $o_pos )/2;
			$used ++;
			$total ++;
			print STDOUT "$n_seq\t$delta_cm\t$pos\t$o_pos\t$n_pos\n";
			print USED "$line\n";
		}
	}else{
		$counter = 1;
		$o_seq = $field[$par[1]-1];
		$o_cm = $field[$par[2]-1];
		$o_pos = $field[$par[3]-1];
#t		print STDERR "to check: $o_seq; line number:$n_line\n";
		$n_chr ++;
		next;
	}

	$o_seq = $n_seq;
	$o_cm = $n_cm;
	$o_pos = $n_pos;	
}

close INFILE;

print STDERR "number of lines: $n_line\n";
print STDERR "\tskipped because of marked: $marked\n";
print STDERR "\tanalyzed: $total\n";
print STDERR "\t\tused: $used\n";
print STDERR "\t\tdiscarded: $discarded.\n\n";
print STDERR "number of chromosomes: $n_chr\n";
# Modification 11/13/14 yj
#
# physical positions that are not sorted alongside with recombination rates will be deleted.  Specifically, the second rows of wherever such cases occur will be deleted.
#
# 4/23/15 yj
# 
# able to print two additinal files: markers incorporated (used) and markers discarded.  Good for sorting it out...
