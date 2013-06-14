#!/usr/bin/env perl

# $Id: treefam_scan.pl 6938 2012-05-30 11:01:00Z jm14 $

use strict;
use warnings;

use Bio::Pfam::Scan::PfamScan;
use Getopt::Long;

my $VERSION = "1.3"; 

#-------------------------------------------------------------------------------

# get the user options
my ( $outfile, $e_seq, $e_dom, $b_seq, $b_dom, $dir, 
     $clan_overlap, $fasta, $align, $help, $as, $pfamB, 
     $json, $only_pfamB, $cpu );

$e_seq = 10;
GetOptions( 'help'         => \$help,
            'outfile=s'    => \$outfile,
            'e_seq=f'      => \$e_seq,
            'e_dom=f'      => \$e_dom,
            'b_seq=f'      => \$b_seq,
            'b_dom=f'      => \$b_dom,
            'dir=s'        => \$dir,
            'clan_overlap' => \$clan_overlap,
            'fasta=s'      => \$fasta,
            'align'        => \$align,
            'h'            => \$help,
            'as'           => \$as,
            'pfamB'        => \$pfamB,
            'only_pfamB'   => \$only_pfamB,
            'json:s'       => \$json,
	    'cpu=i'        => \$cpu
);




help() if $help;
help() unless ( $dir and $fasta ); # required options

my $pfamA;
if ( $only_pfamB ) {
  die qq(FATAL: can't use pfamB and only_pfamB option together) if $pfamB;
  $pfamB=1;
}
else {
  $pfamA=1;
}

my @hmmlib;
push @hmmlib, 'TreeFam.hmm' if $pfamA;
#push @hmmlib, 'TreeFam-B.hmm' if $pfamB;

#-------------------------------------------------------------------------------

# check the input parameters

die qq(FATAL: must specify both "-dir" and "-fasta")
  unless ( defined $dir and defined $fasta );

die qq(FATAL: can't find directory "$dir")
  unless -d $dir;

die qq(FATAL: can't find file "$fasta")
  unless -s $fasta;

#foreach my $hmmlib ( @hmmlib ) {
#  die qq(FATAL: can't find "$hmmlib" and/or "$hmmlib" binaries and/or "$hmmlib.dat" file in "$dir")
#    unless ( -s "$dir/$hmmlib"     and 
#             -s "$dir/$hmmlib.h3f" and
#             -s "$dir/$hmmlib.h3i" and
#             -s "$dir/$hmmlib.h3m" and
#             -s "$dir/$hmmlib.h3p" and
#             -s "$dir/$hmmlib.dat" );
#}

#die qq(FATAL: can't use E-value or bit score threshold with Pfam-B searches; Pfam-B searches use a default cut_off of 0.001)
#  if ( ( $e_seq or $e_dom or $b_seq or $b_dom ) and not $pfamA ); 

die qq(FATAL: can't use E-value and bit score threshold together)
  if ( ( $e_seq and ( $b_seq or $b_dom ) ) or 
       ( $b_seq and ( $e_seq or $e_dom ) ) or 
       ( $b_dom and $e_dom ) );

die qq(FATAL: output file "$outfile" already exists)
  if ( $outfile and -s $outfile );

#if ( $as ) {
#  die qq("FATAL: "-as" option only works on Pfam-A families")
#    unless $pfamA;

  #die qq(FATAL: can't find "active_site.dat" in "$dir")
  #  unless -s "$dir/active_site.dat";
#}

#-------------------------------------------------------------------------------

# build the object
my $ps = Bio::Pfam::Scan::PfamScan->new(
  -e_seq        => $e_seq,
  -e_dom        => $e_dom,
  -b_seq        => $b_seq,
  -b_dom        => $b_dom,
  -dir          => $dir,
  -clan_overlap => $clan_overlap,
  -fasta        => $fasta,
  -align        => $align,
  -as           => $as,
  -hmmlib       => \@hmmlib,
  -version      => $VERSION,
  -cpu          => $cpu
);

# run the search
$ps->search;

# print the results
if ( defined $json ) {

  my $json_object;
  eval {
    require JSON;
    $json_object = new JSON;
  };
  if ( $@ ) {
    die qq(FATAL: can't load JSON module; can't write JSON-format output);
  }

  if ( $json eq 'pretty' ) {
    $json_object->pretty( 1 ) ;
  }
  print $json_object->encode( $ps->results );

}
else {
  $ps->write_results( $outfile, $e_seq, $e_dom, $b_seq, $b_dom );
}

exit;

#-------------------------------------------------------------------------------

sub help {
  print STDERR <<EOF;

treefam_scan.pl: search a FASTA file against a library of TreeFam HMMs

Usage: treefam_scan.pl -fasta <fasta_file> -dir <directory location of TreeFam files>

Additonal options:

  -h              : show this help
  -outfile <file> : output file, otherwise send to STDOUT
  -clan_overlap   : show overlapping hits within clan member families (applies to TreeFam-A families only)
  -align          : show the HMM-sequence alignment for each match
  -e_seq <n>      : specify hmmscan evalue sequence cutoff for TreeFam-A searches (default TreeFam defined)
  -e_dom <n>      : specify hmmscan evalue domain cutoff for TreeFam-A searches (default TreeFam defined)
  -b_seq <n>      : specify hmmscan bit score sequence cutoff for TreeFam-A searches (default TreeFam defined)
  -b_dom <n>      : specify hmmscan bit score domain cutoff for TreeFam-A searches (default TreeFam defined)
  -as             : predict active site residues for TreeFam-A matches
  -json [pretty]  : write results in JSON format. If the optional value "pretty" is given,
                    the JSON output will be formatted using the "pretty" option in the JSON
                    module
  -cpu <n>        : number of parallel CPU workers to use for multithreads (default all)

  For more help, check the perldoc:

      shell\% perldoc treefam_scan.pl

EOF
  exit;

}

#-------------------------------------------------------------------------------

=head1 NAME

treefam_scan.pl -- Search protein sequences against the TreeFam HMM library

=head1 SYNOPSIS

treefam_scan.pl [options] -fasta <fasta_file> -dir <TreeFam_data_file_dir>

=head1 OPTIONS

=over

=item B<-dir> I<TreeFam_data_file_dir>

Directory containing TreeFam data files [required]

=item B<-fasta> I<fasta_file>

Filename of input file containing sequence(s) [required]

=item B<-outfile> I<output_file>

Write output to C<output_file> [default: STDOUT]

=item B<-e_seq>

Sequence E-value cut-off [default: use TreeFam GA cutoff]

=item B<-e_dom> 

Domain E-value cut-off [default: use TreeFam GA cutoff]

=item B<-b_seq>

Sequence bits score cut-off [default: use TreeFam GA cutoff]

=item B<-b_dom>

Domain bits score cut-off [default: use TreeFam GA cutoff]

=item B<-pfamB>

Search against TreeFam-B HMMs [default: false]

=item B<-only_pfamB>

Search against TreeFam-B HMMs only [default: false]

=item B<-clan_overlap>

Allow sequences in different clans to overlap [default: false]

=item B<-align>

Show alignment snippets in results [default: false]

=item B<-as>

Search for active sites on TreeFam-A matches [default: false]

=item B<-json> [I<pretty>]

Write the results in JSON format [default: false]

=item B<-cpu>

Number of parallel CPU workers to use for multithreads [default: all]

=item B<-h>

Display help message

=back

The input must be a FASTA-format file. The C<-fasta> and C<-dir> options are 
mandatory. You cannot specify both an E-value and bits score threshold.  

=head1 OVERVIEW

C<treefam_scan.pl> is a script for searching one or more protein sequences against the
library of HMMs from TreeFam. It requires a local copy of the TreeFam data files, which 
can be obtained from the TreeFam download area:

  http://www.treefam.org/download

You must also have the HMMER3 binaries installed and their locations given by your
C<PATH> environment variable. You can download the HMMER3 package at:

  ftp://selab.janelia.org/pub/software/hmmer3/

=head1 OUTPUT

The output format is:
<seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> <predicted_active_site_residues>

Example output (with -pfamB, -as options):

  Q5NEL3.1      2    224      2    227 PB013481  TreeFam_13481      TreeFam     1   184   226    358.5  1.4e-107  NA NA
  O65039.1     38     93     38     93 PF08246   Inhibitor_I29     TreeFam     1    58    58     45.9   2.8e-12   1 No_clan
  O65039.1    126    342    126    342 PF00112   Peptidase_C1      TreeFam     1   216   216    296.0   1.1e-88   1 CL0125   predicted_active_site[150,285,307]

Most of these values are derived from the output of I<hmmscan> (see HMMER3
documentation for details). The significance value is 1 if the bit score for a
hit is greater than or equal to the curated gathering threshold for the
matching family, 0 otherwise.

=head1 AUTHORS

Fabian Schreiber (fs@ebi.ac.uk), Mateus Patricio (mateus@ebi.ac.uk), Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)

=cut

=head1 COPYRIGHT

Copyright (c) 2009: Genome Research Ltd.

Authors: Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)

This is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 
=cut

