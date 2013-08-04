#
#===============================================================================
#
#         FILE: OtherFunctions.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Fabian Schreiber (), fs9@sanger.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 09/26/2012 16:15:19
#     REVISION: ---
#===============================================================================
package TreeFam::OtherFunctions;

use strict;
use warnings;
 




sub write_to_file{
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	### OPENING FILE
	open my $out, '>>', $file_name or die "Couldn't open '$file_name': $!";
	### Writing file
	print {$out} $text;
	### CLOSING FILE
	close $out or die "Couldn't close '$file_name': $!";
	return 1;
}

sub convert_hmmmer_to_json{
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $infile = $arg_ref->{infile};
	my ( $outfile, $e_seq, $e_dom, $b_seq, $b_dom, $dir, 
     $clan_overlap, $fasta, $align, $help, $as, $pfamB, 
     $json, $only_pfamB, $cpu, @hmmlib );

	my $ps = Bio::Pfam::Scan::PfamScan->new();


#my $infile = "/Users/fs9/bioinformatics/current_projects/d3treeviewer/data/hmmer_output/BRCA.out.txt";
	#my $infile = "treefam.hmm.out";
	$ps->{_hmmresultIO} = Bio::Pfam::HMM::HMMResultsIO->new;
$ps->{_all_results} = $ps->{_hmmresultIO}->parseMultiHMMER3( $infile );
my @search_results = ();


#print Dumper $results;



my $json_object;
  eval {
    require JSON;
    $json_object = new JSON;
  };
  if ( $@ ) {
    die qq(FATAL: can't load JSON module; can't write JSON-format output);
  }
    $json_object->pretty( 1 ) ;
 	return $json_object->encode( $ps->results );


}

1;

