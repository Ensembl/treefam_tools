use Bio::EnsEMBL::Registry;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Getopt::Long;
use Data::Dumper;
use Scalar::Util;
use strict;
use JSON;
use warnings;

use TreeFam::HomologyHelper;
use TreeFam::SearchHelper;
use TreeFam::OtherFunctions;


my %opts;
my $all = 1;
my $hmm = 0;
my $aln = 0;
my $seq = 0;
my $tree = 0;
my $help;

my $families_file = "families.txt";

my $result = GetOptions ("all" => \$all,    # numeric
                               "hmm"   => \$hmm,      
                               "aln"   => \$aln,      
                               "tree"   => \$tree,      
                               "help"   => \$help,      
                               "seq"   => \$seq) ;     
if($help){usage("")}
# decide if we should print all information or not
$all = 0 if ($hmm || $aln || $seq || $tree);
#die "all: $all, hmm: $hmm, aln:$aln, seq:$seq, tree: $tree\n";

my $registry = 'Bio::EnsEMBL::Registry';
Bio::EnsEMBL::Registry->load_all("../registry/production_treefam_reg_conf.pl");
my $homology_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Homology' );
if(!$homology_adaptor){
	die "Could not load homology_adaptor\n";
}
my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
if(!$member_adaptor){
	die "Could not load member_adaptor\n";
}
my $genetree_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTree' );
if(!$genetree_adaptor){
	die "Could not load genetree_adaptor\n";
}


# Output directory
my $output_directory = "treefam_data";
mkdir($output_directory) if !-e $output_directory;

# Read families file
	die "There are no families in  '$families_file' or the file does not exist\n" if(!-e $families_file || !-s $families_file);
my @all_families = `cat $families_file`;
	die "There is no families in 'families.txt'\n" if(!scalar(@all_families));

foreach my $family_id(@all_families){
	chomp($family_id);
	my $family_aln_file = "$output_directory/$family_id.aln";
	my $family_cds_aln_file = "$output_directory/$family_id.cds.aln";
    my $family_hmm_file= "$output_directory/$family_id.hmm";
    my $family_sequences_file= "$output_directory/$family_id.fa";
    my $family_cds_sequences_file= "$output_directory/$family_id.cds.fa";
    my $family_tree_file= "$output_directory/$family_id.nhx";
      
	print "Looking at $family_id\n";
	print "files are $family_aln_file,$family_hmm_file,$family_sequences_file,$family_cds_sequences_file, $family_tree_file\n";
	my $family_tree_adaptor = $genetree_adaptor->fetch_by_stable_id($family_id);
		die "Could not get tree adaptor for $family_id\n" if(!$family_tree_adaptor);
    my $family_root_of_tree = $family_tree_adaptor->root;
	
# Tree
	if($all || $tree){
		my $family_tree_newick = $family_root_of_tree->nhx_format;
			die "Could not get newick tree for $family_id\n" if(!$family_tree_newick);
   		TreeFam::OtherFunctions::write_to_file({file_name => $family_tree_file ,text => $family_tree_newick});
	}
# Sequences
	if($all || $seq){
		my $family_tree_leaves = ( $family_root_of_tree )->get_all_leaves();
			die "Could not get leaves from tree for $family_id\n" if(!scalar(@{$family_tree_leaves}));
   		my ($family_sequences,$family_cds_sequences);
		foreach my $leaf(@{$family_tree_leaves}){
			$family_sequences .= ">".$leaf->stable_id."\n".$leaf->sequence."\n";
			my $member = $member_adaptor->fetch_by_source_stable_id( undef, $leaf->stable_id );
			my $translated_member = $member->get_canonical_Member();
			$family_cds_sequences .= ">".$leaf->stable_id."\n".$translated_member->sequence_cds."\n";
		}
   		TreeFam::OtherFunctions::write_to_file({file_name => $family_sequences_file ,text => $family_sequences});
   		TreeFam::OtherFunctions::write_to_file({file_name => $family_cds_sequences_file ,text => $family_cds_sequences});
	}	
# Protein Alignment
	if($all || $aln){
		my $family_alignment = $family_tree_adaptor->get_SimpleAlign();
			die "Could not get alignment for $family_id\n" if(!$family_alignment);
		my $alignIO = Bio::AlignIO->newFh( -file   => ">$family_aln_file",-format => "fasta" );
		print $alignIO $family_alignment;
	}
	# CDS Sequences 
	#my $family_alignment_cds = $family_tree_adaptor->get_SimpleAlign(-CNDA => 1);
#		die "Could not get alignment for $family_id\n" if(!$family_alignment_cds);
#	my $alignIO = Bio::AlignIO->newFh( -file   => ">$family_cds_aln_file",-format => "fasta" );
#	print $alignIO $family_alignment_cds;
}

sub usage {
   my $message = $_[0];
   if (defined $message && length $message) {
      $message .= "\n"
         unless $message =~ /\n$/;
   }

   my $command = $0;
   $command =~ s#^.*/##;

   print STDERR (
      $message,
      "usage: $command -all -hmm -seq -tree -aln\n" .
      "       -all: gets hmm, alignment, sequences (protein+cds), tree for each family (default option) \n" .
      "       -hmm:  gets hmm for each family\n" .
      "       -seq: gets sequences (protein+cds) for each family\n" .
      "       -tree: gets tree for each family\n" 
   );

   die("\n")
}
