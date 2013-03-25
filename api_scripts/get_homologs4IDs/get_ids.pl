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

my $ids_file = "ids.txt";
my $homologs_file = "homologs.txt";

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

#my $registry = 'Bio::EnsEMBL::Registry';
#Bio::EnsEMBL::Registry->load_all("../production_treefam_reg_conf.pl");
#my $homology_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Homology' );
#if(!$homology_adaptor){
#	die "Could not load homology_adaptor\n";
#}
#my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
#if(!$member_adaptor){
#	die "Could not load member_adaptor\n";
#}
#my $genetree_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTree' );
#if(!$genetree_adaptor){
#	die "Could not load genetree_adaptor\n";
#}


# Output directory
my @all_ids = `cat $ids_file`;
	die "There is no ids in 'ids.txt'\n" if(!scalar(@all_ids));

foreach my $id(@all_ids){
	chomp($id);
	#print "Searching for $id\n";
	my $grepped_line = `grep $id $homologs_file`;
	#print $grepped_line;
	#exit;
	
	chomp($grepped_line);
	
	print "$grepped_line\n";


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
