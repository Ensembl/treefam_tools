#
#===============================================================================
#
#         FILE: MakeSpeciesTree.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Fabian Schreiber (), fs9@sanger.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 02/13/2013 10:10:11 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use Data::Dumper; 
# make a list of clean ids, upload it to NCBI Taxonomy. Save tree in phylip format, remove " ' ", s/ /_/. Swap branches using e.g. FigTree to make the tree look nicer.
#If FigTree can't open the file, try Archaeopteryx.
#Convert newick2Json. Use the script from Pierre Lindebaum (http://www.biostars.org/p/48424/) to do that.
#Should be in the converter directory on Github (https://github.com/fabsta/TreeFam/tree/master/tools/converters)


package TreeFam::Release::SpeciesTree;
  use Moose;
    with 'TreeFam::Release::ReleaseMaster';
    has 'treestring' => (isa => 'string', is => 'ro');
    has 'filename' => (isa => 'string', is => 'ro');
  
# Get the tree from 
# SELECT * FROM method_link_species_set_tag
  
sub new {
   my $obj = shift;
   my $class = ref($obj)?ref($obj) : $obj;
   my $args = shift;

    my $self = bless {}, $class;
   $self->{-db_adaptor} = $args->{db_adaptor} if( $args && ref($args) eq 'HASH' && $args->{db_adaptor});
   $self->{-filename} = $args->{filename} if( $args && ref($args) eq 'HASH' && $args->{filename});
    #print "in constructor: $self->{-db_adaptor}\n";
   return $self;
}



  sub get_db_tree_string {
      my $self = shift;
    
    my $db_adaptor = $self->{-db_adaptor};
    my $db_sth = $db_adaptor->prepare('SELECT * FROM method_link_species_set_tag');
    
	$db_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my %hgnc_hits;
    my $returned_rows = $db_sth->rows;
    die "get_db_tree_string: returned too many rows\n" if $returned_rows > 1;
### Parse
    my $tree_string;
	while ( my @data = $db_sth->fetchrow_array() ){
        $tree_string = $data[2];
    }
    #print $tree_string."\n"; 
    if($tree_string =~ /^\(.*\w+.*\)\d*;$/){
        $self->{-tree_string} = $tree_string;
        return $tree_string;    
    }
    else{
        return undef;
    }
  }

  sub write_to_file {
      my $self = shift;
      my $filename = $self->{-filename};
        die "no filename provided" if (!defined($filename) || $filename eq "");
      my $tree_string = $self->{-tree_string}; 
       open my $FILE_WRITE_FH, '>',$filename or die "\t\t\tCouldn't open $filename\n";
        print {$FILE_WRITE_FH} $tree_string;
        close $FILE_WRITE_FH or die "\t\t\tCould not close fasta file $filename\n"; 
    }
__PACKAGE__->meta->make_immutable;
1; 
