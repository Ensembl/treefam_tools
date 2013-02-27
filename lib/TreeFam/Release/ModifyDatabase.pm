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


package TreeFam::Release::DatabaseHandle;
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



  sub update_gene_tree_root_stable_ids {
      my $self = shift;
        print "TODO: Update stable_id columns in gene_tree_root
                  script in /nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/update_genetrees\n";
    return 1;
    }
 sub  search_replace_columns{
      my $self = shift;
      my $args = shift;
      my $table = $args->{"table"};
      my $column = $args->{"column"};
        # make backup of members table
        print "CREATE TABLE member_backup2 LIKE member;\n";
        print "INSERT INTO member_backup2 SELECT * FROM member;\n";
        # test
        print "update $table set $column = \"TCOGS2_TC002120-PA\" where member_id = 3832490; \n";
        print "update $table set $column = replace(stable_id,':','_')\n";
    return 1;
    }
  sub create_full_search_index {
      my $self = shift;
      my $args = shift;
      my $table = $args->{table};
      print "\tMaking full-text search index for table $table\n";
      print  "SET FOREIGN_KEY_CHECKS=0;\n";
      print  "alter table member engine = MyISAM;\n";
      print "SET FOREIGN_KEY_CHECKS=1;\n";
      print "CREATE FULLTEXT INDEX ft_index_description ON member (description);\n";
                                                     
    return 1;
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
