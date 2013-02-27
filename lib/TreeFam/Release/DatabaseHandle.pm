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
        print "CREATE TABLE ".$table."_backup LIKE member;\n";
        print "INSERT INTO ".$table."_backup SELECT * FROM member;\n";
        # test
        print "update $table set $column = \"TCOGS2_TC002120-PA\" where member_id = 3832490; \n";
        print "update $table set $column = replace(stable_id,':','_')\n";
    return 1;
    }
  sub create_full_search_index {
      my $self = shift;
      my $args = shift;
      my $table = $args->{table};
      my $column = $args->{"column"};
      print "\tMaking full-text search index for table $table\n";
      print  "SET FOREIGN_KEY_CHECKS=0;\n";
      print  "alter table $table engine = MyISAM;\n";
      print "SET FOREIGN_KEY_CHECKS=1;\n";
      print "CREATE FULLTEXT INDEX ft_index_description ON $table ($column);\n";
                                                     
    return 1;
    }
sub backup_table {
      my $self = shift;
      my $args = shift;
      my $table = $args->{table};
      my $backup_table = $table.".backup";
      my $backup_cmd = "CREATE TABLE $backup_table LIKE $table;";
      print $backup_cmd."\n";
      my $backup_copy_cmd = "INSERT INTO $backup_table SELECT * FROM $table;";
      print $backup_copy_cmd."\n";

      print "\tMaking full-text search index for table $table\n";

    }
  sub create_tables_from_file {
      my $self = shift;
      my $args = shift;
      my $create_tables_file = $args->{create_tables_file};

      print "mysql < $create_tables_file\n";

    }
  sub populate_tables_from_file{
      my $self = shift;
      my $args = shift;
      my $populate_tables_file = $args->{populate_tables_file};
    
        open my $POPULATE_TABLES_FH, '<',$populate_tables_file or die "\t\t\tCouldn't open $populate_tables_file\n";
        while(<$POPULATE_TABLES_FH>){
            next if /^$/;
            next if /^#/;
            chomp;
            /LOAD DATA LOCAL INFILE '(.*)' INTO/;
            my $input_file = $1;
            die "empty/no input file: $input_file\n" if(!-e $input_file && -s $input_file);
            print "file is: ".$input_file."\n";
            #print $_."\n";
        }
        close $POPULATE_TABLES_FH or die "Could not close $populate_tables_file\n";
    }


sub write_to_file {
      my $self = shift;
      my $args = shift;
      my $filename = $self->{-filename};
        die "no filename provided" if (!defined($filename) || $filename eq "");
      my $tree_string = $self->{-tree_string}; 
       open my $FILE_WRITE_FH, '>',$filename or die "\t\t\tCouldn't open $filename\n";
        print {$FILE_WRITE_FH} $tree_string;
        close $FILE_WRITE_FH or die "\t\t\tCould not close fasta file $filename\n"; 
}

sub get_all_families {
      my $self = shift;
      my $args = shift;
    my $families_file = $args->{"old_families_file"};
    my $family_file = $args->{"families_file"};
    my $genetree_adaptor = $args->{"genetree_adaptor"};
    
    die "no filename provided" if (!defined($family_file) || $family_file eq "");
    print "read all treefam families from $families_file\n";
    
    my @all_families = `cat $families_file`;
    my %treefam8_hash;
    foreach(@all_families){
        chomp;
        my ($fam,$symbol,$desc,$prevS,$prevF) = split(/\t/);
        $treefam8_hash{$fam}{'symbol'} = $symbol;
        $desc =~ s/'|"//g;
        $treefam8_hash{$fam}{'desc'} = $desc;
        $treefam8_hash{$fam}{'prevS'} = $prevS;
        $treefam8_hash{$fam}{'prevF'} = $prevF;
    }
    if(!keys(%treefam8_hash)){
        die "well, could not read treefam8 families. will stop here\n";
        return undef;
    }
    #my $genetree_adaptor = $db->get_GeneTreeAdaptor;
    #my $tree = $genetree_adaptor->fetch_by_stable_id("TF101003");
    my @all_trees = @{$genetree_adaptor->fetch_all( -CLUSTERSET_ID => "default")};
    my @all_trees; 
    #push(@all_trees,$tree);
    print "Found ".scalar(@all_trees)." trees\n";
    my $familyCount = 0;
    open my $file_out,">", $family_file or die "Could not open file\n";
    print {$file_out} "[\n";
    foreach my $gt(@all_trees){
        my %gt_names;
        my $tagvalue_hashref = $gt->get_tagvalue_hash();
        if(!keys(%$tagvalue_hashref)){
            die "Could not get tagvalue_hashref for tree\n";
        }
	    my ($modelName,$alnPercentIdentity,$alnLength,$geneCount, $pfam) = (
		        $tagvalue_hashref->{model_name},
		    $tagvalue_hashref->{aln_percent_identity},
		    $tagvalue_hashref->{aln_length},
		    $tagvalue_hashref->{gene_count},
		    $tagvalue_hashref->{pfam},
    
	    );
        ## get taxonomic distribution
        ## pfam domains?
        next if $modelName !~ /TF1/;
        #next if $modelName !~ /TF101003/;
	    $alnPercentIdentity = ($alnPercentIdentity eq "")? "NaN": int($alnPercentIdentity);
	    $alnLength = ($alnLength eq "")? "NaN": $alnLength;
	    $geneCount = ($geneCount eq "")? "NaN": $geneCount;
        my @pfams;
        if(defined($pfam) && $pfam ne ''){
            my @array = split(" ",$pfam);         
            foreach(@array){
                my @arr2 = split("->",$_);
                push(@pfams, $arr2[0]);
            }
        }
        my $pfams = (scalar(@pfams))? join(",",@pfams): "No hits";
        $gt_names{$modelName}{geneCount} = ($geneCount eq "")? "NaN": $geneCount;
	    
    my $gt_root = $gt->root;
    
    my $root_taxon = "NaN";
    if($gt_root->has_tag('taxon_name')){
        $root_taxon = $gt_root->get_value_for_tag('taxon_name');
    }

    print {$file_out}  "{\"modelName\":\"$modelName\", \"hgncSymbol\":\"".$treefam8_hash{$modelName}{symbol}."\", \"geneCount\":\"$geneCount\", \"percentIdentity\":\"$alnPercentIdentity\",\"alnLength\":\"$alnLength\",\"rootTaxon\":\"$root_taxon\",\"description\":\"".$treefam8_hash{$modelName}{desc}."\"},\n";
	#last if $familyCount++ > 300;
}

    print {$file_out} "]";
    close $file_out || die "Could not close families file for writing $family_file\n";       

    return (-e $family_file && -s $family_file)? $family_file: undef;    
}

__PACKAGE__->meta->make_immutable;
1; 
