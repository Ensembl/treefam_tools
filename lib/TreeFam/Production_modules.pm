#created
#===============================================================================
#
#         FILE: HomologyHelper.pm
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
package TreeFam::Production_modules;

use strict;
use warnings;
use Bio::Phylo::Factory;
use Bio::Phylo::IO qw'parse unparse';
use Bio::Phylo::Util::CONSTANT qw':objecttypes :namespaces';
use Bio::TreeIO;
use Bio::Phylo::Forest::Tree;
use Data::Dumper;
use Bio::EnsEMBL::Registry; 
use File::Temp qw/ tempdir tempfile /;
use JSON;
sub get_db_connection{
	my ($arg_ref) = @_;
	my $registry = $arg_ref->{registry};
	my $registry_file = $arg_ref->{registry_file};
    
    my $reg = Bio::EnsEMBL::Registry->load_all($registry_file);
    return defined($reg)? $reg : undef;
}
sub get_tree_object{
	my ($arg_ref) = @_;
	my $registry = $arg_ref->{registry};
	my $entry = $arg_ref->{entry};
    
    my $genetree_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTree' );
    if(!$genetree_adaptor){
        die "Could not load genetree_adaptor\n";
    }

    #my $genetree_adaptor = $registry->get_GeneTreeAdaptor;
    my $tree             = $genetree_adaptor->fetch_by_root_id($entry);
    if ( !defined($tree) || $tree eq '' )
    {
        print "Could not get tree\n";
            return;
    }
	return $tree; 
}

sub get_pruned_tree{
	my ($arg_ref) = @_;
	my $superhash= $arg_ref->{superhash};
	my $sequence_href = $arg_ref->{sequences_href};
	my $treefam_name = $arg_ref->{treefam_name};
	my $tree = $arg_ref->{treeObject};

	my @ids2keep;
	my %allowed_species_ids = (
		"10116" => 1,
		"10090" => 1,
		"9606" => 1,
		"9031" => 1,
		"7955" => 1,
		"7227" => 1,
		"6239" => 1,
		"559292" => 1,
		"284812" => 1,
		"3702" => 1,

);
	print "Getting root node: ";
	# get root node:
	my $root = $tree->root;

	my $all_leaves = $root->get_all_leaves();
	print "Found ".scalar(@{$all_leaves})." leaves";
	foreach my $leaf (@{$all_leaves}){
		my ($stable_id, $taxon) = ($leaf->stable_id,$leaf->taxon_id);
		print "$stable_id: $taxon";
		push(@ids2keep, $stable_id) if(exists $allowed_species_ids{$taxon});	
	}
	print "Found ".scalar(@ids2keep)."";
	my $ret_tree = $root->keep_nodes_by_taxon_ids(\@ids2keep);
	print "found tree: $ret_tree\n";
	my $pruned_tree = $ret_tree->nhx_format();
	print Dumper $pruned_tree;
	exit;
}
sub get_alignment_information{
	my ($arg_ref) = @_;
	my $superhash= $arg_ref->{superhash};
	my $sequence_href = $arg_ref->{sequences_href};
	my $treefam_name = $arg_ref->{treefam_name};
	my $tree = $arg_ref->{treeObject};
	my $aln_length = $arg_ref->{aln_length};

	print "Getting root node: ";
	# get root node:
	my $root = $tree->root;
	# get all leaves
	my $nodes = $root->get_all_nodes();
	print "found ".scalar(@{$nodes})." nodes\n";
	foreach my $curr_node(@{$nodes}){
		#print Dumper $curr_node;
		#exit;
		my ($name, $node_id, $cigar_line)  = ($curr_node->name, $curr_node->node_id, $curr_node->consensus_cigar_line);
		#print $name."\t".$node_id."\t".$cigar_line."\n";
		next if not defined($name);
		my $binary_pattern = get_binaryPattern4cigar({ cigar_string => $cigar_line, alignment_length => $aln_length, no_of_splits => 100});
		
		$sequence_href->{$name}{binary_pattern} = $binary_pattern;
		print "saving binary pattern for $name: $binary_pattern\n";#print "binary pattern: $binary_pattern\n";	
		#exit;
	}
	# check if 
	
	#exit;
	return 1;
}

sub get_binaryPattern4cigar{
	my ($arg_ref) = @_;
	my $cigar_string = $arg_ref->{cigar_string};
	my $alignment_length = $arg_ref->{alignment_length};
	my $no_of_splits = $arg_ref->{no_of_splits};
	my $binary_pattern;	
	my $split_size = int($alignment_length / $no_of_splits);
	print "trying to get pattern for length $alignment_length and split size $split_size\n";

	# expand cigar
	
	my @splitted = split(/(\d*\w)/,$cigar_string);
	#my @splitted = split(/[D-Md-m]+/,$cigar_string);
	my @expanded_cigar;
	my $carry_forward;
	print "found ".scalar(@splitted)." splits\n";
	foreach my $split(@splitted){
		print "$split\n";
		next if $split =~ /^$/;
		$split =~ /^(\d*)([\w])/;
		my ($number,$type) = ($1,$2);
			print "$number and $type\n";
			if(!$type){ die "problem with line $split\n";}
			if(!$number){
				push(@expanded_cigar, $type);
			}
			else{
				for(my $i=0;$i<$number;$i++){
					push(@expanded_cigar, $type);
				}
			}		
	}
	my $expanded_string = join("",@expanded_cigar);
	for(my $i=0;$i<$alignment_length; $i+=$split_size){
		print "slicing $i\n";
		my $substring = substr($expanded_string,$i,$split_size);
		print "found subalignment: $substring\n";
		my %letter_hash;
		foreach(split(//, $substring)){
			$letter_hash{$_}++;
		}
	my @sorted = sort { $letter_hash{$a} <=> $letter_hash{$b} } keys %letter_hash; 
	print "result is ".$sorted[0]."\n";
	$binary_pattern .= ($sorted[0] eq "M" || $sorted[0] eq "m")? 1 : 0;
	}
	#print "expanded cigar = ".join("",@expanded_cigar)."\n";
return $binary_pattern;

}

sub clean_old_entries{
	my ($arg_ref) = @_;
	my $entry = $arg_ref->{entry};
	my $db_host = $arg_ref->{db_host};
	my $db_user = $arg_ref->{db_user};
	my $db_pass = $arg_ref->{db_pass};
	my $db_port = $arg_ref->{db_port};
	my $db_database = $arg_ref->{db_database};
 
### do we need to clean old entry?
            #my $treeObject = $superhash{'treeObject'};
            my @tags_to_delete = ("fam_description", "fam_n_full", "hgnc", "fam_n_seed", "fam_full_description","fam_symbol", "fasta_aln", "homologs_array_json","numspecies","pfam","treenhx","treephyloxml","tree_image_png", "sequence_array_json","taxa_count","wikigene", "numSequences", "json_tree", "minimal_species_tree_json");
            #foreach my $tag(@tags_to_delete){
                #print "Deleting tag $tag ...";
                ##print "\n";
                ##exit;
                #my $remove_success = $treeObject->remove_tag($tag);
                #print "".(($remove_success)? "yes":"no")."\n"; 
            #}
            foreach my $tag(@tags_to_delete){
                my $mysql_command = "mysql -h $db_host -P $db_port -u$db_user -p$db_pass -e 'DELETE FROM $db_database.gene_tree_root_tag where root_id = $entry and tag =\"$tag\";'";
                print $mysql_command."\n";
                system($mysql_command);
            }
            # Now delete entries from xrefID2Sequence, 
            my $mysql_command_xrefID2Sequence = "mysql -h $db_host -P $db_port -u$db_user -p$db_pass -e 'DELETE FROM $db_database.xrefID2Sequence where gene_tree_id =  $entry;'";
                print $mysql_command_xrefID2Sequence."\n";
                system($mysql_command_xrefID2Sequence);
            # Now delete entries from xrefID2Sequence, 
            my $mysql_command_xrefID2family = "mysql -h $db_host -P $db_port -u$db_user -p$db_pass -e 'DELETE FROM $db_database.xrefID2Family where gene_tree_id =  $entry;'";
            print $mysql_command_xrefID2family."\n";
            system($mysql_command_xrefID2family);
            
}

sub read_gene_tree_statistics{
	my ($arg_ref) = @_;
	my $gene_tree_string= $arg_ref->{gene_tree_string};
    my $node2gene_count = $arg_ref->{node2gene_count};
    
   my $io = IO::String->new($gene_tree_string);
   my $treeio = Bio::TreeIO->new(-fh => $io,
                                 -format => 'nhx');
  
    #print Dumper $node2gene_count;
    return 1; 
}
sub get_species_tree_statistics{
	my ($arg_ref) = @_;
	my $superhash= $arg_ref->{superhash};
	my $sequence_href = $arg_ref->{sequences_href};
	my $treefam_name = $arg_ref->{treefam_name};

    my $species_count = $arg_ref->{species_count};
    my $no_genes = $arg_ref->{no_genes};

	my $minimal_species_tree = $arg_ref->{minimal_species_tree};
    my $input_tree = "$treefam_name.in.newick";
    my %node2gene_count;
    print "get gene_tree_statistics\n";
    die "get_species_tree_statistics: no species counts available. Stopping\n" if !keys(%{$species_count}); 
    
    #my interesting_taxa = ("Mammalia","Arthropoda", "Primates","Nematoda", );

   my $tree = Bio::Phylo::IO->parse(
    '-file' => $minimal_species_tree,
    '-format' => 'newick',
 )->first;
       #print Dumper $tree; 
   #my $proj = parse(
    #'-format' => 'newick',
    #'-file' => $minimal_species_tree,
    #'-as_project' => 1,
    #);
#my ($forest) = @{ $proj->get_items(_FOREST_) };

#my $tree=$forest->next;
#print "here hash has ".keys(%{$sequence_href})." members \n";

my %species_tree_counts = (
                "Metazoa" => 104,
                "Primates" => 10,
                "Eukaryota" =>109,
                "Mammalia" => 40,
                "Nematoda" =>11,
                "Euarchontoglires" => 18,
                "Laurasiatheria" => 13,
                "Ecdysozoa" => 38,
                "Arthropoda" => 27,
                "Glires" => 7,
                "non-Mammals" => 8,
                "Vertebrata" => 56,
                "non-Vertebrata" => 2,
                "Chordata" => 58,
                "Outgroup" => 4,
);
my $result = traverse_species_tree({
                    "node" => $tree->get_root,
                    "species_count" => $species_count,
                    "no_genes" => $no_genes,
                    "species_tree_counts" => \%species_tree_counts,
            });
    print  JSON->new->encode($result);
    $superhash->{"minimal_species_tree_json"} = JSON->new->encode($result);
    return 1;
 

sub traverse_species_tree {
    my ($arg_ref) = @_;
	#print "started traverse_species_tree with ";
    #print Dumper $arg_ref->{node};
    my $no_genes = $arg_ref->{no_genes};
    my $node = $arg_ref->{node};
	my $species_count = $arg_ref->{species_count};
    my $species_tree_counts = $arg_ref->{species_tree_counts};
    my $node_name = $node->get_name;
    my $no_species = keys(%{$species_count});
    my $node = $arg_ref->{node};
    print "looking at node: $node_name\n";
   

    # we need to map the following taxa
    die "no species counts available. Stopping\n" if !keys(%{$species_count}); 
    #count number of species that match node
    my ($found_species,$found_genes) = (0,0); 
    foreach my $species(keys(%{$species_count})){
        my $classification = $species_count->{$species}{'classification'};
        print "has classification: $classification\n";
        my @matches = grep(/$node_name/,split(" ",$classification));
        if(scalar(@matches)){
            #print "Found matches: for $node_name \n";
            $found_genes += $species_count->{$species}{'count'};
            $found_species += 1;
            
        }
    }
    $no_species = $species_tree_counts->{$node_name};
    $no_species = (exists $species_tree_counts->{$node_name})?$species_tree_counts->{$node_name}: 0;
    print "$node_name: species (".$found_species." / $no_species), genes: (".$found_genes."/$no_genes)\n";  

    my $missing_species = ($found_species)?$no_species - $found_species:0;
    my $missing_genes = ($found_genes)?$no_genes - $found_genes:0;
	
	# map node names? 
	$node_name = ($node_name eq 'non-Mammals')?'Frogs/Lizards/Birds':$node_name;
	$node_name = ($node_name eq 'non-Vertebrata')?'Tunicates':$node_name;

   	my $result = { name=>$node_name , 
                    species_total=>int($no_species),
                    species_presence=>($found_species)?[int($missing_species),int($found_species)]: [$no_species,0],
                    gene_presence=> ($found_genes)?[int($missing_genes), int($found_genes)]: [$no_genes,0],
			};		
	if ( my @children = @{ $node->get_children } ) {
        print "traversing child\n";
        $result->{'children'} = [ map { 
           traverse_species_tree({
                    "node" => $_,
                    "species_count" => $species_count,
                    "no_genes" => $no_genes,
                    "species_tree_counts" => $species_tree_counts,
        }) 
            } @children ];
	}

	return $result;

}
}



sub newickspeciestree2json{
	my ($arg_ref) = @_;
	my $superhash= $arg_ref->{superhash};
	my $sequence_href = $arg_ref->{sequences_href};
	my $treefam_name = $arg_ref->{treefam_name};
	my $registry = $arg_ref->{registry};
	my $newick_tree = $arg_ref->{treeNhx};
	my $species_tree = $arg_ref->{species_tree};
	my $species2gene  = $arg_ref->{species2gene};    	

    #open my $nw_tree_out, ">", $input_tree or die "Could not open $input_tree\n";
    #print {$nw_tree_out} $newick_tree."\n";
    #close $nw_tree_out || die "Could not close $input_tree\n";
    #if(!-e $input_tree || ! -s $input_tree){
    #    warn "[newick2json] Problem saving tree to file $input_tree\n";
    #    return 0;
    #} 


my $tree = Bio::Phylo::IO->parse(
    '-file' => $species_tree,
    '-format' => 'newick',
 )->first;
print "species tree here hash has ".keys(%{$sequence_href})." members \n";

my $result = species_tree_traverse({
                    "node" => $tree->get_root,
                    "sequence_href" => $sequence_href,
                    "species2gene" => $species2gene,
            });

#my $result = traverse({"node"=> $tree->get_root, "sequence_href" => $sequence_href});
    #$superhash->{"json_tree"} = JSON->new->pretty->encode($result);
    
	$superhash->{"json_species_tree"} = JSON->new->encode($result);
    #unlink($output_tree) if -e $output_tree;
    return 1;

sub species_tree_traverse {
    my ($arg_ref) = @_;
	my $node = $arg_ref->{node};
	my $sequence_href = $arg_ref->{sequence_href};
	my $species2gene  = $arg_ref->{species2gene};    	
	my $node_name = $node->get_name;
    my @domains_for_seq;
 	print "checking species name: : for $node_name: ";
				#print Dumper $sequence_href->{$node_name};
			#print Dumper $allNodesMapping{$node_name};
 	if(exists $species2gene->{$node_name}){
		print "present!!!\n";	
	}
	else{
		print "not\n";
	}
		my @empty_array;
		my $result = { 		
					name=>($node_name)? $node_name:"NaN", 
					#duplication=> (exists $allNodesMapping{$node_name}{D} && $allNodesMapping{$node_name}{D} eq "Y")?"Y":"N", 
					taxon=> ($node_name)? $node_name:"NaN", 
					sequences=>(exists $species2gene->{$node_name})?$species2gene->{$node_name}:\@empty_array,
		};
        if ( my @children = @{ $node->get_children } ) {
			$result->{'children'} = [ map { 
           		species_tree_traverse({
                    "node" => $_,
                    "sequence_href" => $sequence_href,
                    "species2gene" => $species2gene,
        	}) 
            } @children ];
	}
	#print Dumper $result;

	return $result;
}
}
sub newick2json{
	my ($arg_ref) = @_;
	my $superhash= $arg_ref->{superhash};
	my $sequence_href = $arg_ref->{sequences_href};
	my $treefam_name = $arg_ref->{treefam_name};
	my $registry = $arg_ref->{registry};
	my $newick_tree = $arg_ref->{treeNhx};
    my $input_tree = "$treefam_name.in.newick";
	my $species2gene  = $arg_ref->{species2gene};    	

    open my $nw_tree_out, ">", $input_tree or die "Could not open $input_tree\n";
    print {$nw_tree_out} $newick_tree."\n";
    close $nw_tree_out || die "Could not close $input_tree\n";
    if(!-e $input_tree || ! -s $input_tree){
        warn "[newick2json] Problem saving tree to file $input_tree\n";
        return 0;
    } 
my $output_tree = "$treefam_name.output.tre";
my $taxon2species_file = "species_info.table";
my $taxon2color_file = "species_info.table";
my $debug = 0;
my %leafNodesMapping;
my %internalNodesMapping;
my %allNodesMapping;
my %taxid2species_hash;
my %taxid2color_hash;
my %mapBack2ID;
# read taxon2species

unlink($output_tree) if -e $output_tree;

# add colors
my @taxids2colors = `cat $taxon2color_file`;
foreach (@taxids2colors){
    chomp;
    #print "Have: $_\n";
    my ($taxid,$species_name,$common_name,$color) = split("\t",$_);
    print "Found $taxid,$species_name,$common_name,$color\n";
	$taxid2color_hash{$taxid} = ucfirst($color);
}

my @taxids2species = `cat $taxon2species_file`;
foreach (@taxids2species){
	next if /^taxon_id/;
	#/(\d+)\s+(\w+)/;
	chomp;
	#$taxid2species_hash{$1} = ucfirst($2);
	my ($id,$species_name,$common_name,$color) = split("\t",$_);
	$species_name =~ s/ /_/g;
	$taxid2species_hash{$id}{species_name} = ucfirst($species_name);
	$taxid2species_hash{$id}{common_name} = ucfirst($common_name);
	print "TAXA:saving $id: $species_name, $common_name\n";
}
if(!keys(%taxid2species_hash)){
	die "problem reading taxid_species information\n";
}
if(!read_tree_io({"input_tree"=> $input_tree, 
					"output_tree" => $output_tree, 
					"leafNodesMapping" => \%leafNodesMapping, 
					"internalNodesMapping" => \%internalNodesMapping,
					"allNodesMapping" => \%allNodesMapping,
					"mapBack2ID" => \%mapBack2ID
		})){
	die "could not read tree $input_tree\n";
}

unlink($input_tree);
#if(-e $input_tree){ warn "INPUT TREE STILL THERE\n"}
#if(-e $output_tree){ warn "OUTPUT TREE STILL THERE\n"}
my @terminal_ids = keys(%leafNodesMapping);
my %seqIDLength;
#unlink($input_tree) if -e $input_tree;

# get color codes
# have a predefined set of internal nodes

my %ID2Color = (
  "13735" => "green"  
);
# get sequence length
read_sequence_length({"registry" => $registry, "terminal_ids" => \@terminal_ids,  "seqIDLength" => \%seqIDLength});

################################################################################
########     Parse as Bio::Phylo::Project
################################################################################
# we parse the newick as a project so that we end
# up with an object that has both the tree and the
# annotated OTUs
my $proj = parse(
    '-format' => 'newick',
    '-file' => $output_tree,
    '-as_project' => 1,
);
unlink($output_tree) if -e $output_tree;
my ($forest) = @{ $proj->get_items(_FOREST_) };

my $tree=$forest->next;
print "here hash has ".keys(%{$sequence_href})." members \n";

my $result = traverse({
                    "node" => $tree->get_root,
                    "sequence_href" => $sequence_href,
                    "seqIDLength" => \%seqIDLength,
                    "species2gene" => $species2gene,
                    "taxid2species_hash" => \%taxid2species_hash,
                    "ID2color" => \%ID2Color,
                    "taxid2color_href" => \%taxid2color_hash
            });

#my $result = traverse({"node"=> $tree->get_root, "sequence_href" => $sequence_href});
    #$superhash->{"json_tree"} = JSON->new->pretty->encode($result);
    $superhash->{"json_tree"} = JSON->new->encode($result);
    unlink($output_tree) if -e $output_tree;
	unlink($input_tree) if -e $input_tree;
    return 1;

sub traverse {
    my ($arg_ref) = @_;
	my $node = $arg_ref->{node};
	my $sequence_href = $arg_ref->{sequence_href};
	my $seqIDLength= $arg_ref->{seqIDLength};
	my $ID2Color= $arg_ref->{ID2Color};
    my $taxid2color_href  = $arg_ref->{taxid2color_href};
    my $taxid2species_href  = $arg_ref->{taxid2species_hash};
	my $species2gene  = $arg_ref->{species2gene};    	
	my $node_name = $node->get_name;
    my @domains_for_seq;
   if(exists $leafNodesMapping{$node_name}){
   # add domain information
    print "domain? $node_name hash has ".keys(%{$sequence_href})."\n";

    foreach my $pfam_hit(keys(%{$sequence_href->{$node_name}{pfam_hits}})){
            #my $ratio = $ext_counts->{"pfam_counts"}{$pfam_hit} / keys(%{$sequence_href});
           		foreach my $domain_hit(keys(%{$sequence_href->{$node_name}{pfam_hits}{$pfam_hit}})){
                	my ($id,$name, $a_start,$a_end,$evalue)  =  ($sequence_href->{$node_name}{pfam_hits}{$pfam_hit}{$domain_hit}{hmm_name},
                                                            $sequence_href->{$node_name}{pfam_hits}{$pfam_hit}{$domain_hit}{description},
                                                            $sequence_href->{$node_name}{pfam_hits}{$pfam_hit}{$domain_hit}{alignment_start},
                                                            $sequence_href->{$node_name}{pfam_hits}{$pfam_hit}{$domain_hit}{alignment_end},
                                                            $sequence_href->{$node_name}{pfam_hits}{$pfam_hit}{$domain_hit}{evalue} );
            print "Hit: $id,$name, $a_start,$a_end,$evalue\n";
            my $temp_hash = {domain_start=>$a_start, domain_stop=>$a_end, evalue=> $evalue, name=>$name, id => $id};
               push(@domains_for_seq, $temp_hash);
        	}
		} 
    }
				print "checking SWISS: for $node_name: ";
				#print Dumper $sequence_href->{$node_name};
				#print Dumper $allNodesMapping{$node_name};
 	my $result = { name=>(exists $mapBack2ID{$node_name})? $mapBack2ID{$node_name}:$node_name , 
					duplication=> (exists $allNodesMapping{$node_name}{D} && $allNodesMapping{$node_name}{D} eq "Y")?"Y":"N", 
					taxon=> (exists $allNodesMapping{$node_name}{T})? $taxid2species_hash{$allNodesMapping{$node_name}{T}}{species_name}:"N/A", 
					bootstrap=>(exists $allNodesMapping{$node_name}{B})?$allNodesMapping{$node_name}{B}:"N/A",
					
					common_name=>(exists $taxid2species_href->{$allNodesMapping{$node_name}{T}}{common_name})?$taxid2species_href->{$allNodesMapping{$node_name}{T}}{common_name}:"no common name",
					
					binary_alignment=>(exists $sequence_href->{$node_name}{binary_pattern})?$sequence_href->{$node_name}{binary_pattern}:"N/A",
					uniprot_name=>(exists $sequence_href->{$node_name}{uniprot_hits}{'UniProtKB-ID'})?$sequence_href->{$node_name}{uniprot_hits}{'UniProtKB-ID'}:"N/A",
					display_label=>(exists $sequence_href->{$node_name}{display_label})?$sequence_href->{$node_name}{display_label}:"N/A",
					swissprot_gene_name=>(exists $sequence_href->{$node_name}{swissprot_hits}{'gene_name'})?$sequence_href->{$node_name}{swissprot_hits}{'gene_name'}:"N/A",
					swissprot_protein_name=>(exists $sequence_href->{$node_name}{swissprot_hits}{'protein_name'})?$sequence_href->{$node_name}{swissprot_hits}{'protein_name'}:"N/A",
					swissprot_lca_gene_name=>(exists $allNodesMapping{$node_name}{'SwissProtGene'})?$allNodesMapping{$node_name}{'SwissProtGene'}:"N/A",
			};
                   #print "adding display_label: ".$sequence_href->{$node_name}{display_label}."\n"; 

					$result->{type} = "leaf" if exists $leafNodesMapping{$node_name};
                    $result->{type} = "node" if not exists $leafNodesMapping{$node_name};
                    #$result->{domains} = \@domains_for_seq if exists $leafNodesMapping{$node_name}; 
		            $result->{domains} = [];
                    if( exists $leafNodesMapping{$node_name}){
                        foreach my $hashref (@domains_for_seq){
                                push(@{$result->{domains}}, $hashref);
                        }
                    }
		            #print "DOMAINS: ".join(",",@{$result->{domains}})."\n";

                    #print "looking for length of $node_name is ".$seqIDLength->{$node_name}."\n";	
                    $result->{seq_length} = $seqIDLength->{$node_name} if exists $seqIDLength->{$node_name}; 
		            #print "looking for color for ".$allNodesMapping{$node_name}{T}." is ".$ID2Color->{$allNodesMapping{$node_name}{T}}."\n";	
                    $result->{color} = $taxid2color_href->{$allNodesMapping{$node_name}{T}} if exists $taxid2color_href->{$allNodesMapping{$node_name}{T}}; 
                    #$result->{go} = $sequence_href->{$node_name}{go_hits} if exists $sequence_href->{$node_name}{go_hits};
	
	## save for building species tree
		my $species_name = $taxid2species_href->{$allNodesMapping{$node_name}{T}}{species_name};
		if($species_name && $species_name ne ""){
			print "SPECIES2GENE for ".$species_name."\n";
			push(@{$species2gene->{$species_name}}, $result);
			print "currently have ".$species2gene->{$species_name}."\n";
			print "SPECIES2GENE has: ".keys(%{$species2gene})." elements\n";
		}
	if ( my @children = @{ $node->get_children } ) {
		$result->{'children'} = [ map { 
           traverse({
                    "node" => $_,
                    "sequence_href" => $sequence_href,
                    "seqIDLength" => \%seqIDLength,
                    "taxid2species_hash" => $taxid2species_href,
                    "ID2color" => \%ID2Color,
                    "species2gene" => $species2gene,
                    "taxid2color_href" => $taxid2color_href
        }) 
            } @children ];
	}

	return $result;
}

sub read_tree_io{
	my ($arg_ref) = @_;
    my $input_tree = $arg_ref->{input_tree};
	my $output_tree = $arg_ref->{output_tree};
	my $leafNodesMapping = $arg_ref->{leafNodesMapping};
	my $internalNodesMapping = $arg_ref->{internalNodesMapping};
	my $allNodesMapping = $arg_ref->{allNodesMapping};
	my $mapBack2ID   = $arg_ref->{mapBack2ID};
	my %alreadyMapped;
	#my $temp_tree = $treefam_id."temp_tree.nhx";
	print "Reading tree from $input_tree\n";
    my $in = Bio::TreeIO->new( '-format' => 'nhx', '-file' => $input_tree );
	while( my $tree = $in->next_tree ) {
		my @nodes = $tree->get_nodes();
		print "Found ".scalar(@nodes)." nodes in total\n";
		die "Could not get nodes from tree" if !scalar(@nodes);
		foreach my $node(@nodes){
		    print "node_name is ".$node->id."  ";
            if($node->is_Leaf){
                print "\tleaf\n";
				my $name = $node->id;
				my @tags = $node->get_all_tags;
				foreach(@tags){
					#print "get value for $_\n";
					my @values = $node->get_tag_values($_);
					$leafNodesMapping{$name}{$_}= $values[0];	
					$allNodesMapping{$name}{$_}= $values[0];	
				}		
			}
			# internal node
			# need to replace ids
			else{
                print "\tinner node\n";
				my $name = $node->id;
				next if !defined $name || $name eq "";
				my @tags = $node->get_all_tags;
				my $newID = $name."_".$alreadyMapped{$name}++;
				$mapBack2ID{$newID} = $name;
				print "name is $name ($newID)\t" if $debug;
				foreach(@tags){
					my @values = $node->get_tag_values($_);
					$internalNodesMapping{$newID}{$_}= $values[0];	
					$allNodesMapping{$newID}{$_}= $values[0];	
					print "\t\tvalue for $_ is ".$values[0]."\n" if $debug;
				}	
				#print "BOOTSTRAP saved: ".$node->bootstrap."\n";
				$internalNodesMapping{$newID}{'bootstrap'} = $node->bootstrap;
				$allNodesMapping{$newID}{'bootstrap'} = $node->bootstrap;
				#$node->id($name);
				#$internalNodesMapping{$newID}{'bootstrap'} = $node->bootstrap;
				$node->id($newID);
			}
		}
	my $out = new Bio::TreeIO(-file => ">$output_tree", -format => 'nhx');
    $out->write_tree($tree);
}
	#unlink($input_tree) if -e $input_tree;
	return (-e $output_tree && -s $output_tree)? 1:undef;
}

}




sub get_summary{
	my ($arg_ref) = @_;
	my $tree = $arg_ref->{treeObject};
	my $c = $arg_ref->{superhash};
	my $entry = $arg_ref->{entry};

	#get tag-values for this tree
    my $tagvalue_hashref = $tree->get_tagvalue_hash();
    if ( !keys(%$tagvalue_hashref) )
    {
        die "Could not get tagvalue_hashref for tree\n";
    }
	$c->{numSequences}              = $tagvalue_hashref->{gene_count};
    $c->{tree_max_branch}           = $tagvalue_hashref->{tree_max_branch};
    $c->{tree_num_human_peps}       = $tagvalue_hashref->{tree_num_human_peps};
    $c->{tree_num_dup_nodes}        = $tagvalue_hashref->{tree_num_dup_nodes};
    $c->{aln_method}                = $tagvalue_hashref->{aln_method};
    $c->{aln_percent_identity}      = $tagvalue_hashref->{aln_percent_identity};
    $c->{aln_num_residues}          = $tagvalue_hashref->{aln_num_residues};
    $c->{tree_num_spec_nodes}       = $tagvalue_hashref->{tree_num_spec_node};
    $c->{aln_length}                = $tagvalue_hashref->{aln_length};
    $c->{aln_runtime}               = $tagvalue_hashref->{aln_runtime};
    $c->{tree_max_length}           = $tagvalue_hashref->{tree_max_length};
    $c->{buildhmm_runtime_msec}     = $tagvalue_hashref->{buildhmm_runtime_msec};
    $c->{njtree_phyml_runtime_msec} = $tagvalue_hashref->{njtree_phyml_runtime_msec};
    $c->{orthotree_runtime_msec}    = $tagvalue_hashref->{orthotree_runtime_msec};
    $c->{tree_num_leaves}           = $tagvalue_hashref->{tree_num_leaves};
    return 1;
}
sub get_sequences{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $sequence_hash = $arg_ref->{sequences_hash};
	my $tree = $arg_ref->{treeObject};
	my $species_count = $arg_ref->{species_count};
	
# Get all sequences
    my $all_leaves = ( $tree->root )->get_all_leaves();
    my @leaf_array;    
    my %leaf_hash;
    my $counter = 0;
    my @array_of_arrays;
    my $symbol4family = "NaN";
   # get data in array format
    foreach my $leaf(@{$all_leaves}){
        my @tmp_array;
        my $prot_seq_id = $leaf->stable_id;
        ### Save each sequence with 
        # use member id to search hmmer_scores 
        $sequence_hash->{$prot_seq_id}{"member_id"} = $leaf->dbID;
        # gene member id
        my $gene4prot = $leaf->gene_member;
        $sequence_hash->{$prot_seq_id}{"geneID"} = $gene4prot->stable_id;
        # sequence object
        #$sequence_hash->{$prot_seq_id}{"sequence_object"} = $leaf;
        my $species_name = ($leaf->taxon)->binomial;
        if(!defined($species_name) || $species_name eq ""){
            warn "Problem with ".$leaf->taxon->name."\n";
            $species_name = $leaf->taxon->name;
        }
        # taxon ID
        $sequence_hash->{$prot_seq_id}{"taxonID"} = $leaf->taxon_id;
        $sequence_hash->{$prot_seq_id}{"display_label"} = uc($gene4prot->display_label);
        # Taxon name
        $sequence_hash->{$prot_seq_id}{"taxon_name"} = $species_name;
        $species_name =~ s/\s/_/;
        #my $species_image_file = "http://localhost:3000/static/images/species_pictures/species_files/thumb_".$species_name.".png";
        my $species_image_file = "thumb_".$species_name.".png";
        # Taxon image
        $sequence_hash->{$prot_seq_id}{"taxon_image"} = $species_image_file;
	    # Gene description
        $sequence_hash->{$prot_seq_id}{"description"} = $gene4prot->description ;

        #push(@tmp_array, "<img src=''$species_image_file'' /> ".($leaf->taxon)->binomial);
        #push(@tmp_array,$leaf->stable_id );
        #push(@tmp_array,(defined $leaf->display_label ? $leaf->display_label: 'NaN'));
		#if(defined $leaf->display_label){ 
		#$symbol4family = $leaf->display_label;
		#	# grep description
		#	$symbol4family =~ s/-00\d//;
		#	my $grep_cmd = "grep -w \"$symbol4family\" /Users/fs9/Downloads/mart_export.txt";
		#	#print "$grep_cmd\n";
		#	my $grepped_line;# = `$grep_cmd`;
		#	chomp($grepped_line);
		#	my ($EnsemblGeneID,$EnsemblProteinID,$HGNCsymbol,$PFAMID,$WikiGeneDescription,$GOTermAccession,$GOTermName,$GOTermDefinition,$GOdomain,$GOTermEvidenceCode,$RfamID) = split("\t",$grepped_line);
		#	#print "grepped_line: $grepped_line";
		#	#$c->stash->{treefam}{pfamID} = (defined $PFAMID ? $PFAMID : "No Pfam ID");
		#	#$c->stash->{treefam}{Wikigenedescription} = (defined $WikiGeneDescription ? $WikiGeneDescription : "No WikiGene description");
		#	#$c->stash->{treefam}{GOTermAccession} = (defined $GOTermAccession ? $GOTermAccession : "No GOTermAccession");
		#	#
		#}
		#push(@tmp_array,$leaf->description );
        #$c->log->debug( "Found line:  ".join(",",@tmp_array)."  members" ) if $c->debug;
        #push(@array_of_arrays,\@tmp_array);
    }
    #$c->{'sequence_array_json'} =  encode_json \@array_of_arrays;
   
    #$c->log->debug('Family::Tree::get_summary decoded all sequences') if $c->debug; 
    # Count species
    #$c->{'numSpecies'} = keys(%species_count);
    foreach(@$all_leaves){ 
        $species_count->{$_->taxon->name}{'count'} += 1;
		# grep entry from file
		my $grep_line = "grep ".$_->taxon->taxon_id." id2classification.txt";
		my $grepped_classification = `$grep_line`;
		#print "grep ".$_->taxon->taxon_id." id2classification.txt\n";
		chomp($grepped_classification);
        $species_count->{$_->taxon->name}{'classification'} = $grepped_classification;
    }
    $c->{'numSpecies'} = keys(%{$species_count});
    $c->{'numSequences'} = scalar(@{$all_leaves});
    $c->{'tree_num_leaves'} = scalar(@{$all_leaves});
    
    #$c->log->debug( "Found " . keys(%species_count) . " species" ) if $c->debug;
    return 1;
}


sub transfer_annotations{
	my ($arg_ref) = @_;
    my $treeNhx= $arg_ref->{treeNhx};
    my $sequences_hash= $arg_ref->{sequences_hash};
    my $member_adaptor= $arg_ref->{member_adaptor};
    my $homology_adaptor= $arg_ref->{homology_adaptor};
	my $treefam_name = $arg_ref->{treefam_name};
	my %annotation_counter; 
	my %swissprot_sequences;
	my $transfer_type = "swissprot-not";	
	
    #if($transfer_type eq "swissprot"){
    	
		if($transfer_type eq "swissprot"){
		# collect all sequences with swissprot annotation
    	foreach my $ensembl_prot_id(keys(%{$sequences_hash})){
				if(exists $sequences_hash->{$ensembl_prot_id}{'uniprot_hits'}{'UniProtKB-AC'}){
					my $uniprot_ac = $sequences_hash->{$ensembl_prot_id}{'uniprot_hits'}{'UniProtKB-AC'};
					#print "have uniprot entry '$uniprot_ac'\n";
					my $swissprot_hits = &get_swissprot_hits({"id" => $uniprot_ac,  "db_adaptor" => $member_adaptor});	
					next if !keys(%{$swissprot_hits});
					$swissprot_sequences{$ensembl_prot_id} = $swissprot_hits;
					#print "swissprot hits: \n";
					#print Dumper $swissprot_hits;
				#	die "Lets stop here\n";
			}	
	
		}
		}
		else {
    		foreach my $ensembl_prot_id(keys(%{$sequences_hash})){
					my %seq_hash;
					print "getting display label for $ensembl_prot_id: ".$sequences_hash->{$ensembl_prot_id}{display_label}."\n";
					my $display_label = $sequences_hash->{$ensembl_prot_id}{display_label};
					$swissprot_sequences{$ensembl_prot_id}{protein_names} = $display_label; 
					$swissprot_sequences{$ensembl_prot_id}{gene_names} =  $display_label;
					$annotation_counter{$display_label}++; 
					#$swissprot_sequences{$ensembl_prot_id} = \%seq_hash;
			}
		}
		# remove bad annotations. e.g. ones that occur only once
		foreach my $ensembl_prot_id(keys(%swissprot_sequences)){
			if($annotation_counter{$swissprot_sequences{$ensembl_prot_id}{protein_names}} <= 1){
				$swissprot_sequences{$ensembl_prot_id}{protein_names} = ""; 
				$swissprot_sequences{$ensembl_prot_id}{gene_names} =  "";
				print "removing annotation for $ensembl_prot_id\n";
			}
		}
		print Dumper %swissprot_sequences;
		print "Found ".keys(%swissprot_sequences)." sequences with swissprot annotation\n";
		my %already_annotated;
		foreach my $ensembl_id(keys(%swissprot_sequences)){
				# get orthologs
				print "swissprot sequence: $ensembl_id\n";
				my ($gene_name, $protein_name) = ($swissprot_sequences{$ensembl_id}{'gene_names'},
														  $swissprot_sequences{$ensembl_id}{'protein_names'});
				next if $protein_name eq "";
				next if $already_annotated{$ensembl_id};
				my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLPEP',$ensembl_id);
				if(!defined($member) || $member eq ''){
					warn "Could not get member object for $ensembl_id\n";
				}
				my $gene_member = $member->gene_member();
				# then you get the homologies where the member is involved
				my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);
				warn "Could not get homologies for ".$gene_member->stable_id."\n" if !scalar(@{$homologies});
				my ($genename4group,$proteinname4group);
				foreach my $homology (@{$homologies}) {
				print "HOMMO: ".$homology->description," ", $homology->subtype,"\n";	
					next if $homology->description =~ /paralog/;
					foreach my $hom_member (@{$homology->get_all_Members}) {
						my $orth_prot_id = $hom_member->stable_id;
						next if $already_annotated{$orth_prot_id};
						next if $hom_member->stable_id eq $ensembl_id;
						if($swissprot_sequences{$ensembl_id}{'protein_names'} ne ""){
							print $hom_member->stable_id."is an ortholog, but already has annotation: ".$swissprot_sequences{$ensembl_id}{'protein_names'}."\n";
							next;
						}	
						# transfer annotation here
						print "\ttransfer annotation to ".$hom_member->stable_id." with ".$protein_name."\n";
							$sequences_hash->{$orth_prot_id}{'swissprot_hits'}{'protein_name'} = ($protein_name eq '')? 'NaN': $protein_name;
							if($gene_name =~ / /){
								my @split_gene_name = split(/ /,$gene_name);
								print "only taking first one: ".$split_gene_name[0]."\n";
								$gene_name = uc($split_gene_name[0]);
								$genename4group = uc($split_gene_name[0]);
							}
							else{ 
								$gene_name = ($gene_name eq '')? 'NaN': uc($gene_name);
							}
							$sequences_hash->{$orth_prot_id}{'swissprot_hits'}{'gene_name'} = $gene_name;
							$already_annotated{$orth_prot_id}	 = 1; # mark as annotated
							$proteinname4group = $protein_name;
						
						print "collect swiss entries: SWISS: $orth_prot_id:";
						#print Dumper $sequences_hash->{$orth_prot_id};
					}
				}
			# add annotation of sequence itself
				$sequences_hash->{$ensembl_id}{'swissprot_hits'}{'protein_name'} = $swissprot_sequences{$ensembl_id}{'protein_name'};
				#$sequences_hash->{$ensembl_id}{'swissprot_hits'}{'gene_name'} = $genename4group;
				my $gene_name = $swissprot_sequences{$ensembl_id}{'gene_names'};
				if($gene_name =~ / /){
								my @split_gene_name = split(/ /,$gene_name);
								print "only taking first one: ".$split_gene_name[0]."\n";
								$gene_name = uc($split_gene_name[0]);
								$genename4group = uc($split_gene_name[0]);
							}
							else{ 
								$gene_name = ($gene_name eq '')? 'NaN': uc($gene_name);
							}
				$sequences_hash->{$ensembl_id}{'swissprot_hits'}{'gene_name'} = $gene_name;
				#$swissprot_sequences{$ensembl_id}{'gene_names'} = $gene_name;
				$already_annotated{$ensembl_id}	 = 1; # mark as annotated
			}


	#}
	#die "for testing we die here\n";	

	# now we could label mrca of 
  	#print "trying to read tree: $treeNhx\n"; 
	#my $tree = Bio::Phylo::IO->parse(
   # 	'-string' => $treeNhx,
   # 	'-format' => 'Newick',
 	#)->first;
    my $input_tree = "$treefam_name.in.newick";
	open my $nw_tree_out, ">", $input_tree or die "Could not open $input_tree\n";
    print {$nw_tree_out} $$treeNhx."\n";
    close $nw_tree_out || die "Could not close $input_tree\n";
    if(!-e $input_tree || ! -s $input_tree){
        warn "[transfer annotation] Problem saving tree to file $input_tree\n";
        return 0;
    } 	
	my $treeio = Bio::TreeIO->new(-format => 'nhx',
			      -file => $input_tree);
	my $tree = $treeio->next_tree; 
	#print Dumper $tree;
	my %same_annotation_hash;
	my @terminals = $tree->get_leaf_nodes;
	#my @terminals = @{ $tree->get_terminals };
	foreach my $terminal(@terminals){
		print "Looking at ".$terminal->id."\n";
		my $gene_name = $sequences_hash->{$terminal->id}{'swissprot_hits'}{'gene_name'};
		push(@{$same_annotation_hash{$gene_name}}, $terminal);
		print "added to array: $gene_name, has now ".scalar(@{$same_annotation_hash{$gene_name}})." entries\n";
	}
	foreach my $same_annotation(keys(%same_annotation_hash)){
		my @nodes = @{$same_annotation_hash{$same_annotation}};
		next if scalar(@nodes) < 2 || $same_annotation eq "";		
		print "checking annotation $same_annotation with ".scalar(@nodes)." entries \n";
		my $lca = $tree->get_lca(-nodes => \@nodes);
		#my $mrca = $tree->get_mrca(\@nodes);
		#$mrca->set_generic( 'SwissProtGene' => $same_annotation );
		$lca->set_tag_value('SwissProtGene' , $same_annotation);
		print "found mrca with ".$lca->id."\n";
	}
	my $write2tree = "$treefam_name.nhx";
    my $out = new Bio::TreeIO(-file => ">$write2tree",
                          -format => 'nhx');
    $out->write_tree($tree);
	my $new_tree = `cat $treefam_name.nhx`;
	chomp($new_tree);
	$$treeNhx = $new_tree;
	unlink("$write2tree");
	#my $new_tree = Bio::Phylo::IO->unparse(
    #'-phylo' => $tree,                         
    #'-format' => 'Newick',
	#-nhxstyle => 'nhx'
 #);
#print $new_tree."\n";
#die "stopping here\n";
unlink($input_tree);	

	return 1;
}


sub get_sequence_annotations{
	my ($arg_ref) = @_;
	my $sequences_hash = $arg_ref->{sequences_hash};
	my $json_entry = $arg_ref->{sequence_array_json};
    my $hgnc_file= $arg_ref->{hgnc_file};
    my $ext_counts = $arg_ref->{ext_counts};
    my $wikigene_file= $arg_ref->{wikigene_file};
    my $pfam_file = $arg_ref->{pfam_file};
    my $hmmer_scores_file= $arg_ref->{hmmer_scores_file};
    my $uniprot_file= $arg_ref->{uniprot_file};
    my $single_gene_list= $arg_ref->{single_gene_list};
    my $db_adaptor= $arg_ref->{db_adaptor};
    
    my $counter = 0;
    my @all_sequences_array;
    print "Getting annotations for ".keys(%{$sequences_hash})." sequences\n";
    foreach my $ensembl_prot_id(keys(%{$sequences_hash})){
        #print "looking at $ensembl_prot_id\n";
        #######################################
        # HGNC
        #######################################
        my $gene4prot = $sequences_hash->{$ensembl_prot_id}{"geneID"};
        $sequences_hash->{$ensembl_prot_id}{"hgnc_hits"} = &get_hgnc_hits({"id" => $gene4prot, "file_to_search" => $hgnc_file, "db_adaptor" => $db_adaptor});
        $sequences_hash->{$ensembl_prot_id}{"wikigene_hits"} = &get_wikigene_hits({"id" => $gene4prot, "file_to_search" => $wikigene_file, "db_adaptor" => $db_adaptor });
        #print "WIKIGENE all: ".$sequences_hash->{$ensembl_prot_id}{"wikigene_hits"}."\n"; 
        #my $description_hits = &get_description_hits($leaf);
        #my $go_hits = &get_go_hits($leaf);
        $sequences_hash->{$ensembl_prot_id}{"pfam_hits"} = &get_pfam_hits({"id" => $ensembl_prot_id, "file_to_search" => $pfam_file, "pfam_counts" => \%{$ext_counts->{"pfam_counts"}}, "db_adaptor" => $db_adaptor});
        #print Dumper $sequences_hash->{$ensembl_prot_id}{"go_hits"} ;
        my $member_id4prot = $sequences_hash->{$ensembl_prot_id}{"member_id"};
        #$sequences_hash->{$ensembl_prot_id}{"hmmer_hits"} = &get_hmmer_hits({"id" => $member_id4prot, "file_to_search" => $hmmer_scores_file, "db_adaptor" => $db_adaptor});
        $sequences_hash->{$ensembl_prot_id}{"uniprot_hits"} = &get_uniprot_hits({"id" => $ensembl_prot_id, "file_to_search" => $uniprot_file, "db_adaptor" => $db_adaptor});
        if($ensembl_prot_id =~ /ENSP\d+/){
            #print "checking SINGLE: $ensembl_prot_id\n";
            my @single_gene_hit = `grep -w '$ensembl_prot_id' $single_gene_list`;
            #print "SINGLE: ".join(@single_gene_hit)."\n";
            $sequences_hash->{$ensembl_prot_id}{"single_copy_gene"} = scalar(@single_gene_hit)? 1:0;
            $sequences_hash->{$ensembl_prot_id}{"go_hits"} = &get_go_hits({"id" => $ensembl_prot_id, "file_to_search" => $pfam_file, "pfam_counts" => \%{$ext_counts->{"pfam_counts"}}, "db_adaptor" => $db_adaptor});
        }
    } 
    foreach my $ensembl_prot_id(keys(%{$sequences_hash})){
        my @pfam_hits;
        my @seq_array;
        #print Dumper $sequences_hash;
        #exit;
        # Seq order is: taxon, prot_id,
        push(@seq_array, (exists $sequences_hash->{$ensembl_prot_id}{"taxon_name"} && $sequences_hash->{$ensembl_prot_id}{"taxon_name"} ne '')?$sequences_hash->{$ensembl_prot_id}{"taxon_name"}: "NaN");
        push(@seq_array, $ensembl_prot_id);
        push(@seq_array, (exists $sequences_hash->{$ensembl_prot_id}{"hmmer_hits"} && $sequences_hash->{$ensembl_prot_id}{"hmmer_hits"} ne '')?$sequences_hash->{$ensembl_prot_id}{"hmmer_hits"}: "No hits" );
        if(exists $sequences_hash->{$ensembl_prot_id}{"pfam_hits"}){
                #print "$ensembl_prot_id has ".keys(%{$sequences_hash->{$ensembl_prot_id}{"pfam_hits"}})." Pfam hits. \n";
                foreach my $acc(keys(%{$sequences_hash->{$ensembl_prot_id}{"pfam_hits"}})){
                	foreach my $domain_hit(keys(%{$sequences_hash->{$ensembl_prot_id}{"pfam_hits"}{$acc}})){
                    	#print "$acc hit has ".$ext_counts->{'pfam_counts'}{$acc}." other sequences (total: ".keys(%{$sequences_hash}).")....";
                    	my $ratio = $ext_counts->{"pfam_counts"}{$acc} / keys(%{$sequences_hash});
                    	if ( $ext_counts->{"pfam_counts"}{$acc} && $ratio < 0.05){
                #       		print "SKIP  this entry (ratio: $ratio) \n";
                 #       	print "delete this entry ";
                        	#delete($sequences_hash->{$ensembl_prot_id}{"pfam_hits"}{$acc});
                    	}
                    	else{
                  #      	print "\n";
                        	push(@pfam_hits, $sequences_hash->{$ensembl_prot_id}{"pfam_hits"}{$acc}{$domain_hit}{"description"});
                    	}
					}
                }
        }
        
        if(!scalar(@pfam_hits)){
            push(@seq_array, "No hits");
        }
        else{
            push(@seq_array, join(",",@pfam_hits));
        }
        #my $hgnc_hits_string;
        #if(keys(%{$sequences_hash->{$ensembl_prot_id}->{"hgnc_hits"}})){ 
        #    $hgnc_hits_string =  join(",",keys(%{$sequences_hash->{$ensembl_prot_id}{"hgnc_hits"}}));
        #}
        #else{
        #    $hgnc_hits_string = "NaN";
        #}
        #push(@seq_array, $hgnc_hits_string);
        if(!exists($sequences_hash->{$ensembl_prot_id}{"description"}) ||  !defined($sequences_hash->{$ensembl_prot_id}{"description"}) || $sequences_hash->{$ensembl_prot_id}{"description"} eq ''){
            #print "no description for $ensembl_prot_id\n";
            push(@seq_array, "No description");
        }
        else{
            push(@seq_array, $sequences_hash->{$ensembl_prot_id}{"description"});
        }
        push( @all_sequences_array,\@seq_array);
       #last;
    }
       #die "we stop here, ok?"; 
    #$c->{"hgnc"}  = join(",",@hgnc_hits_array);
    #$c->{"pfam"}  = join(",",keys(%pfam_hits));
    $$json_entry = encode_json \@all_sequences_array;

    return 1;
    }

sub get_family_annotations{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $tree = $arg_ref->{treeObject};
	my $species_count = $arg_ref->{species_count};
	my $sequences_hash = $arg_ref->{sequences_hash};
	my $previous_treefam_file = $arg_ref->{previous_treefam_file};
  ## get previous treefam annotations;
    my ($acc,$symbol,$desc,$seed,$full,$date1,$date2,$type,$full_desc,$comment) = &get_previous_treefam_annotation({"id" => $tree->stable_id, "file_to_search" => $previous_treefam_file});
    $c->{fam_symbol} = $symbol || "NaN";
    $c->{fam_description} = $desc || "NaN";
    $c->{fam_n_seed} = $seed;
    $c->{fam_n_full} = $full;
    #$c->{fam_full_description} = $full_desc;
    print "Symbol: $symbol, desc: $desc, n_seed: $seed, n_full: $full\n";
    #return 1;
    ### will do some summarizing things here as well
    # get hgnc symbols
    my $hgnc_string;
    my %pfam_hits;
    my %hgnc_hits;
    my $total_number_of_seqs =  keys(%{$sequences_hash});
    foreach my $seq(%{$sequences_hash}){
       # wikigene
       my %wikigene_hits;
        use Data::Dumper;
        #print Dumper $sequences_hash->{$seq}{'wikigene_hits'};
        print "WIKIGENE: ".keys(%{$sequences_hash->{$seq}{'wikigene_hits'}})." entries\n";
        if(keys(%{$sequences_hash->{$seq}{'wikigene_hits'}})){
            foreach my $wikigene_hit(keys(%{$sequences_hash->{$seq}->{wikigene_hits}})){
                $wikigene_hits{$wikigene_hit} = 1;
            print "WIKIGENE: add $wikigene_hit\n";
            }
        }
        # HGNC
        if(keys(%{$sequences_hash->{$seq}{'hgnc_hits'}})){
            foreach my $hgnc_hit(keys(%{$sequences_hash->{$seq}{hgnc_hits}})){
                $hgnc_hits{$hgnc_hit} = $sequences_hash->{$seq}{hgnc_hits}{$hgnc_hit};
            }
        }
        # PFAM
        # to avoid counting multiple domains for the same sequence
        my %domain_per_sequence;
        if(keys(%{$sequences_hash->{$seq}{'pfam_hits'}})){
            foreach my $pfam_hit(keys(%{$sequences_hash->{$seq}{pfam_hits}})){
            	foreach my $domain_hit(keys(%{$sequences_hash->{$seq}{pfam_hits}{$pfam_hit}})){
                	$pfam_hits{$pfam_hit}{id} = $sequences_hash->{$seq}{pfam_hits}{$pfam_hit}{$domain_hit}{hmm_name};
                	$pfam_hits{$pfam_hit}{name} = $sequences_hash->{$seq}{pfam_hits}{$pfam_hit}{$domain_hit}{description};
                	$pfam_hits{$pfam_hit}{count}++ if !exists $domain_per_sequence{$seq}{$pfam_hits{$pfam_hit}{name}};
            		print "PFAM: count: ".$pfam_hits{$pfam_hit}{count}." id: ".$pfam_hits{$pfam_hit}{id}." name: ".$pfam_hits{$pfam_hit}{name}."\n"  ;
					$domain_per_sequence{$seq}{$pfam_hits{$pfam_hit}{name}} = 1;
				}
			}
        }
        my @temp_hgnc_array;
        foreach my $id(keys(%hgnc_hits)){
            if(!exists $hgnc_hits{$id}){
                print "PARSE_ERROR: HGNC no id for $id\n";
            }
            else{
                push(@temp_hgnc_array, $id."->NaN");
            }
        }
        $c->{'hgnc'} = join(" ", @temp_hgnc_array);
        # PFAM
        my @temp_pfam_array;
        my @sorted_pfam_domains = sort { $pfam_hits{$a} <=> $pfam_hits{$b} } keys %pfam_hits;
        foreach my $id(@sorted_pfam_domains){
            my $percent = (int(($pfam_hits{$id}{count} * 100)/$total_number_of_seqs)) || 1;   
            push(@temp_pfam_array, $pfam_hits{$id}{name}."->".$pfam_hits{$id}{id}."->".$percent);
        }

        $c->{'pfam'} = join(" ",@temp_pfam_array);
        $c->{'wikigene'} = join(" ",keys(%wikigene_hits));
    }
   
    ## get overrepresented taxa
    my $total_no_seq_species = 0;

    my @species_array;
    my $top = 3;
    my $count_top = 1;
    foreach my $key (sort { $species_count->{$b}{'count'} <=> $species_count->{$a}{'count'}} keys %$species_count ){
           push(@species_array, $key."->".$species_count->{$key});
           print "TAXA_COUNT: ".$key."->".$species_count->{$key}."\n";
           last if $count_top++ >= $top;
    }
    $c->{'taxa_count'} = join(" ",@species_array);
    
    print "TAXA_COUNT: ".$c->{'taxa_count'};
    print "HGNC: ".$c->{'hgnc'};
    print "PFAM ".$c->{'pfam'};
   #exit; 
    return 1;
}


sub get_tree_nhx{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $tree = $arg_ref->{treeObject};
	my $entry = $arg_ref->{entry};
    my $root_of_tree = $tree->root;
    $c->{'treeNhx'} = $root_of_tree->nhx_format;
  
    return defined($c->{'treeNhx'})? 1 : undef ;
}

sub get_alignment{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $tree = $arg_ref->{treeObject};
	my $entry = $arg_ref->{entry};

# Get Alignment
    my $simpleAlign = $tree->get_SimpleAlign();
    my $fasta_string;
    foreach my $seq ( $simpleAlign->each_seq ){
                $fasta_string .= ">" . $seq->id . "\n" . $seq->seq . "\n";
    }
    $c->{'fasta_aln'} = $fasta_string;
    
  return 1;
  
}
sub get_hgnc_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
	my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select hgnc_id,app_symbol from hgnc_mapping where ensembl_gene_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
 
    my %hgnc_hits;
### Parse
	while ( my ($hgnc_id,$app_symbol) = $extID2seq_sth->fetchrow_array() ){
        $hgnc_hits{$app_symbol} = $hgnc_id;
    }
    if(keys(%hgnc_hits)){
         #print Dumper %hgnc_hits;
        #die "finished for now";
    } 
    return  \%hgnc_hits;
}

sub get_wikigene_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select wikigene_name,wikigene_id,wikigene_description from wikigene where ensembl_gene_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
 
    my %wikigene_hits;
    ### Parse
	while ( my ($wikigene_name,$wikigene_id,$wikigene_description) = $extID2seq_sth->fetchrow_array() ){
            $wikigene_hits{$wikigene_name}{id} = $wikigene_id;
            $wikigene_hits{$wikigene_name}{description} = $wikigene_description;
    } 
    return  \%wikigene_hits;
}

sub get_go_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
	#my $go_counts = $arg_ref->{go_counts};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select * from go_mapping where seq_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
 
    ### Parse
    my %go_counts;
    my %go_hits;
    my $no_seq = 0;
	while ( my ($acc,$go_id,$go_name,$go_definition,$go_evidence,$go_namespace) = $extID2seq_sth->fetchrow_array() ){
### Parse
    #foreach (@result_lines){
        #print "fetched $_\n";    
        #my @array = split(/\t/,$_);
        #my ($ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan) = split(/\s+/,$_);
        print "$acc,$go_id,$go_name,$go_definition,$go_evidence,$go_namespace\n";
        #my $no_seq = exists $go_counts{$acc}? $go_counts{$acc}:1;
        $go_hits{$no_seq}{go_id} = $go_id;
        $go_hits{$no_seq}{go_name} = $go_name;
        $go_hits{$no_seq}{go_definition} = $go_definition;
        $go_hits{$no_seq}{go_evidence} = $go_evidence;
        $go_hits{$no_seq}{go_namespace} = $go_namespace;
        $no_seq++;
        $go_counts{$acc}++;
    } 
    #print Dumper $pfam_hits_href;
    #exit;
    return \%go_hits;
}

sub get_pfam_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
	my $pfam_counts = $arg_ref->{pfam_counts};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select * from pfam_hits where seq_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
 
    my %wikigene_hits;
    ### Parse
    my %pfam_hits;
	 
	while ( my ($ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan) = $extID2seq_sth->fetchrow_array() ){
    #my $grep_command = "grep \"$id\" $file_to_search";
    #my @result_lines = `$grep_command`;
    #my @hgnc_hits;
### Parse
    #foreach (@result_lines){
        #print "fetched $_\n";    
        #my @array = split(/\t/,$_);
        #my ($ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan) = split(/\s+/,$_);
        #print "$ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan\n";
        my $counter = (exists $pfam_counts->{$acc})?$pfam_counts->{$acc}: 1;   
		$pfam_hits{$acc}{$counter}{alignment_start} = $startA;
        $pfam_hits{$acc}{$counter}{alignment_end} = $endA;
        $pfam_hits{$acc}{$counter}{count}++;
        $pfam_hits{$acc}{$counter}{description} = $name;
        $pfam_hits{$acc}{$counter}{hmm_name} = $acc;
        $pfam_hits{$acc}{$counter}{bitscore} = $bitscore;
        $pfam_hits{$acc}{$counter}{confidence} = $evalue;
        $pfam_hits{$acc}{$counter}{evalue} = $evalue;
        $pfam_counts->{$acc}++;
    } 
    #print Dumper %pfam_hits;
    #exit;
    return \%pfam_hits;
}

sub get_uniprot_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select * from uniprot_mapping where seq_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my %uniprot_hits;
### Parse
	while ( my ($seq_id,$db,$ext_seq_id) = $extID2seq_sth->fetchrow_array() ){
        $uniprot_hits{$db} = $ext_seq_id;
    } 
    return \%uniprot_hits;
}
sub get_swissprot_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select * from swissprot_mapping where entry = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my %swissprot_hits;
### Parse
	while ( my ($entry,$entry_name,$protein_names,$gene_names) = $extID2seq_sth->fetchrow_array() ){
        $swissprot_hits{"protein_names"} = $protein_names;
        $swissprot_hits{"gene_names"} = $gene_names;
        $swissprot_hits{"entry_name"} = $entry_name;
    } 
    return \%swissprot_hits;
}

sub get_hmmer_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select evalue from hmmer_scores where member_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    my $hmm_score; 
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my %uniprot_hits;
### Parse
	while ( my ($evalue) = $extID2seq_sth->fetchrow_array() ){
        $hmm_score = $evalue;
        last;
    } 
    return $hmm_score;
}

sub get_previous_treefam_annotation {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $grep_command = "grep \"^$id\" $file_to_search";
    my @result_lines = `$grep_command`;
    #print "grep $grep_command\n";
### Parse
    my ($acc,$symbol,$desc,$empty,$seed,$full,$empty2,$date1,$date2,$type,$full_desc,$comment) = split(/\t/,$result_lines[0]);
    #print "$acc,$symbol,$desc,$empty,$seed,$full,$empty2,$date1,$date2,$type,full:$full_desc,$comment\n";
    #exit;
    return ($acc,$symbol,$desc,$seed,$full,$date1,$date2,$type,$full_desc,$comment);
}

sub read_species_info_table{
	my ($arg_ref) = @_;
	my $species_file = $arg_ref->{species_file};
	my $taxonIDMappings_href = $arg_ref->{taxonIDMappings_href};
	my @lines_of_file = `cat $species_file`;
	foreach my $species(@lines_of_file){
		chomp($species);
		my ($taxID,$sciName,$comName) = split(/\t/,$species);
		$taxonIDMappings_href->{$taxID}{'scientific_name'} = $sciName;
		$taxonIDMappings_href->{$taxID}{'common_name'} = $comName;
	}
	print "Found ".scalar(@lines_of_file)." species infos\n";
}
sub read_domain_info_file{
	my ($domain_file,$seqIDMappings_href) = (@_);
	my @lines_of_file = `cat $domain_file`;
	my %seq_counter;
	foreach my $line(@lines_of_file){
		chomp($line);
		my ($seqID,$alignment_start,$alignment_end,$envelope_start,$envelope_end,$hmm_acc,$hmm_name,$type,$hmm_start,$hmm_end,$hmm_length,$bit_score,$E_value,$significance,$clan,$predicted_active_site_residues) = split(/\s+/,$line);
		my $counter = (exists $seq_counter{$seqID})? $seq_counter{$seqID}+1 : 1;
		$seqIDMappings_href->{$seqID}{$counter}{"alignment_start"} = $alignment_start;
		$seqIDMappings_href->{$seqID}{$counter}{"alignment_end"} = $alignment_end;
		$seqIDMappings_href->{$seqID}{$counter}{"hmm_acc"} = $hmm_acc;
		$seqIDMappings_href->{$seqID}{$counter}{"hmm_name"} = $hmm_name;
		$seqIDMappings_href->{$seqID}{$counter}{"confidence"} = $E_value;
		$seq_counter{$seqID}++;
	}
	#print "Found ".scalar(@lines_of_file)." species infos\n";
	#print Dumper $seqIDMappings_href;
#exit;
}

sub read_sequence_length{
	my ($arg_ref) = @_;
	my $registry = $arg_ref->{registry};
	my $id_aref = $arg_ref->{terminal_ids};
	my $seqIDLengthMappings_href = $arg_ref->{seqIDLength};
    
    my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
	foreach my $id(@{$id_aref}){
    	print "READ_SEQUENCE_LENGTH: ENSEMBLPEP -> $id\n";
	    my $member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLPEP",$id);
		my $member_length = $member->seq_length;
		$seqIDLengthMappings_href->{$id} = $member_length;
	}
	print "Found seq length for ".keys(%{$seqIDLengthMappings_href})." ids\n";
}

sub get_tree_phyloxml{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $sequences_hash = $arg_ref->{sequences_hash};
	my $tree = $arg_ref->{treeObject};
	my $newickTree = $arg_ref->{treeNhx};
	my $registry = $arg_ref->{registry};
	my $entry = $arg_ref->{entry};
	my $treefam_id = $arg_ref->{treefam_name};
    my $phyloxml_file = $entry."test.phyloxml";
    #my $image_prefix = "http://static.ensembl.org/i/species/48/";
    my $image_prefix = "/static/images/species_pictures/species_files/";
    my %internalNodesMapping;
    my %leafNodesMapping;
    my %alreadyMapped;
    my $debug = 1;
    my ($fh, $tree_file) = tempfile( );
   
    # save newick string as file
    open my $nw_tree_out, ">", $tree_file or die "Could not open $tree_file\n";
    print {$nw_tree_out} $newickTree."\n";
    close $nw_tree_out || die "Could not close $tree_file\n";
    if(!-e $tree_file || ! -s $tree_file){
        warn "Problem saving tree to file $tree_file\n";
        return 0;
    } 

my $temp_tree = $treefam_id."temp_tree.nhx";
my $in = Bio::TreeIO->new( '-format' => 'nhx', '-file' => $tree_file );
while( my $tree = $in->next_tree ) {
	my @nodes = $tree->get_nodes();
	#print "Found ".scalar(@nodes)." nodes in total\n";
	die "Could not get nodes from tree" if !scalar(@nodes);
	foreach my $node(@nodes){
		
		if($node->is_Leaf){
			my $name = $node->id;
			my @tags = $node->get_all_tags;
			foreach(@tags){
				#print "get value for $_\n";
				my @values = $node->get_tag_values($_);
				$leafNodesMapping{$name}{$_}= $values[0];	
			}		
			#print "Found leaf: $name with \n";
			#print Dumper @tags;
			#exit;
		}
		# internal node
		# need to replace ids
		else{
			my $name = $node->id;
			next if !defined $name || $name eq "";
			my @tags = $node->get_all_tags;
			my $newID = $name."_".$alreadyMapped{$name}++;
			foreach(@tags){
				#print "get value for $_\n";
				my @values = $node->get_tag_values($_);
				$internalNodesMapping{$newID}{$_}= $values[0];	
			}	
            print "BOOTSTRAP saved: ".$node->bootstrap."\n";
			$internalNodesMapping{$newID}{'bootstrap'} = $node->bootstrap;
			$node->id($newID);
				
	
		}
		
	}
   #my $bio_phylo_tree = Bio::Phylo::Forest::Tree->new_from_bioperl($tr);

 #print Dumper %leafNodesMapping;
 #print Dumper %internalNodesMapping;
	my $out = new Bio::TreeIO(-file => ">$temp_tree", -format => 'nhx');
    $out->write_tree($tree);
#print $tree->to_string;
}
################################################################################
########     Parse as Bio::Phylo::Project
################################################################################
# we parse the newick as a project so that we end
# up with an object that has both the tree and the
# annotated OTUs
my $proj = parse(
    '-format' => 'newick',
    '-file' => $temp_tree,
    '-as_project' => 1,
);
#unlink($temp_tree);
# here we make the OTUs
my ($forest) = @{ $proj->get_items(_FOREST_) };
#print "Found ".scalar(@{ $proj->get_items(_FOREST_) })." nodes\n";
#my @nodes = @{ $proj->get_items(_NODE_) };
# it's easier to make a factory object for creating the annotations
my $fac = Bio::Phylo::Factory->new;
my $tree = $forest->first;

################################################################################
########     Species info
################################################################################

my $species_info_file = "species_info.table";
my %taxonIDMappings;
&read_species_info_table({"species_file" => $species_info_file, "taxonIDMappings_href" => \%taxonIDMappings});
die "problem reading species infos\n" if !keys(%taxonIDMappings);

################################################################################
########     Sequence length info
################################################################################
# Collect terminal taxa and get sequence length
my @all_terminal_taxa = @{$tree->get_terminals()};
my @all_terminal_ids = map {$_->get_name } @all_terminal_taxa;
my %seqIDLength;

&read_sequence_length({"registry" => $registry, "terminal_ids" => \@all_terminal_ids, "seqIDLength" => \%seqIDLength});

die "problem reading seq_length infos\n" if !keys(%seqIDLength);
#print Dumper %seqIDLength;
$tree->visit_depth_first(
		    '-in' => sub {
				    my $current_node = shift;			
				    return if !defined($current_node->get_name) || $current_node->get_name eq "";
				    print "looking at node ".$current_node->get_name."\n" if $debug;
				    ## INNER NODES
				    if($current_node->is_internal){ 
							    my $node_name = $current_node->id;
							    print "\tis an internal node\n" if $debug;
							    if(!exists $internalNodesMapping{$node_name}){
									    warn "Could not find ".$node_name." \n";
							    }
							    if(exists $internalNodesMapping{$node_name}{'bootstrap'}){
									    #print "setting bootstrap for $node_name\n";
									    #print "before: ".$current_node->get_generic('bootstrap');
									    #$current_node->set_generic( 'bootstrap' => $internalNodesMapping{$node_name}{'bootstrap'} );
									    #$current_node->set_score( $internalNodesMapping{$node_name}{'bootstrap'} );
									    #print "after: ".$current_node->get_generic('bootstrap')."\n";
									    #$current_node->score($internalNodesMapping{$node_name}{'bootstrap'});
									    my $bootstrap =$internalNodesMapping{$node_name}{'bootstrap'} ;
                                        if(!$bootstrap || $bootstrap eq ""){
                                            $bootstrap = 0;
                                        }
                                        print "setting bootstrap: $bootstrap\n"; 
                                        #$current_node->add_meta(
									        #$fac->create_meta(
									    	    #'-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
									    	    #'-triple' => { 'pxml:confidence' => $bootstrap },
									        #)
        							    #);
							    }
							    #return;
						        ## label events	
							    if(exists $internalNodesMapping{$node_name}{'D'}){
								    my $type_of_event =$internalNodesMapping{$node_name}{'D'} ;
								    #$current_node->score($internalNodesMapping{$node_name}{'bootstrap'});
								    my $event = _create_dummy_event(($type_of_event eq 'N')?1:0);
        						    print "EVENT: $event\n";
                                    #$event = '<speci><>';
                                    $current_node->add_meta(
                					    $fac->create_meta(
                        				    '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
                        				    '-triple' => { 'pxml:events' => $event },
                					    )
        						    );
							    }
							    # replace name
							    $node_name =~ s/_\d+//g;
							    $current_node->set_name($node_name);
				    }
				    ## TERMINAL
				    if($current_node->is_terminal){
					    my $node_name = $current_node->id;
                        # we could rename leaves to a uniprot-style: CCNB1_HUMAN
					    print "\tterminal node: ".$node_name."\n";
                        if(exists $sequences_hash->{$node_name}{'uniprot_hits'}->{'UniProtKB-ID'}){
					        #if(exists $ens2uniprot{$node_name}){
                            print "\twe can replace it with ".$sequences_hash->{$node_name}{'uniprot_hits'}->{'UniProtKB-ID'}."!\n";
                            $current_node->id($sequences_hash->{$node_name}{'uniprot_hits'}->{'UniProtKB-ID'});
					    }
				    #print "\tis an terminal node\n";
							    if(exists $leafNodesMapping{$node_name}{'T'}){
								        my $tax_id = $leafNodesMapping{$node_name}{'T'};
								        my $scientific_name = $taxonIDMappings{$tax_id}{'scientific_name'};
								        my $common_name = $taxonIDMappings{$tax_id}{'common_name'};
	 							        if($tax_id eq "" || $scientific_name eq "" || $common_name eq ""){
										        warn "\tCould not find tax_id: $tax_id and/or scientific_name: $scientific_name and/or common_name: $common_name\n";
								        }
								        $scientific_name = "NaN" if $scientific_name eq "";
								        $tax_id = "0" if $tax_id eq "";
								        $common_name = "NaN" if $common_name eq "";
								        print "taxonomy with $tax_id,$scientific_name,$common_name\n";
								        my $image_file = "thumb_";
								        my $abused_scientific_name = $scientific_name;
								        $abused_scientific_name =~ s/ /_/g; 
								        my $image_path = "http://web-treefamdev.internal.sanger.ac.uk:3000/".$image_prefix."/thumb_".$abused_scientific_name.".png";
								        print "image_path: $image_path\n";
								        my $tax = _create_dummy_taxonomy($tax_id,$scientific_name,$common_name, $image_path);
								        $current_node->add_meta(
                				        $fac->create_meta(
                        			        '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
                        			        '-triple' => { 'pxml:taxonomy' => $tax },
                					        )
        						        );
							    }
						        print "checking if $node_name is in seq hash ....";
                                if(exists $sequences_hash->{$node_name}{'pfam_hits'}){
									        print "yes\n";
                                            my $number_of_predictions = keys(%{$sequences_hash->{$node_name}{'pfam_hits'}});
										        print "\tfound $number_of_predictions domain predictions for $node_name\n";
									        if($number_of_predictions){
										        if(!exists $seqIDLength{$node_name}){ die "no seq length for $node_name\n";}
										        my $sequence_length = $seqIDLength{$node_name};
										        if($sequence_length eq ""){ die "no given seq length for $node_name\n";}
                                                my $arch = _create_dummy_architecture($sequence_length,\%{$sequences_hash->{$node_name}{'pfam_hits'}});
											        $current_node->add_meta(
											        $fac->create_meta(
												        '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
												        '-triple' => { 'pxml:sequence' => $arch },
											        )
										        );
									    }
                                        else{
                                           my $sequence_length = $seqIDLength{$node_name};
                                                print "Writing empty architecture for $node_name with seq length: $sequence_length\n";
                                                my $arch = _create_dummy_empty_architecture($sequence_length);
											        $current_node->add_meta(
											        $fac->create_meta(
												        '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
												        '-triple' => { 'pxml:sequence' => $arch },
											        )
										        );
 
 
                                            
                                            }
						    }
                            else{
                                                                                print "no\n";
                                }
    
				    }
		    }
    );
    # now write the output
    #print unparse( '-format' => 'phyloxml', '-phylo' => $proj );
    
    open my $phyloxmlOUT_FH, ">", $phyloxml_file or die "Could not open $phyloxml_file";
    print {$phyloxmlOUT_FH} unparse( '-format' => 'phyloxml', '-phylo' => $proj );
    close $phyloxmlOUT_FH or die "Could not close $phyloxml_file\n";

    ### HACK!
    # unable to set type of confidence
    my $in_place_replacement = "perl -i -pe 's/<confidence/<confidence type=\"bootstrap\"/g' $phyloxml_file";
    print "REPLACE $in_place_replacement\n";
    `$in_place_replacement`;
    # perl  -pe 's/<confidence/<confidence type="bootstrap"/g' 911644test.phyloxml
    $c->{'treephyloxml'} =  `cat $phyloxml_file`;
    $c->{'treephyloxml'} =~ s/\n//g ;
    #unlink($phyloxml_file);
	return 1; 
}

sub _create_dummy_architecture {
	my ($sequence_length,$domain_info_href) = (@_);
	my $domain_string = "<domain_architecture length=\"$sequence_length\">";
	foreach my $domain_number(keys(%$domain_info_href)){
		my($from,$to,$confidence, $hmm_name) = ($domain_info_href->{$domain_number}{'alignment_start'},
												$domain_info_href->{$domain_number}{'alignment_end'},
												$domain_info_href->{$domain_number}{'confidence'}, 
												$domain_info_href->{$domain_number}{'hmm_name'} );
		$domain_string .= "<domain from=\"$from\" to=\"$to\" confidence=\"$confidence\">$hmm_name</domain>";
	}
		$domain_string .= "</domain_architecture>";
}
sub _create_dummy_empty_architecture {
	my ($sequence_length) = (@_);
	my $domain_string = "<domain_architecture length=\"$sequence_length\">";
		$domain_string .= "<domain from=\"1\" to=\"2\" confidence=\"1\">No domains</domain>";
	$domain_string .= "</domain_architecture>";
}
sub _create_dummy_event {
	my $speciation = shift;
	return ($speciation == 1 )? '<speciations>1</speciations>' :  '<duplications>1</duplications>' ;
}

sub _create_dummy_bootstrap {
	my $bootstrap = shift;
	return "<confidence type=\"bootstrap\">$bootstrap</confidence>";
}
sub _create_dummy_taxonomy {
	my ($tax_id,$scientific_name,$common_name,$image_path) = (@_);
	## make it nice
	if($common_name =~ m/kangaroo/){
		print "before $common_name\n";
		
	}
	$scientific_name =~ s/\'//;
	$common_name =~ s/\'//;
	print "after $common_name\n";
	return  
	"<scientific_name>$scientific_name</scientific_name>
	<common_name>$common_name</common_name>";
    #<uri>$image_path</uri>";
} 
1;
