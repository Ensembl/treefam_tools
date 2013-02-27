#
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
package TreeFam::HomologyHelper;

use strict;
use warnings;
 
sub get_member_object{
	my ($arg_ref) = @_;
	my $member_object = $arg_ref->{'member_adaptor'};
	my $source_id = $arg_ref->{'to_search'};
	my $member = $member_object->fetch_by_source_stable_id(undef,$source_id);
	print "getting member object\n";
	return defined($member)?$member:undef;
}
sub get_genetree_object{
	my ($arg_ref) = @_;
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $genetree = $genetree_adaptor->fetch_by_stable_id($to_search);
	return defined($genetree)?$genetree:undef;
}
sub get_all_species_for_family{
	my ($arg_ref) = @_;
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my %have_species;
	my @all_species;
	my $gt = &get_genetree_object({"genetree_adaptor" => $genetree_adaptor,"to_search" => $to_search});
	if(!$gt){ print "Could not find genetree for $to_search\n"; return undef;}
	my $all_members = ($gt->root)->get_all_leaves();
	foreach my $member(@{$all_members}){
	my $species_name = $member->taxon->name;
		push(@all_species, $species_name)  if !exists $have_species{$species_name};
		$have_species{$species_name} = 1;
	}

    return scalar(@all_species)?\@all_species:undef;
}
sub get_all_homologs_for_family{
	my ($arg_ref) = @_;
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $node_id = $arg_ref->{'to_search'};
	my $type = $arg_ref->{'type'};
	my $homologies = $homology_adaptor->fetch_all_by_tree_node_id($node_id);
	my $encoded_homologies = encode_homologies($homologies,$type);
	return defined($encoded_homologies)?$encoded_homologies:undef
}
sub get_member_pair_species{
	my ($arg_ref) = @_;
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $second_species = $arg_ref->{'second_species'};
	my $source_member_object = $arg_ref->{'source_member_object'};
	my $homologies = $homology_adaptor->fetch_all_by_Member_paired_species($source_member_object, $second_species);
    return defined($homologies)?$homologies:undef
}
sub get_homologs_for_gene{
	my ($arg_ref) = @_;
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $source_member_object = $arg_ref->{'source_object'};
	my $type = $arg_ref->{'type'};
	my $homology_type= $arg_ref->{'homology_type'};
	my $memberID = $source_member_object->dbID;
	# the api requires member ids (int) rather than member objects
	print "have homology adaptor";
	my $homologies = $homology_adaptor->fetch_all_by_Member($source_member_object);
	#my $homologies = $homology_adaptor->fetch_all_by_Member($source_member_object);
	#return undef if !scalar(@$homologies);
	my $genetree = &get_genetree_for_member({"genetree_adaptor" => $genetree_adaptor, "source_object" => $source_member_object});
	my $tf_family = $genetree->[0]->stable_id."\n";
	
	my $encoded_homologies = encode_homologies({ "homologies" =>  $homologies, "type" => $type, "homology_type" => $homology_type,"tf_family" => $tf_family});
	#print Dumper $encoded_homologies;
	return defined($encoded_homologies)?$encoded_homologies:undef
}
sub get_homology_relation_for_gene_pair{
	my ($arg_ref) = @_;
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $source_member_object = $arg_ref->{'source_object'};
	my $source_member2_object = $arg_ref->{'source_object2'};
	my ($memberID1,$memberID2) = ($source_member_object->dbID, $source_member2_object->dbID);
	# the api requires member ids (int) rather than member objects
	my $homologies = $homology_adaptor->fetch_by_Member_id_Member_id($memberID1,$memberID2);
	my $genetree = &get_genetree_for_member({"genetree_adaptor" => $genetree_adaptor, "source_object" => $source_member_object});
	my $tf_family = $genetree->[0]->stable_id."\n";
	my $encoded_homologies = encode_homologies($homologies,"JSON", "",$tf_family);
    return defined($homologies)?$homologies:undef;
}
sub get_genetree_for_member{
	my ($arg_ref) = @_;
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $member = $arg_ref->{'source_object'};
	print "\tsearching for genetree \n";
	my $genetrees = $genetree_adaptor->fetch_all_by_Member($member, -CLUSTERSET_ID => "default");
	return scalar(@{$genetrees})? $genetrees : undef;
}	
sub get_homologs_for_family_level{
	my ($arg_ref) = @_;
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $gt_object = $arg_ref->{'genetree_object'};
	my $type = $arg_ref->{'type'};
	my $taxon_level = $arg_ref->{'taxon_level'};
	
	print "\tsearching for $taxon_level with ".$gt_object->stable_id."\n";
	my $root = $gt_object->root;
	## get all nodes that match taxon level 
	# this can be multiple
	#my $homologies = $source_genetree_object->get_all_nodes_by_tagvalue('taxon_name'=>$taxon_level);
	my $taxon_level_nodes= $root->get_all_nodes_by_tag_value('taxon_name'=>$taxon_level);
	my @all_homologies;
	foreach my $node(@{$taxon_level_nodes}){
		my $homologies4node = $homology_adaptor->fetch_all_by_tree_node_id($root->node_id);
		push(@all_homologies, @{$homologies4node});
	}
    return defined(\@all_homologies)?\@all_homologies:undef
}
sub encode_homologies{
	my ($arg_ref) = @_;
	my $homologies = $arg_ref->{'homologies'};
	my $type = $arg_ref->{'type'};
	my $homology_type = $arg_ref->{'homology_type'};
	my $tf_family = $arg_ref->{'tf_family'};
	my @results_array;
	my @sequences_array;
	my %sequences_remember_hash;
	foreach my $this_homology (@$homologies) {
		if(defined($homology_type)  && $homology_type ne "" && $this_homology->description ne $homology_type){
			#print "skipping type: ".$this_homology->description."\n";
			next;
		}
		my $homologue_genes = $this_homology->gene_list;
    	if(scalar(@{$homologue_genes}) > 2){
			die "too many genes: ".join(" and ", @$homologue_genes), " are ", $this_homology->description, "\n";
		}
		else{
			my ($source,$target) = ($homologue_genes->[0],$homologue_genes->[1]);
			my %pair_hash;
			my ($source_protein_id,$target_protein_id) = ($source->get_canonical_Member->stable_id,$target->get_canonical_Member->stable_id);
			my ($source_protein_sequence,$target_protein_sequence) = ($source->get_canonical_Member->sequence(),$target->get_canonical_Member->sequence());
			# source
					$pair_hash{"source"} = {
							"protein_id" => $source_protein_id,
							"species" => $source->taxon->name,
							"id" => $source->stable_id,
							"sequence" => $source_protein_sequence
						};
			# target
					$pair_hash{"target"} = {
							"protein_id" => $target_protein_id,
							"species" => $target->taxon->name,
							"id" => $target->stable_id,
							"sequence" => $target_protein_sequence
						};
			# type
				$pair_hash{"tf_family"} = $tf_family;
				$pair_hash{"type"} = $this_homology->description();
				$pair_hash{"subtype"} = $this_homology->subtype();
				
				#my $target_fasta = ">".$target_protein_id."\n".$target_protein_sequence;	
				#my $source_fasta = ">".$source_protein_id."\n".$source_protein_sequence;
				#push(@sequences_array, \$target_fasta) if !exists $sequences_remember_hash{$target_protein_id};
				#push(@sequences_array, \$source_fasta) if !exists $sequences_remember_hash{$source_protein_id};
				push(@results_array, \%pair_hash);		
			# Remember we saved the sequences
			#$sequences_remember_hash{$target_protein_id} = 1;
			#$sequences_remember_hash{$source_protein_id} = 1;
		}
  	}
	#print "found ".scalar(@results_array)." homologies and ".scalar(@sequences_array)." sequences\n";
  #return \@results_array,\@sequences_array;
  return \@results_array;
}


sub get_homologs_for_family_orthoxml{
	my ($arg_ref) = @_;
	my $genetree_object = $arg_ref->{'genetree_object'};
	my $write_type = $arg_ref->{'write_type'};
	
	my ($w,$string_handle,$file_handle);	
	my $file_name = "output.xml";
	if($write_type eq "file"){
		$file_handle = IO::File->new($file_name, 'w');
  		$w = Bio::EnsEMBL::Compara::Graph::OrthoXMLWriter->new(	-SOURCE => 'Ensembl', -SOURCE_VERSION => 67, -HANDLE => $file_handle);
	}
	elsif($write_type eq "string"){
		$string_handle = IO::String->new();
  		$w = Bio::EnsEMBL::Compara::Graph::OrthoXMLWriter->new(
    		-SOURCE => 'Ensembl', -SOURCE_VERSION => 67, -HANDLE => $string_handle
  		);
	}
	else{ die "Wrong type for writing given for get_homology_for_family_orthoxml: $write_type\n";}
  	$w->write_trees($genetree_object);
  	$w->finish(); #YOU MUST CALL THIS TO WRITE THE FINAL TAG
	if($write_type eq "string"){
  		my $xml_scalar_ref = $string_handle->string_ref();
		return defined($xml_scalar_ref)?$xml_scalar_ref:undef;	
	}
	elsif($write_type eq "file"){
		$file_handle->close();
		return (-e $file_name && -s $file_name)? $file_name:undef;
	}
	else{
		die "wrong option\n";
	}
}

sub check_existence {
	my ($object,$type) = shift;
	if(!$object){ die "Could not get member object for $type. Stopping\n"}
}
1;
