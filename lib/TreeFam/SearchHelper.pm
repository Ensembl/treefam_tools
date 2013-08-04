#
#===============================================================================
#
#         FILE: SearchHelper.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Fabian Schreiber (), fs9@sanger.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 09/26/2012 16:13:16
#     REVISION: ---
#===============================================================================
package TreeFam::SearchHelper;

use strict;
use warnings;
use Data::Dumper;
use TreeFam::HomologyHelper;


sub search_family_with_id{
	my ($arg_ref) = @_;
	my $genetree_adaptor= $arg_ref->{'genetree_adaptor'};
	my $to_search= $arg_ref->{'to_search'};
	if (!defined($to_search) || $to_search eq '' || !defined($genetree_adaptor) || $genetree_adaptor eq ''){
		warn "missing parameter (to_search: $to_search, genetree_adaptor: $genetree_adaptor)\n";
		return undef;
	}
 	my $tree             = $genetree_adaptor->fetch_by_stable_id($to_search);
	if ( !defined($tree) || $tree eq '' ){        
		  #   $c->log->debug('Family::Tree::get_summary Could not find tree using stable_id') if $c->debug;
	    #try with root_id
	   # here we can only use numbers, no characters
	 	if($to_search =~ m/^\d+$/){
				print "searching root with $to_search\n";
			$tree = $genetree_adaptor->fetch_by_root_id($to_search);
		}
	}
	return (defined($tree)?$tree:undef);
}
sub search_gene_members_stable_id{
	my ($db,$to_search) = (@_);
	my $member_adaptor = $db->get_MemberAdaptor;
	my $member = $member_adaptor->fetch_by_source_stable_id( "ENSEMBLGENE", $to_search );
	return $member;
}
sub search_members_stable_id{
	my ($db,$to_search) = (@_);
	my $member_adaptor = $db->get_MemberAdaptor;
	my $member = $member_adaptor->fetch_by_source_stable_id( "ENSEMBLPEP", $to_search );
	return $member;
}
sub search_members_member_id{
	my ($db,$to_search) = (@_);
	my $member_adaptor = $db->get_MemberAdaptor;
    print "searching members with $to_search\n";	
	my $member = $member_adaptor->fetch_by_source_stable_id("undef", $to_search );
	return $member;
}
sub search_members_by_undef_id{
	my ($arg_ref) = @_;
	my $member_adaptor= $arg_ref->{'member_adaptor'};
	my $to_search= $arg_ref->{'to_search'};
    print "searching members with $to_search\n";	
	my $member = $member_adaptor->fetch_by_source_stable_id("undef", $to_search );
	return $member;
}
sub get_all_homologs_for_gene{
	my ($arg_ref) = @_;
	my $member= $arg_ref->{'member'};
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $homologies = $homology_adaptor->fetch_all_by_Member($member);
    return defined($homologies)?$homologies:undef
}


sub search_ext_ids{
	my ($file,$id) = (@_);
	my $search_cmd = "grep $id $file";
	my %hits4id;
	my ($further_search_id,$db_type);
	my @hits = `$search_cmd`;
    foreach (@hits){
		chomp;
		my ($hit_id,$db,$value) = split;
		$further_search_id = $hit_id;
		$hits4id{$hit_id}{$db} = $value;
		## determine accession type
		if($value eq $id){ $db_type = $db;}
		elsif($hit_id eq $id){	$db_type = "ENSEMBL_PRO";}
		else{$db_type = "FP"; $further_search_id = ''}
	}	
	return (\%hits4id,$further_search_id, $db_type);
}


#------------------------------------------------------------------------------------
#------------------------------ MATEUS ----------------------------------------------
#------------------------------------------------------------------------------------
sub get_member_by_xref
{
	my ($arg_ref)		= @_;
	my $member_adaptor	= $arg_ref->{'member_adaptor'};
	my $to_search		= $arg_ref->{'to_search'};
	my $limit 			= $arg_ref->{'limit'};
	my $extID2seq_sth_member;
	my $extID2seq_sth;

	$extID2seq_sth_member = $member_adaptor->prepare('select * from xrefID2Sequence where external_db_id = ? or external_db_id_name = ? limit ?');
	$extID2seq_sth_member->bind_param(1,$to_search);
	$extID2seq_sth_member->bind_param(2,$to_search);
	$extID2seq_sth_member->bind_param(3,$limit);
	$extID2seq_sth = $extID2seq_sth_member ;

	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
	#my $sequences = $extID2seq_sth->fetchall_arrayref();
	while ( my ($extID,$extName,$dbname,$gtID,$gtName,$description) = $extID2seq_sth->fetchrow_array())
	{
		print "$extID\t$extName\t$dbname\t$gtID\t$gtName\t$description\n";
	}
}

sub is_superTree
{
	my ($arg_ref)		= @_;
	my $member_adaptor	= $arg_ref->{'member_adaptor'};
	my $TF				= $arg_ref->{'TF'};
	my $supertree		= 0;

	my $superTree_sth_member = $member_adaptor->prepare('SELECT tree_type FROM gene_tree_root where stable_id = ?');
	$superTree_sth_member->bind_param(1,$TF);
	my $superTree_sth = $superTree_sth_member ;

	$superTree_sth->execute() or die "SQL Error: $DBI::errstr\n";
	#my $sequences = $extID2seq_sth->fetchall_arrayref();
	while ( my ($tree_type) = $superTree_sth->fetchrow_array())
	{
		if ($tree_type eq "supertree")
		{
			$supertree = 1;
		}
	}
	return $supertree;
}
=begin
	#my ($dbh,$db,$to_search, $column) = (@_);

	my $extID2seq_sth_member;
	my $extID2seq_sth_ext;
	my $extID2seq_sth;
		print "search using external: '$to_search'\n";
		$extID2seq_sth_ext = $dbh->prepare('select * from xrefID2Sequence where external_db_id= ? or external_db_id_name = ?');
		$extID2seq_sth_ext->bind_param(1,$to_search);
		$extID2seq_sth_ext->bind_param(2,$to_search);
		$extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
	while ( my ($extID,$extName,$dbname,$memberID,$gtID,$gtName,$description) = $extID2seq_sth->fetchrow_array() ){
		print "$extID,$extName,$dbname,$memberID,$gtID,$gtName,$description\n";
	   	my @line_array;
		my %line_hash = ("ExtID"=>$extID,"extName"=>$extName,"db"=>$dbname,"memberID"=>$memberID,"gtID"=>$gtID,"gtName"=>$gtName, "description" =>$description);	
		# try to see how long it takes to get member attributes as well
		#my $member_adaptor = $db->get_MemberAdaptor;
		#my $member = $member_adaptor->fetch_by_dbID($memberID);
		#return defined($member)? $member : undef;

	}
	return undef;
=cut
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
sub get_family_by_member_sequence{
	my ($arg_ref) = @_;
	my $to_search = $arg_ref->{'to_search'};
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	print "get family by member for $to_search\n";
	#print Dumper $member_adaptor;
	my $member = &check_valid_sequence({member_adaptor => $member_adaptor, to_search => $to_search});


	#my $member = &search_members_by_undef_id({member_adaptor => $member_adaptor, to_search => $to_search});
	#my $member = &search_members_member_id($genetree_adaptor, $to_search);
	#print Dumper $member;
	if(!$member){
		print "search by undef member found nothing\n";
		return undef;
	}
	else{
				print "search by undef member found something!!!! ".$member->taxon_id."\n";
				my $all_trees = $genetree_adaptor->fetch_all_by_Member($member);
				my $tree_id = (scalar(@{$all_trees}))? $all_trees->[0]->stable_id: "not in a family";
		# get gt object
			return $tree_id;
	}

	return undef;
}	
sub get_wikipedia4species{
	my ($arg_ref) = @_;
	my $db_adaptor = $arg_ref->{'db_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $extID2seq_sth_family;
	my $extID2seq_sth_ext;
	my $extID2seq_sth;
	print "search using external: '$to_search'\n";
	$extID2seq_sth_ext = $db_adaptor->prepare('select * from species2wikipedia where ncbi_id= ?');
	$extID2seq_sth_ext->bind_param(1,$to_search);
	$extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
	while ( my ($ncbi_id,$wikipedia_id,) = $extID2seq_sth->fetchrow_array() ){
		#my $genetree = $genetree_adaptor->fetch_by_dbID($gtID);
		return (defined($wikipedia_id)? $wikipedia_id:undef);
		#return $genetree;
	}	
	return undef;
}
sub get_family_by_xref_sequence{
	my ($arg_ref) = @_;
	my $db_adaptor = $arg_ref->{'db_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $column = $arg_ref->{'column'};
	my $extID2seq_sth_family;
	my $extID2seq_sth_ext;
	my $extID2seq_sth;
	print "search using external: '$to_search'\n";
	$extID2seq_sth_ext = $db_adaptor->prepare('select * from xrefID2Sequence where external_db_id= ? or external_db_id_name = ?');
	$extID2seq_sth_ext->bind_param(1,$to_search);
	$extID2seq_sth_ext->bind_param(2,$to_search);
	$extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
	while ( my ($extID,$extName,$dbname,$memberid,$gtID,$gtName,$description) = $extID2seq_sth->fetchrow_array() ){
		print "$extID,$extName,$dbname,$gtID,$gtName,$description\n";
	   	my @line_array;
		my %line_hash = ("ExtID"=>$extID,"extName"=>$extName,"db"=>$dbname,"gtID"=>$gtID,"gtName"=>$gtName, "description" =>$description);	
		# try to see how long it takes to get member attributes as well
		#my $genetree_adaptor = $db->get_GeneTreeAdaptor;
		#my $genetree = $genetree_adaptor->fetch_by_dbID($gtID);
		return (defined($gtName)? $gtName:undef);
		#return $genetree;
	}	
	return undef;
}
sub get_family_by_xref{
	my ($arg_ref) = @_;
	my $db_adaptor = $arg_ref->{'db_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $column = $arg_ref->{'column'};
	my $extID2seq_sth_family;
	my $extID2seq_sth_ext;
	my $extID2seq_sth;
	print "search using external: '$to_search'\n";
	$extID2seq_sth_ext = $db_adaptor->prepare('select * from xrefID2Family where external_db_id= ? or external_db_id_name = ?');
	$extID2seq_sth_ext->bind_param(1,$to_search);
	$extID2seq_sth_ext->bind_param(2,$to_search);
	$extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
	while ( my ($extID,$extName,$dbname,$gtID,$gtName,$description) = $extID2seq_sth->fetchrow_array() ){
		print "$extID,$extName,$dbname,$gtID,$gtName,$description\n";
	   	my @line_array;
		my %line_hash = ("ExtID"=>$extID,"extName"=>$extName,"db"=>$dbname,"gtID"=>$gtID,"gtName"=>$gtName, "description" =>$description);	
		# try to see how long it takes to get member attributes as well
		#my $genetree_adaptor = $db->get_GeneTreeAdaptor;
		#my $genetree = $genetree_adaptor->fetch_by_dbID($gtID);
		return (defined($gtName)? $gtName:undef);
		#return $genetree;
	}	
	return undef;
}
# Requires a genetree adaptor
sub check_valid_family{
		my ($arg_ref) = @_;
		my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
		my $member_adaptor = $arg_ref->{'member_adaptor'};
		my $to_search = $arg_ref->{'to_search'};
		my $tree;
		if($to_search =~ /^TF\d+$/){
			return $to_search;
		}
		
		if($to_search =~ /^\d+$/){
			$tree = &search_family_with_id({"genetree_adaptor" => $genetree_adaptor, "to_search" => $to_search});
		}
		if(defined($tree)){
			return $tree->stable_id;
		}
		else{
			# is it a member?
			print "\tsearching in member db with $to_search ...\n";
			$tree = &get_family_by_member_sequence({member_adaptor => $member_adaptor,  "genetree_adaptor" => $genetree_adaptor, "to_search" => $to_search});
			if($tree){
				print "Found something, we can stop here\n";
				return $tree;
			}	
			print "\tsearching in ext references with $to_search ...\n";
			$tree = &get_family_by_xref({"db_adaptor" => $genetree_adaptor, "to_search" => $to_search});
			if(defined($tree)){
				print "\tfound in xref table. redirect to /family/$tree\n";
				return $tree;
			}
			else{
				# last chance: we assume a sequence id was entered
				print "not found!!! now trying xref_sequence\n";
				$tree = &get_family_by_xref_sequence({"db_adaptor" => $genetree_adaptor, "to_search" => $to_search});
				if(defined($tree)){
					print "\tfound in xref sequence table. redirect to /family/$tree\n";
					return $tree;
				}
				return 0;
			}
		}
	return undef;
}

sub check_valid_sequence{
	my ($arg_ref) = @_;
	my $to_search = $arg_ref->{'to_search'};
	my $member_adaptor = $arg_ref->{'member_adaptor'};
   	my ($found_member, $prot_member);	
	my $member = search_members({"member_adaptor" => $member_adaptor, "to_search" => $to_search});
	if($member){
		print "done (found) extract information\n";
		$found_member = 1;
		if($member->get_canonical_Member){   # is a protein member
				$prot_member = $member->get_canonical_Member();
				print "switched to protein member\n";
		}
		return $member;
	}
	else{ 
		warn  " trying to search xrefs\n";
		my $extIDs = search_xref_hits({"member_adaptor" => $member_adaptor, "to_search" => $to_search, "type" => "external"});
		if(!scalar(@$extIDs)){
					return undef;
			}
			else{
				my $to_search = $extIDs->[0]{"memberId"};
				print "getting member object for $to_search\n";
				$member = $member_adaptor->fetch_by_dbID(  $to_search );
				#print Dumper $member;
				if(!$member){
					warn "Could not get member object\n";
				}
				else{
					return $member;
					$found_member = 1;
				}
			}
	}
	return ($found_member)?$member : undef;
}


sub get_MemberInformation4seq{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	print "\tAssumning ensemblID $to_search \n"; 
	my ($member,$sequence, $sequence_cds,$resultset,$found_member,$prot_member);
	$member = search_members({"member_adaptor" => $member_adaptor, "to_search" => $to_search});
	if($member){
		print "done (found) extract information\n";
		$found_member = 1;
		if($member->get_canonical_Member){   # is a protein member
				$prot_member = $member->get_canonical_Member();
				print "switched to protein member\n";
		}
		#print "get external references now, searching for ".$prot_member->dbID."\n";
		#my $extIDs = search_xref_hits({"member_adaptor" => $member_adaptor, "to_search" => $prot_member->dbID, "type" => "member"});
		#push(@{$resultset}, {"xref" => $extIDs});
		#$resultset->{"xref"} = $extIDs;
	}
	else{ 
		warn  " trying to search xrefs\n";
		my $extIDs = search_xref_hits({"member_adaptor" => $member_adaptor, "to_search" => $to_search, "type" => "external"});
		#push(@{$resultset}, {"xref" => $extIDs});
		$resultset->{"xrefs"} = $extIDs;
		#print Dumper $extIDs;
		if(!scalar(@$extIDs)){
					print "done (Could not find id $to_search in external ids)\n";
		}
		else{
				my $to_search = $extIDs->[0]{"memberId"};
				print "getting member object for $to_search\n";
				$member = $member_adaptor->fetch_by_dbID(  $to_search );
				#print Dumper $member;
				if(!$member){
					warn "Could not get member object\n";
				}
				else{
					$found_member = 1;
				}
		}
	}
	if($found_member){
		my @members;
		#push(@members,$member);
		#my $encoded_members = encode_members({"members" => \@members,"genetree_adaptor" => $genetree_adaptor,"limit" =>20});
	#push(@{$resultset}, {"member" => $encoded_members});
		#$resultset->{"member"} = $encoded_members;
		$resultset->{"member"}{"stable_id"} = $member->stable_id;
		$resultset->{"member"}{"taxon_name"} = ($member->taxon)->name;
		$resultset->{"member"}{"taxon_id"} = $member->taxon_id;
		$resultset->{"member"}{"description"} = defined($member->description)? $member->description: "NaN";
		my $all_trees = $genetree_adaptor->fetch_all_by_Member($member);
		$resultset->{"member"}{"family"} = (scalar(@{$all_trees}))? $all_trees->[0]->stable_id: "not in a family";
	
		#print "get homologs\n";
		#my $homologies = TreeFam::HomologyHelper::get_homologs_for_gene({"homology_adaptor" => $homology_adaptor, "genetree_adaptor" => $genetree_adaptor, "source_object" => $member, "type" =>  $type, "homology_type" => $homology_type});
		#if($homologies){
		#	#push(@{$resultset}, {"homologies" => $homologies});
		#	$resultset->{"homologies"} = $homologies;
		#}
		if($member->get_canonical_Member){   # is a protein member
				$prot_member = $member->get_canonical_Member();
				print "switched to protein member\n";
			$sequence = $prot_member->sequence();
			$sequence_cds = $prot_member->sequence_cds();
		}

	
	print "results\n";
	$sequence = ($sequence eq '') ? "NaN": $sequence;
	$sequence_cds = ($sequence_cds eq '')? "NaN": $sequence_cds;
	my %seq_hash = ("seq" => $sequence,"molecule"=>"protein","description"=> $member->description );
    my %seq_cds_hash = ("seq"=>$sequence_cds,"molecule"=>"cds","description"=> $member->description );
	my @seq_array = [\%seq_hash,\%seq_cds_hash];
	my %sequences_hash =("sequences" => @seq_array) ;
	#print Dumper %sequences_hash;
	#push(@{$resultset}, \%sequences_hash);
	$resultset->{"sequences"} = {"protein" => \%seq_hash, "cds"=> \%seq_cds_hash};
	
	}	

	return $resultset;
}



sub get_all_for_sequence_id{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $homology_adaptor = $arg_ref->{'homology_adaptor'};
	my $type = $arg_ref->{'type'};
	my $homology_type = $arg_ref->{'homology_type'};
	my $to_search = $arg_ref->{'to_search'};
	#print "\tAssumning ensemblID ... \n"; 
	my ($member,$sequence, $sequence_cds,$resultset,$found_member,$prot_member);
	$member = search_members({"member_adaptor" => $member_adaptor, "to_search" => $to_search});
	if($member){
		print "done (found) extract information\n";
		$found_member = 1;
		if($member->get_canonical_Member){   # is a protein member
				$prot_member = $member->get_canonical_Member();
				print "switched to protein member\n";
		}
		print "get external references now, searching for ".$prot_member->dbID."\n";
		my $extIDs = search_xref_hits({"member_adaptor" => $member_adaptor, "to_search" => $prot_member->dbID, "type" => "member"});
		#push(@{$resultset}, {"xref" => $extIDs});
		$resultset->{"xref"} = $extIDs;
	}
	else{ 
		warn  " trying to search xrefs\n";
		my $extIDs = search_xref_hits({"member_adaptor" => $member_adaptor, "to_search" => $to_search, "type" => "external"});
		#push(@{$resultset}, {"xref" => $extIDs});
		$resultset->{"xrefs"} = $extIDs;
		#print Dumper $extIDs;
		if(!scalar(@$extIDs)){
					print "done (Could not find id $to_search in external ids)\n";
			}
			else{
				my $to_search = $extIDs->[0]{"memberId"};
				print "getting member object for $to_search\n";
				$member = $member_adaptor->fetch_by_dbID(  $to_search );
				#print Dumper $member;
				if(!$member){
					warn "Could not get member object\n";
				}
				else{
					$found_member = 1;
				}
			}
	}
	if($found_member){
		my @members;
		push(@members,$member);
		my $encoded_members = encode_members({"members" => \@members,"genetree_adaptor" => $genetree_adaptor,"limit" =>20});
	#push(@{$resultset}, {"member" => $encoded_members});
		$resultset->{"member"} = $encoded_members;
		print "get homologs\n";
		my $homologies = TreeFam::HomologyHelper::get_homologs_for_gene({"homology_adaptor" => $homology_adaptor, "genetree_adaptor" => $genetree_adaptor, "source_object" => $member, "type" =>  $type, "homology_type" => $homology_type});
		if($homologies){
			#push(@{$resultset}, {"homologies" => $homologies});
			$resultset->{"homologies"} = $homologies;
		}
		if($member->get_canonical_Member){   # is a protein member
				$prot_member = $member->get_canonical_Member();
				print "switched to protein member\n";
			$sequence = $prot_member->sequence();
			$sequence_cds = $prot_member->sequence_cds();
		}

	}
	print "results\n";
	$sequence = ($sequence eq '') ? "NaN": $sequence;
	$sequence_cds = ($sequence_cds eq '')? "NaN": $sequence_cds;
	my %seq_hash = ("seq" => $sequence,"molecule"=>"protein","description"=> $member->description );
    my %seq_cds_hash = ("seq"=>$sequence_cds,"molecule"=>"cds","description"=> $member->description );
	my @seq_array = [\%seq_hash,\%seq_cds_hash];
	my %sequences_hash =("sequences" => @seq_array) ;
	#print Dumper %sequences_hash;
	#push(@{$resultset}, \%sequences_hash);
	$resultset->{"sequences"} = [\%seq_hash,\%seq_cds_hash];
	
	

	return $resultset;
}










sub search_xref_hits{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $column = $arg_ref->{'type'};
	my $extID2seq_sth_member;
	my $extID2seq_sth_ext;
	my $extID2seq_sth;
	print "type: $column\n";
    if($column eq "external"){
		print "search using external: '$to_search'\n";
		$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Sequence where external_db_id= ? or external_db_id_name = ?');
		$extID2seq_sth_ext->bind_param(1,$to_search);
		$extID2seq_sth_ext->bind_param(2,$to_search);
		$extID2seq_sth = $extID2seq_sth_ext ;
	}
	elsif($column eq "external_sequence"){
		print "search using external_sequence '$to_search'\n";
		$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Sequence where (external_db_id= ? or external_db_id_name = ?) and not db ="Pfam" and not db="HGNC" and not db="Wikigene"');
		$extID2seq_sth_ext->bind_param(1,$to_search);
		#$extID2seq_sth_ext->bind_param(2,$to_search);
		$extID2seq_sth = $extID2seq_sth_ext ;
	}elsif($column eq "autocomplete"){
		print "search using autocomplete: '$to_search'\n";
		$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Sequence where MATCH(external_db_id, external_db_id_name) AGAINST (?) and db = "HGNC"');
		$extID2seq_sth_ext->bind_param(1,$to_search);
		#$extID2seq_sth_ext->bind_param(2,$to_search);
		$extID2seq_sth = $extID2seq_sth_ext ;
	}elsif($column eq "member"){
		print "search using member $to_search\n";
		$extID2seq_sth_member = $member_adaptor->prepare('select * from xrefID2Sequence where member_id = ?');
		$extID2seq_sth_member->bind_param(1,$to_search);
		$extID2seq_sth = $extID2seq_sth_member ;
	}
	elsif($column eq "tree"){
		print "search using member $to_search\n";
		$extID2seq_sth_member = $member_adaptor->prepare('select * from xrefID2Sequence where gene_tree_stable_id = ?');
		$extID2seq_sth_member->bind_param(1,$to_search);
		$extID2seq_sth = $extID2seq_sth_member ;
	}
else{
		die "need to know which column to search in xref\n";
	}
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
	my $sequences = $extID2seq_sth->fetchall_arrayref();
	my $encoded_sequences = encode_xrefs({"sequences" => $sequences, "limit" => 100, "type" => "sequences"});
	return $encoded_sequences;
}
sub search_xref_families_hits{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $column = $arg_ref->{'type'};
	my $extID2seq_sth_member;
	my $extID2seq_sth_ext;
	my $extID2seq_sth;
	if($column eq "external"){
		print "search using external: '$to_search'\n";
		$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Family where external_db_id= ? or external_db_id_name = ?');
		$extID2seq_sth_ext->bind_param(1,$to_search);
		$extID2seq_sth_ext->bind_param(2,$to_search);
		$extID2seq_sth = $extID2seq_sth_ext ;
	}
	elsif($column eq "autocomplete"){
		print "search using autocomplete: '$to_search'\n";
		#$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Family where MATCH(external_db_id, external_db_id_name) AGAINST (? IN BOOLEAN MODE) and db = "hgnc"');
		my $t = "*".$to_search."*";
		print "search with $t\n";
		$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Family where MATCH(external_db_id, external_db_id_name) AGAINST (? IN BOOLEAN MODE) and db = "hgnc"');
		$extID2seq_sth_ext->bind_param(1,$t);
		$extID2seq_sth = $extID2seq_sth_ext ;
	}
	elsif($column eq "search_all"){
		print "search using search_all: '$to_search'\n";
		my $t = "*".$to_search."*";
		print "search with $t\n";
		$extID2seq_sth_ext = $member_adaptor->prepare('select * from xrefID2Family where MATCH(external_db_id, external_db_id_name, description) AGAINST (? IN BOOLEAN MODE) and db = "hgnc"');
		$extID2seq_sth_ext->bind_param(1,$t);
		$extID2seq_sth = $extID2seq_sth_ext ;
	}elsif($column eq "tree"){
		print "search using member $to_search\n";
		$extID2seq_sth_member = $member_adaptor->prepare('select * from xrefID2Family where gene_tree_stable_id = ?');
		$extID2seq_sth_member->bind_param(1,$to_search);
		$extID2seq_sth = $extID2seq_sth_member ;
	}else{
		die "need to know which column to search in xref\n";
	}
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
	my $sequences = $extID2seq_sth->fetchall_arrayref();
	my $encoded_sequences;
	if($column eq "autocomplete" || $column eq 'search_all'){
		$encoded_sequences = encode_xrefs({"sequences" => $sequences, "limit" => 100, "type" => "autocomplete"});
	}
	else{
		$encoded_sequences = encode_xrefs({"sequences" => $sequences, "limit" => 1000, "type" => "normal"});
	}
	return $encoded_sequences;
}


# expects a arraref

sub encode_xrefs{
	my ($arg_ref) = @_;
	my $sequences = $arg_ref->{'sequences'};
	my $limit = $arg_ref->{'limit'};
	my $type = $arg_ref->{'type'};
	my @results_array;
	my @sequences_array;
	my $counter = 0;
	my %redundancy_hash;
	my %sequences_remember_hash;
	foreach my $seq_array(@{$sequences}){	
		my %pair_hash;
		if($type eq "normal"){
			my ($extID,$extName,$dbname,$memberID,$gtID,$gtName,$description) = @{$seq_array};
			$pair_hash{"extID"} = defined($extID)? $extID : "NaN";
			$pair_hash{"extName"} = defined($extName)? $extName: "NaN";
			$pair_hash{"dbName"} = defined($dbname)? $dbname : "NaN";
			$pair_hash{"gtName"} = defined($gtName)? $gtName : "NaN";
			$pair_hash{"description"} = defined($description)? $description : "NaN";
			my $key =$pair_hash{"extID"}."".$pair_hash{"extName"}."".$pair_hash{"dbName"}."".$pair_hash{"gtName"};
			next if exists($redundancy_hash{$key});
			next if ($extID eq "NaN" );
			$redundancy_hash{ $key} = 1;
		}
		elsif($type eq "sequences"){
			my ($extID,$extName,$dbname,$memberID,$gtID,$gtName,$description) = @{$seq_array};
			$pair_hash{"extID"} = defined($extID)? $extID : "NaN";
			$pair_hash{"extName"} = defined($extName)? $extName: "NaN";
			$pair_hash{"dbName"} = defined($dbname)? $dbname : "NaN";
			$pair_hash{"gtName"} = defined($gtName)? $gtName : "NaN";
			$pair_hash{"memberId"} = defined($memberID)? $memberID : "NaN";
			$pair_hash{"description"} = defined($description)? $description : "NaN";
			my $key =$pair_hash{"extID"}."".$pair_hash{"extName"}."".$pair_hash{"dbName"}."".$pair_hash{"gtName"};
			next if exists($redundancy_hash{$key});
			next if ($extID eq "NaN" );
			$redundancy_hash{ $key} = 1;
		}else{
			my ($extID,$extName,$dbname,$gtID,$gtName,$description) = @{$seq_array};
			$pair_hash{"name"} = defined($description)? $description : "NaN";
			$pair_hash{"id"} = defined($extName)? $extName: "NaN";
			$pair_hash{"symbol"} = defined($extID)? $extID : "NaN";
			$pair_hash{"gtName"} = defined($gtName)? $gtName : "NaN";
			my $key =$pair_hash{"name"}."".$pair_hash{"id"}."".$pair_hash{"symbol"};
			next if exists($redundancy_hash{$key});
			next if ($extID eq "NaN" );
			$redundancy_hash{ $key} = 1;
		}
		push(@results_array, \%pair_hash);		
		last if(defined($limit) && $counter++ > $limit);
	}
# Remember we saved the sequences
  return \@results_array;
}

sub search_description{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $member_sth_member;
	my $member_sth_ext;
	my $member_sth;
	print "search using members '$to_search'\n";
	my $members = $member_adaptor->generic_fetch("MATCH(description) AGAINST (\'+".$to_search."\' IN BOOLEAN MODE)" );
    my @hits;
	print "\tdone, found ".scalar(@$members)." hits in members. doing some formatting now\n";
	my $encoded_members = encode_members({"members" => $members,"genetree_adaptor" => $genetree_adaptor,"limit" =>20});
	# can then be converted to json
	return $encoded_members;
}

sub encode_members{
	my ($arg_ref) = @_;
	my $members = $arg_ref->{'members'};
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $limit = $arg_ref->{'limit'};
	my @results_array;
	my @sequences_array;
	my $counter = 0;
	my %sequences_remember_hash;
	foreach my $member (@{$members}){
		my %pair_hash;
				$pair_hash{"member_id"} = $member->member_id;
				$pair_hash{"stable_id"} = $member->stable_id;
				$pair_hash{"version"} = $member->version;
				$pair_hash{"source_name"} = $member->source_name;
				$pair_hash{"taxon_id"} = $member->taxon_id;
				$pair_hash{"taxon_name"} = ($member->taxon)->name;
				$pair_hash{"description"} = defined($member->description)? $member->description: "NaN";
				$pair_hash{"display_label"} = defined($member->display_label)? $member->display_label : "NaN";
				my $all_trees = $genetree_adaptor->fetch_all_by_Member($member);
				$pair_hash{"family"} = (scalar(@{$all_trees}))? $all_trees->[0]->stable_id: "not in a family";
				push(@results_array, \%pair_hash);		
				last if(defined($limit) && $counter++ > $limit);
		}
  return \@results_array;
}
sub search_members{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	print "getting member $to_search\n";
	my $member = $member_adaptor->fetch_by_source_stable_id( undef, $to_search );
	return (defined $member)? $member: undef ;
}
sub search_members_memberId{
	my ($arg_ref) = @_;
	my $member_adaptor = $arg_ref->{'member_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $member = $member_adaptor->fetch_by_dbID( "ENSEMBLPEP", $to_search );
	return $member;
}

####### searching for gene trees using TreeFam id (e.g. TF101001)
sub search_trees{
	my ($arg_ref) = @_;
	my $genetree_adaptor = $arg_ref->{'genetree_adaptor'};
	my $to_search = $arg_ref->{'to_search'};
	my $tree = $genetree_adaptor->fetch_by_stable_id($to_search);
	my @found_tree;  
	my $id = $tree->stable_id;
	my $gene_count = ($tree->get_value_for_tag('gene_count'))? $tree->get_value_for_tag('gene_count'): 'NaN';
	my $aln_percent_identity = ($tree->get_value_for_tag('aln_percent_identity'))? $tree->get_value_for_tag('aln_percent_identity'): 'NaN';
	my $aln_length = ($tree->get_value_for_tag('aln_length'))? $tree->get_value_for_tag('aln_length'): 'NaN';
	@found_tree = [$id,$gene_count,$aln_percent_identity, $aln_length];
	return \@found_tree;
}
sub search_families{
	my ($arg_ref) = @_;
	my $genetree_adaptor= $arg_ref->{'genetree_adaptor'};
	my $member = $arg_ref->{'member'};
	my $all_trees = $genetree_adaptor->fetch_all_by_Member($member, "default");
	return $all_trees;
}
1;
