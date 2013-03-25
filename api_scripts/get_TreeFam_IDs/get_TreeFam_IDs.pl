use Data::Dumper;
#use JSON;
use strict;
use warnings;
#use TreeFam::HomologyHelper;
use TreeFam::SearchHelper;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';
Bio::EnsEMBL::Registry->load_all("../../registry/production_treefam_reg_conf.pl");

my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
if(!$member_adaptor){
	die "Could not load member_adaptor\n";
}

my $IDs_file	= "input_IDs.txt";
my $limit		= 200;
my @IDs_vec;

open(IDs,$IDs_file) || die "Could not open file:input_IDs.txt";
while(<IDs>)
{
	chomp($_);
	push (@IDs_vec,"$_");
}
close(IDs);

foreach my $id (@IDs_vec)
{
	my $xref_families_hits = TreeFam::SearchHelper::get_member_by_xref({"member_adaptor" => $member_adaptor, "to_search" => $id, "type" => 'external_db_id', "limit" => $limit});
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
