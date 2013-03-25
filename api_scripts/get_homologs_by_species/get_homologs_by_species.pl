use strict;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Bio::EnsEMBL::Registry;

my $numArgs     = $#ARGV + 1;

if ($numArgs != 1)
{
    print "Use: perl tool_box_get_species_homologs.pl species_list_file\n";
    exit();
}
my $species_list_file= $ARGV[0];

#Temporary directory.
#Used to unzip.
#--------------------
my $tmp_dir = "./tmp";

#Output file.
#------------
my $OUT_FILE = "$tmp_dir/homology_list.txt";

#Get Registry.
#-------------
my $registry = 'Bio::EnsEMBL::Registry';
Bio::EnsEMBL::Registry->load_all("../../registry/production_treefam_reg_conf.pl");

#Get Adaptor.
#------------
my $genomedb_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GenomeDB' );
if(!$genomedb_adaptor)
{
	die "Could not load genomedb_adaptor\n";
}

#Lists used to store the pairs.
#------------------------------
my @species_list_vec;
my %species_list;

#Open file list.
#---------------
open(IN,$species_list_file);
while(<IN>)
{
	chop($_);
#	$_ ~= s/\s+/_/;:q
	push(@species_list_vec,$_);
}
close(IN);

#Get the pairs and sort them in order to get the correct entries in th DB.
#-------------------------------------------------------------------------
my $pair = 0;
for (my $i=0; $i < (scalar(@species_list_vec)-1); $i++)
{
	for (my $j=$i+1; $j < scalar(@species_list_vec); $j++)
	{
		my @tmp_vec;
		push(@tmp_vec,$species_list_vec[$i]);
		push(@tmp_vec,$species_list_vec[$j]);
		@tmp_vec = sort(@tmp_vec);
		@{$species_list{$pair}} = @tmp_vec;
		$pair++;
	}
}

#open Output file.
#-----------------
open(OUT,">$OUT_FILE");
foreach my $pair (sort num keys %species_list)
{
	my ($species_A,$species_B) = @{$species_list{$pair}};

	#DB query.
	#---------	
	print "$pair:\tspecies_A=$species_A\tspecies_B=$species_B\n";
	my $pairwise = $genomedb_adaptor->prepare('SELECT file FROM pairwise_homology WHERE species_A=? and species_B=?');
	$pairwise->bind_param(1,$species_A);
	$pairwise->bind_param(2,$species_B);

	$pairwise->execute() or die "SQL Error: $DBI::errstr\n";

	#Get all the gzip files stored as BLOBs, uncompress and concatenate.
	#-------------------------------------------------------------------
	my $gzip_file = "$tmp_dir/$species_A-$species_B.gz";
	my $buffer_gzip;
	my $buffer_txt;
	while ( my ($BLOB) = $pairwise->fetchrow_array())
	{
		open(ZIP,">$gzip_file");
		print ZIP $BLOB;
		close(ZIP);
	}
	gunzip $gzip_file => \$buffer_txt;
	print OUT $buffer_txt;
	unlink($gzip_file);
}
close(OUT);

sub num { $a <=> $b }
