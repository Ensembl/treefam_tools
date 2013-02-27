#
#===============================================================================
#
#         FILE: ReleaseMaster.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Fabian Schreiber (), fs9@sanger.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 02/13/2013 10:39:09 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
 
package TreeFam::Release::ReleaseMaster;
use Moose::Role;
use Bio::EnsEMBL::Registry; 

  has 'registry' => (isa => 'string', is => 'ro', clearer => 'clear_x' );
  has 'registry_file' => (isa => 'Int', is => 'rw', clearer => 'clear_y');

sub new {
   my $obj = shift;
   my $class = ref($obj)?ref($obj) : $obj;
   my $args = shift;

    my $self = bless {}, $class;
   $self->{-registry_file} = $args->{registry_file} if( $args && ref($args) eq 'HASH' && $args->{registry_file});
   $self->{-registry} = $args->{registry} if( $args && ref($args) eq 'HASH' && $args->{registry});
    #print "in releasemaster constructor: ".$args->{registry}." and ".$args->{registry_file}."\n";
   return $self;
}

sub get_db_connection{
	my ($self) = @_;
	#my $registry = $arg_ref->{registry};
	#my $registry_file = $arg_ref->{registry_file};
    my $reg = $self->{-registry}->load_all($self->{-registry_file});
    #my $reg = Bio::EnsEMBL::Registry->load_all($registry_file);
    return defined($reg)? $reg : undef;
}

1;
