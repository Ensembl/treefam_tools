#
#===============================================================================
#
#         FILE: OtherFunctions.pm
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
package TreeFam::OtherFunctions;

use strict;
use warnings;
 




sub write_to_file{
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	### OPENING FILE
	open my $out, '>>', $file_name or die "Couldn't open '$file_name': $!";
	### Writing file
	print {$out} $text;
	### CLOSING FILE
	close $out or die "Couldn't close '$file_name': $!";
	return 1;
}

1;

