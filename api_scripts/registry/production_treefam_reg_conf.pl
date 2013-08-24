
# Release Coordinator, please update this file before starting every release
# and check the changes back into CVS for everyone's benefit.

# Things that normally need updating are:
#
# 1. Release number
# 2. Check the name prefix of all databases
# 3. Possibly add entries for core databases that are still on genebuilders' servers

use strict;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


   
# Extra core databases that live on genebuilders' servers:

# Compara databases used during the release (master, previous and current compara dbs)


# Master
Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'mysql-treefam-public.ebi.ac.uk',
    -user => 'treefam_ro',
    -pass => '',
    -port => 4418,
    -species => 'TreeFam',
    -dbname => 'treefam_production_9_69',
);
1;
