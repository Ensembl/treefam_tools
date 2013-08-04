
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


# The majority of core databases live on two staging servers:
#Bio::EnsEMBL::Registry->load_registry_from_url('mysql://ensro@ens-livemirror/67');
#Bio::EnsEMBL::Registry->load_registry_from_url('mysql://anonymous@mysql.ebi.ac.uk:4157');

#Bio::EnsEMBL::Registry->load_registry_from_url(
  #'mysql://ensro@ens-staging1/68');

  #Bio::EnsEMBL::Registry->load_registry_from_url(
    #'mysql://ensro@ens-staging2/68');

#Bio::EnsEMBL::Registry->load_registry_from_db(
#    -host => 'mysql.ebi.ac.uk',
#        -port => 4157,
#            -user => 'anonymous',
#                -verbose => 1 );

#Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
    #{
            #-host => 'mysql.ebi.ac.uk',
                #-port => 4157,
                    #-user => 'anonymous',
                        #-verbose => 1,
            #-db_version => 68, 
   #},
    #{
            #-host => 'ens-livemirror',
            #-user => 'ensro',
            #-port => '3306',
            #-verbose => 1
    #}
    #);
    
# Extra core databases that live on genebuilders' servers:

# Compara databases used during the release (master, previous and current compara dbs)


# Master
Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'mysql-treefam-public',
    -user => 'treefam_ro',
    -pass => '',
    -port => 4418,
    -species => 'TreeFam',
    -dbname => 'treefam_production_9_69',
    #-dbname => 'treefam_homology_67hmm',
);
#Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    #-host => 'web-mei-treefam',
    #-user => 'treefam_admin',
    #-pass => 'treefam_king1982',
    #-port => 3365,
    #-species => 'compara_master',
    #-dbname => 'treefam_master9',
#);
# Prev release
#Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    #-host => 'mysql.ebi.ac.uk',
    #-user => 'anonymous',
    #-port => 4157,
    #-species => 'compara_prev',
    #-dbname => 'ensembl_compara_68',
#);

1;
