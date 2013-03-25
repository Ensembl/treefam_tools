use strict;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


# Connection parameters to the TreeFam db
Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'mysql-treefam',
    -user => 'treefam_user',
    -pass => '',
    -port => 4401,
    -species => 'TreeFam',
    -dbname => 'treefam_production_9_69',
    #-dbname => 'treefam_homology_67hmm',
);


1;
