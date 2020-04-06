#!/bin/bash
module remove openmpi/4.0.1
module load openmpi/3.1.2

#lscpu
#mpirun -n $procs ./spmv_random "input/"$input"_random_"$procs".mtx" $nodes
of="../GraphConversionAndPartition/partition/2D_partition/output/"
mpirun -n $procs ./spmv_random $of"delaunay_n24_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"hugetrace-00020_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"inf-europe_osm_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"inf-italy_osm_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"soc-livejourna_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"soc-BlogCatalog_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"soc-flickr_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"oc-flixster_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"soc-lastfm_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"soc-youtube_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"great-britain_osm_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"hugetrace-00010_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"it-2004_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"netherlands_osm_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random $of"rel9_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random $of"relat9_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random $of"sk-2005_random_"$procs".mtx" $nodes


mpirun -n $procs ./spmv_random "input/cnr_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/eu_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/in_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/NACA0015_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/road_usa_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/road_central_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/nlpkkt120_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/nlpkkt160_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/nlpkkt200_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/nlpkkt240_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/cage15_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/G3_circuit_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/kkt_power_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/ecology1_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/uk2002_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/uk2007_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/asia_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/NLR_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/coPapersDBLP_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/citationCiteseer_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/coPapersCiteseer_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/copter2_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/333SP_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/af_shell10_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/AS365_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/europe_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/M6_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/thermal2_random_"$procs".mtx" $nodes
