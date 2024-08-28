#!/bin/bash

start=$1
end=$2

for cn in $(seq ${start} 1 ${end}); do
	Rscript findAllMarkers_louvain_parallel.R wilcox ${cn}
        #job2pbs.py -l "Rscript findAllMarkers_louvain_parallel.R negbinom ${cn}"
done
