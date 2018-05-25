#!/bin/bash
for name in `ls -f /work/daphne/24genomes2017/5-sweeps/sweed/GrpT/filled/*.txt`
do
    filename=`basename $name`
    python /work/daphne/24genomes2017/5-sweeps/filter_windows/filter_windows.py $name /work/daphne/24genomes2017/BcinB0510_woFASTA.sorted.gff /work/daphne/24genomes2017/5-sweeps/sweed/GrpT/filled_wf/
done