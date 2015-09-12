#!/bin/bash

TEST=45_0_0_10_1 

cd $TEST
KT3D_LVA kt3d.par > output.txt 2> error.txt
RES="`perl absDiff.pl kriging.out kriging.out_ok | wc -l`"
if [ $RES -eq 0 ]; then
	echo "[kt3d-lva] $TEST: PASSED"
else
	echo "[kt3d-lva] $TEST: NOT PASSED. Compare output file ${TEST}/kriging.out with baseline file ${TEST}/kriging.out_ok"
fi
cd ..	
