#!/bin/bash

TEST=45_0_0_10_1 

cd $TEST
GAMV_LVA gamv.par > output.txt 2> error.txt
RES="`diff gamv.out gamv.out_ok | wc -l`"
if [ $RES -eq 0 ]; then
	echo "[gamv-lva] $TEST: PASSED"
else
	echo "[gamv-lva] $TEST: NOT PASSED. Compare output file ${TEST}/gamv.out with baseline file ${TEST}/gamv.out_ok"
fi
cd ..	
