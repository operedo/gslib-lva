#!/bin/bash

source /Soft/intelf/composer_xe_2013.3.163/bin/compilervars.sh ia32

TEST=45_0_0_10_1 

cp ../bin/SGS_LVA $TEST/
cp ../../Boost_dijkstra/Boost_dijkstra $TEST/
cd $TEST
./SGS_LVA sgs_lva.par > output.txt 2> error.txt
RES="`perl absDiff.pl sgs.out sgs.out_ok | wc -l`"
if [ $RES -eq 0 ]; then
	echo "[sgs-lva] $TEST: PASSED"
else
	echo "[sgs-lva] $TEST: NOT PASSED. Compare output file ${TEST}/sgs.out with baseline file ${TEST}/sgs.out_ok"
fi
cd ..	
