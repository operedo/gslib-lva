#!/bin/sh
### Directivas para el gestor de colas (modificar los valores NAMEOFJOB y USERNAME, y mantener la opción "-S")
# Cambiar el nombre del trabajo
#$ -N sgs-lva
# Especificar un shell
#$ -S /bin/sh
# Enviame un correo cuando empiece el trabajo y cuando acabe...
#$ -m be
# ... a esta dirección de correo
#$ -M operedo@ac.upc.edu



source /Soft/intelf/composer_xe_2013.3.163/bin/compilervars.sh ia32

TEST=45_0_0_10_1 
cd $TEST
#cp /homeA/o/operedo/gslib-lva/sgs-lva/bin/SGS_LVA /homeA/o/operedo/gslib-lva/sgs-lva/test/$TEST/
#cp /homeA/o/operedo/gslib-lva/Boost_dijkstra/Boost_dijkstra /homeA/o/operedo/gslib-lva/sgs-lva/test/$TEST/
#
#rm -rf /scratch/nas/4/operedo/$TEST
#cp -R /homeA/o/operedo/gslib-lva/sgs-lva/test/$TEST/ /scratch/nas/4/operedo/

#cd /scratch/nas/4/operedo/$TEST
./SGS_LVA sgs_lva.par > output.txt 2> error.txt
RES="`perl absDiff.pl sgs.out sgs.out_ok | wc -l`"
if [ $RES -eq 0 ]; then
	echo "[sgs-lva] $TEST: PASSED"
else
	echo "[sgs-lva] $TEST: NOT PASSED. Compare output file ${TEST}/sgs.out with baseline file ${TEST}/sgs.out_ok"
fi
cd -	
