#/bin/bash
for ht in 0pt1 0pt3 0pt4 
do
	echo 'Running:'$ht
	##Ejemplo del nombre del fichero: holland_250_h0pt1.ini
	mpirun -np 2 python tcrm.py -v -c example/holland_250_h$ht.ini
done
