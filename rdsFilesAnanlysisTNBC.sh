#!/bin/bash

for f in *.rds; 
do
    mkdir "${f%.rds}"
    cp "${f%.rds}".rds "${f%.rds}"
    cp *.py "${f%.rds}"
    cp *.R "${f%.rds}"
    cd "${f%.rds}"
        echo "${f%.rds}".rds > fileName.txt
        filename="${f%.rds}".rds
	suffix="${filename##*[0-9]}"
	number="${filename%"$suffix"}"
	number="${number##*[!-0-9]}"
	echo $number > number.txt
	echo "["$number"]" > numberFile.txt

	Rscript PxToCoordTNBC.R
	Rscript IHC_tumor_total_PolygonsTissue.R
	python CoordsLists.py
	Rscript RanalysisWholeTissue.R
	python PolygonsLists.py
	python PtoRRandom.py
	Rscript RanalysisWholeTissueRand.R
    cd ..
	

    
done
