#!/bin/bash
if test "$#" -ne 2; then
    echo "Usage: getfilenames.bash [inputfile] [pathtofile]"
    exit -1
fi

stilts tpipe in=$1 cmd="keepcols run" out=tmp.csv

tail -n +2 tmp.csv > tmp2.csv

for f in $(cat tmp2.csv); do
	subfolder=r${f:0:3}
	echo \/$subfolder\/r$f.fits.fz >> tmp3.csv
done

sed -e "s/^/$2/" tmp3.csv > tmp.csv

cat tmp.csv

rm tmp2.csv 
rm tmp3.csv

# sed 's/$/.fits.fz/' tmp.csv

