#!/bin/bash
if test "$#" -ne 2; then
    echo "Usage: getfilenames.bash [inputfile] [pathtofile]"
    exit -1
fi

stilts tpipe in=$1 cmd="keepcols url" out=tmp.csv

tail -n +2 tmp.csv > tmp2.csv

grep -o 'r[0-9][0-9][0-9]/r[0-9]*-[1-4].fits.fz' tmp2.csv > tmp.csv

rm tmp2.csv 

sed -e "s/^/$2\//" tmp.csv 

rm tmp.csv

