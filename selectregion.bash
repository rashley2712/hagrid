#!/bin/bash
if test "$#" -lt 3; then
    echo "Usage: selectregion.bash [inputfile] [ra] [dec] optional: [radius] [output]"
	echo output will default to a file called out.fits
    exit -1
fi
if test "$#" -gt 3; then
    radius=$4
    echo using radius of $4
  else
    radius=1
fi
if test "$#" -gt 4; then
    output=$5
else
    output=out.fits
fi
ra=$2
dec=$3
rd=12
echo Selecting: RA: $ra DEC: $dec Radius: $radius
echo Writing to: $output

# echo ‘select skyDistanceDegrees(ra, dec, $ra, $dec)< $rd’ > cmds.list

stilts -verbose tpipe in=$1 out=$output cmd="select 'skyDistanceDegrees(ra, dec, $ra, $dec)< $4' "


