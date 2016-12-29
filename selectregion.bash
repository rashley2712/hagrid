#!/bin/bash
if test "$#" -ne 3; then
    echo "Usage: selectfilter.bash [inputfile] [ra] [dec] optional: [radius] [output]"
    exit -1
fi

stilts -verbose tpipe in=$1 out=out.fits cmd="select 'skyDistanceDegrees(ra, dec, $2, $3)<1' "


