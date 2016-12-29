#!/bin/bash
if test "$#" -ne 2; then
    echo "Usage: selectfilter.bash [inputfile] [filter]"
    exit -1
fi

stilts -verbose tpipe in=$1 cmd="select 'equals(band, \"$2\")'" out=out.fits
