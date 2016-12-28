#!/bin/bash
if test "$#" -ne 2; then
    echo "Usage: selectregion.bash [ra] [dec] optional: [radius] [output]"
    exit -1
fi
echo $1 $2

echo stilts tpipe in=iphas-images.fits.gz out=out.fits cmd=" select 'skyDistanceDegrees(ra, dec, $1, $2)<1' "

stilts tpipe in=iphas-images.fits.gz out=out.fits cmd="select 'skyDistanceDegrees(ra, dec, $1, $2)<1' "


