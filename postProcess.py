#!/usr/bin/env python

import os, sys, argparse, shutil, re, subprocess
import sourceClasses
from astropy.coordinates import SkyCoord
from astropy import units as u


def calcDistance(s1, s2):
	separation = s1['coord'].separation(s2['coord'])
	# print separation.arcminute
	return separation.arcminute

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='Perform some post-processing of a merged catalogue of sources.')
	parser.add_argument('catalogue', type=str, help='The name of the FITS file containing the input catalogue.')

	arg = parser.parse_args()
	print arg

	# First, check if the source data is there
	if not os.path.exists(arg.catalogue):
		print "File for source data %s could not be found. Exiting."%arg.catalogue
		sys.exit()
	
	sources = sourceClasses.sourcesClass()

	new, total = sources.addSourcesFromFITS(arg.catalogue)
	
	
	
	print "Loaded %d sources."%new
	print "Columns are: %s"%sources.getColumnNames()
	
	sources.sortSources("mean")
	sources.addSkyCoords()
	
	metaSources = []
	gaussianRadius = 0.5   			# Gaussian radius in arcmin
	ras = [s['ra'] for s in sources.sources]
	decs = [s['dec'] for s in sources.sources]
	# allSources = [s['coord'] for s in sources.sources]
	allSources = SkyCoord(ras*u.degree, decs*u.degree)
	print len(allSources)
	topSource = sources.sources[0]
	topSourceCoords = SkyCoord(ra = topSource['ra']* u.degree, dec = topSource['dec']*u.degree)
	print "Searching around:", topSourceCoords
	for index in range(1, len(sources.sources)):
		distance = calcDistance(topSource, sources.sources[index])
		if distance < 5:
			print "Close match: %f [%d]"%(distance, index)
	print 
	idxc, idxcatalog, d2d, d3d = allSources.search_around_sky(topSourceCoords, 1*u.deg)
	# idx, d2d, d3d = topSourceCoords.match_to_catalog_sky(allSources) 
	print idxc
		