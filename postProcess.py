#!/usr/bin/env python

import os, sys, argparse, shutil, re, subprocess
import sourceClasses
import astropy, numpy
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot 


def calcDistance(s1, s2):
	separation = s1['coord'].separation(s2['coord'])
	# print separation.arcminute
	return separation.arcminute
	
def gaussian(x, width):
	return 1.0 * numpy.exp((-4.0 * numpy.log(2) * x * x)/(width * width))

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='Perform some post-processing of a merged catalogue of sources.')
	parser.add_argument('catalogue', type=str, help='The name of the FITS file containing the input catalogue.')
	parser.add_argument('--limit', type=int, help='Limit the number of sources to process. For debugging purposes only.')
	parser.add_argument('-o', '--output', type=str, help='Output filename (.fits) is automatically added. Default is "out.fits".')
	parser.add_argument('-p', '--parameter', type=str, help='Which column to use a ranking parameter. Sensible options are "peak" or "mean". Default is "mean".')

	parser.add_argument('-dr', '--debugrow', type=int, help='Debug the contributors to a certain source matching debugrow.')

	arg = parser.parse_args()
	
	if arg.limit is None:
		limit = False
		limitNumber = 0
	else:
		limit = True
		limitNumber = arg.limit
		print "Warning: Limiting number of sources processed to %d."%limitNumber
	
	if arg.debugrow is None:
		debug = False
		debugRow = 0
	else:
		debug = True
		debugRow = arg.debugRow
		print "Warning: Will debug source number %d."%debugRow
	
	if arg.output is None:
		outputFilename = "out.fits"
	else:
		outputFilename = arg.output + ".fits"
	
	if arg.parameter is None:
		rankColumn = "mean"
	else:
		rankColumn = arg.parameter
	
	# First, check if the source data is there
	if not os.path.exists(arg.catalogue):
		print "File for source data %s could not be found. Exiting."%arg.catalogue
		sys.exit()
	
	sources = sourceClasses.sourcesClass()

	new, total = sources.addSourcesFromFITS(arg.catalogue)
	
	print "Astropy version: %s"%astropy.__version__
	
	print "Loaded %d sources."%new
	print "Columns are: %s"%sources.getColumnNames()
	
	sources.sortSources(rankColumn)
	# sources.addSkyCoords()
	
	metaSources = []
	gaussianRadius = 0.2  			# Gaussian radius in arcmin
	ras = [s['ra'] for s in sources.sources]
	decs = [s['dec'] for s in sources.sources]
	# allSources = [s['coord'] for s in sources.sources]
	allSources = SkyCoord(ras*u.degree, decs*u.degree)
	
	for sourceNumber, chosenSource in enumerate(sources.sources):
		print sourceNumber, "Original %s: %f"%(rankColumn, chosenSource[rankColumn]), 
		sourceCoords = SkyCoord(ra = chosenSource['ra']* u.degree, dec = chosenSource['dec']*u.degree)
		separationFromChosenSource = sourceCoords.separation(allSources)
		matchArray = []
		for index, s in enumerate(separationFromChosenSource):
			source = {'n': index, 'separation': s.arcminute, 'coord':s}
			matchArray.append(source)
		sortedMatches = sorted(matchArray, key=lambda object: object['separation'], reverse = False)
		for index in range(1, len(sortedMatches)):
			s = sortedMatches[index]
			contribution =  gaussian(s['separation'], gaussianRadius) * sources.sources[s['n']][rankColumn]
			if debug and sourceNumber==debugRow:
				 print index, s['separation'], s['n'], sources.sources[s['n']][rankColumn], gaussian(s['separation'], gaussianRadius), contribution
			if s['separation'] > 5*gaussianRadius: break
			chosenSource[rankColumn]+= contribution
		print "Final %s: %f"%(rankColumn, chosenSource[rankColumn])
		if limit and sourceNumber>=limitNumber-1: break
	
	sources.sortSources(rankColumn)
	
	minimumSeparation = 1.0
	firstSource = sources.sources[0]	
	finalSources = [firstSource]
	
	finalRAs = [firstSource['ra']]
	finalDECs = [firstSource['dec']]
	finalCoords = SkyCoord(finalRAs*u.degree, finalDECs*u.degree)
	print "Applying minimum separation criterion: %f arcminutes"%minimumSeparation 
	print
	for index, s in enumerate(sources.sources):
		sourceCoords = SkyCoord(ra = s['ra']* u.degree, dec = s['dec']*u.degree)
		separations = sourceCoords.separation(finalCoords).arcminute
		# if type(separations) == float: separations = [separations]
		# print index, separations
		add = True
		for sep in separations:
			if sep < minimumSeparation: add = False
		if add:
			finalSources.append(s)
			finalRAs.append(s['ra'])
			finalDECs.append(s['dec'])
			finalCoords = SkyCoord(finalRAs*u.degree, finalDECs*u.degree)
		if (index%10) == 0: 
			sys.stdout.write("\rChecking: %d of %d."%(index, len(sources.sources)))
			sys.stdout.flush()
		if limit and index >=limitNumber-1: break
	sys.stdout.write("\rChecking: %d of %d.\n"%(index, len(sources.sources)))
	sys.stdout.flush()
	
		# if index>4: break
		
	print "Final list of sources is %d long"%len(finalSources)
	outputSources = sourceClasses.sourcesClass()
	outputSources.sources = finalSources
	outputSources.writeToFile(outputFilename)
	
	""" # For sanity sake, plot the Gaussian
	gaussianPlot = matplotlib.pyplot.figure("Gaussian function", figsize=(10, 10))
	gaussianPlot.set_tight_layout(True)
	xValues = numpy.arange(-3, 3, 0.1)
	yValues = gaussian(xValues, gaussianRadius)
	matplotlib.pyplot.plot( xValues, yValues, color = 'r')
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show()
	"""		 
		
	
	sys.exit()
	
	print "Searching around:", topSourceCoords
	for index in range(1, len(sources.sources)):
		distance = calcDistance(topSource, sources.sources[index])
		if distance < 5:
			print "Close match: %f [%d]"%(distance, index)
	print 
	print topSourceCoords.separation(allSources)
	idxc, idxcatalog, d2d, d3d = allSources.search_around_sky(topSourceCoords, 1*u.deg)
	# idx, d2d, d3d = topSourceCoords.match_to_catalog_sky(allSources) 
	print idxc
		
