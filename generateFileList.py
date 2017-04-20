#!/usr/bin/env python

import os, sys, json, argparse, shutil, re, subprocess
from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.vo.client.conesearch import list_catalogs
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u

class locationsClass:
	def __init__(self, filename = None):
		self.locationList = []
		self.filename = filename
		defaultLocation = { 'name': 'cygnus_ha', 
								'criteria': [ 	{'type' : 'coord',  'ra': 304.1963, 'dec': 37.1692, 'radius': 0.25} , 
					      	  					{'type' : 'filter', 'band': 'halpha' } 
					    					] }
		self.appendLocation(defaultLocation)
		if filename is not None:
			self.loadFromFile(filename)
		
	def locationList(self):
		return [ l['name'] for l in self.locationList ]
		
	def loadFromFile(self, filename):
		print "Loading locations from file:", filename
		jsonData = open(filename).read()
		jsonObject = json.loads(jsonData)
		print jsonObject
		
	def appendLocation(self, location):
		self.locationList.append(location)
		
	def getLocationByName(self, name):
		for l in self.locationList:
			print l['name']
			if l['name'] == name.lower(): return l
		return None
		
	def saveLocations(self):
		print json.dumps(self.locationList)
		


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A simple Python tool build a file list of IPHAS images for input into hagrid (batch mode).')
	parser.add_argument('-a', '--archive', type=str, help='The root folder in which the IPHAS images are contained (ie the IPHAS image archive). ')
	parser.add_argument('-db', '--database', type=str, help='The IPHAS database file (usually "iphas-images.fits.gz", stored in the archive folder).')
	parser.add_argument('-w', '--workingpath', type=str, help='A root folder for the temporary and output files.')
	parser.add_argument('-o', '--output', type=str, help='The output file. Otherwise results will do to stdout and to a file called <criterion name>.list.')
	arg = parser.parse_args()
	#print arg
	
	# Could get the db file from: http://www.iphas.org/data/images/iphas-images.fits.gz
	

	if arg.archive is None: archivePath="."
	else: archivePath = arg.archive
	
	if arg.database is None: 
		dbFilename = os.path.join(archivePath, "iphas-images.fits.gz")
	else:
		dbFilename = arg.database
	
	if arg.workingpath is None: workingPath="."
	else: workingPath = arg.workingpath

	name = 'Cygnus'
	locations = locationsClass(filename = 'temp.json')
	location = locations.getLocationByName(name)
	if location is None:
		print "No location to match."
		sys.exit()
	
	try:
		IPHASdb = Table.read(dbFilename)
	except IOError as e:
		print "Could not find the IPHAS database file: %s"%dbFilename
		sys.exit()
	
	
	
	for c in location['criteria']:
		if c['type'] == 'coord':
			print "Applying a coord criterion..."
			center = SkyCoord(ra = c['ra'] * u.degree, dec = c['dec'] * u.degree)
			print "Looking for matches to position %f, %f with a radius of %2.2f degrees."%(c['ra'], c['dec'], c['radius'])
			catalog_RAs = IPHASdb['ra']
			catalog_DECs = IPHASdb['dec']
			catalog = SkyCoord(ra = catalog_RAs * u.degree, dec = catalog_DECs * u.degree )
	
			results = center.separation(catalog).degree
			matches = []
			for idx,r in enumerate(results):
				if r < c['radius'] : matches.append(idx)
			print "%d images have centres within %2.2f degrees of %s."%(len(matches), c['radius'], center.to_string('hmsdms'))
			IPHASdb = IPHASdb[matches]
			
		if c['type'] =='filter':
			print "Applying a filter criterion..."
			filters = IPHASdb['band']
			bands = c['band']
			if isinstance(bands, list): 
				print "More than one filter"
			else:
				bands = [c['band']]
			
			matches = []
			for b in bands:
				print "Checking for matches in the %s band."%b
				for idx, f in enumerate(filters):
					if b == f: matches.append(idx)
			
			print "In total, %s images match the filter criterion."%(len(matches)) 
			IPHASdb = IPHASdb[matches]
			
	
	
	# Now create the filename list
	fileSubFolders = [ u.split('/')[-2] + '/' + u.split('/')[-1] for u in IPHASdb['url'] ]
	
	if arg.output is None:
		for f in fileSubFolders:
			sys.stdout.write("%s\n"%(os.path.join(archivePath, f)))
		sys.stdout.flush()
		outFilename = os.path.join(workingPath, location['name'] + '.list')
	else:
		outFilename = os.path.join(workingPath, arg.output)
		
	print "Writing results to: %s"%outFilename
	outputFile = open(outFilename, 'wt')
	for f in fileSubFolders:
		outputFile.write(os.path.join(archivePath, f) + '\n')
	outputFile.close()
	
	
	sys.exit()
