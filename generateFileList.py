#!/usr/bin/env python

import os, sys, json, argparse, shutil, re, subprocess
from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.vo.client.conesearch import list_catalogs
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy

class locationsClass:
	def __init__(self, filename = None):
		self.locationList = []
		self.filename = filename
		if filename is not None:
			self.loadFromFileCSV(filename)
		else:
			defaultFilename = 'locations.csv'
			defaultLocations = 	[ 	['cygnus',  304.1963, 37.1692, 0.25] , 
									['perseus', 39.6075, 62.3381, 1.0] ]
			self.appendLocations(defaultLocations)
			self.saveToFileCSV(defaultFilename)
		
	def locationNames(self):
		return [ l['name'] for l in self.locationList ]
		
	def loadFromFileCSV(self, filename):
		print "Loading locations from file:", filename
		inputFile = open(filename, 'rt')
		for f in inputFile:
			items = f.strip().split(',')
			if items[0][0] == '#': continue
			try:
                                if "x" in items[3]:
                                    sqsize=items[3].split("x")
                                    location = { 'name': items[0], 'ra': float(items[1]), 'dec': float(items[2]), 'ra_size': float(sqsize[0]), 'dec_size': float(sqsize[1])}
                                else:
				    location = { 'name': items[0], 'ra': float(items[1]), 'dec': float(items[2]), 'radius': float(items[3])}
				self.locationList.append(location)
			except (ValueError, IndexError):
				print "Could not interpret the line:", f
		inputFile.close()
		
		
	def appendLocations(self, locations):
		for l in locations:
			location = {'name': l[0], 'ra': l[1], 'dec': l[2], 'radius': l[3]}
			self.locationList.append(location)
					
	def getLocationByName(self, name):
		for l in self.locationList:
			if l['name'] == name.lower(): return l
		return None
			
	def saveToFileCSV(self, filename):
		outfile = open(filename, 'wt')
		outfile.write("# name, ra, dec, radius\n")
		for l in self.locationList:
			outfile.write("%s, %f, %f, %f\n"%(l['name'], l['ra'], l['dec'], l['radius']))
		outfile.close()	
		


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A simple Python tool build a file list of IPHAS images for input into hagrid (batch mode).')
	parser.add_argument('location', type=str, help='A string specifying the location you want.')
	parser.add_argument('-a', '--archive', type=str, help='The root folder in which the IPHAS images are contained (ie the IPHAS image archive). ')
	parser.add_argument('-db', '--database', type=str, help='The IPHAS database file (usually "iphas-images.fits.gz", stored in the archive folder or installed with the hagrid app).')
	parser.add_argument('-w', '--workingpath', type=str, help='A root folder for the temporary and output files.')
	parser.add_argument('-o', '--output', type=str, help='The output file. Otherwise results will dump to stdout and to a file called <criterion name>.list.')
	parser.add_argument('-locationlist', type=str, help='Specify an alternative list of locations. By default this is ''location.csv'' and found in the install directory of ''hagrid''.')
	parser.add_argument('--ra', type=float, help='The RA in degrees')
	arg = parser.parse_args()
	#print arg
	installPath = os.path.realpath(__file__).rsplit('/',1)[0]
	
	if arg.archive is None: archivePath="."
	else: archivePath = arg.archive
	
	if arg.database is None: 
		print "The program file is installed at:", installPath
		dbFilename = os.path.join(installPath, "iphas-images.fits.gz")
		# Could get the db file from: http://www.iphas.org/data/images/iphas-images.fits.gz
	else:
		dbFilename = arg.database
	
	if arg.workingpath is None: workingPath="."
	else: workingPath = arg.workingpath

	name = arg.location
	
	if os.path.exists(os.path.join(installPath, 'locations.csv')):
		locations = locationsClass(filename = os.path.join(installPath, 'locations.csv'))
	else:
		locations = locationsClass()
	location = locations.getLocationByName(name)
	if location is None:
		print "No location to match."
		sys.exit()
	
	try:
		IPHASdb = Table.read(dbFilename) 
		print "Loaded %d rows from %s."%(len(IPHASdb), dbFilename)
	except IOError as e:
		print "Could not find the IPHAS database file: %s"%dbFilename
		sys.exit()
	
	print "Applying a coord criterion... for %s"%location['name']
	center = SkyCoord(ra = location['ra'] * u.degree, dec = location['dec'] * u.degree)
	catalog_RAs = IPHASdb['ra']
	catalog_DECs = IPHASdb['dec']
	catalog = SkyCoord(ra = catalog_RAs * u.degree, dec = catalog_DECs * u.degree )
	
	matches = []
        if 'radius' in location:
	    print "Looking for matches to position %f, %f with a radius of %2.2f degrees."%(location['ra'], location['dec'], location['radius'])
	    results = center.separation(catalog).degree
	    for idx,r in enumerate(results):
                if r < location['radius'] : matches.append(idx)
	    print "%d images have centres within %2.2f degrees of %s."%(len(matches), location['radius'], center.to_string('hmsdms'))
        else:
	    print "Looking for matches to position %f, %f with a box of %2.2fx%2.2f degrees."%(location['ra'], location['dec'], location['ra_size'], location['dec_size'])
	    dra,ddec = center.spherical_offsets_to(catalog)
	    for idx,r in enumerate(dra):
                if r.degree < location['ra_size']/2. and ddec[idx].degree < location['dec_size']/2.: matches.append(idx)
	    print "%d images have centres within a box of %2.2fx%2.2f degrees of %s."%(len(matches), location['ra_size'], location['dec_size'], center.to_string('hmsdms'))
	IPHASdb = IPHASdb[matches]
			
	print "Applying a filter criterion..."
	filters = IPHASdb['band']
	band = 'halpha'
	matches = []
	for idx, f in enumerate(filters):
            grade=IPHASdb["qcgrade"][idx]
            if grade.startswith("C") or grade.startswith("D"):
              if IPHASdb["in_dr2"][idx]=="false": continue
	    if band == f: matches.append(idx)
			
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
