#!/usr/bin/env python

import os, sys, json, argparse, shutil, re, subprocess, datetime
from astropy.io import fits


class FITScollection:
	def __init__(self):
		self.objectList = []
			
	def additem(self, filename):
		fitsObject = {}
		fitsObject['filename'] = filename
		fitsObject['processed'] = False
		self.objectList.append(fitsObject)
	
	def getFilenames(self):
		return [f['filename'] for f in self.objectList]
		
	def sort(self):
		self.objectList = sorted(self.objectList, key=lambda object: object['filename'], reverse = False)
		
	def writeToCSV(self, filename):
		outfile = open(filename, 'wt')
		for o in self.objectList:
			outfile.write("%s, %s\n"%(o['filename'], o['processed']))
		outfile.close()

	def loadFromCSV(self, filename):
		infile = open(filename, 'rt')
		
		for l in infile:
			fitsObject = {}
			l = l.strip()
			items = l.split(',')
			filename = items[0].strip()
			if items[1].strip() == 'True':
				status = True
			else:
				status = False
			fitsObject['filename'] = filename
			fitsObject['processed'] = status
			self.objectList.append(fitsObject)
		infile.close()
		
	def updateStatus(self, newStatus):
		for o in self.objectList:
			filename = o['filename']
			found = False
			changed = False
			oldStatus = o['processed']
			for c in newStatus.objectList:
				if c['filename'] == filename:
					found = True
					if oldStatus != c['processed']:
						changed = True
						o['processed'] = c['processed']
						print "changed status of", filename, "to", o['processed']

				
	def __str__(self):
		return "%d objects in the list."%len(self.objectList)
	
class FITSdata:
	def __init__(self):
		self.sources = []
		
	def appendFromFile(self, filename):
		hdulist = fits.open(filename)
		header = hdulist[0].header
		# print(repr(header))
		tableData = hdulist[1].data
		#print tableData
		cols = hdulist[1].columns
		added = 0
		for d in tableData:
			rowObject = {}
			for c in cols.names:
				rowObject[c] = d[c]
			# print rowObject
			added+= 1
			self.sources.append(rowObject)
		return (added, len(self.sources))
		
	def sort(self):
		self.sources = sorted(self.sources, key=lambda object: object['mean'], reverse = True)
		
	def writeToFile(self, filename):
		objects = self.sources
		hdu = fits.PrimaryHDU()
		cols = []
		cols.append(fits.Column(name='id', format='16A', array = [o['id'] for o in objects]))
		cols.append(fits.Column(name='ra', format='E', array = [o['ra'] for o in objects]))
		cols.append(fits.Column(name='dec', format = 'E', array = [o['dec'] for o in objects]))
		cols.append(fits.Column(name='xmax', format = 'E', array = [o['xmax'] for o in objects]))
		cols.append(fits.Column(name='ymax', format = 'E', array = [o['ymax'] for o in objects]))
		cols.append(fits.Column(name='mean', format = 'E', array = [o['mean'] for o in objects]))
		cols.append(fits.Column(name='peak', format = 'E', array = [o['peak'] for o in objects]))
		cols.append(fits.Column(name='variance', format = 'E', array = [o['variance'] for o in objects]))
		cols.append(fits.Column(name='type', format = '8A', array = [o['type'] for o in objects]))
		cols.append(fits.Column(name='CCD', format = '4A', array = [o['CCD'] for o in objects]))
		cols = fits.ColDefs(cols)
		tbhdu = fits.BinTableHDU.from_columns(cols)
			
		prihdr = fits.Header()
		prihdr['COMMENT'] = "Created by Hagrid (mergeOutput) on %s."%( datetime.datetime.ctime(datetime.datetime.now()))
			
		prihdu = fits.PrimaryHDU(header=prihdr)
		thdulist = fits.HDUList([prihdu, tbhdu])
		thdulist.writeto(filename, clobber=True)
	


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='Merge the output files produced by Hagrid into one big FITS file.')
	parser.add_argument('output', type=str, help = "The output filename. It will be a FITS file.")
	parser.add_argument('-p', '--path', type=str, help='The folder in which the IPHAS images are contained.')
	parser.add_argument('-w', '--workingpath', type=str, help='A root folder for the temporary and output files.')
	parser.add_argument('-s', '--sky', action='store_true', help='Look for sky pointings rather than bright sources.')
	arg = parser.parse_args()
	print arg

	if arg.path is None:
		dataPath="."
	else: 
		dataPath= arg.path
	if arg.workingpath is None:
		workingPath="."
	else: 
		workingPath= arg.workingpath
		
	
	# Get a list of files in the folder
	# First, check if the source data is there
	if not os.path.exists(dataPath):
		print "The folder for the source data %s could not be found. Exiting."%dataPath
		sys.exit()
	
	searchString = "r[0-9]{6}-[1-4].sources.fits"
	if arg.sky:
		searchString = "r[0-9]{6}-[1-4].sky.fits"
	search_re = re.compile(searchString)
		 
	(_, _, filenames) = os.walk(dataPath).next()
	allObjects = FITScollection()
	print "looking for matches to:", searchString
	
	for file in filenames:
		m = search_re.match(file)
		if (m): 
			allObjects.additem(file)

	print allObjects
	allObjects.sort()
	
	
	dataObject = FITSdata()

	for index, f in enumerate(allObjects.objectList):
		numAdded, total = dataObject.appendFromFile(f['filename'])
		print "Added %d new rows from data in %s. Total: %d"%(numAdded, f['filename'], total)
	
	dataObject.sort()
	
	print "Writing %d rows to %s."%(len(dataObject.sources), arg.output)
	dataObject.writeToFile(arg.output)
		

		
	
