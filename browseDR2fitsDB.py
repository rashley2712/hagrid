#!/usr/bin/env python

import argparse, sys
import datetime, time
import numpy
import ppgplot
import generalUtils, configHelper, os
import curses
from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.table import Table, vstack
from astropy.utils import data
import astropy.table	



if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Loads an IFAS DR2 catalogue FITS file and inspects the contents.')
	parser.add_argument('filename', type=str, help='The DR2 FITS file.')
	parser.add_argument('--save', action="store_true", help='Write the input parameters to the config file as default values.')
	parser.add_argument('--ignorecache', action="store_true", help='Ignore the cached DR2 catalogue. Will overwrite one if it already exists.')
	args = parser.parse_args()
	print args
	
	config = configHelper.configClass("browseDRfitsDB")
	
	if args.save:
		config.save()
	
	hdulist = fits.open(args.filename)
	
	print hdulist.info()
	
	metadata = hdulist[0].data
	
	for index, card in enumerate(hdulist):
		print "Card number:%d"%index
		print card.header.keys()
		print repr(card.header)
		
	dataKeys = hdulist[1].header.keys()
	
	# Get a list of the data columns
	columns = []
	columnIndex = []
	for key in dataKeys:
		if "TTYPE" in key:
			index = int(key[5:])
			label = hdulist[1].header[key]
			columns.append(label)
			columnIndex.append(index)
	
	# For each column, get the corresponding datatype and column description
	columnDescriptions = []
	columnFormats = []
	for i, c in enumerate(columns):
		index = columnIndex[i]
		commentKey = "TCOMM%d"%index
		formatKey = "TFORM%d"%index
		print index, commentKey, formatKey, c
		foundDescription = False
		for key in dataKeys:
			if commentKey == key:
				columnDescriptions.append(hdulist[1].header[key])
				foundDescription = True
			if formatKey == key:
				columnFormats.append(hdulist[1].header[key])
		if not foundDescription:
			columnDescriptions.append(c)
		
	print "Table description:"			
	for i, index in enumerate(columnIndex):
		print "[%d] \t%s \t%s \t(%s)"%(index, columns[i], columnDescriptions[i], columnFormats[i])
	
	dataTable = hdulist[1].data
	
	table = []
	counter = 0;
	totalExpectedRows = len(dataTable)
	for data in dataTable:
		counter+=1
		tableRow = {}
		for i, c in enumerate(columns):
			tableRow[c] = data[i]
		if (counter%10)==0:
			percent = (float(counter)+1.)/float(totalExpectedRows) * 100.0
			sys.stdout.write("\r%d rows parsed of a total of %d or %2.2f%%"%(counter, totalExpectedRows, percent))
			sys.stdout.flush()
		# if counter>1000: break
		table.append(tableRow)
	
	sys.stdout.write("\n")
	sys.stdout.flush()
	
	print "Parsed %d rows into dict object."%(counter+1)