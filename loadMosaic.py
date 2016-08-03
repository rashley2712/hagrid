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

	
def redraw():
	ppgplot.pggray(boostedImage, xlimits[0], xlimits[1]-1, ylimits[0], ylimits[1]-1, imageMinMax[0], imageMinMax[1], imagePlot['pgPlotTransform'])
	
		
def checkMargins(margins):
	if margins[0][0] < margins[1][0]:
		temp = margins[0][0]
		margins[0][0] = margins[1][0]
		margins[1][0] = temp
	if margins[0][1] < margins[1][1]:
		temp = margins[0][1]
		margins[0][1] = margins[1][1]
		margins[1][1] = temp
	return margins
		
def withinMargins(table, namedColumns):
	print namedColumns
	return True
	

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Loads an IFAS mosaic image. Displays it with PGPLOT.')
	parser.add_argument('filename', type=str, help='The FITS image file.')
	parser.add_argument('--save', action="store_true", help='Write the input parameters to the config file as default values.')
	parser.add_argument('--ignorecache', action="store_true", help='Ignore the cached DR2 catalogue. Will overwrite one if it already exists.')
	args = parser.parse_args()
	print args
	
	config = configHelper.configClass("inspectIPHASImage")

	"""exposureTime = config.getProperty("ExposureTime")
	if args.exposureTime!=None:
		exposureTime = args.exposureTime
		config.ExposureTime = exposureTime
	if exposureTime == None:
		print "Please specify an exposure time or save one in the config file."
		sys.exit()
	"""
	
	paperSize = 6  # Paper size in inches
	imageMinMax = (0, 255)
	plotSources = True
	invertedColours = True
	
	if args.save:
		config.save()
	
	hdulist = fits.open(args.filename)
	
	print hdulist.info()
	
	
	for card in hdulist:
		print card.header.keys()
		print repr(card.header)
	
	imageData =  hdulist[0].data
	
	wcsSolution = WCS(hdulist[0].header)
	
	hdulist.close()
	
	(height, width) = numpy.shape(imageData)
	
	aspectRatio = float(height)/float(width)
	print aspectRatio
	
	""" Set up the PGPLOT windows """
	imagePlot = {}
	imagePlot['pgplotHandle'] = ppgplot.pgopen('/xs')
	ppgplot.pgpap(paperSize, aspectRatio)
	ppgplot.pgsvp(0.0, 1.0, 0.0, 1.0)
	ppgplot.pgswin(0, width, 0, height)
	
	# ppgplot.pgenv(0., width,0., height, 1, -2)
	imagePlot['pgPlotTransform'] = [0, 1, 0, 0, 0, 1]
	
	boostedImage = generalUtils.percentiles(imageData, 20, 99)
	ppgplot.pggray(boostedImage, 0, width-1, 0, height-1, 0, 255, imagePlot['pgPlotTransform'])
	
	# Determine the RA, DEC of the centre of the image, using the WCS solution found in the FITS header
	imageCentre = [ width/2, height/2]
	
	
	ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
	
	
	positionString = generalUtils.toSexagesimal((ra, dec))
	print "RA, DEC of image centre is: ", positionString, ra, dec
	margins = wcsSolution.all_pix2world([[0, 0], [width, height]], 1)
	margins = checkMargins(margins)
	print "ra, dec limits:", margins
	
	xlimits = (0, width)
	ylimits = (0, height)
		
		
	# try: 
	x=width/2
	y=height/2
	newWidth = width
	newHeight = height
	ppgplot.pgsci(3)
	keyPressed = None
	while keyPressed != 'q':
		ch = ppgplot.pgcurs(x, y)
		x=ch[0]
		y=ch[1]
		keyPressed = ch[2]
		# print "Key pressed:", ch[2]
			
		if keyPressed== 'w':
			imageMinMax = (imageMinMax[1], imageMinMax[0])
			if invertedColours: invertedColours= False
			else: invertedColours=True
			redraw()
			
		if keyPressed=='i':
			print "Zoom requested at (%0.0f, %0.0f)"%(x, y)
			zoomFactor = 1.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			if (newWidth<100) or (newHeight<100):
				print "Maximum zoom reached."
				continue
			else:
				xlimits = (int(x - newWidth/2), int(x + newWidth/2)) 
				if xlimits[0] < 0:
					xlimits = (0, xlimits[1] + abs(xlimits[0]))
				if xlimits[1] > width:
					xlimits = (width - newWidth, width)
				ylimits = (int(y - newHeight/2), int(y + newHeight/2)) 
				if ylimits[0] < 0:
					ylimits = (0, ylimits[1] + abs(ylimits[0]))
				if ylimits[1] > height:
					ylimits = (height - newHeight, height)
				
			xlimits = (int(xlimits[0]), int(xlimits[1]))
			ylimits = (int(ylimits[0]), int(ylimits[1]))
			print "new limits:", xlimits, ylimits
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			print "new limits (world)", ra_limits, dec_limits
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]
			redraw()
			
		if keyPressed=='o':
			if newWidth == width: continue
			print "Zoom out requested at (%0.0f, %0.0f)"%(x, y)
			zoomFactor = 0.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			
			if (newWidth >= width) or (newHeight>=height):
				print "Back to 1:1 scale."
				newWidth = width
				newHeight = height
				
			xlimits = (int(x - newWidth/2), int(x + newWidth/2)) 
			if xlimits[0] < 0:
				xlimits = (0, xlimits[1] + abs(xlimits[0]))
			if xlimits[1] > width:
				xlimits = (width - newWidth, width)
			ylimits = (int(y - newHeight/2), int(y + newHeight/2)) 
			if ylimits[0] < 0:
				ylimits = (0, ylimits[1] + abs(ylimits[0]))
			if ylimits[1] > height:
				ylimits = (height - newHeight, height)
				
			xlimits = (int(xlimits[0]), int(xlimits[1]))
			ylimits = (int(ylimits[0]), int(ylimits[1]))
			print "new limits:", xlimits, ylimits
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]	
			print "new limits (world)", ra_limits, dec_limits
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			redraw()
		
		if plotSources: sourceStatus = "ON"
		else: sourceStatus = "OFF" 
		if invertedColours: invertStatus = "ON"
		else: invertStatus = "OFF" 
		print "PlotSources[%s] Invert Greyscale[%s]"%(sourceStatus, invertStatus)
			
			
	# except KeyboardInterrupt:
	#	print "Ctrl-C pressed, but I dealt with it. "


	ppgplot.pgclos()
	
	
	
