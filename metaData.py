#!/usr/bin/env python

import argparse, sys, os, re
import datetime, time
import configHelper, numpy
from astropy.io import fits
from astropy.wcs import WCS

class fitsObject:
	def __init__(self):
		self.ra = 0;
		self.dec = 0;
		self.filename = ""
		self.path = ""
		self.filter = None
		self.CCD = None
		
	def setCentre(self, ra, dec):
		self.ra = ra
		self.dec = dec

def joinPaths(path1, path2):
	if path1[-1] == '/': return path1 + path2
	return path1 + '/' + path2

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Crawls through the directories and gathers file and FITS header meta-data.')
	parser.add_argument('--datapath', type=str, help='The root path for the IPHAS images.')
	parser.add_argument('--save', action="store_true", help='Write the input parameters to the config file as default values.')
	args = parser.parse_args()
	print(args)

	if args.save:
		config.save()
	
	config = configHelper.configClass("metaData")
	configDefaults  = {
		"RootPath": "/Users/rashley/astro/IPHAS",
		"SearchPath": ".*(r[0-9]{3})"
	}
	config.setDefaults(configDefaults)
	rootPath = config.assertProperty("RootPath", args.datapath)
	searchPath = config.assertProperty("SearchPath", None)
	search_re = re.compile(searchPath)
	
	fitsObjects = []
	# Find all folders in data path
	folders = os.walk(rootPath)
	IPHASFolders = []
	for f in folders:
		m = search_re.match(f[0])
		if (m): 
			IPHASFolders.append(f)
			print(m.group(0))
	
	outCSV = open("data.csv", 'wt')
	outCSV.write("filename, CCD, filter, ra, dec, cached\n")
	for f in IPHASFolders:
		print("Looking at the files in %s."%(f[0]))
		fitsFiles = []
		for file in f[2]:
			if "fits.fz" in file:
				sys.stdout.write("\r%s    "%file)
				sys.stdout.flush()
		
				newfitsObject = fitsObject()
				fitsFiles.append(file)
				hdulist = fits.open(joinPaths(f[0], file))
				newfitsObject.filter = hdulist[1].header['WFFBAND']
				newfitsObject.CCD = hdulist[1].header['DASCHAN']
				wcsSolution = WCS(hdulist[1].header)
				(height, width) = numpy.shape(hdulist[1].data)
				imageCentre = [ width/2, height/2]
				ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
				newfitsObject.setCentre(ra, dec)
				newfitsObject.path = f[0]
				newfitsObject.filename = file
				cacheFilename = file[:9] + "_dr2_cache.fits"
				if os.path.exists(joinPaths(f[0], cacheFilename)):
					sys.stdout.write(" found cached dr2 data in:" + str(cacheFilename))
					sys.stdout.flush()
					newfitsObject.cached = True
				else:
					newfitsObject.cached = False
					
				hdulist.close()
				fitsObjects.append(newfitsObject)
				outCSV.write("%s, %s, %s, %f, %f, %s\n"%(newfitsObject.filename, newfitsObject.CCD, newfitsObject.filter, newfitsObject.ra, newfitsObject.dec, newfitsObject.cached))
				outCSV.flush()
				
		sys.stdout.write("\n")
		sys.stdout.flush()
		print("%d files in this folder"%len(fitsFiles))

	
	outCSV.close()
	sys.exit()
	
	
	
	for card in hdulist:
		print(card.header.keys())
		print(repr(card.header))
		for key in card.header.keys():
			if 'WFFBAND' in key:
				filter = card.header[key]
	
	imageData =  hdulist[1].data

	
	wcsSolution = WCS(hdulist[1].header)
	
	hdulist.close()
	
	(height, width) = numpy.shape(imageData)
	
	aspectRatio = float(height)/float(width)
	print(aspectRatio)
	
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
	print("RA, DEC of image centre is: ", positionString, ra, dec)
	margins = wcsSolution.all_pix2world([[0, 0], [width, height]], 1)
	margins = checkMargins(margins)
	print("ra, dec limits:", margins)
	
	
	filenameParts = args.filename.split('.')
	dr2Filename = filenameParts[0] + "_dr2_cache.fits"
	cached = False
	if not args.ignorecache:
		print("Looking for a cached copy of the DR2 catalogue:", dr2Filename)
		if os.path.exists(dr2Filename):
			cached = True



	if not cached:
		for index, yCentre in enumerate(numpy.arange(height/4, height, height/4)):
			print(yCentre)
			imageCentre = [ width/2, yCentre]
			ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
			print(index, imageCentre, ra, dec)
			with data.conf.set_temp('remote_timeout', 60):
				try: 
					search = conesearch(center=(ra, dec),
	                    radius=0.1,
	                    verb=3,
	                    catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
				except: 
					print("Failed to retrieve any results from Vizier.")
			dr2nearbyTemp = search.to_table()
			if index==0: 
				dr2nearby = dr2nearbyTemp
			else:
				print("Table was %d rows. Additional data is %d rows."%(len(dr2nearby), len(dr2nearbyTemp)))
				dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
				dr2nearby = dr2nearbyWhole
				print("New table is %d rows."%len(dr2nearby))
			print(dr2nearby)
		
		radius = 0.05	
		# Get bottom left corner from vizier
		ra, dec = wcsSolution.all_pix2world([[0, 0]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		# Get bottom right corner from vizier
		ra, dec = wcsSolution.all_pix2world([[width, 0]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		# Get top left corner from vizier
		ra, dec = wcsSolution.all_pix2world([[0, height]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		# Get top right corner from vizier
		ra, dec = wcsSolution.all_pix2world([[width, height]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		print("Checking for duplicate objects.")
		IPHASnames = []
		for index, row in enumerate(dr2nearby):
			IPHASname = row['IPHAS2']
			if IPHASname in IPHASnames:
				dr2nearby.remove_row(index)
			else:
				IPHASnames.append(IPHASname)
			if (index%100) == 0:
				sys.stdout.write("\rChecking row: %d    "%index)
				sys.stdout.flush()
		
		sys.stdout.write("\n")
		sys.stdout.flush()
		
		print("Trimming out objects that lie outside the limits of this image.")
		margins = checkMargins(margins)
		print("Margins: " + str(margins))
	
		for index, row in enumerate(dr2nearby):
			ra = row['RAJ2000']
			dec = row['DEJ2000']
			if ra>margins[0][0] or ra<margins[1][0] or dec>margins[0][1] or dec<margins[1][1]:
				sys.stdout.write("\rRejecting this row: %d %f, %f is outside image borders."%(index, ra, dec))
				sys.stdout.flush()
				dr2nearby.remove_row(index)
		
		sys.stdout.write("\n")
		sys.stdout.flush()
				
		dr2nearby.write(dr2Filename, format='fits', overwrite=True)
	
	else:
		dr2nearby = Table.read(dr2Filename)
	
	print("Data columns found in the DR2 catalogue: " + str(dr2nearby.colnames))
	
	dr2nearby.remove_column('errBits2')
	
	print("Length of the dr2 table: %d"%len(dr2nearby))
	
	# Move table into a dictionary object
	dr2Objects = []
	for row in dr2nearby:
		IPHAS2name = row['IPHAS2']
		dr2Object={}
		dr2Object['name'] = IPHAS2name
		dr2Object['ra'] = row['RAJ2000']
		dr2Object['dec'] = row['DEJ2000']
		dr2Object['class'] = row['mergedClass']
		dr2Object['pStar'] = row['pStar']
		dr2Object['iClass'] = row['iClass']
		dr2Object['haClass'] = row['haClass']
		x, y = wcsSolution.all_world2pix([dr2Object['ra']], [dr2Object['dec']], 1)
		dr2Object['x'] = x[0]
		dr2Object['y'] = y[0]
		dr2Objects.append(dr2Object)
		
	
	
	# Run through all objects and find the ones that are haClass = "+1" but iClass != "+1" and overall class = "+1"
	extendedHaSources = []
	for index, d in enumerate(dr2Objects):
		if (d['class'] ==1) and (d['haClass'] == 1) and (d['iClass'] != 1):
			extendedHaSources.append(index)
	
	print("Found %d extended Ha sources out of %d total objects in DR2."%(len(extendedHaSources), len(dr2Objects)))
	
	
	xlimits = (0, width)
	ylimits = (0, height)
	
	plotCircles(dr2Objects, margins)

	help = []
	helpItem = {'key': "p", 'text': "Toggle the plotting of DR2 sources on/off."}
	help.append(helpItem)
	helpItem = {'key': "h", 'text': "Plot only the Ha sources."}
	help.append(helpItem)
	helpItem = {'key': "w", 'text': "Invert the gray-scale."}
	help.append(helpItem)
	helpItem = {'key': "i/o", 'text': "Zoom in/out."}
	help.append(helpItem)
	helpItem = {'key': "q", 'text': "Quit."}
	help.append(helpItem)
		
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
		if keyPressed== '?':
			print("Help:")
			for helpItem in help:
				print("\t" + helpItem['key'] + " : " +  helpItem['text'])
			print()
		if keyPressed== 'p':
			if plotSources == True:
				plotSources = False
			else:
				plotSources = True
				plotHa = False
			redraw()
		
		if keyPressed== 'h':
			plotHa = True
			plotSources = False
			redraw() 
		
			
		if keyPressed== 'w':
			imageMinMax = (imageMinMax[1], imageMinMax[0])
			if invertedColours: invertedColours= False
			else: invertedColours=True
			redraw()
			
		if keyPressed=='i':
			print("Zoom requested at (%0.0f, %0.0f)"%(x, y))
			zoomFactor = 1.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			if (newWidth<100) or (newHeight<100):
				print("Maximum zoom reached.")
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
			print("new limits:", xlimits, ylimits)
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			print("new limits (world)", ra_limits, dec_limits)
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]
			redraw()
			
		if keyPressed=='o':
			if newWidth == width: continue
			print("Zoom out requested at (%0.0f, %0.0f)"%(x, y))
			zoomFactor = 0.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			
			if (newWidth >= width) or (newHeight>=height):
				print("Back to 1:1 scale.")
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
			print("new limits:", xlimits, ylimits)
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]	
			print("new limits (world)", ra_limits, dec_limits)
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			redraw()
		
		if plotSources: sourceStatus = "ON"
		else: sourceStatus = "OFF" 
		if plotHa: HaStatus = "ON"
		else: HaStatus = "OFF" 
		if invertedColours: invertStatus = "ON"
		else: invertStatus = "OFF" 
		print("Filter: %s PlotSources[%s]  Plot Ha extended[%s]  Invert Greyscale[%s]"%(filter, sourceStatus, HaStatus, invertStatus))
			
			
	# except KeyboardInterrupt:
	#	print "Ctrl-C pressed, but I dealt with it. "


	ppgplot.pgclos()
	
	
	
