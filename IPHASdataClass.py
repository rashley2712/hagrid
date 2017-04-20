from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.vo.client.conesearch import list_catalogs
from astropy.table import Table, vstack
from astropy.utils import data
from matplotlib.path import Path

import numpy, math, os, sys, json
import generalUtils
import astroquery
import matplotlib.pyplot
import datetime

def distance(p1, p2):

	return math.sqrt( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 )
	
def distanceP(p1, p2):
	return math.sqrt( (p1.x-p2.x)**2 + (p1.y-p2.y)**2)
		
class Pointing:
	def __init__(self):
		self.x1 = 0
		self.y1 = 0
		self.x = 0    # Position of the centre of the superpixel
		self.y = 0    # ....
		self.mean = 0
		self.ra = 0
		self.dec = 0
		self.data = None
		self.maxPosition = (0, 0)
		self.peak = 0
		self.length = 0
		self.type = "Maximum"
		
	def __str__(self):
		return "mean: %3.2f  pos: (%d, %d) masked: %d"%(self.mean, self.x, self.y, numpy.ma.count_masked(self.data))
	
	def computeAbsoluteLocation(self, wcsSolution):
		xoffset = self.maxPosition[1]
		yoffset = self.length - 2 - self.maxPosition[0]
		self.AbsoluteLocationPixels = (self.x1 + xoffset, 4096 - (self.y1 + yoffset))
		self.ra, self.dec = wcsSolution.all_pix2world([self.AbsoluteLocationPixels[0]], [self.AbsoluteLocationPixels[1]], 1)
		print "ra, dec", self.ra, self.dec
		
	def computeMax(self):
		""" Finds the max pixel in the data and saves the position as (xmax, ymax) """
		print "Number of masked pixels in this data:", numpy.ma.count_masked(self.data)
		if self.type=="Maximum":
			maxPixel = numpy.ma.max(self.data)
			position = numpy.unravel_index(numpy.ma.argmax(self.data), self.data.shape)
		else:
			maxPixel = numpy.ma.min(self.data)
			position = numpy.unravel_index(numpy.ma.argmin(self.data), self.data.shape)
		print "max: %4.2f pos: (%3.2f, %3.2f)"%(maxPixel, position[0], position[1])
		self.peak = maxPixel
		self.maxPosition = position
		
	def getPixelPosition(self):
		# return (self.y, self.x)
		return ( self.y1 + self.maxPosition[0], self.x1 + self.maxPosition[1])
		
	def toJSON(self):
		jsonObject = {}
		jsonObject['xc'] = float(self.x)
		jsonObject['yc'] = float(self.y)
		jsonObject['xmax'] = float(self.x + self.maxPosition[1])
		jsonObject['ymax'] = float(self.y + self.maxPosition[0])
		jsonObject['ra'] = float(self.ra) 
		jsonObject['dec'] = float(self.dec) 
		jsonObject['mean'] = float(self.mean)
		jsonObject['peak'] = float(self.peak)
		return json.dumps(jsonObject)
	
chipNames = {
	'CCD1': 'A5506-4',
	'CCD2': 'A5383-17-7',
	'CCD3': 'A5530-3',
	'CCD4': 'A5382-1-7'
}
		
catalogMetadata = {
	'tycho': {
		'columns': {
			'ra': 'RA_ICRS_',
			'dec': 'DE_ICRS_',
			'B': 'BTmag',
			'V': 'VTmag',
			'mag': 'VTmag' },
		'catalog_db': "http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/259/tyc2&-out.all&",
		'colour': 'blue',
		'VizierName': 'I/259/tyc2',
		'VizierLookup': 'tycho'
		}, 
	'usno': {
		'columns': {
					'ra': 'RAJ2000',
					'dec': 'DEJ2000',
					'B': 'Bmag',
					'R': 'Rmag',
					'mag': 'Rmag' },
		'colour': 'blue',
		'VizierName': 'I/252/out',
		'VizierLookup': 'usno'
		},
	'dr2': {
		'columns': {
			'ra': 'RAJ2000', 
			'dec': 'DEJ2000',
			'i': 'i',
			'r': 'r',
			'H': 'ha',
			'mag': 'r',
		    'class': 'mergedClass',
			'pStar': 'pStar', 
		    'iclass': 'iClass', 
		    'haClass': 'haClass',
		    'pixelFWHM': 'haSeeing'},
		'catalog_db': "http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&",
		'VizierLookup': 'dr2',
		'VizierName': 'II/321/iphas2',
		'colour': 'green'
	}
}

def maskRadius(object, catalogName, CCDseeing):
	r = 0
	if catalogName=='dr2':
		if CCDseeing!=0:
		    r = object['pixelFWHM'] * 2.8 * CCDseeing
		else:
		    r = 4.0 * object['pixelFWHM'] 
	else:
		if object['mag']>12:
			r = 40*math.exp((-object['mag']+12)/4)
		elif object['mag']<8:
			r = 250
			if object['mag']<7:
				r = 350
		else: 
			r = 50
						
	return r
	


class IPHASdataClass:
	def __init__(self):
		print "Initialising an empty IPHAS data class"
		self.originalImageData = None
		self.boostedImage = None
		self.FITSHeaders = {}
		self.filter = None
		self.pixelScale = None
		self.centre = None
		self.filename = None
		self.rootname = "unknown"
		self.ignorecache = False
		self.catalogs = {}
		self.figSize = 8.
		self.previewSize = 4.
		self.magLimit = 18
		self.mask = None
		self.borderSize = 50
		self.superPixelSize = 50
		self.spacingLimit = 60./60.  # Minimum spacing of pointings in arcminutes
		self.rejectTooManyMaskedPixels = 0.70
		self.varianceThreshold = 5
		self.fullDebug = False
		self.objectStore = {}
		self.activeColour = 'r'
		self.CCD = "unknown"
		self.cachedir = "/tmp/hagrid/"
		self.CCDseeing = 0
		return None
		
	
		
	def setProperty(self, property, value):
		truths = ["true", "yes", "on", "1", "Y", "y", "True"]
		falses = ["false", "no", "off", "0", "N", "n", "False"]
		if property=='maglimit':
			self.__dict__['magLimit'] = float(value)
		if property=="ignorecache":
			if value in truths:
				self.ignorecache = True
			if value in falses:
				self.ignorecache = False
		if property=='superpixelsize':
			self.__dict__['superPixelSize'] = int(value)
		if property=='spacinglimit':
			self.__dict__['spacingLimit'] = float(value)
		if property=='plotwindowsize':
			self.__dict__['figSize'] = float(value)
		if property=='debug':
			if value in truths:
				self.fullDebug = True
			if value in falses:
				self.fullDebug = False
		if property=='colour' or property=='color':
			self.__dict__['activeColour'] = str(value)
		if property=='cachedir':
			self.__dict__['cachedir'] = str(value)
			

			
	def getStoredObject(self, name):
		try:
			return self.objectStore[name]
		except KeyError:
			print "Could not find an object called %s in internal object storage."%name
		return None	
		
	def loadFITSFile(self, filename):
		hdulist = fits.open(filename)
		self.filename = filename
		self.baseFilename = os.path.basename(filename)
		self.rootname = self.baseFilename.split(".")[0]
		FITSHeaders = []
		self.rawHeaders = hdulist[0].header
		for card in hdulist:
			# print(card.header.keys())
			# print(repr(card.header))
			for key in card.header.keys():
				self.FITSHeaders[key] = card.header[key]
				if 'WFFBAND' in key:
					self.filter = card.header[key]
		import astropy.io.fits as pf
		self.originalImageData = pf.getdata(filename, uint=False, do_not_scale_image_data=False)
		# self.originalImageData =  hdulist[1].data
		self.height, self.width = numpy.shape(self.originalImageData)
		self.wcsSolution = WCS(hdulist[1].header)
		print "width, height", self.width, self.height, "shape:", numpy.shape(self.originalImageData)
		self.getRADECmargins()
		imageCentre = (self.width/2, self.height/2)
		ra, dec = self.wcsSolution.all_pix2world([imageCentre], 1)[0]
		self.centre = (ra, dec)
		positionString = generalUtils.toSexagesimal((ra, dec))
		print "RA, DEC of image centre is: ", positionString, ra, dec
		try:
			ccdidentifier = self.FITSHeaders['CCDNAME']
			for k in chipNames.keys():
				if ccdidentifier==chipNames[k]:
					self.CCD = k
			print "CCD is:", self.CCD
		except Exception as e:
			print "WARNING: Could not identify CCD number (1, 2, 3 or 4)"
			print e
		
		try:
			self.CCDseeing = self.FITSHeaders['SEEING']
		except Exception as e:
			print "WARNING: Could not find the 'seeing' value for this image in the FITS headers."
			
		hdulist.close()
		
	def showVizierCatalogs(self):
		(ra, dec) = self.centre
		from astroquery.vizier import Vizier
		Vizier.ROW_LIMIT = 50
		from astropy import coordinates
		from astropy import units as u
		c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
		skyHeight= coordinates.Angle(self.raRange, unit = u.deg)
		results = Vizier.query_region(coordinates = c, radius= 1.0 * u.deg)
		print results
		
		
	def getVizierObjects(self, catalogName):
		""" Make a request to Vizier to get an Astropy Table of catalog object for this field. """
		(ra, dec) = self.centre
		
		availableCatalogs = catalogMetadata.keys()
		if catalogName not in availableCatalogs:
			print "The definitions for this catalogue are unknown. Available catalogues are:", availableCatalogs
			return
		
		# First look for a cached copy of this data
		filenameParts = self.filename.split('.')
		catalogCache = filenameParts[0] + "_" + catalogName + "_cache.fits"
		cached = False
		if not self.ignorecache:
			localCacheFolder = self.cachedir
			localCacheFile = os.path.join(localCacheFolder, filenameParts[0][0:4], catalogCache)
			localCacheSubFolder = os.path.join(localCacheFolder, filenameParts[0][0:4])
			print "Looking for a local cached copy of the catalogue:", localCacheFile, 
			if os.path.exists(localCacheFile):
				print "FOUND"
				newCatalog = Table.read(localCacheFile)
				cached = True
			
			if not cached: 
				print "NOT FOUND"
				print "Going online to fetch %s results from Vizier with mag limit %f."%(catalogName, self.magLimit)
				from astroquery.vizier import Vizier
				Vizier.ROW_LIMIT = 1E5
				Vizier.column_filters={"r":"<%d"%self.magLimit}
				from astropy import coordinates
				from astropy import units as u
				c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
				skyRA  = coordinates.Angle(self.raRange, unit = u.deg)
				skyDEC = coordinates.Angle(self.decRange, unit = u.deg)
				print "Sky RA, DEC range:", skyRA, skyDEC
				print "going to Astroquery for:", catalogMetadata[catalogName]['VizierLookup']
				result = Vizier.query_region(coordinates = c, width = skyRA, height = skyDEC, catalog = catalogMetadata[catalogName]['VizierName'], verbose=True)
				print result
				newCatalog = result[catalogMetadata[catalogName]['VizierName']]
				newCatalog.pprint()
			
			
				# Write the new catalog to the cache file
				# newCatalog.write(catalogCache, format='fits', overwrite=True)
				# We might need to create the folder 
				print "Checking for existance of cache folder:", localCacheSubFolder
				if not os.path.exists(self.cachedir):
					os.mkdir(self.cachedir)
				if not os.path.exists(localCacheSubFolder):
					 os.mkdir(localCacheSubFolder)    
				newCatalog.write(localCacheFile)
		
		self.addCatalog(newCatalog, catalogName)
		
		return
		
		
	def printCatalog(self, catalogName):
		catalog = self.catalogs[catalogName]
		for b in catalog:
			print b
		print "%d rows printed."%len(catalog)
		
	def typeObject(self, objectName):
		try:
			objects = self.objectStore[objectName]
			for index, o in enumerate(objects):
				print index, ":", o
		except KeyError:
			print "Could not find an object called %s stored internally."%objectName
	
	def addCatalog(self, catTable, catalogName):
		newCatalog = []
		columnMapper = catalogMetadata[catalogName]['columns']
		
		
		
		skippedRowCount = 0
		for index, row in enumerate(catTable):
			object={}
			skipRow = False
			for key in columnMapper.keys():
				object[key] = row[columnMapper[key]]
				if numpy.isnan(float(row[columnMapper[key]])): skipRow = True
			if skipRow: 
				skippedRowCount+= 1
				continue		
			x, y = self.wcsSolution.all_world2pix([object['ra']], [object['dec']], 1)
			object['x'] = x[0]
			object['y'] = y[0]
			
			newCatalog.append(object)
			if  ((index+1)%100) == 0:
				sys.stdout.write("\rCopying: %d of %d."%(index+1, len(catTable)))
				sys.stdout.flush()
		sys.stdout.write("\rCopying: %d of %d.\n"%(index+1, len(catTable)))
		sys.stdout.flush()
		
		print "Skipped %d rows because they contained some null data."%skippedRowCount		
		trimmedCatalog = []
		for row in newCatalog:
			if row['x']<0: continue
			if row['x']>self.width: continue
			if row['y']<0: continue
			if row['y']>self.height: continue
			trimmedCatalog.append(row)
		print "Rejected %d points for being outside of the CCD x, y pixel boundaries."%(len(newCatalog)-len(trimmedCatalog))
		newCatalog = trimmedCatalog

		print "Adding catalog %s to list of stored catalogs."%catalogName
		self.catalogs[catalogName] =  newCatalog
		
		return
				
	def getRADECmargins(self):
		boundingBox = self.wcsSolution.all_pix2world([[0, 0], [0, self.width], [self.height, self.width], [self.height, 0]], 1, ra_dec_order = True)
		# boundingBox = self.wcsSolution.all_pix2world([[0, 0], [0, self.height], [self.width, self.height], [self.width, 0]], 1, ra_dec_order = True)
		print "Bounding box:", boundingBox
		pixelDiagonal = math.sqrt(self.height**2 + self.width**2)
		pixel1 = boundingBox[0]
		pixel2 = boundingBox[2]
		skyDiagonal = distance(pixel1, pixel2)
		print "Diagonal size:", pixelDiagonal, skyDiagonal
		self.pixelScale = (skyDiagonal / pixelDiagonal) * 3600.
		self.pixelScale = 0.333
		raMin = numpy.min([r[0] for r in boundingBox])
		raMax = numpy.max([r[0] for r in boundingBox])
		decMin = numpy.min([r[1] for r in boundingBox])
		decMax = numpy.max([r[1] for r in boundingBox])
		print "RA, DEC min/max:", raMin, raMax, decMin, decMax
		raRange = raMax - raMin
		decRange = decMax - decMin
		print "RA range, DEC range", raRange, decRange, raRange*60, decRange*60
		self.raRange = raRange
		self.decRange = decRange
		print "Pixel scale: %6.4f \"/pixel"%self.pixelScale
		self.boundingBox = boundingBox
		
	def showFITSHeaders(self):
		headersString = ""
		for key in self.FITSHeaders.keys():
			print key + " : " + str(self.FITSHeaders[key])
			headersString+= str(key) + " : " + str(self.FITSHeaders[key]) + "\n"
		return headersString
			
	def getFITSHeader(self, key):
		try:
			print key + " : " + str(self.FITSHeaders[key]) 
			return self.FITSHeaders[key]
		except KeyError:
			print "Could not find a header with the name:", key
			return None 
			
	def plotCatalog(self, catalogName):
		try:
			catalog = self.catalogs[catalogName]
		except KeyError:
			print "Could not find a catalog called %s."%catalogName
			return
		
		catalogColour = catalogMetadata[catalogName]['colour']
		try:
			xArray = []
			yArray = []
			rArray = []
			for o in catalog:
				# Check that the catalog has a class flag
				if 'class' in o.keys():
					if o['class'] != -1: continue   # Skip objects that are not stars  
				xArray.append(o['x'] - 1)
				yArray.append(self.height - 1 - o['y'] )
				rArray.append(maskRadius(o, catalogName, self.CCDseeing))
	
			# Nick Wright 
			# R / pixels = 8192/M^2 + 1000/M + 100 
			matplotlib.pyplot.figure(self.figure.number)
			patches = [matplotlib.pyplot.Circle((x_, y_), s_, fill=False, linewidth=1) for x_, y_, s_ in numpy.broadcast(xArray, yArray, rArray)]
			collection = matplotlib.collections.PatchCollection(patches, alpha = 0.25, color = catalogColour)
			ax = matplotlib.pyplot.gca()
			ax.add_collection(collection)
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			matplotlib.pyplot.pause(0.01)
			# matplotlib.pyplot.savefig("test.png", bbox_inches='tight')
		except AttributeError as e:
			print "There is no drawing surface defined yet. Please use the 'draw' command first."
			print e
		except Exception as e:
			print e
			
			
	def drawMask(self):
		if self.mask is None:
			print "There is no mask defined yet."
			return
		self.maskFigure = matplotlib.pyplot.figure(self.filename + " mask", figsize=(self.figSize/1.618, self.figSize))
		self.maskFigure.frameon = False
		self.maskFigure.set_tight_layout(True)
		axes = matplotlib.pyplot.gca()
		axes.set_axis_off()
		self.maskFigure.add_axes(axes)
		imgplot = matplotlib.pyplot.imshow(numpy.flipud(self.mask), cmap="gray_r", interpolation='nearest')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
		matplotlib.pyplot.pause(0.01)
		
		return
	
	def maskBadPixels(self):
		print "About to mask out the bad pixels for %s"%self.CCD
		
		installPath = os.path.realpath(__file__).rsplit('/',1)[0]
		print "Looking for a local copy in", installPath
		bpmName = "bpm_" + self.CCD + ".fits.fz"
		url = "http://www.devicecloud.co.uk/WFC/" + bpmName
		bpmFilename = installPath + "/" + bpmName
		if not os.path.exists(bpmFilename):
			generalUtils.downloadFITS(url, bpmFilename)
		
		hdulist = fits.open(bpmFilename)
		FITSHeaders = []
		self.badPixelMask = hdulist[0].data
		# print "old shape:", numpy.shape(self.badPixelMask)
		
		# Trim the bad pixel mask to match observed image size
		startX = 47
		startY = 1
		self.badPixelMask = self.badPixelMask[startY:startY+4096, startX:startX+2048]
		# print "new shape:", numpy.shape(self.badPixelMask)
		# print "reqd shape:", numpy.shape(self.originalImageData)
		
		if self.mask is None:
			self.mask = numpy.zeros(numpy.shape(self.originalImageData))
			print "Creating a new blank mask of size:", numpy.shape(self.mask)

		if self.CCD == "CCD3":
			# Add a clipped area to the vignetted part of CCD3
			vignetteMask = numpy.zeros(numpy.shape(self.originalImageData))
			for x in range(530):
				ylim = int(x*(-1.207) + 640)
				for y in range(ylim):
					vignetteMask[y,x] = 132
			cornerMask = vignetteMask == 132
			self.mask[cornerMask] = 132
	
	
		isMasked = self.badPixelMask!=0
		self.mask[isMasked] = 132
		self.drawMask()
		return
			
			
	def maskCatalog(self, catalogName):
		if self.mask is None:
			self.mask = numpy.zeros(numpy.shape(self.originalImageData))
			print "Creating a new blank mask of size:", numpy.shape(self.mask)

		# Mask the border areas
		if catalogName == 'border':
			border = self.borderSize
			self.mask[0:border, 0:self.width] = 132
			self.mask[self.height-border:self.height, 0:self.width] = 132
			self.mask[0:self.height, 0:border] = 132
			self.mask[0:self.height, self.width-border:self.width] = 132
			self.drawMask()
			return
			
		# Retrieve the catalogue
		try:
			catalog = self.catalogs[catalogName]
		except KeyError:
			print "Could not find a catalog called %s."%catalogName
			return
		

		xArray = []
		yArray = []
		rArray = []
		for o in catalog:
			xArray.append(o['x'])
			yArray.append(o['y'])
			rArray.append(maskRadius(o, catalogName, self.CCDseeing))
			
		index = 1	
		for x, y, r in zip(xArray, yArray, rArray):
			sys.stdout.write("\rMasking: %d of %d."%(index, len(catalog)))
			sys.stdout.write(" %f %f %f "%(x, y, r))
			sys.stdout.flush()
			self.mask = generalUtils.gridCircle(y, x, r, self.mask)
			index+= 1
		sys.stdout.write("\n")
		sys.stdout.flush()
	
		self.drawMask()
		
	def plotObject(self, objectName):
		objects = self.getStoredObject(objectName)
		
		colour = self.activeColour
		
		# Get the main plotting figure
		matplotlib.pyplot.figure(self.figure.number)
		
		for index, o in enumerate(objects):
			position = o.getPixelPosition()
			# print position
			# matplotlib.pyplot.plot(o.x, o.y, color = 'r', marker='o', markersize=25, lw=4, fillstyle='none')
			xoffset = o.maxPosition[1]
			yoffset = self.superPixelSize - 2 - o.maxPosition[0]
			# print "offsets", xoffset, yoffset
			matplotlib.pyplot.plot(o.x1 + xoffset, o.y1 + yoffset , color = colour, marker='o', markersize=15, mew=3, fillstyle='none')
			matplotlib.pyplot.annotate(str(index), (o.x1+xoffset+20, o.y1+yoffset), color=colour, fontweight='bold', fontsize=15)
			# if index==2: break
			
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
		matplotlib.pyplot.pause(0.01)
		return
		
		
	def drawPreview(self, pointingsName, index, title=None):
		if title is None:
			title = "Preview of pointing number %d in %s"%(index, pointingsName)
		print "Creating preview: %s"%title
		
		
		objectList = self.getStoredObject(pointingsName)
		if objectList is None: return
		
		pointingObject = objectList[index]
		
		print "mean: %f"%pointingObject.mean
		self.previewFigure = matplotlib.pyplot.figure(title, figsize=(self.previewSize, self.previewSize))
		self.previewFigure.frameon = False
		self.previewFigure.set_tight_layout(True)
		axes = matplotlib.pyplot.gca()
		axes.cla()
		axes.set_axis_off()
		self.previewFigure.add_axes(axes)
		imgplot = matplotlib.pyplot.imshow(numpy.flipud(pointingObject.data), cmap="hsv", interpolation='nearest')
		matplotlib.pyplot.plot(pointingObject.maxPosition[1], self.superPixelSize - 2 - pointingObject.maxPosition[0], color = 'r', marker='o', markersize=25, lw=4, fillstyle='none')
		matplotlib.pyplot.plot(10, 10, color = 'g', marker='x')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
		matplotlib.pyplot.pause(0.01)
		
		print pointingObject.data
		print pointingObject.ra, pointingObject.dec, generalUtils.toSexagesimal((pointingObject.ra, pointingObject.dec))
		return
		
			
	def drawBitmap(self):
		if self.boostedImage is None:
			print "Boosting the image"
			self.boostedImage = generalUtils.percentiles(numpy.copy(self.originalImageData), 20, 99)
		matplotlib.pyplot.ion()
		# mplFrame = numpy.rot90(self.boostedImage)
		mplFrame = self.boostedImage
		mplFrame = numpy.flipud(mplFrame)
		self.figure = matplotlib.pyplot.figure(self.filename, figsize=(self.figSize/1.618, self.figSize))
		self.figure.frameon = False
		self.figure.set_tight_layout(True)
		axes = matplotlib.pyplot.gca()
		axes.set_axis_off()
		self.figure.add_axes(axes)
		imgplot = matplotlib.pyplot.imshow(mplFrame, cmap="gray_r", interpolation='nearest')
		
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.draw()
		matplotlib.pyplot.pause(0.01)
		
		# matplotlib.pyplot.savefig("test.png",bbox_inches='tight')
		
	def applyMask(self):
		if self.mask is None:
			print "There is no mask defined. Define one with the 'mask' command."
			return
		
		if self.originalImageData is None:
			print "There is no source bitmap defined. Load one with the 'load' command."
			return
			
			
		booleanMask = numpy.ma.make_mask(numpy.flipud(self.mask))
		maskedImageData = numpy.ma.masked_array(self.originalImageData,  numpy.logical_not(booleanMask))
		
		self.maskedImage = maskedImageData
		
		matplotlib.pyplot.figure(self.figure.number)
		axes = matplotlib.pyplot.gca()
		imgplot = matplotlib.pyplot.imshow(self.maskedImage, cmap="gray_r", interpolation='nearest')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
		matplotlib.pyplot.pause(0.01)

	def makeSuperPixels(self):
		superPixelList = []
		superPixelSize = self.superPixelSize
		borderMask = self.borderSize	
		width = self.width
		height = self.height
		
		# Draw the grid on the matplotlib panel
		matplotlib.pyplot.figure(self.figure.number)
		# axes = matplotlib.pyplot.gca()
		for yStep in range(borderMask, self.height-borderMask, superPixelSize):
			matplotlib.pyplot.plot([borderMask, self.width - borderMask], [yStep, yStep], ls=':', color='g', lw=2)
		for xStep in range(borderMask, self.width-borderMask, superPixelSize):
			matplotlib.pyplot.plot([xStep, xStep], [borderMask, self.height - borderMask], ls=':', color='r', lw=2)
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
		matplotlib.pyplot.pause(0.01)
		# End of drawing
		
		imageCopy = numpy.copy(self.originalImageData)
		booleanMask = numpy.ma.make_mask(self.mask)
		maskedImageCopy = numpy.ma.masked_array(imageCopy, booleanMask)
			
		pixelBitmapWidth = int((width - 2.*borderMask) / superPixelSize) + 1
		pixelBitmapHeight = int((height - 2.*borderMask) / superPixelSize) + 1
		pixelBitmap = numpy.zeros((pixelBitmapHeight, pixelBitmapWidth))
		pixelBitmap.fill(99E9) 
		
		rejectMaskCount = 0
		rejectVarCount = 0
		index = 0
		for yStep in range(borderMask, self.height-borderMask, superPixelSize):
			# matplotlib.pyplot.plot([borderMask, self.width - borderMask], [yStep, yStep], ls=':', color='g')
			for xStep in range(borderMask, self.width-borderMask, superPixelSize):
				"""index+=1
				if index>30: return
				"""
				x1 = xStep
				x2 = xStep + superPixelSize - 1
				y1 = yStep
				y2 = yStep + superPixelSize - 1
				# print xStep, yStep, x1, x2, y1, y2, self.height-y2, self.height-y1
				superPixel = maskedImageCopy[self.height-y2:self.height-y1, x1:x2]
				superPixelObject = {}
				mean = float(numpy.ma.mean(superPixel))
				if math.isnan(mean): continue;
				superPixelObject['mean'] = mean
				superPixelObject['median'] = numpy.ma.median(superPixel)
				superPixelObject['min'] = numpy.ma.min(superPixel)
				superPixelObject['max'] = numpy.ma.max(superPixel)
				superPixelObject['x1'] = x1
				superPixelObject['y1'] = y1
				superPixelObject['x2'] = x2
				superPixelObject['y2'] = y2
				superPixelObject['xc'] = x1 + superPixelSize/2.
				superPixelObject['yc'] = y1 + superPixelSize/2.
				superPixelObject['data'] = superPixel
				
				
				bitmapX = (x1-borderMask)/superPixelSize
				bitmapY = (y1-borderMask)/superPixelSize
				
				if self.fullDebug:
					matplotlib.pyplot.figure(self.figure.number)
					matplotlib.pyplot.plot(superPixelObject['xc'], superPixelObject['yc'], color = 'r', marker='x', markersize=25, lw=4, fillstyle='none')
					matplotlib.pyplot.draw()
					matplotlib.pyplot.show()
					matplotlib.pyplot.pause(0.001)
				
					self.previewFigure = matplotlib.pyplot.figure("Superpixel", figsize=(self.previewSize, self.previewSize))
					self.previewFigure.frameon = False
					self.previewFigure.set_tight_layout(True)
					axes = matplotlib.pyplot.gca()
					axes.cla()
					axes.set_axis_off()
					self.previewFigure.add_axes(axes)
					imgplot = matplotlib.pyplot.imshow(numpy.flipud(superPixelObject['data']), cmap="gray_r", interpolation='nearest')
					matplotlib.pyplot.draw()
					matplotlib.pyplot.show()
					matplotlib.pyplot.pause(0.001)
					print x1, x2, y1, y2, superPixelObject['mean'], superPixelObject['median']
					raw_input("Press Enter to continue...")
					

				variance = numpy.ma.var(superPixel)
				numPixels= numpy.ma.count(superPixel)
				superPixelObject['varppixel'] = variance/numPixels
				if superPixelObject['varppixel']>self.varianceThreshold: 
					rejectVarCount+= 1
					continue
				
				numMaskedPixels = numpy.ma.count_masked(superPixel)
				superPixelObject['maskedpixels'] = numMaskedPixels
				maskedRatio = float(numMaskedPixels)/float(numPixels)
				if maskedRatio>self.rejectTooManyMaskedPixels: 
					# print "too many masked pixels here. Rejecting."
					rejectMaskCount+=1
					continue;
				superPixelList.append(superPixelObject)
				pixelBitmap[bitmapY, bitmapX] = mean
				
		print "%d pixels rejected for having too many masked pixels. Masked pixel ratio > %2.2f%%"%(rejectMaskCount, self.rejectTooManyMaskedPixels)
		print "%d pixels rejected for having too large variance. Variance per pixel > %2.2f"%(rejectVarCount, self.varianceThreshold)
		
		self.sampledImageFigure = matplotlib.pyplot.figure("Sampled Image", figsize=(self.figSize/1.618, self.figSize))
		self.sampledImageFigure.frameon = False
		self.sampledImageFigure.set_tight_layout(True)
		axes = matplotlib.pyplot.gca()
		axes.cla()
		axes.set_axis_off()
		self.sampledImageFigure.add_axes(axes)
		
		maskedPixelImage = numpy.ma.masked_equal(pixelBitmap, 99E9)
		
		
		# minimumPixel = numpy.min(pixelBitmap)
		# pixelBitmap[pixelBitmap==99E9] = minimumPixel
		
		imgplot = matplotlib.pyplot.imshow(maskedPixelImage, cmap="hsv", interpolation='nearest')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
				
		
		self.superPixelList = superPixelList
		return 
		
	def getRankedPixels(self, number=50):
		# Top sources
		top = True
		if number<0:
			top = False
			number = abs(number)
			
		# Sort superpixels
		if top: self.superPixelList.sort(key=lambda x: x['mean'], reverse=True)
		else: self.superPixelList.sort(key=lambda x: x['mean'], reverse=False)
		
		
		
		pointings = []
		distanceLimitPixels = self.spacingLimit*60/self.pixelScale
		
		for index, s in enumerate(self.superPixelList):
			print index, s['mean'], s['varppixel'], s['xc'], s['yc']
			pointingObject = Pointing()
			pointingObject.length = self.superPixelSize
			pointingObject.x1 = s['x1']
			pointingObject.y1 = s['y1']
			pointingObject.x2 = s['x2']
			pointingObject.y2 = s['y2']
			
			pointingObject.x = s['xc']
			pointingObject.y = s['yc']
			pointingObject.mean = s['mean']
			pointingObject.varppixel = s['varppixel']
			pointingObject.data = s['data']
			if top: pointingObject.type = "Maximum"
			else: pointingObject.type = "Minimum"
			# Check if this is not near to an existing pointing
			reject = False
			for p in pointings:
				if distanceP(p, pointingObject) < distanceLimitPixels: 
					reject=True
					break
			if not reject: pointings.append(pointingObject)
			if len(pointings)>=number: break;
		
		# Compute the position of the max for each pointing and store it internally
		for p in pointings:
			p.computeMax()	
			p.computeAbsoluteLocation(self.wcsSolution)
		return pointings
		
	def clearFigure(self):
		""" Clears the main drawing window """
		print "Clearing the current figure."
		matplotlib.pyplot.figure(self.figure.number)
		matplotlib.pyplot.clf()
		return
		
	def dumpImage(self, filename):
		matplotlib.pyplot.figure(self.figure.number)
		filename = filename.format(root = self.rootname)
		extension = os.path.splitext(filename)[1]
		if not extension==".png":
			filename+= ".png" 
		matplotlib.pyplot.savefig(filename,bbox_inches='tight')
		return
		
	
	def dumpObject(self, objectName, filename, outputFormat):
		print "About to dump %s"%objectName
		objects = self.getStoredObject(objectName)
			
		filename = filename.format(root = self.rootname)
			
		if outputFormat=="json":
			objectList = [o.toJSON() for o in objects]
			outputFile = open(filename, "wt")
			outputFile.write(json.dumps(objectList))
			outputFile.close()
			return
		
		if outputFormat=="fits" or outputFormat=="votable":
			objectTable = Table()
			ids = []
			for index, o in enumerate(objects):
				id = self.rootname + "-%02d"%index
				if o.type=="Minimum": id = "sky-" + self.rootname + "-%02d"%index
				ids.append(id)
			
			hdu = fits.PrimaryHDU()
			cols = []
			cols.append(fits.Column(name='id', format='16A', array = ids))
			cols.append(fits.Column(name='ra', format='E', array = [o.ra for o in objects]))
			cols.append(fits.Column(name='dec', format = 'E', array = [o.dec for o in objects]))
			cols.append(fits.Column(name='xmax', format = 'E', array = [o.AbsoluteLocationPixels[0] for o in objects]))
			cols.append(fits.Column(name='ymax', format = 'E', array = [o.AbsoluteLocationPixels[1] for o in objects]))
			cols.append(fits.Column(name='mean', format = 'E', array = [o.mean for o in objects]))
			cols.append(fits.Column(name='peak', format = 'E', array = [o.peak for o in objects]))
			cols.append(fits.Column(name='variance', format = 'E', array = [o.varppixel for o in objects]))
			cols.append(fits.Column(name='type', format = '8A', array = [o.type for o in objects]))
			cols.append(fits.Column(name='CCD', format = '4A', array = [self.CCD for i in objects]))
			cols = fits.ColDefs(cols)
			tbhdu = fits.BinTableHDU.from_columns(cols)
			
			prihdr = fits.Header()
			prihdr['COMMENT'] = "Created by Hagrid on %s."%( datetime.datetime.ctime(datetime.datetime.now()))
			
			for key in self.FITSHeaders.keys():
				try:
					# prihdr[key] = str(self.FITSHeaders[key]).encode("ascii")
					prihdr[key] = self.FITSHeaders[key]
				except ValueError:
					print "Warning... could not transfer the header %s to the new FITS file."%(key)
			
			prihdu = fits.PrimaryHDU(header=prihdr)
			thdulist = fits.HDUList([prihdu, tbhdu])
			try:
				thdulist.writeto(filename, overwrite=True)
			except:
				thdulist.writeto(filename, clobber=True)
			
			return
		
	
	def trim(self, topObject, bottomObject):
		
		import scipy.optimize

		def exponentialDecay(x, a0, a1, a2):
			y = a0 * numpy.exp(a1 *x) + a2
			return y
		
		def exponentialRise(x, a0, a1, a2):
			y = a2 - (a2-a0) * numpy.exp(a1 *x) 
			return y
		
		
		print "Going to trim from %s and %s"%(topObject, bottomObject)
		topSources = self.getStoredObject(topObject)
		bottomSources = self.getStoredObject(bottomObject)
		topMeans = [t.mean for t in topSources]
		bottomMeans = [b.mean for b in bottomSources]
		#print topMeans
		#print bottomMeans
		topRange = numpy.max(topMeans) - numpy.min(topMeans)
		bottomRange = numpy.max(bottomMeans) - numpy.min(bottomMeans)
		print "Top:", numpy.max(topMeans), numpy.min(topMeans), topRange
		print "Bottom:", numpy.max(bottomMeans), numpy.min(bottomMeans), bottomRange
		
		# Plot top means
		chart = matplotlib.pyplot.figure("Mean values", figsize=(10, 8))
		axes = matplotlib.pyplot.gca()
		axes.set_xlabel("Superpixel rank")
		axes.set_ylabel("Mean counts")
		chart.add_axes(axes)
		matplotlib.pyplot.plot(topMeans, color='r')
		xfit = numpy.arange(0, len(topSources))
		a0 = numpy.max(topMeans) - numpy.min(topMeans)
		a1 = -0.05
		a2 = numpy.min(topMeans)
		guess = [a0, a1, a2]
		results, covariance = scipy.optimize.curve_fit(exponentialDecay, xfit, topMeans, guess)
		print "fit result:", results
		errors = numpy.sqrt(numpy.diag(covariance))
		print "fit errors:", errors
		a0 = results[0]
		a1 = results[1]
		a2 = results[2]
		yfit = exponentialDecay(xfit, a0, a1, a2)
		matplotlib.pyplot.plot(yfit, color='r', ls='dashed')
		topFolding = -1 / a1
		
		xfit = numpy.arange(0, len(bottomSources))
		a0 = numpy.min(bottomMeans)
		a1 = -0.05
		a2 = numpy.max(bottomMeans)
		guess = [a0, a1, a2]
		results, covariance = scipy.optimize.curve_fit(exponentialRise, xfit, bottomMeans, guess)
		print "fit result:", results
		errors = numpy.sqrt(numpy.diag(covariance))
		print "fit errors:", errors
		a0 = results[0]
		a1 = results[1]
		a2 = results[2]
		yfit = exponentialRise(xfit, a0, a1, a2)
		matplotlib.pyplot.plot(yfit, color='b', ls='dashed')
		bottomFolding = -1 / a1
		if topFolding<len(topSources):
			matplotlib.pyplot.plot([topFolding, topFolding],  [axes.get_ylim()[0], axes.get_ylim()[1]],  color='r', ls='dotted')
		if bottomFolding<len(bottomSources):
			matplotlib.pyplot.plot([bottomFolding, bottomFolding],  [axes.get_ylim()[0], axes.get_ylim()[1]],  color='b', ls='dotted')
		
		matplotlib.pyplot.plot(bottomMeans, color='b')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(0.01)
		
		print "Top e-folding after %d pixels."%topFolding
		print "Bottom e-folding after %d pixels."%bottomFolding
		
		if topFolding<len(topSources):
			trimmedTop = []
			for index in range(int(topFolding)): trimmedTop.append(topSources[index])
			self.objectStore[topObject] = trimmedTop
			
		if bottomFolding<len(bottomSources):
			trimmedBottom = []
			for index in range(int(bottomFolding)): trimmedBottom.append(bottomSources[index])
			self.objectStore[bottomObject] = trimmedBottom
		
		return
		
	def listPixels(self, number=0):
		for index, s in enumerate(self.superPixelList):
			print s['mean'], s['xc'], s['yc']
			print s
			if number!=0 and index==number:
				return
		return
		
		"""print "Original range:", numpy.min(self.originalImageData), numpy.max(self.originalImageData)
		print "Masked range:", numpy.min(self.maskedImage), numpy.max(self.maskedImage)
		"""
		
		
