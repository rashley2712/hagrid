from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.vo.client.conesearch import list_catalogs
from astropy.table import Table, vstack
from astropy.utils import data
import datetime

import numpy, math, os, sys, json
import generalUtils
import astroquery
import matplotlib
import matplotlib.pyplot
from matplotlib.path import Path

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
		self.AbsoluteLocationPixels = (self.x1 + xoffset, self.y1 + yoffset)
		self.ra, self.dec = wcsSolution.all_pix2world([self.AbsoluteLocationPixels[0]], [self.AbsoluteLocationPixels[1]], 1)
		
	def computeMax(self):
		""" Finds the max pixel in the data and saves the position as (xmax, ymax) """
		if self.type=="Maximum":
			maxPixel = numpy.ma.max(self.data)
			position = numpy.unravel_index(numpy.ma.argmax(self.data), self.data.shape)
		else:
			maxPixel = numpy.ma.min(self.data)
			position = numpy.unravel_index(numpy.ma.argmin(self.data), self.data.shape)
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
        'A5506-4'   : 'CCD1',
        'A5383-17-7': 'CCD2',
        'A5530-3'   : 'CCD3',
        'A5382-1-7' : 'CCD4'
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
		    'iclass': 'iClass', 
		    'haClass': 'haClass',
			'pixelFWHM': 'haSeeing',
		    'class': 'mergedClass',
			'pStar': 'pStar'
		    },
		'catalog_db': "http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&",
		'VizierLookup': 'dr2',
		'VizierName': 'II/321/iphas2',
		'colour': 'green'
	}
}

# 		'VizierName': 'II/321/iphas2',


def maskRadius(object, catalogName, CCDseeing):
	r = 0
	if catalogName=='dr2':
		if CCDseeing!=0:
		    r = 3.0 * 1.5 * CCDseeing
		else:
		    r = 3.0 * 1.5 * object['pixelFWHM'] 
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
	
def maskRadiusArray(objects, catalogName, CCDseeing):
	r = numpy.zeros(len(objects),dtype=numpy.float)
	if catalogName=='dr2':
		if CCDseeing!=0:
		    r = r + 3.0 * 1.5 * CCDseeing
		else:
		    r = 3.0 * 1.5 * objects['pixelFWHM'] 
	else:
                r = r + 50 # base value
                r = numpy.where(objects['mag']>12, 40*numpy.exp((-objects['mag']+12)/4.), r)
                r = numpy.where(objects['mag']<8, 250, r)
                r = numpy.where(objects['mag']<7, 350, r)
						
	return r
	


class IPHASdataClass:
	def __init__(self):
		print "Initialising an empty IPHAS data class"
		self.originalIPHASdb = None
		self.originalImageData = None
		self.archivePath = "."
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
		self.rejectTooManyMaskedPixels = 0.7/1.7 #0.70
                self.rBandRatioLimit = 0.7
		self.varianceThreshold = 5
		self.fullDebug = False
		self.objectStore = {}
		self.activeColour = 'r'
		self.CCD = "unknown"
		self.cachedir = "/tmp/hagrid/"
		self.cacheImages = False
		self.CCDseeing = 0
		self.rBandImageData = None
		self.rBandFITSHeaders = {}
		self.figure = None
		self.autoplot = True
		return None
		
	
		
	def setProperty(self, property, value):
		truths = ["true", "yes", "on", "1", "y"]
		falses = ["false", "no", "off", "0", "n"]
		if property=='archive':
			self.__dict__['archivePath'] = str(value)
		if property=="autoplot":
			if value.lower() in truths:
				self.autoplot = True
			if value.lower() in falses:
				self.autoplot = False
		if property=='cachedir':
			self.__dict__['cachedir'] = str(value)
		if property=="cacheimages":
			if value.lower() in truths:
				self.cacheImages = True
			if value.lower() in falses:
				self.cacheImages = False
		if property=='colour' or property=='color':
			self.__dict__['activeColour'] = str(value)
		if property=='debug':
			if value.lower() in truths:
				self.fullDebug = True
			if value.lower() in falses:
				self.fullDebug = False
		if property=="ignorecache":
			if value.lower() in truths:
				self.ignorecache = True
			if value.lower() in falses:
				self.ignorecache = False
		if property=='maglimit':
			self.__dict__['magLimit'] = float(value)
		if property=='plotwindowsize':
			self.__dict__['figSize'] = float(value)
		if property=='rbandratiolimit':
			self.__dict__['rBandRatioLimit'] = float(value)
		if property=='superpixelsize':
			self.__dict__['superPixelSize'] = int(value)
		if property=='spacinglimit':
			self.__dict__['spacingLimit'] = float(value)
			
	def getStoredObject(self, name):
		try:
			return self.objectStore[name]
		except KeyError:
			print "Could not find an object called %s in internal object storage."%name
		return None	
		
	def loadIPHASdb(self):
	        installPath = os.path.dirname(os.path.realpath(__file__))
		dbFilename = os.path.join(installPath, "iphas-images.fits.gz")
	
		try:
			self.originalIPHASdb = Table.read(dbFilename) 
			print "Loaded the %d rows from %s."%(len(self.originalIPHASdb), dbFilename)
		except IOError as e:
			print "Could not find the IPHAS database file: %s"%dbFilename
			return False
                return True

	def loadFITSFile(self, filename):
		
		self.filename = generalUtils.getFITSfilename(filename, self.archivePath, self.cachedir)
		print "Debugging self.filename", self.filename
		if self.filename == -1: return -1
		self.baseFilename = os.path.basename(self.filename)
		self.rootname = self.baseFilename.split(".")[0]
		
		self.originalImageData, self.FITSHeaders = fits.getdata(self.filename, header=True, uint=False, do_not_scale_image_data=False)
                if 'WFFBAND' in self.FITSHeaders:
			self.filter = self.FITSHeaders["WFFBAND"]
		self.height, self.width = numpy.shape(self.originalImageData)
		self.wcsSolution = WCS(self.FITSHeaders)
		print "width, height", self.width, self.height, "shape:", numpy.shape(self.originalImageData)
		self.getRADECmargins()
		imageCentre = (self.width/2, self.height/2)
		ra, dec = self.wcsSolution.all_pix2world([imageCentre], 1)[0]
		self.centre = (ra, dec)
		positionString = generalUtils.toSexagesimal((ra, dec))
		print "RA, DEC of image centre is: ", positionString, ra, dec
                try:
			self.CCD = chipNames[self.FITSHeaders['CCDNAME']]
		except Exception as e:
			print "WARNING: Could not identify CCD number (1, 2, 3 or 4)"
			print e
		
                if "SEEING" in self.FITSHeaders:
			self.CCDseeing = self.FITSHeaders['SEEING']
                else:
			print "WARNING: Could not find the 'seeing' value for this image in the FITS headers."
			
		return
		
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
		print "Debugging: filenameParts", filenameParts
		catalogCache = filenameParts[-3].split('/')[-1] + "_" + catalogName + "_cache.fits"
		cached = False
		if not self.ignorecache:
			localCacheFolder = self.cachedir
			cacheSubFolder = filenameParts[0].split('/')[-1][0:4]
			localCacheFile = os.path.join(localCacheFolder, cacheSubFolder, catalogCache)
			localCacheSubFolder = os.path.join(localCacheFolder, cacheSubFolder)
			print "Looking for a local cached copy of the catalogue:", localCacheFile, 
			if os.path.exists(localCacheFile):
				print "FOUND"
				newCatalog = Table.read(localCacheFile)
				cached = True
			
			if not cached: 
				print "NOT FOUND"
				print "Going online to fetch %s results from Vizier with mag limit %f."%(catalogName, self.magLimit)
				from astroquery.vizier import Vizier
				# Vizier.ROW_LIMIT = 25000
				from astropy import coordinates
				from astropy import units as u
				c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
				skyRA  = coordinates.Angle(self.raRange, unit = u.deg)
				skyDEC = coordinates.Angle(self.decRange, unit = u.deg)
				if self.CCD == "CCD2":
					width = coordinates.Angle(23, unit = u.arcmin)
					height = coordinates.Angle(12, unit = u.arcmin)
				else:
					height = coordinates.Angle(23, unit = u.arcmin)
					width = coordinates.Angle(12, unit = u.arcmin)

				print "Sky RA, DEC range:", skyRA, skyDEC
				print "Width, Height:", width, height

				print "going to Astroquery for:", catalogMetadata[catalogName]['VizierLookup']
				v = Vizier(columns=["all"], catalog= catalogMetadata[catalogName]['VizierName'])				
				v.ROW_LIMIT = 25000
				v.column_filters={"r":"<%d"%self.magLimit}
				# result = Vizier.query_region(coordinates = c, width = skyRA, height = skyDEC, catalog = catalogMetadata[catalogName]['VizierName'], verbose=True)
				result = v.query_region(coordinates = c, width = width, height = height, verbose=True)
				print result
				newCatalog = result[catalogMetadata[catalogName]['VizierName']]
				newCatalog.pprint()
			
			
				# Write the new catalog to the cache file
				# We might need to create the folder 
				print "Checking for existence of cache folder:", localCacheSubFolder
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

		newCatalog = Table()
		columnMapper = catalogMetadata[catalogName]['columns']
		
                # map columns to new table
                skipRow = numpy.zeros(len(catTable),dtype=numpy.bool)
		for key in columnMapper.keys():
                        newCatalog[key] = catTable[columnMapper[key]]
                        skipRow = numpy.logical_or(skipRow, numpy.isnan(newCatalog[key].astype(float)))

                # skip rows with missing data
		skippedRowCount = numpy.sum(skipRow)
                if skippedRowCount>0:
                        newCatalog.remove_rows(numpy.where(skipRow)[0])

		newCatalog['x'], newCatalog['y'] = self.wcsSolution.wcs_world2pix(newCatalog['ra'], newCatalog['dec'], 1)
			
		print "Skipped %d rows because they contained some null data."%skippedRowCount		
                # remove entries outside CCD
                trimRow = numpy.zeros(len(newCatalog),dtype=numpy.bool)
                trimRow = numpy.logical_or(trimRow, newCatalog['x']<0)
                trimRow = numpy.logical_or(trimRow, newCatalog['x']>self.width)
                trimRow = numpy.logical_or(trimRow, newCatalog['y']<0)
                trimRow = numpy.logical_or(trimRow, newCatalog['y']>self.height)

                trimmedRowCount = numpy.sum(trimRow)
		print "Rejected %d points for being outside of the CCD x, y pixel boundaries."%trimmedRowCount
                if trimmedRowCount>0:
                        newCatalog.remove_rows(numpy.where(trimRow)[0])

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
		# matplotlib.pyplot.pause(0.01)
		
		return
	
	def maskBadPixels(self):
		print "About to mask out the bad pixels for %s"%self.CCD
		
                if self.originalIPHASdb is None:
                        self.loadIPHASdb()

                # get used confidence map
                cpmFilename = os.path.join(self.archivePath,"confmaps",self.FITSHeaders["CONFMAP"]+".fz")
		print "Looking for a local copy at", cpmFilename
		if not os.path.exists(cpmFilename):
		        url = "http://www.iphas.org/data/images/confmaps/" + self.FITSHeaders["CONFMAP"]+".fz"
                        cpmPath=os.path.split(cpmFilename)[0]
                        if not os.path.exists(cpmPath):
                                os.makedirs(os.path.normpath(cpmPath))
			generalUtils.downloadFITS(url, cpmFilename)
		
		hdulist = fits.open(cpmFilename)
                self.badPixelMask = hdulist[int(self.CCD[3:4])].data
		hdulist.close()
		
		if self.mask is None:
			self.mask = numpy.zeros(numpy.shape(self.originalImageData))
			print "Creating a new blank mask of size:", numpy.shape(self.mask)

		isMasked = self.badPixelMask<90
		self.mask[isMasked] = 132
                if self.autoplot: self.drawMask()
		return
			
	def old_maskBadPixels(self):
		print "About to mask out the bad pixels for %s"%self.CCD
		
		installPath = os.path.dirname(os.path.realpath(__file__))
		print "Looking for a local copy in", installPath
		bpmName = "bpm_" + self.CCD + ".fits.fz"
		bpmFilename = installPath + "/" + bpmName
		if not os.path.exists(bpmFilename):
		        url = "http://www.devicecloud.co.uk/WFC/" + bpmName
			generalUtils.downloadFITS(url, bpmFilename)
		
		hdulist = fits.open(bpmFilename)
		self.badPixelMask = hdulist[0].data
		# print "old shape:", numpy.shape(self.badPixelMask)
		
		# Trim the bad pixel mask to match observed image size
                # use the information from trimsec
		startX = 54
		startY = 0
		self.badPixelMask = self.badPixelMask[startY:startY+self.height, startX:startX+self.width]
		
		if self.mask is None:
			self.mask = numpy.zeros(numpy.shape(self.originalImageData))
			print "Creating a new blank mask of size:", numpy.shape(self.mask)

		if self.CCD == "CCD3":
			# Add a clipped area to the vignetted part of CCD3
			vignetteMask = numpy.zeros(numpy.shape(self.originalImageData))
                        xcorner=1200 #530
                        ycorner=1200 #640
			for x in range(xcorner):
				ylim = ycorner - (x*ycorner)/xcorner
				for y in range(ylim):
					vignetteMask[y,x] = 132
			cornerMask = vignetteMask == 132
			self.mask[cornerMask] = 132
	
                # augment bad pixel masks (access array with y,x)
		if self.CCD == "CCD1":
		        self.badPixelMask[3575,  531] = 1
                        self.badPixelMask[3639, 1560:1569] = 1
                        self.badPixelMask[2185:2191, 1555:1560] = 1
                        self.badPixelMask[ 995:1005, 1762:1769] = 1
		if self.CCD == "CCD2":
                        self.badPixelMask[1350:4096, 837] = 1
                        self.badPixelMask[1063:1068, 1983:1988] = 1
		if self.CCD == "CCD3":
                        self.badPixelMask[2801:2880, 1250] = 1
                        self.badPixelMask[1723:4096, 1224] = 1
		if self.CCD == "CCD4":
                    self.badPixelMask[ 900:4096, 548:558] = 1
                    self.badPixelMask[ 990:1600, 464] = 1
                    self.badPixelMask[ 980:4096, 388:391] = 1
	
		isMasked = self.badPixelMask!=0
		self.mask[isMasked] = 132
                if self.autoplot: self.drawMask()
		return
			
	def maskrBand(self):
		print "About to mask out the pixels based on rband ratio"
		if self.mask is None:
			self.mask = numpy.zeros(numpy.shape(self.originalImageData))
			print "Creating a new blank mask of size:", numpy.shape(self.mask)
                # currently assumes that Halpha and r align
		f=12./(120./float(self.rBandFITSHeaders["EXPTIME"]))
                print "exposure time factor is",f
                a = self.originalImageData / self.rBandImageData
                print numpy.mean(self.originalImageData),numpy.median(self.originalImageData),numpy.min(self.originalImageData),numpy.max(self.originalImageData)
                print numpy.mean(self.rBandImageData),numpy.median(self.rBandImageData),numpy.min(self.rBandImageData),numpy.max(self.rBandImageData)
                print numpy.mean(a),numpy.median(a),numpy.min(a),numpy.max(a)
                isMasked = a < 0.7*f
                #isMasked = a < 0.7*numpy.median(a)
		self.mask[isMasked] = 132
                if self.autoplot: self.drawMask()
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
                        if self.autoplot: self.drawMask()
			return
			
		# Retrieve the catalogue
		try:
			catalog = self.catalogs[catalogName]
		except KeyError:
			print "Could not find a catalog called %s."%catalogName
			return
		
                rArray = maskRadiusArray(catalog, catalogName, self.CCDseeing)
		index = 1	
		for x, y, r in zip(catalog['x'], catalog['y'], rArray):
			#sys.stdout.write("\rMasking: %d of %d."%(index, len(catalog)))
			#sys.stdout.write(" %f %f %f "%(x, y, r))
			#sys.stdout.flush()
			self.mask = generalUtils.gridCircle(y, x, r, self.mask)
			index+= 1
		#sys.stdout.write("\n")
		#sys.stdout.flush()
	
                if self.autoplot: self.drawMask()
		
	def plotObject(self, objectName):
		objects = self.getStoredObject(objectName)
		
		colour = self.activeColour
		
		# Get the main plotting figure
		matplotlib.pyplot.figure(self.figure.number)
		print "Setting figure number to: ", self.figure.number
		
		for index, o in enumerate(objects):
			position = o.getPixelPosition()
			xoffset = o.maxPosition[1]
			yoffset = self.superPixelSize - 2 - o.maxPosition[0]
                        x=o.x1 + xoffset
                        y=self.height-(o.y1 + yoffset)
			matplotlib.pyplot.plot(x, y, color = colour, marker='o', markersize=15, mew=3, fillstyle='none')
			matplotlib.pyplot.annotate(str(index), (x+20, y), color=colour, fontweight='bold', fontsize=15)
			
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
		# matplotlib.pyplot.pause(0.01)
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
		if self.figure is None:
			self.figure = matplotlib.pyplot.figure(self.filename, figsize=(self.figSize/1.618, self.figSize))
		else:
			matplotlib.pyplot.figure(self.figure.number)
		self.figure.frameon = False
		self.figure.set_tight_layout(True)
		axes = matplotlib.pyplot.gca()
		axes.set_axis_off()
		self.figure.add_axes(axes)
		imgplot = matplotlib.pyplot.imshow(mplFrame, cmap="Oranges", interpolation='nearest')
		
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=False)
		# matplotlib.pyplot.pause(0.01)
		
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
		
	def applyMedianFilter(self, size):
		""" Applies a median filter of sizexsize to image
		"""
		import scipy.signal
		print "Applying a median filter of %d x %d to the image."%(size, size)
		
		originalImage = self.originalImageData
		filteredImage = scipy.signal.medfilt2d(originalImage, size)
		self.originalImageData = filteredImage
		return

	def makeSuperPixels(self):

		superPixelSize = self.superPixelSize
		borderMask = self.borderSize	
		width = self.width
		height = self.height
		
		# Draw the grid on the matplotlib panel
                if self.autoplot:
                    matplotlib.pyplot.figure(self.figure.number)
		    # axes = matplotlib.pyplot.gca()
		    for yStep in range(borderMask, self.height-borderMask, superPixelSize):
			matplotlib.pyplot.plot([borderMask, self.width - borderMask], [yStep, yStep], ls=':', color='g', lw=2)
		    for xStep in range(borderMask, self.width-borderMask, superPixelSize):
			matplotlib.pyplot.plot([xStep, xStep], [borderMask, self.height - borderMask], ls=':', color='r', lw=2)
		    matplotlib.pyplot.draw()
		    matplotlib.pyplot.show()
		    # matplotlib.pyplot.pause(0.01)
		# End of drawing
		
		pixelBitmapWidth = int((width - 2.*borderMask) / superPixelSize) + 1
		pixelBitmapHeight = int((height - 2.*borderMask) / superPixelSize) + 1
		pixelBitmap = numpy.zeros((pixelBitmapHeight, pixelBitmapWidth))
		pixelBitmap.fill(99E9) 
		
                # cut image area used for superpixels
                x1=borderMask
                x2=borderMask+pixelBitmapWidth*superPixelSize
                y1=borderMask
                y2=borderMask+pixelBitmapHeight*superPixelSize
                # which image data
                if False:
                    # original image
                    imageCopy = numpy.copy(self.originalImageData[y1:y2,x1:x2])
                    maskCopy = numpy.copy(self.mask[y1:y2,x1:x2])
                else:
                    # median filtered
                    from scipy.signal import medfilt2d
                    imageCopy = medfilt2d(self.originalImageData)[y1:y2,x1:x2]
                    maskCopy = numpy.copy(self.mask[y1:y2,x1:x2])
                    #from scipy.ndimage.filters import convolve
                    #maskCopy = convolve(self.mask,numpy.array([[1,1,1],[1,1,1],[1,1,1]]),mode="constant",cval=False)
			
                # create 4D view into 2D image
                from numpy.lib.stride_tricks import as_strided as ast
                block=(superPixelSize,superPixelSize)
                shape= (imageCopy.shape[0]/ block[0], imageCopy.shape[1]/ block[1])+ block
                strides= (block[0]* imageCopy.strides[0], block[1]* imageCopy.strides[1])+ imageCopy.strides
                blockImage = ast(imageCopy, shape= shape, strides= strides)
                strides= (block[0]* maskCopy.strides[0], block[1]* maskCopy.strides[1])+ maskCopy.strides
                blockMask= ast(maskCopy, shape= shape, strides= strides)
                blockBooleanMask = numpy.ma.make_mask(blockMask)
                blockMaskImage = numpy.ma.masked_array(blockImage, blockBooleanMask)
                # Create Superpixel table
		superPixelList=Table()
                superPixelList['x1'] = range(x1,x2,superPixelSize)*pixelBitmapHeight
                superPixelList['x2'] = superPixelList['x1'] + superPixelSize - 1
                superPixelList['xc'] = superPixelList['x1'] + superPixelSize/2.
                superPixelList['y1'] = numpy.repeat(range(y1,y2,superPixelSize),pixelBitmapWidth)
                superPixelList['y2'] = superPixelList['y1'] + superPixelSize - 1
                superPixelList['yc'] = superPixelList['y1'] + superPixelSize/2.
		superPixelList['mean'] = numpy.ma.mean(blockMaskImage,axis=(2,3)).flatten()
		superPixelList['median'] = numpy.ma.median(blockMaskImage,axis=(2,3)).flatten()
		superPixelList['min'] = numpy.ma.min(blockMaskImage,axis=(2,3)).flatten()
		superPixelList['max'] = numpy.ma.max(blockMaskImage,axis=(2,3)).flatten()
                # have not found a good solution here yet (this costs TIME)
                dataList=[]
                for y in range(blockMaskImage.shape[0]):
                    for x in range(blockMaskImage.shape[1]):
                        dataList.append(blockMaskImage[y,x,...])
		superPixelList['data'] = dataList

		variance = numpy.ma.var(blockMaskImage,axis=(2,3)).flatten()
		numPixels= numpy.ma.count(blockMaskImage,axis=(2,3)).flatten()
		superPixelList['varppixel'] = variance/numPixels
		numMaskedPixels = numpy.ma.count_masked(blockMaskImage,axis=(2,3)).flatten()
		superPixelList['maskedpixels'] = numMaskedPixels
                #superPixelList.write("test_spl.fits",overwrite=True)

                # reject "bad" rows
		rejectNanCount = len(superPixelList)
                superPixelList.remove_rows(numpy.where(numpy.isnan(superPixelList['mean'].astype(float)))[0])
		rejectNanCount -= len(superPixelList)
				
		rejectVarCount = len(superPixelList)
                superPixelList.remove_rows(numpy.where(superPixelList['varppixel']>self.varianceThreshold)[0])
		rejectVarCount -= len(superPixelList)

		rejectMaskCount = len(superPixelList)
		maskedRatio = superPixelList['maskedpixels'].astype(float)/(superPixelSize*superPixelSize)
                superPixelList.remove_rows(numpy.where(maskedRatio>self.rejectTooManyMaskedPixels)[0])
		rejectMaskCount -= len(superPixelList)

		print "%d pixels rejected for being NaN."%(rejectNanCount)
		print "%d pixels rejected for having too many masked pixels. Masked pixels > %2.2f%%"%(rejectMaskCount, self.rejectTooManyMaskedPixels*100)
		print "%d pixels rejected for having too large variance. Variance per pixel > %2.2f"%(rejectVarCount, self.varianceThreshold)
		
		bitmapX = (superPixelList['x1']-borderMask)/superPixelSize
		bitmapY = -1-(superPixelList['y1']-borderMask)/superPixelSize
		pixelBitmap[bitmapY, bitmapX] = superPixelList['mean']

                if self.autoplot:
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
                # smallest value comes first - hence we need to flip for top
                sortedIndex = numpy.argsort(self.superPixelList['mean'])
                if top:
                        sortedIndex=numpy.flipud(sortedIndex)
		
		
		pointings = []
		distanceLimitPixels = self.spacingLimit*60/self.pixelScale
		
                for index in sortedIndex:
                        s=self.superPixelList[index]
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
                        # Additional checks if we can choose this pointing
			reject = False
			# Check if this is not near to an existing pointing
			for p in pointings:
				if distanceP(p, pointingObject) < distanceLimitPixels: 
					reject = True
					break
		        # Compute the position of max and rband flux
			if not reject:
			        pointingObject.computeMax()	
			        pointingObject.computeAbsoluteLocation(self.wcsSolution)
                        # check on mean to rband ratio (avoids stars)
			if not reject and self.rBandImageData is not None:
                                self.attachRBand(pointingObject,single=True)
                                if pointingObject.rBandValue is None:
                                        reject = True
                                else:
                                        ratio = pointingObject.mean/pointingObject.rBandValue
                                        if ratio<self.rBandRatioLimit: reject = True
                        # end of checks
			if not reject: pointings.append(pointingObject)
			if len(pointings)>=number: break;
		
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
                        if self.originalIPHASdb is None:
                                self.loadIPHASdb()
		        # Find our field information from run number
		        run = int(self.FITSHeaders['RUN'])
                        runIndex = numpy.where(self.originalIPHASdb["run"]==run)[0][0]
		        # Find rband field information from run number
		        run = int(self.rBandFITSHeaders['RUN'])
                        rBandrunIndex = numpy.where(self.originalIPHASdb["run"]==run)[0][0]
                        # number of rows
                        nr = len(objects)

			objectTable = Table()
			ids = []
			for index, o in enumerate(objects):
				oid = self.rootname + "-%02d"%index
				if o.type=="Minimum": oid = "sky-" + oid
				ids.append(oid)
			
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
			cols.append(fits.Column(name='CCD', format = '4A', array = [self.CCD]*nr))
			cols.append(fits.Column(name='Ha_sky', format = 'E', array = [self.originalIPHASdb["skylevel"][runIndex]]*nr))
			cols.append(fits.Column(name='r', format = 'E', array = [o.rBandValue for o in objects]))
			cols.append(fits.Column(name='r_sky', format = 'E', array = [self.originalIPHASdb["skylevel"][rBandrunIndex]]*nr))
			cols = fits.ColDefs(cols)
			tbhdu = fits.BinTableHDU.from_columns(cols)
			
			prihdr = self.FITSHeaders.copy()
			prihdr['COMMENT'] = "Created by Hagrid on %s."%( datetime.datetime.ctime(datetime.datetime.now()))
			
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
		
	def attachRBand(self, pointings, single=False):
		""" Attach r-band values to the Halpha pointings
		"""
		if self.rBandImageData is None:
			print "There is no r-band image loaded."
			return

                if single:
                        sources = [pointings]
                else:
		        sources = self.getStoredObject(pointings)

		for s in sources:
			x, y = self.rBandwcsSolution.all_world2pix(s.ra, s.dec, 1)
			print s.AbsoluteLocationPixels[0], s.AbsoluteLocationPixels[1], "in Halpha translates to", x[0], y[0], "in r-band."
			index_x = int(round(x[0]))
			index_y = int(round(y[0]))
			try: 
                                r_value = numpy.mean(self.rBandImageData[index_y-1:index_y+1, index_x-1:index_x+1])
                                # adjust to exposure time difference
				r_value/=12./(120./float(self.rBandFITSHeaders["EXPTIME"]))
				s.rBandValue = r_value
			except IndexError:
				s.rBandValue = None
				print "Outside of CCD limits!"
		return
		
		

	def findMatch(self):
		""" Find a match to the current CCD
		"""
		
		# First check that we have a valid image loaded
		if self.originalImageData is None:
			print "There is no image loaded. Nothing to try to match to. Load one with the 'load' command."
			return -1
		
		JD = self.FITSHeaders['JD']
		ra, dec = self.centre
		radius = 5./60.
		print "Looking for a match to %s centred at RA: %f, DEC: %f taken on %f"%(self.filename, ra, dec, JD)
		
                if self.originalIPHASdb is None:
                        self.loadIPHASdb()

		# Find our field information from run number
		run = int(self.FITSHeaders['RUN'])
                runIndex = numpy.where(self.originalIPHASdb["run"]==run)[0][0]
                fid = self.originalIPHASdb["fieldid"][runIndex]

		# Filter out all but the r-band images for this field and CCD
                matches = numpy.where((self.originalIPHASdb["band"]=="r")&(self.originalIPHASdb["fieldid"]==fid)&(self.originalIPHASdb["ccd"]==int(self.CCD[3:4])))[0]
			
		print "%d images match the filter criterion."%(len(matches))
		IPHASdb = self.originalIPHASdb[matches]

		dates = IPHASdb['utstart']
		from astropy.time import Time
		t = Time(dates, format='isot', scale='utc')
		print "Taken on:", t.jd

		timeDifference = [(JD - j) * 24 * 60 for j in t.jd]
		minimum = 1E8
		closest = -1
		for index, t in enumerate(timeDifference):	
			print "%d: %s was taken %.1f minutes"%(index, IPHASdb[index]['url'], abs(t)), 
			if t<0: print "after",
			else:
				print "before",
			print "the original" 
			if abs(t) < minimum:
				minimum = abs(t)
				closest = index

		closestImage = IPHASdb[closest] 
		lastSlash =  self.filename.rfind('/')
		secondLastSlash = self.filename[:lastSlash].rfind('/')
		archivePath = self.filename[:secondLastSlash]
		
		url = closestImage['url']
		filenameParts = url.split('/')[-2:]
		filename = generalUtils.getFITSfilename(filenameParts[1], self.archivePath, self.cachedir)
		print "r band image filename:", filename
		print "Loading the image:", filename
		

		self.rBandImageData, self.rBandFITSHeaders = fits.getdata(filename, header=True, uint=False, do_not_scale_image_data=False)
		self.rBandwcsSolution = WCS(self.rBandFITSHeaders)

                if self.autoplot:
		        print "Boosting the image"
		        rBandBoostedImage = generalUtils.percentiles(numpy.copy(self.rBandImageData), 20, 99)
		        matplotlib.pyplot.ion()
		        # mplFrame = numpy.rot90(self.boostedImage)
		        mplFrame = rBandBoostedImage
		        mplFrame = numpy.flipud(mplFrame)
		        self.figure = matplotlib.pyplot.figure(filename, figsize=(self.figSize/1.618, self.figSize))
		        self.figure.frameon = False
		        self.figure.set_tight_layout(True)
		        axes = matplotlib.pyplot.gca()
		        axes.set_axis_off()
		        self.figure.add_axes(axes)
		        imgplot = matplotlib.pyplot.imshow(mplFrame, cmap="Reds", interpolation='nearest')
		
		        matplotlib.pyplot.draw()
		        matplotlib.pyplot.show(block=False)
		        matplotlib.pyplot.draw()
		        # matplotlib.pyplot.pause(0.01)
		

	def old_findMatch(self):
		""" Find a match to the current CCD
		"""
		
		# First check that we have a valid image loaded
		if self.originalImageData is None:
			print "There is no image loaded. Nothing to try to match to. Load one with the 'load' command."
			return -1
		
		JD = self.FITSHeaders['JD']
		ra, dec = self.centre
		radius = 5./60.
		print "Looking for a match to %s centred at RA: %f, DEC: %f taken on %f"%(self.filename, ra, dec, JD)
		
		installPath = os.path.dirname(os.path.realpath(__file__))
		dbFilename = os.path.join(installPath, "iphas-images.fits.gz")
	
		try:
			IPHASdb = Table.read(dbFilename) 
			print "Loaded the %d rows from %s."%(len(IPHASdb), dbFilename)
		except IOError as e:
			print "Could not find the IPHAS database file: %s"%dbFilename
			return -1

		# Filter out all but the r-band images
		band = 'r'
		filters = IPHASdb['band']
		matches = []
		for idx, f in enumerate(filters):
			if band == f: matches.append(idx)
			
		print "%d images match the filter criterion: band = %s."%(len(matches), band) 
		IPHASdb = IPHASdb[matches]

		# Find the nearest image by position
		from astropy.coordinates import SkyCoord
		from astropy import units as u
		center = SkyCoord(ra = ra * u.degree, dec = dec * u.degree)
		catalog_RAs = IPHASdb['ra']
		catalog_DECs = IPHASdb['dec']
		catalog = SkyCoord(ra = catalog_RAs * u.degree, dec = catalog_DECs * u.degree )
	
		results = center.separation(catalog).degree
		# results = sorted(results)
		matches = []
		for idx,r in enumerate(results):
			if r < radius : matches.append(idx)
		print "%d images have centres within %2.2f degrees of %s."%(len(matches), radius, center.to_string('hmsdms'))
		IPHASdb = IPHASdb[matches]

		dates = IPHASdb['utstart']
		from astropy.time import Time
		t = Time(dates, format='isot', scale='utc')
		print "Taken on:", t.jd

		timeDifference = [(JD - j) * 24 * 60 for j in t.jd]
		minimum = 1E8
		closest = -1
		for index, t in enumerate(timeDifference):	
			print "%d: %s was taken %.1f minutes"%(index, IPHASdb[index]['url'], abs(t)), 
			if t<0: print "after",
			else:
				print "before",
			print "the original" 
			if abs(t) < minimum:
				minimum = abs(t)
				closest = index

		closestImage = IPHASdb[closest] 
		lastSlash =  self.filename.rfind('/')
		secondLastSlash = self.filename[:lastSlash].rfind('/')
		archivePath = self.filename[:secondLastSlash]
		
		url = closestImage['url']
		filenameParts = url.split('/')[-2:]
		filename = generalUtils.getFITSfilename(filenameParts[1], self.archivePath, self.cachedir)
		print "r band image filename:", filename
		print "Loading the image:", filename
		
		self.rBandImageData, self.rBandFITSHeaders = fits.getdata(filename, header=True, uint=False, do_not_scale_image_data=False)
		self.rBandwcsSolution = WCS(self.rBandFITSHeaders)

                if self.autoplot:
		        print "Boosting the image"
		        rBandBoostedImage = generalUtils.percentiles(numpy.copy(self.rBandImageData), 20, 99)
		        matplotlib.pyplot.ion()
		        # mplFrame = numpy.rot90(self.boostedImage)
		        mplFrame = rBandBoostedImage
		        mplFrame = numpy.flipud(mplFrame)
		        self.figure = matplotlib.pyplot.figure(filename, figsize=(self.figSize/1.618, self.figSize))
		        self.figure.frameon = False
		        self.figure.set_tight_layout(True)
		        axes = matplotlib.pyplot.gca()
		        axes.set_axis_off()
		        self.figure.add_axes(axes)
		        imgplot = matplotlib.pyplot.imshow(mplFrame, cmap="Reds", interpolation='nearest')
		
		        matplotlib.pyplot.draw()
		        matplotlib.pyplot.show(block=False)
		        matplotlib.pyplot.draw()
		        # matplotlib.pyplot.pause(0.01)
		

