from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import datetime

class sourcesClass:
	def __init__(self):
		self.sources = []
		
	def addSourcesFromFITS(self, filename):
		try:
			hdulist = fits.open(filename)
			header = hdulist[0].header
			tableData = hdulist[1].data
			cols = hdulist[1].columns
			added = 0
			for d in tableData:
				rowObject = {}
				for c in cols.names:
					rowObject[c] = d[c]
				added+= 1
				self.sources.append(rowObject)
				totalRows = len(self.sources)
		except IOError as e:
			print "Could not load: %s"%filename
			return (-1, -1)
		
		hdulist.close()	
		return (added, totalRows)
		  
	def printSources(self, n=-1):
		if n==-1:
			for h in self.sources:
				print h
		else:
			for index in range(n):
				print self.sources[index]
				
	def addSkyCoords(self):
		for s in self.sources:
			ra = s['ra']
			dec = s['dec']
			c = SkyCoord(ra= ra *u.degree, dec = dec * u.degree)
			s['coord'] = c
			
	
	def getColumnNames(self):
		retStr = ""
		for k in self.sources[0].keys():
			retStr+=str(k) + ", "
		return retStr[:-2]
		
	def addTempColumn(self, columnName):
		""" Adds a duplicate column to the source list
		"""
		for s in self.sources:
			s['temp'] = s[columnName]
			
			
	def sortSources(self, column):
		sampleObject = self.sources[0]
		if column not in sampleObject.keys():
			print "Sort: Could not find the column called %s"%column
		self.sources = sorted(self.sources, key=lambda object: object[column], reverse = True)
		
						
	def plotSourceHistogram(self):
		topSourceArray = [ s['mean'] for s in self.HaSources ]
		topSourceArray = sorted(topSourceArray, reverse = True)
		skyArray = [s['mean'] for s in self.skySources ]
		skyArray = sorted(skyArray)
		sourceMedian = numpy.median(topSourceArray)
		sourceMean = numpy.mean(topSourceArray[3:])
			
		skyMedian = numpy.median(skyArray)
		skyMean = numpy.mean(skyArray[3:])
		chart = matplotlib.pyplot.figure("Mean values", figsize=(10, 8))
		axes = matplotlib.pyplot.gca()
		matplotlib.pyplot.clf()
		axes.set_xlabel("Superpixel rank")
		axes.set_ylabel("Mean counts")
		chart.add_axes(axes)
		matplotlib.pyplot.scatter(range(0, len(topSourceArray)), topSourceArray, color='r')
		trimmedSources = [t['mean'] for t in self.trimmedSources]
		trimmedSky = [t['mean'] for t in self.trimmedSky]
		matplotlib.pyplot.scatter(range(0, len(trimmedSources)), trimmedSources, color='g')
		matplotlib.pyplot.scatter(range(0, len(skyArray)), skyArray, color='b')
		matplotlib.pyplot.scatter(range(0, len(trimmedSky)), trimmedSky, color='y')
		matplotlib.pyplot.plot([0, len(topSourceArray)], [sourceMedian, sourceMedian], color='r', ls='dashed')
		matplotlib.pyplot.plot([0, len(topSourceArray)], [sourceMean, sourceMean], color='r', ls='dotted')
		matplotlib.pyplot.plot([0, len(skyArray)], [skyMedian, skyMedian], color='b', ls='dashed')
		matplotlib.pyplot.plot([0, len(skyArray)], [skyMean, skyMean], color='b', ls='dotted')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.pause(10)
		
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
		cols.append(fits.Column(name='sky_mean', format = 'E', array = [o['sky_mean'] for o in objects]))
		cols.append(fits.Column(name='r', format = 'E', array = [o['r'] for o in objects]))
		cols = fits.ColDefs(cols)
		tbhdu = fits.BinTableHDU.from_columns(cols)
			
		prihdr = fits.Header()
		prihdr['COMMENT'] = "Created by Hagrid (postProcess.py) on %s."%( datetime.datetime.ctime(datetime.datetime.now()))
			
		prihdu = fits.PrimaryHDU(header=prihdr)
		thdulist = fits.HDUList([prihdu, tbhdu])
		thdulist.writeto(filename, overwrite=True)
		
	
