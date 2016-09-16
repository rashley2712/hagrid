#!/usr/bin/env python

import os, sys, json, argparse, shutil, re

class FITSfileObject:
	def __init__(self):
			filename = "unknown"
			processed = False
			
	def __str__(self):
		return "filename: %s"%self.filename
		
	

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='A simple Python tool to run hagrid in a batch mode over a folder full of files.')
	parser.add_argument('-p', '--path', type=str, help='The folder in which the IPHAS images are contained.')
	parser.add_argument('-w', '--workingpath', type=str, help='A root folder for the temporary and output files.')
	arg = parser.parse_args()
	print arg

	if arg.path is None:
		dataPath="."
	else: 
		dataPath= arg.path
	if arg.workingpath is None:
		workingPath="."
	else: 
		workingPata= arg.workingpath
		
	
	# Get a list of files in the folder
	# First, check if the source data is there
	if not os.path.exists(dataPath):
		print "The folder for the source data %s could not be found. Exiting."%dataPath
		sys.exit()
	
	searchString = ".*.(fits|fits.gz|fits.fz|fit)"
	search_re = re.compile(searchString)
		 
	(_, _, filenames) = os.walk(dataPath).next()
	FITSFilenames = []
	for file in filenames:
		m = search_re.match(file)
		if (m): 
			FITSFilenames.append(file)

	FITSFilenames = sorted(FITSFilenames)
	print FITSFilenames
	
	# Second, check to see if the working directory already exists
	if not os.path.exists(workingPath):
		debug("Creating folder %s"%workingPath)
		os.makedirs(workingPath)
		
	
