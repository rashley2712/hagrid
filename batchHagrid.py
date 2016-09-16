import os, sys, json

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='A simple Python tool to run hagrid in a batch mode over a folder full of files.')
	parser.add_argument('-p', '--path', type=str, help='The folder in which the IPHAS images are contained.')
	arg = parser.parse_args()
	
	print arg
