#!/usr/bin/env python

import commandsIPHAS
import argparse
import sys


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='[Ha grid] A tool for loading, manipulating and plotting IPHAS reduced images in order to find and export Halpha bright regions.')
	parser.add_argument('script', type=str, nargs='?', help='The name of a script file containing commands to execute.')
	arg = parser.parse_args()
	
	commands = commandsIPHAS.commandClass
	
	if arg.script != None:
		if arg.script=="restore":
			commands().do_restore("")
		else: 
			input = open(arg.script, 'rt')
			print "Running the commands found in :", arg.script
			try:
				commands(stdin=input).cmdloop()
			finally:
				input.close()
	commands.prompt = 'hagrid> '
	commands.use_rawinput = True
	commands().cmdloop()

	sys.exit()
