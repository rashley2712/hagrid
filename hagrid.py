#!/usr/bin/env python


import argparse
import sys

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='[Ha grid] A tool for loading, manipulating and plotting IPHAS reduced images in order to find and export Halpha bright regions.')
	parser.add_argument('script', type=str, nargs='?', help='The name of a script file containing commands to execute.')
	parser.add_argument('--noplot', action='store_true', help='Set mathplotlib backend to Agg.')
	arg = parser.parse_args()

        if arg.noplot:
          print "noplot",arg.noplot
          # Force matplotlib to not use any Xwindows backend.
          import matplotlib
          matplotlib.use('Agg')
	
        # this import has to happen AFTER the noplot check
        import commandsIPHAS
	commands = commandsIPHAS.commandClass

	if arg.script is not None:
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
