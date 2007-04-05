
srcfile="../msc/GPLv2.txt"
import sys

def printfile():
	try:
		file = open(srcfile, 'r')
	except IOError:
		sys.exit(0)
	print file
	for line in file.readlines():
		if line.endswith("\n"): 
			line=line.replace("\n","")
			print line

