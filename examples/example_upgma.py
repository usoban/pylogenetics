import sys
sys.path.append('../')

from pylogen.treebuilder.upgma import upgma

if len(sys.argv) > 1:
	print upgma(sys.argv[1]) 
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)
