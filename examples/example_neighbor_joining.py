import sys
sys.path.append('../')

from pylogen.treebuilder.nj import neighbor_joining

if len(sys.argv) > 1:
	print neighbor_joining(sys.argv[1]) 
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)
