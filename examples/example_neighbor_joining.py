import sys
sys.path.append('../')

from pylogen.treebuilder.nj import neighbor_joining
from pylogen.util.data      import import_mega_csv

if len(sys.argv) > 1:

	names, distanceMatrix = import_mega_csv(sys.argv[1])
	
	print neighbor_joining(distanceMatrix, names) 
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)
