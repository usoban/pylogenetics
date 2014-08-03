import sys
sys.path.append('../')

from pylogen.treebuilder.nj import neighbor_joining
from pylogen.util.data      import import_mega_csv
from pylogen.util.visualize import visualize_tree

if len(sys.argv) > 1:

	names, distanceMatrix = import_mega_csv(sys.argv[1])
	
	root = neighbor_joining(distanceMatrix, names) 
	visualize_tree(root)
	
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)

