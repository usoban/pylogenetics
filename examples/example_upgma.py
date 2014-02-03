import sys
sys.path.append('../')

from pylogen.treebuilder.upgma import upgma
from pylogen.util.data         import import_mega_csv

if len(sys.argv) > 1:
	names, distMtx = import_mega_csv(sys.argv[1])
	
	print upgma(distMtx, names) 
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)
