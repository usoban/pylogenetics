import sys
sys.path.append('../')

from pylogen.util.data             import import_fasta_sequences
from pylogen.distance.jukes_cantor import jc_distance
from pylogen.treebuilder.upgma     import upgma

if len(sys.argv) > 1:
	
	names, sequences = import_fasta_sequences(sys.argv[1])
	distanceMatrix   = jc_distance(sequences)
	
	print upgma(distanceMatrix, names) 
	
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)

