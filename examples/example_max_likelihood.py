import sys
sys.path.append('../')

from pylogen.util.data             import import_fasta_sequences
#from pylogen.distance.jukes_cantor import jc_distance
from pylogen.treebuilder.ml        import max_likelihood

if len(sys.argv) > 1:
	
	names, sequences = import_fasta_sequences(sys.argv[1])

	print max_likelihood(sequences, names) 
	
else:
	print "Please supply distance matrix CSV file as first argument"
	exit(1)

