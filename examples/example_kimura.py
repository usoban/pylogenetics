import sys
sys.path.append('../')

from pylogen.distance.kimura import kimura_distance
from pylogen.util.data       import import_fasta_sequences

if len(sys.argv) > 1:
	names, sequences = import_fasta_sequences(sys.argv[1])
	
	print ','.join(names)
	print kimura_distance(sequences)
else:
	print "Please supply sequences file in FASTA format as first argument"
	exit(1)


