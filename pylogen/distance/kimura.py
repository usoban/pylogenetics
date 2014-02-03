from numpy import zeros, log

def _compare(seq1, seq2):
	"""
	Compares two sequences and counts transitions and transversions
	"""
	purines     = ["A", "G"]
	pyrimidines = ["C", "T"]
	bases       = purines + pyrimidines
		
	compared      = 0
	transitions   = 0
	transversions = 0
	for i in range(0, len(seq1)):		
		if seq1[i] in bases and seq2[i] in bases:
			compared = compared + 1			
			if seq1[i] != seq2[i]:
				if seq1[i] in purines and seq2[i] in purines or seq1[i] in pyrimidines and seq2[i] in pyrimidines: 
					transitions = transitions + 1
				else:
					transversions = transversions + 1
	
	return (float(transitions)/float(compared), float(transversions)/float(compared))
		
def kimura_distance(sequences):
	"""
	Computes Kimura (2 parameter) distance between given sequences and returns
	estimates in 2D array of lower-triangular form
	"""
	otus          = len(sequences)
	mtx           = zeros((otus, otus))
	mtx[mtx == 0] = float("inf")

	for i in range(0, otus):
		for j in range(0, i):
			p, q = _compare(sequences[i], sequences[j])
			mtx[i][j] = -0.5 * log(1 - 2 * p - q) - 0.25 * log(1 - 2 * q)
			
	return mtx
