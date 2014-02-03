from numpy import zeros, log

def _p_distance(seq1, seq2):
	bases    = ["A", "G", "T", "C"]
	compared = 0
	matching = 0
	
	for i in range(0, len(seq1)):
		if seq1[i] in bases and seq2[i] in bases:
			compared = compared + 1
			if seq1[i] == seq2[i]: matching = matching + 1
	
	return 1 - float(matching)/float(compared)
	
def jc_distance(sequences):
	"""
	Computes Jukes-Cantor distance between given sequences and returns
	estimates in 2D array of lower-triangular form
	"""
	
	otus          = len(sequences)
	mtx           = zeros((otus, otus))
	mtx[mtx == 0] = float("inf")
	
	for i in range(0, otus):
		for j in range(0, i):
			pDistance = _p_distance(sequences[i], sequences[j])
			mtx[i][j] = -(3.0/4.0) * log(1 - (4.0/3.0)*pDistance)
			
	return mtx
