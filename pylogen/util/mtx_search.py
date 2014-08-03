from numpy import argmin, unravel_index

def find_min(mtx):
	"""
	Finds first minimim element of matrix and returns its index in form of
	tuple (row_idx, column_idx)
	"""
	minIdx  = argmin(mtx)
	return unravel_index(minIdx, mtx.shape)

def relative_frequencies(sequences, chars):
	"""
	Counts number of char occurances in given sequences and returns their
	relative frequencies in a dict 
	"""
	total = 0 
	freq  = {}
	# initialize dictionary
	for c in chars: freq[c] = 0
	
	for seq in sequences:
		for i in range(0, len(seq)):
			if seq[i] in chars: 
				freq[seq[i]] = freq[seq[i]] + 1
				total = total + 1
	
	# make relative
	for c in chars: freq[c] = freq[c]/float(total)
	
	return freq
	
	
