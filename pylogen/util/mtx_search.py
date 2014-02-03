from numpy import argmin, unravel_index

def find_min(mtx):
	"""
	Finds first minimim element of matrix and returns its index in form of
	tuple (row_idx, column_idx)
	"""
	minIdx  = argmin(mtx)
	return unravel_index(minIdx, mtx.shape)
