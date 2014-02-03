"""
 Implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean) 
 method for building phylogenetic tree 

 (c) 2014 Urban Soban <u.soban@gmail.com>
 
 Contributors:
	Primoz Turnsek <primoz.turnsek@gmail.com>
"""
from pylogen.util.data       import import_mega_csv
from pylogen.util.mtx_search import find_min 
from pylogen.tree            import Edge, UPGMANode
from numpy                   import *
import sys

def upgma(matrixFilePath):
	
	def recompute_matrix(dimMtx, fst, snd):
		"""
		well, could (should, must) be done better.
		"""
		w, h   = dimMtx.shape
		newMtx = zeros((w-1, h-1))
		newMtx[newMtx == 0] = float("inf")
		
		newMtxRowIdx = 0
		
		for i in range(0, h):
			newMtxColIdx = 0

			if i == fst: 
				# skip row of first matching item (this item is eliminated)
				continue
			elif i == snd and i < 2:
				# skip for of second matching item, and increse row counter of new matrix (this item represents combined item)
				# no need to compute distances (lower triangular form)			
				newMtxRowIdx = newMtxRowIdx + 1
				continue
				
			combined = i == snd
			
			for j in range(0, i):
				if j == fst:
					# skip column of first matching item (this item is eliminated)
					continue
					
				elif j == snd:	
					# compute average - dimMtx not symmetric, so watch for Infs!
					if (isfinite(dimMtx[i][fst])):
						fstVal = dimMtx[i][fst]
					else:
						fstVal = dimMtx[fst][i]
						
					if (isfinite(dimMtx[i][snd])):
						sndVal = dimMtx[i][snd]
					else:
						sndVal = dimMtx[snd][i]
						 
					newMtx[newMtxRowIdx][newMtxColIdx] = (fstVal + sndVal) / 2
				else:
					newMtx[newMtxRowIdx][newMtxColIdx] = dimMtx[i][j]
				
				newMtxColIdx = newMtxColIdx + 1
					
			newMtxRowIdx = newMtxRowIdx + 1
				 
		return newMtx
	
	names, distMtx = import_mega_csv(matrixFilePath) 
	nodes          = {}
	root           = None

	while any(isfinite(distMtx)): 
		match   = find_min(distMtx)
		dist    = distMtx[match[0]][match[1]]
		distMtx = recompute_matrix(distMtx, match[1], match[0]) 

		matchNames      = [names[match[1]], names[match[0]]]
		ancestorName    = "[%s + %s]" % (matchNames[0], matchNames[1])
		names[match[0]] = ancestorName
		names.pop(match[1]) 
		
		commonAncestor = UPGMANode(ancestorName)
		if matchNames[0] in nodes:
			fstNode = nodes[matchNames[0]]
		else:
			fstNode = UPGMANode(matchNames[0])
		
		if matchNames[1] in nodes:
			sndNode = nodes[matchNames[1]]
		else:
			sndNode = UPGMANode(matchNames[1])
			
		fstEdge = Edge(commonAncestor, fstNode, dist/2)
		sndEdge = Edge(commonAncestor, sndNode, dist/2)
		
		nodes[ancestorName] = commonAncestor
		root = commonAncestor

	return root	
