from numpy                   import *
from pylogen.tree            import Edge, NeighborJoiningNode
from pylogen.util.data       import import_mega_csv
from pylogen.util.mtx_search import find_min

def neighbor_joining(matrixFilePath):

	def compute_s_measure(dMtx):
		"""
		Computes 'S' measure matrix. 'S' measures distance from a given 
		OTU (operational taxonomic unit) to all other OTUs.
		"""
		w, h     = dMtx.shape
		nspecies = w
		measures = []
		
		for i in range(0, nspecies):
			# matrix not symmetric, gather OTUs row and column
			#                             row      U     column
			distances = concatenate((dMtx[i].ravel(), dMtx[:, [i]].ravel()))
			distances = [e for e in distances if isfinite(e)]
			measures.append( sum(distances)/(nspecies-2) )
			
		return measures
		
	def compute_m_matrix(dMtx, sMeasure):
		"""
		Computes 'M' measure matrix. 'M' measures distance between pairs
		using formula M(i, j) = D(i,j) - S(i) - S(j)
		"""
		w, h = dMtx.shape
		mMtx = zeros(dMtx.shape)
		mMtx[mMtx == 0] = float("inf")
		
		for i in range(0, h):
			for j in range(0, w):
				if isfinite(dMtx[i][j]):
					mMtx[i][j] = dMtx[i][j] - sMeasure[i] - sMeasure[j]
		
		return mMtx

	def recompute_d_matrix(dMtx, bestPair):
		"""
		Recomputes a distance matrix
		"""
		w, h   = dMtx.shape
		newMtx = delete(dMtx, bestPair[0], 0)      # delete a row
		newMtx = delete(newMtx, bestPair[0], 1)    # delete a col
		
		# correct column
		for i in range(bestPair[1]+1, h-1):
			d_idx = i+1 if bestPair[0] <= i else i # correct OTU index
			
			d_ik = dMtx[d_idx][bestPair[0]] if isfinite(dMtx[d_idx][bestPair[0]]) else dMtx[bestPair[0]][d_idx] 
			d_jk = dMtx[d_idx][bestPair[1]] if isfinite(dMtx[d_idx][bestPair[1]]) else dMtx[bestPair[1]][d_idx]
			d_ij = dMtx[bestPair[0]][bestPair[1]] if isfinite(dMtx[bestPair[0]][bestPair[1]]) else dMtx[bestPair[1]][bestPair[0]]
			
			newMtx[i][bestPair[1]] = (d_ik + d_jk - d_ij) / 2
		
		# correct row
		for j in range(0, bestPair[1]+1):
			d_idx = j+1 if bestPair[0] <= j else j # correct OTU index
			
			d_ik = dMtx[d_idx][bestPair[0]] if isfinite(dMtx[d_idx][bestPair[0]]) else dMtx[bestPair[0]][d_idx] 
			d_jk = dMtx[d_idx][bestPair[1]] if isfinite(dMtx[d_idx][bestPair[1]]) else dMtx[bestPair[1]][d_idx]
			d_ij = dMtx[bestPair[0]][bestPair[1]] if isfinite(dMtx[bestPair[0]][bestPair[1]]) else dMtx[bestPair[1]][bestPair[0]]
			
			newMtx[bestPair[1]][j] = (d_ik + d_jk - d_ij) / 2
		
		return newMtx


	names, dMtx = import_mega_csv(matrixFilePath)
	nodes = {}
	root  = None

	while any(isfinite(dMtx)):
		if dMtx.shape[0] > 2:
			sMeasure = compute_s_measure(dMtx)
			mMtx     = compute_m_matrix(dMtx, sMeasure)
			minPair  = find_min(mMtx)
			
			"""
			match[0] is max. index, match[1] is min. index.
			recomputation of distance matrix will remove match[0], that is 
			row and col with max. index of the best pair. To retain node 
			names in correct order, we remove the node name in max. index, 
			and rename the node name in min. index to ancestor's name 
			"""
			matchNames      = (names[minPair[0]], names[minPair[1]])
			ancestorName    = "[%s + %s]" % matchNames
			names[minPair[1]] = ancestorName
			names.pop(minPair[0]) 
			
			commonAncestor = NeighborJoiningNode(ancestorName)
			if matchNames[0] in nodes:
				sndNode = nodes[matchNames[0]]
			else:
				sndNode = NeighborJoiningNode(matchNames[0])
			
			if matchNames[1] in nodes:
				fstNode = nodes[matchNames[1]]
			else:
				fstNode = NeighborJoiningNode(matchNames[1])

			d_ij    = dMtx[minPair[0]][minPair[1]]
			s_i     = sMeasure[minPair[1]]
			s_j     = sMeasure[minPair[0]]
			
			fstEdge = Edge(commonAncestor, fstNode, 0.5*d_ij + 0.5 * (s_i - s_j))
			sndEdge = Edge(commonAncestor, sndNode, 0.5*d_ij + 0.5 * (s_j - s_i))
			nodes[ancestorName] = commonAncestor
			root = commonAncestor
			
			dMtx = recompute_d_matrix(dMtx, minPair)
		else:
			d              = dMtx[1][0]
			matchNames     = (names[0], names[1])
			ancestorName   = "[%s + %s]" % matchNames
			commonAncestor = NeighborJoiningNode(ancestorName)  
			
			if matchNames[0] in nodes:
				sndNode = nodes[matchNames[0]]
			else:
				sndNode = NeighborJoiningNode(matchNames[0])
				
			if matchNames[1] in nodes:
				fstNode = nodes[matchNames[1]]
			else:
				fstNode = NeighborJoiningNode(matchNames[1])
			
			f = Edge(commonAncestor, fstNode, d)
			s = Edge(commonAncestor, sndNode, d)
			root = commonAncestor
			break

	return root
