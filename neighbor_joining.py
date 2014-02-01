from numpy import *
import math, sys, csv

class Edge:
	
	def __init__(self, startNode, endNode, distance):
		self.distance  = distance

		self.setEndNode(endNode)
		self.setStartNode(startNode)

	def setStartNode(self, node):
		self.startNode = node
		self.startNode.addChildEdge(self)
		
	def setEndNode(self, node):
		self.endNode = node
		self.endNode.setParentEdge(self)

class PhyloNode:
	
	def __init__(self, label):
		self.label      = label
		self.children   = []
		self.parentEdge = None
		
	def addChildEdge(self, childEdge):
		self.children.append(childEdge)
		
	def setParentEdge(self, parentEdge):
		self.parentEdge = parentEdge

	def __str__(self):
		"""
		Ouptuts phylogenetic tree in Newick format
		"""
		if len(self.children) > 0:
			childrenReps = []
			
			for e in self.children:
				childrenReps.append(str(e.endNode))
			
			childrenReps = ','.join(childrenReps)
			
			if self.parentEdge is None:
				return '(%s);' % childrenReps
			else:
				return '(%s):%f' % (childrenReps, self.parentEdge.distance)
		else:
			return "%s:%f" % (self.label, self.parentEdge.distance)


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
	
def find_best_match(mMtx):
	"""
	Finds a pair with minimum M-measure value
	"""
	minIdx  = argmin(mMtx)
	return unravel_index(minIdx, mMtx.shape)

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
		#print "correcting column %d, in row %d with value %f  [d_idx was %d]" % (bestPair[1], i, newMtx[i][bestPair[1]], d_idx)
	
	# correct row
	for j in range(0, bestPair[1]+1):
		d_idx = j+1 if bestPair[0] <= j else j # correct OTU index
		
		d_ik = dMtx[d_idx][bestPair[0]] if isfinite(dMtx[d_idx][bestPair[0]]) else dMtx[bestPair[0]][d_idx] 
		d_jk = dMtx[d_idx][bestPair[1]] if isfinite(dMtx[d_idx][bestPair[1]]) else dMtx[bestPair[1]][d_idx]
		d_ij = dMtx[bestPair[0]][bestPair[1]] if isfinite(dMtx[bestPair[0]][bestPair[1]]) else dMtx[bestPair[1]][bestPair[0]]
		
		newMtx[bestPair[1]][j] = (d_ik + d_jk - d_ij) / 2
	
	return newMtx

###
# Data preparation
###
if len(sys.argv) > 1:
	csvFile = open(sys.argv[1])
	data    = csv.reader(csvFile, delimiter = ',')
	names   = []
	lines   = []
	
	for row in data:
		if len(row) > 0:
			names.append(row[0])
			lines.append(row[1:])
		else:
			break

	# First col is specie name
	nspecies = len(lines[-1])
	distMtx  = zeros((nspecies, nspecies))
	distMtx[distMtx==0] = float('inf')
	
	for i in range(0, len(lines)):
		
		for j in range(0, len(lines[i])):
			try:
				distMtx[i][j] = float(lines[i][j])
			except ValueError:
				break

else:
	print "Please supply distance matrix as CSV file"
	exit(1)

dMtx  = distMtx
nodes = {}
root  = None

while any(isfinite(dMtx)):
	if dMtx.shape[0] > 2:
		sMeasure = compute_s_measure(dMtx)
		mMtx     = compute_m_matrix(dMtx, sMeasure)
		minPair  = find_best_match(mMtx)
		
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
		
		commonAncestor = PhyloNode(ancestorName)
		if matchNames[0] in nodes:
			sndNode = nodes[matchNames[0]]
		else:
			sndNode = PhyloNode(matchNames[0])
		
		if matchNames[1] in nodes:
			fstNode = nodes[matchNames[1]]
		else:
			fstNode = PhyloNode(matchNames[1])

		d_ij    = dMtx[minPair[0]][minPair[1]]
		s_i     = sMeasure[minPair[1]]
		s_j     = sMeasure[minPair[0]]
		
		fstEdge = Edge(commonAncestor, fstNode, 0.5*d_ij + 0.5 * (s_i - s_j))
		sndEdge = Edge(commonAncestor, sndNode, 0.5*d_ij + 0.5 * (s_j - s_i))
		nodes[ancestorName] = commonAncestor
		root = commonAncestor
		
		dMtx     = recompute_d_matrix(dMtx, minPair)
	else:
		print dMtx
		print(names)
		d              = dMtx[1][0]
		matchNames     = (names[0], names[1])
		ancestorName   = "[%s + %s]" % matchNames
		commonAncestor = PhyloNode(ancestorName)  
		
		if matchNames[0] in nodes:
			sndNode = nodes[matchNames[0]]
		else:
			sndNode = PhyloNode(matchNames[0])
			
		if matchNames[1] in nodes:
			fstNode = nodes[matchNames[1]]
		else:
			fstNode = PhyloNode(matchNames[1])
		
		f = Edge(commonAncestor, fstNode, d)
		s = Edge(commonAncestor, sndNode, d)
		root = commonAncestor
		break

print root
