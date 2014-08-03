from random                   import randint
from numpy                    import *
from numpy.random             import permutation
from pylogen.util.mtx_search  import relative_frequencies
from pylogen.tree             import MLNode, Edge, UniEdge
import itertools

class Graph:
	
	def __init__(self, nodeNames):
		nNames     = len(nodeNames)
		self.names = nodeNames 
		self.mtx   = numpy.zeros((nNames, nNames), dtype=float)

	def getAdjecents(self, nodeName):
		idx    = getNodeIndex(nodeName)
		colAdj = self.mtx[idx, range(0, idx)]
		rowAdj = self.mtx[(idx+1):, idx]
		adj    = colAdj + rowAdj
		
		return 
		
	def getNodeIndex(self, nodeName):
		return self.mtx.index(nodeName)
		
	def getPairIndices(self, nodeA, nodeB):
		idxA = self.getNodeIndex(nodeA)
		idxB = self.getNodeIndex(nodeB)
		
		# lower triangular form - unidirected graph
		r = max(idxA, idxB)
		c = min(idxA, idxB)
		
		return (r, c)
	
	def getDistance(self, nodeA, nodeB):
		r, c = self.getPairIndices(nodeA, nodeB)
		return self.mtx[r][c]

	def setDistance(self, nodeA, nodeB, d)
		r, c = self.getPairIndices(nodeA, nodeB)
		self.mtx[r][c] = d
		
	def getNodeStateLikelihood(self, nodeName, state):
		adjecents = self.getAdjecents(nodeName)
		nodeIdx   = self.names[self.getNodeIndex(nodeName)]
			
		if len(adjecents) == 1: # a leaf
			return self.states[nodeIdx] == state ? 1.0 : 0.0
		else:
			for b in range(0, 3): # over all bases
				
			# for all possible states s_i, s_j:
				# adj_i 
				# adj_j
				# P(state, s_i), P(state, s_j)
				# L_si(i), L_sj(j)
		
	def getGraphLikelihood(self):
		# choose an edge to break. let it be the longest one THAT DOES NOT CONNECT A LEAF @TODO
		lowerTriangular = numpy.tril(self.mtx, -1)
		aIdx, bIdx = numpy.unravel_index(numpy.argmax(self.mtx), self.mtx.shape)

			
		

def max_likelihood(sequences, names):
	
	nOTU      = len(names)
	K         = len(sequences[0])
	base_freq = relative_frequencies(sequences, ["A", "C", "G", "T"])
	
	q   = lambda v: exp(-v)
	p   = lambda v: 1 - q(v)
	v_i = lambda q: -log(q)
	
	def build_trees(perms):
		# all possible trees
		perms     = list(itertools.permutations(range(0, nOTU), nOTU))
		curr_perm = 0
		
		while curr_perm < len(perms):			
			edges = []
			otu1  = MLNode(names[perms[curr_perm][0]])
			otu2  = MLNode(names[perms[curr_perm][1]])
			
			#p0 = p(0)
			q0 = q(0)						
			edges.append(UniEdge(otu1, otu2, v_i(q0)))
			
			for i in range(2, nOTU):
				# find an edge to break - the one, whose resulting tree gives the highest likelihood
				for e in edges:
					t_edges = 

			
			curr_perm = curr_perm + 1
		
	#tnames = names
	# initial tree
	# first node...
	#otu1Idx = randint(0, nOTU)   
	#otu1    = MLNode(names[otu1Idx])
	#tnames.pop(otu1Idx)
	
	# second node...
	#otu2Idx = randint(0, len(tnames))
	#otu2    = MLNode(names[otu2Idx])
	#tnames.pop(otu2Idx)
	
	# connect them.
	#p0  = p(0)
	#q0  = q(0) 
	#e  = UniEdge(otu1, otu2, v_i(q0))
	#edges = [e]
	
	while len(tnames) > 0:
		# pick an OTU
		otuIdx = randint(0, len(tnames))
		
		# choose an edge to break
		for e in edges:
			
		
