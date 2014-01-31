"""
 Implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean) 
 method for building phylogenetic tree 

 (c) 2014 Urban Soban <u.soban@gmail.com>
 
 Contributors:
	Primoz Turnsek <primoz.turnsek@gmail.com>
"""
from numpy  import *
import csv, sys

###
# Tree classes
###
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
			sharedDist   = self.children[0].distance
			
			for e in self.children:
				childrenReps.append(str(e.endNode))
			
			childrenReps = ','.join(childrenReps)
			
			if self.parentEdge is None:
				return '(%s);' % childrenReps
			else:
				return '(%s):%f' % (childrenReps, self.parentEdge.distance - sharedDist)
		else:
			return "%s:%f" % (self.label, self.parentEdge.distance)

###
# UPGMA procedures
###
def find_best_match(dimMtx):
	flatMtx = dimMtx.ravel() 
	minIdx  = argmin(dimMtx)
		
	return unravel_index(minIdx, dimMtx.shape)
	
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

###
# UPGMA & Tree building
###
nodes   = {}
root    = PhyloNode(None)
root    = None

while any(isfinite(distMtx)): 
	match   = find_best_match(distMtx)
	dist    = distMtx[match[0]][match[1]]
	distMtx = recompute_matrix(distMtx, match[1], match[0]) 

	matchNames      = [names[match[1]], names[match[0]]]
	ancestorName    = "[%s + %s]" % (matchNames[0], matchNames[1])
	names[match[0]] = ancestorName
	names.pop(match[1]) 
	
	commonAncestor = PhyloNode(ancestorName)
	if matchNames[0] in nodes:
		fstNode = nodes[matchNames[0]]
	else:
		fstNode = PhyloNode(matchNames[0])
	
	if matchNames[1] in nodes:
		sndNode = nodes[matchNames[1]]
	else:
		sndNode = PhyloNode(matchNames[1])
		
	fstEdge = Edge(commonAncestor, fstNode, dist/2)
	sndEdge = Edge(commonAncestor, sndNode, dist/2)
	
	nodes[ancestorName] = commonAncestor
	root = commonAncestor

print root	
