from abc import ABCMeta, abstractmethod

class Edge:
	"""
	Phylogenetic tree edge. Connects start and end node with some given
	distance measure 
	"""
	
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

class Node:
	"""
	Generic phylogenetic tree node. 
	
	Has a label, list of edges that connect it to its children, and a parent
	edge that connects it to its parent.
	
	Priting a node outputs subtree with root in the given node.
	"""
	
	__metaclass__ = ABCMeta
	
	def __init__(self, label):
		self.label      = label
		self.children   = []
		self.parentEdge = None
		
	def addChildEdge(self, childEdge):
		self.children.append(childEdge)
		
	def setParentEdge(self, parentEdge):
		self.parentEdge = parentEdge

	@abstractmethod
	def output(self):
		pass

	def __str__(self):
		return self.output()

class NeighborJoiningNode(Node):
	"""
	Neighbor joining tree node
	"""
	
	def output(self):
		"""
		Ouptuts phylogenetic subtree in NEWICK format
		"""
		if len(self.children) > 0:
			childrenReps = ','.join([str(e.endNode) for e in self.children])
						
			if self.parentEdge is None:
				return '(%s);' % childrenReps
			else:
				return '(%s):%f' % (childrenReps, self.parentEdge.distance)
		else:
			return "%s:%f" % (self.label, self.parentEdge.distance)
			
class UPGMANode(Node):
	"""
	UPGMA tree node
	"""
	
	def output(self):
		"""
		Ouptuts phylogenetic subtree in NEWICK format
		"""
		if len(self.children) > 0:
			childrenReps = ','.join([str(e.endNode) for e in self.children])
			sharedDist   = self.children[0].distance
			
			if self.parentEdge is None:
				return '(%s);' % childrenReps
			else:
				return '(%s):%f' % (childrenReps, self.parentEdge.distance - sharedDist)
		else:
			return "%s:%f" % (self.label, self.parentEdge.distance)
