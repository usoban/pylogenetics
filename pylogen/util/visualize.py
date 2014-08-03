import matplotlib.pyplot as plt
import networkx as nx

def visualize_graph(graphMatrix, verticeNames):
	"""
	Visualizes undirected graph (unrooted phylogenetic tree).
	"""
	
	n_tips = len(verticeNames)
	G      = nx.Graph()
	h, w   = graphMatrix.shape
	edge_labels = {}
	node_labels = {}	
	
	for ri in range(0, h):
		if ri < n_tips: node_labels[ri] = verticeNames[ri]
		else:           node_labels[ri] = ri+1
		
		for ci in range(0, ri):
			if graphMatrix[ri][ci] > 0:
				G.add_edge(ri, ci, distance = graphMatrix[ri][ci])
				edge_labels[(ri, ci)] = graphMatrix[ri][ci]

	pos = nx.graphviz_layout(G, prog='dot')

	# tip nodes
	nx.draw_networkx_nodes(G,pos, weight = 'distance', nodelist=range(0, len(verticeNames)), node_color='b', node_size = 900)
	nx.draw_networkx_nodes(G,pos, weight = 'distance', nodelist=range(len(verticeNames), w), node_color='r', node_size = 450)

	# edges
	nx.draw_networkx_edges(G, pos)

	# labels
	nx.draw_networkx_labels(G, pos, font_size=15, font_family='sans-serif', labels = node_labels)
	nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels)

	plt.axis('off')
	plt.savefig("weighted_graph.png") # save as png

def visualize_tree(root):
	G = nx.Graph()
	root.xnet(G)

	pos = nx.graphviz_layout(G, prog='dot', root=root.label)
	
	nx.draw_networkx_nodes(G, pos)
	nx.draw_networkx_edges(G, pos)
	
	plt.axis('off')
	plt.savefig('neighbor_joining_full.png')
