import sys
sys.path.append('../')

from pylogen.util.visualize import visualize_graph
from numpy                  import *

mtx = array([
[0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 1],
[0, 0, 0, 0, 0, 0, 1, 0],
[0, 0, 0, 0, 0, 0, 1, 0],
[1, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 1, 0, 0, 5],
[0, 0, 1, 0, 0, 1, 5, 0]
])

visualize_graph(mtx, ["ab", "cd", "ef", "gg", "hf"])
