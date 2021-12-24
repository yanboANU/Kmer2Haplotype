#########################################################################
# File Name: read_matrix_draw_graph.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Thu 20 May 2021 10:09:10 PM AEST
#########################################################################
#!/bin/bash
import os
import sys

import matrix
import myGraph
import networkx as nx

from networkx.drawing.nx_pydot import write_dot

edges  = matrix.load_all_edges_from_matrix2( sys.argv[1] )
edgesList = matrix.filter_edges_according_weight(edges, 0)
originalG = myGraph.myGraph()
originalG.add_weighted_edges_from(edgesList)
totalScore = originalG.displayGraph()
#originalG.add_type_A_edges()
print ("After add A edges, number of components: ", nx.number_connected_components(originalG))
    
write_dot(originalG, "chr1_chr2" + ".dot")
