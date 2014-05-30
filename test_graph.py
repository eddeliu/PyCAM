#!/usr/bin/env python
import graph

g = graph.Graph()
g.read_file('6.txt')
#edges = g.goemans(lambda f1,f2 : 1) # spanning tree
edges = g.goemans(lambda f1,f2 : (f1+f2)%2) # perfect matching
g.display_selected_edges(edges)
