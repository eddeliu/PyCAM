#!/usr/bin/env python
import graph
from sys import argv

g = graph.Graph()
g.read_file(argv[1])
edges = g.solve_spanning_tree()
g.display_selected_edges(edges, "spanning")
print (("*"*60)+"\n")*15
edges = g.solve_perfect_matching()
g.display_selected_edges(edges, "matching")
