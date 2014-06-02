#!/usr/bin/env python
import graph

g = graph.Graph()
g.read_file('6.txt')
edges = g.solve_spanning_tree()
g.display_selected_edges(edges, "spanning")
print (("*"*60)+"\n")*15
edges = g.solve_perfect_matching()
g.display_selected_edges(edges, "matching")
