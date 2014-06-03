#!/usr/bin/env python
import graph
from sys import argv

g = graph.Graph()
g.read_file(argv[1])
hamiltonian_path = g.christofides(0)
print hamiltonian_path
