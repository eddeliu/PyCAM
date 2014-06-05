#!/usr/bin/env python
import solvinggraph
from sys import argv

g = solvinggraph.SolvingGraph()
g.read_file(argv[1])
hamiltonian_path = g.christofides()
print hamiltonian_path
