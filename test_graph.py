#!/usr/bin/env python
import graph
from sys import argv

g = graph.Graph()
g.read_file(argv[1])
g.christofides()
