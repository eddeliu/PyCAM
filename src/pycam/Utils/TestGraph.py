#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of PyCAM.

PyCAM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PyCAM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PyCAM.  If not, see <http://www.gnu.org/licenses/>.
"""

if __name__ == "__main__" :
    """A regression test for Christofides implementation

    given a graph (here passed as arg by line command)
    as a matrix of distances between vertexes
    it should output the hamiltonian path corresponding
    to the graph, which should be a 2-approximation"""
    import solvinggraph
    from sys import argv

    g = solvinggraph.SolvingGraph()
    g.read_file(argv[1])
    hamiltonian_path = g.christofides()
    print hamiltonian_path
