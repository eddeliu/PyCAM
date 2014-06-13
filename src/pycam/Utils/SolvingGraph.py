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

from numpy import empty as numpy_empty
from numpy import int as numpy_int
try :
    from pycam.Utils.VariablePriorityQueue import VariablePriorityQueue
    from pycam.Utils.DisjointSet import DisjointSet
    import pycam.Utils.log
    log = pycam.Utils.log.get_logger()
except ImportError :
    # we are not inside of PyCAM : it used run as a test
    class Foo(object) :
        def info(self, s) :
            print s
    log = Foo()
    import VariablePriorityQueue
    import DisjointSet

class SolvingGraph(object) :
    """A class that holds an implementation of Christofides' algorithm

    Given a graph (its matrix N*N holding distances between vertices)
    this class output an hamiltonian path through this graph.
    Christofides' algorithm produce a 3/2-approximation,
    but it is actually a 2-approximation cause the perfect matching
    sub-problem is itself a 2-approximation.
    For both perfect matching and spanning tree problems we use the
    general algorithm from Goemans et Williamson (see references)."""

    def __init__(self, matrix=None) :
        """Initialize the Solver

        If no matrix is yield, it should be next read from file
        using SolvingGraph().read_file"""
        if matrix != None :
            self.size = len(matrix)
            self.matrix = matrix

    def christofides(self) :
        """Computes the hailtonian path for the graph

        It computes the minimum weighted spanning tree for the graph,
        then a mimimum weighted perfect matching with 2-approximation
        amongst odd degree vertices of the spanning tree.
        Finally, it computes an eulerian path from sum of the two
        previous results and make it eulerian."""
        if self.size == 1 :
            return [[0,0]] # else algorithm would not find path
        log.info("spanning tree computation begin ...")
        spanning_tree = self.solve_spanning_tree(range(self.size))
        log.info("spanning tree computation ended")
        # we have now to find odd degree vertices
        edges_by_vertex = {}
        for edge in spanning_tree :
            for vertex in edge :
                if vertex not in edges_by_vertex :
                    edges_by_vertex[vertex] = [edge]
                else : edges_by_vertex[vertex].append(edge)
        odd_degree_vertexes = \
            filter(lambda v : len(edges_by_vertex[v])%2, edges_by_vertex.keys())
        log.info("perfect matching computation begin ...")
        perfect_matching = self.solve_perfect_matching(odd_degree_vertexes)
        log.info("perfect matching computation ended")
        for edge in perfect_matching :
            for vertex in edge :
                edges_by_vertex[vertex].append(edge)
        # finally computes the eulerian path
        def find_path(starting_from) :
            path = []
            start = starting_from
            vertexes = edges_by_vertex
            while True :
                edge = vertexes[start][0]
                vertexes[start].remove(edge)
                end = edge[0] if start == edge[1] else edge[1]
                vertexes[end].remove(edge)
                path.append([start, end])
                if end == starting_from : break
                else : start = end
            return path
        eulerian_path = find_path(edges_by_vertex.keys()[0])
        edges_left = len(spanning_tree) + len(perfect_matching) \
                                        - len(eulerian_path)
        while edges_left :
            for edge in eulerian_path :
                if edges_by_vertex[edge[0]] :
                    additional_path = find_path(edge[0])
                    edges_left -= len(additional_path)
                    index_start = eulerian_path.index(edge)
                    eulerian_path = eulerian_path[:index_start] + \
                                additional_path + eulerian_path[index_start:]
        # this final step is to make the eulerian path hamiltonian
        # by removing the vertices visited multiple times,
        # and thus gain time (according to triangle inequality)
        seen = {}
        i = 0
        while i < len(eulerian_path) :
            edge = eulerian_path[i]
            if edge[1] not in seen :
                seen[edge[1]] = edge[0]
                i += 1
            else :
                eulerian_path[i][1] = eulerian_path[(i+1)%len(eulerian_path)][1]
                eulerian_path.pop((i+1)%len(eulerian_path))
        return eulerian_path # which is now hamiltonian

    def goemans(self, cf, vertexes) :
        """A general algorith for soving graph problems

        See 'a general approximation technique for constrained
        forest problems' by Goemans et Williamson.
        We use it both for min spanning tree and min perfect matching.
        It needs a cf function that f(C1 U C2) = cf(f(C1),f(C2)).
        This is explained in the paper (read it, really !)
        We also need to know the vertices we are working on : these are
        the indexes of the vertices in the matrix."""
        F = [] # holds the solution
        class QueueElement(object) :
            """Utility class for our VariablePriorityQueue

            It represents an edge between two vertices (i and j) and
            have a priority (it is aimed to be an element of
            the VariablePriorityQueue).
            The methods are very hacky, and does not make any sense together,
            they are, i repeat, only meant for utility, not meaning."""
            __slots__ = ['i', 'j', 'priority'] # hack to gain speed
            def __init__(self, i, j, priority) :
                self.i = i
                self.j = j
                self.priority = priority
            def __hash__(self) :
                return hash(str(self.i) + ' ' + str(self.j))
            def __eq__(self, other) :
                return hash(self) == hash(other)
            def __ne__(self, other) :
                return (self.i != other.i) or (self.j != other.i)
            def __lt__(self, other) :
                return self.priority < other.priority
            def __str__(self) :
                return "QElem (%s,%s) [%s]" % (self.i, self.j, self.priority)
            def __repr__(self) :
                return "QueueElement(%s,%s, %s)"%(self.i, self.j, self.priority)
        # declare edges priority queue and initialize it
        queue = VariablePriorityQueue()
        for line in vertexes :
            for column in vertexes :
                if line > column :
                    weight = self.matrix.item(line, column)
                    edge = QueueElement(line,column,weight/2)
                    queue.add(edge)

        #define union-find structure and initialize it
        C = DisjointSet()
        f = {} # also store the values of f function
                # (should access it only through nodes representatives)
        d = {} # store nodes weights
        for vertex in vertexes :
            C.add(vertex)
            f[vertex] = 1
            d[vertex] = 0

        # store the number of trees C_r such that f(C_r) = 1
        # we need this info to check quickly for then end of
        # computation (when it reaches 1)
        number_of_f1_trees = len(vertexes)
        
        # this is the function that computes the priority of an edge
        def epsilon(i, j, d, f) :
            sum_f = f[C.find(i)] + f[C.find(j)]
            if sum_f == 0 :
               sum_f = 0.00001 # 1/100.000
            return (self.matrix.item(i, j) -d[i] -d[j])/(sum_f)
        
        # we also need a matrix that holds the best edge between two trees
        # can only be accessed through nodes representatives
        # we encode the edge to be able to store it in an integer
        best_edges = numpy_empty(shape=(self.size,self.size), dtype=numpy_int)
        for line in vertexes :
            for column in vertexes :
                if line != column :
                    best_edges[line][column] = line*self.size + column
                else :
                    best_edges[line][column] = 0

        # initialize lower bound
        LB = 0

        # initialize time
        # in priority queue the real priority we want
        # is epsilon(i,j)
        # however, it is too costly to store it since any event
        # would require modifying all priorities
        # so we instead use a global time T
        # and define priorities as epsilon(i,j) + T
        # see paper for more details !
        T = 0

        print "iteration starting ..."
        # main loop
        while number_of_f1_trees > 1 :
            # take best edge
            best_edge = queue.get()
            if C.find(best_edge.i) == C.find(best_edge.j) :
                continue # skip edges in same tree

            if best_edges.item(best_edge.i, best_edge.j) != \
                            best_edge.i*self.size+best_edge.j :
                continue # skip edges we do not update because not best ones
            # else add it to the solution
            F.append([best_edge.i, best_edge.j])
            e = epsilon(best_edge.i, best_edge.j, d, f)
            # update time
            T += e

            # update d
            for v in vertexes :
                d[v] += e * f[C.find(v)]

            # update lower bound
            for r in C.retrieve() :
                LB += e * f[r]

            # now last and complex part of the algorithm :
            # doing the fusion and weights update
            # starts by updating number of trees with f=1
            r_i = C.find(best_edge.i)
            r_j = C.find(best_edge.j)
            v_i = f[r_i]
            v_j = f[r_j]

            new_f_value = cf(v_i, v_j)

            # compute new amount of trees with f value equel to 1
            if v_i + v_j > new_f_value :
                number_of_f1_trees -= 2 - cf(1, 1) # one or two less
            # now, fuse
            C.union(r_i, r_j)
            # update f
            r_union = C.find(r_i)

            f[r_union] = new_f_value

            # modify best edges going out of r_union tree
            # start with line
            for column in vertexes :
                if (column != r_i) and (column != r_j) :
                    # we have to use min and max to ensure line > column
                    # because we only use the bottom half of the matrix
                    codedge1 = best_edges.item(max(r_i,column), min(r_i,column))
                    codedge2 = best_edges.item(max(r_j,column), min(r_j,column))
                    edge1 = [codedge1 // self.size, codedge1 % self.size]
                    edge2 = [codedge2 // self.size, codedge2 % self.size]
                    e1 = epsilon(edge1[0], edge1[1], d, f)
                    e2 = epsilon(edge2[0], edge2[1], d, f)
                    best_edges[r_union][column] = codedge1 if e1 < e2 \
                                            else codedge2

            # continue with column
            for line in vertexes :
                if (line != r_i) and (line != r_j) :
                    # and same for min and max
                    coded_edge1 = best_edges.item(max(line,r_i), min(line,r_i))
                    coded_edge2 = best_edges.item(max(line,r_j), min(line,r_j))
                    edge1 = [coded_edge1 // self.size, coded_edge1 % self.size]
                    edge2 = [coded_edge2 // self.size, coded_edge2 % self.size]
                    e1 = epsilon(edge1[0], edge1[1], d, f)
                    e2 = epsilon(edge2[0], edge2[1], d, f)
                    best_edges[line][r_union] = coded_edge1 if e1 < e2 \
                                            else coded_edge2

            # final step : modifying priority queue
            # loop on all best edges leaving r_union tree
            # all edges leaving r_union tree which are not best
            # are therefore NOT UPDATED and will have ERRONEOUS priorities
            # this is why we check edges at beginning of the loop
            for column in vertexes :
                if C.find(column) < r_union and \
                        best_edges.item(r_union, column) == column :
                    coded_edge = best_edges.item(r_union, column)
                    edge = [coded_edge // self.size, coded_edge % self.size]
                    priority = epsilon(r_union, column, d, f) + T
                    changing_edge = QueueElement(r_union,column,priority)
                    queue.update(changing_edge)
        return F

    def solve_spanning_tree(self, vertexes) :
        """Solves the spanning tree problem

        It uses Goemans' algorithm to do so,
        and that's pretty all."""
        # the lamba function is easy : whatever the parity
        # of both trees, keep priority=1 for the new one
        return self.goemans(lambda f1,f2 : 1, vertexes)

    def solve_perfect_matching(self, vertexes) :
        """Solves the perfect matching problem

        It uses Goemans' algortihm but do a processing
        on the results as needed (see the paper for why)"""
        # here the lamba function prioritizes fusions
        # between an odd and an even tree, and make
        # fusion between two odd or even ones, having less priority
        F = self.goemans(lambda f1,f2 : (f1+f2)%2, vertexes)

        # create a map that holds edges for each vertex
        edges_by_vertex = {}
        for edge in F :
            for vertex in edge :
                if vertex not in edges_by_vertex :
                    edges_by_vertex[vertex] = [edge]
                else : edges_by_vertex[vertex].append(edge)

        # returns from a graph (represented as the dic created above)
        # the number of components containing an even number of vertices
        def even_sized_connected_component_number(edges_by_vertex, removed_edge):
            if (removed_edge != None) and (edges_by_vertex[removed_edge[0]] == 1 \
                                        or edges_by_vertex[removed_edge[1]] == 1) :
                return None
            seen = {}
            even_sized = 0
            current_component = 1
            for vertex in edges_by_vertex.keys() :
                if vertex not in seen :
                    vertices_left = [vertex]
                    card = 0
                    while vertices_left :
                        v = vertices_left.pop()
                        if v not in seen :
                            seen[v] = current_component
                            card += 1
                            for edge in edges_by_vertex[v] :
                                if edge != removed_edge :
                                    vertices_left.append(edge[0] if edge[1] == v else edge[1])
                    if card%2 == 0 : even_sized += 1
                    current_component += 1
            return even_sized

        # now try to remove one by one each edge
        # for each edge we simply remove it
        # and see on remaining graph if there exists
        # more connected component of odd size
        # if yes, delete this edge else keep it
        current_edge_index = 0
        previous_number = even_sized_connected_component_number(edges_by_vertex, None)
        while current_edge_index < len(F):
            edge = F[current_edge_index]
            new_number = even_sized_connected_component_number(edges_by_vertex, edge)
            if (new_number is None) or (new_number <= previous_number) :
                current_edge_index += 1
            else:
                for vertex in edge :
                    edges_by_vertex[vertex].remove(edge)
                F.remove(edge)
                previous_number = new_number

        # and still one iteration for creating shortcuts
        for vertex in edges_by_vertex.keys() :
            while len(edges_by_vertex[vertex]) > 2 :
                len1_vertexes = []
                len_sup_2_vertexes = []
                # analyze connex vertexes
                for edge in edges_by_vertex[vertex] :
                    num_of_edges = len(edges_by_vertex[edge[0]+edge[1]-vertex]) # math hack : a+b-a = b
                    if num_of_edges == 1 : len1_vertexes.append(edge)
                    elif num_of_edges > 2 : len_sup_2_vertexes.append(edge)
                # find good pair
                if len(len1_vertexes) >= 2 :
                    edge1,edge2 = len1_vertexes[0], len1_vertexes[1]
                elif len(len_sup_2_vertexes) >= 2 :
                    edge1,edge2 = len_sup_2_vertexes[0], len_sup_2_vertexes[1]
                # create shortcut and update edges
                shortcut = [edge1[edge1.index(edge1[0]+edge1[1]-vertex)],
                            edge2[edge2.index(edge2[0]+edge2[1]-vertex)]] # same math hack
                for vertex in shortcut :
                    edges_by_vertex[vertex].append(shortcut)
                for vertex1,vertex2 in zip(edge1,edge2) :
                    edges_by_vertex[vertex1].remove(edge1)
                    edges_by_vertex[vertex2].remove(edge2)
                # and finally update solution
                F.remove(edge1)
                F.remove(edge2)
                F.append(shortcut)

        return F # is now a real perfect matching

    def read_file(self, filename) :
        """Load the matrix from a file

        The file should :
        - exists
        - have on the first line the integer number N on vertices
        - have on the lines 2:N+1 the float distances with all N
          vertices (including itself), each separated by a space
        so it must look like a matrix.
        It is used mainly for debugging or regression test."""
        matrix_file = open(filename, 'r')
        self.size = int(matrix_file.readline().rstrip())
        matrix = numpy_empty(shape=(self.size,self.size))
        current_line_number = 0
        for line in matrix_file :
            strings = line.rstrip().split(' ')
            floats = map(float, strings)
            matrix[current_line_number] = floats
            current_line_number += 1
        self.matrix = matrix

    def display_selected_edges(self, edges, name="goemans") :
        """Create a beatiful DOT file to see what's going on

        It is used for debugging because it is very handy.
        It will output 2 files :
        - name.dot which will be a complete graph,
          but it starts to look awful for N > 10
        - name_sol.dot which will contains just the solution
          and make it simple to check for correctness
        The DOT file can be further processed to produce
        a PNG file which make everything clear, to do so :
        dot -Tpng name.dot -o name.png"""
        with open(name+".dot", 'w') as edges_file :
            with open(name+"_sol.dot", 'w') as sol_file :
                h = {}
                for edge in edges :
                    h[edge[0]*self.size+edge[1]] = 1
                edges_file.write("graph {\n")
                sol_file.write("graph {\n")
                for line in xrange(self.size) :
                    for column in xrange(self.size) :
                        if line != column :
                            if line*self.size+column in h :
                                edges_file.write("\t%s -- %s [label=%s, color=red];\n" % (line, column, self.matrix.item(line, column)))
                                sol_file.write("\t%s -- %s;\n" % (line, column))
                            elif column > line :
                                edges_file.write("\t%s -- %s [label=%s];\n" % (line, column, self.matrix.item(line, column)))
                edges_file.write("}\n")
                sol_file.write("}\n")
