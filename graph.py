import numpy
import variablepriorityqueue
import disjointset


class Graph(object) :

    def read_file(self, filename) :
        matrix_file = open(filename, 'r')
        self.size = int(matrix_file.readline().rstrip())
        self.matrix = numpy.zeros(shape=(self.size,self.size))
        current_line_number = 0
        for line in matrix_file :
            strings = line.rstrip().split(' ')
            floats = map(float, strings)
            self.matrix[current_line_number] = floats
            current_line_number += 1

    def christofides(self) :
        spanning_tree = self.solve_spanning_tree(range(self.size))
        self.display_selected_edges(spanning_tree, "spanning")
        edges_by_vertex = {}
        for edge in spanning_tree :
            for vertex in edge :
                if vertex not in edges_by_vertex :
                    edges_by_vertex[vertex] = [edge]
                else : edges_by_vertex[vertex].append(edge)
        odd_degree_vertexes = filter(lambda v : len(edges_by_vertex[v])%2, edges_by_vertex.keys())
        perfect_matching = self.solve_perfect_matching(odd_degree_vertexes)
        self.display_selected_edges(perfect_matching, "matching")
        for edge in perfect_matching :
            for vertex in edge :
                edges_by_vertex[vertex].append(edge)
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
        edges_left = len(spanning_tree) + len(perfect_matching) - len(eulerian_path)
        while edges_left :
            for edge in eulerian_path :
                if edges_by_vertex[edge[0]] :
                    additional_path = find_path(edge[0])
                    edges_left -= len(additional_path)
                    index_start = eulerian_path.index(edge)
                    eulerian_path = eulerian_path[:index_start] + additional_path + eulerian_path[index_start:]
        seen = {}
        i = 0
        while i < len(eulerian_path) :
            edge = eulerian_path[i]
            if edge[1] not in seen :
                seen[edge[1]] = edge[0]
                i += 1
            else :
                eulerian_path[i][1] = eulerian_path[i+1][1]
                eulerian_path.pop(i+1)
        return eulerian_path # which is now hamiltonian

    #see 'a general approximation technique for constrained forest problems' by goemans et williamson
    #we use it both for min spanning tree and min perfect matching
    #we need a cf function such that f(C1 U C2) = cf(f(C1),f(C2)) #see paper for what is f
    # we also need the vertexes that the algorithm should work with, these are indexes in self.matrix as well
    # because goemans is applied first to find spanning tree of the initial graph, then
    # to find perfect matching amongst odd degree vertexes of the spanning tree # see Christofides
    # both graphs are assumed to be complete (degree of each vertex equals number of vertexes minus one)
    def goemans(self, cf, vertexes) :
        F = [] #holding solution
        class QueueElement(object) :
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
                return "QueueElement(%s,%s, %s)" % (self.i, self.j, self.priority)
        #declare edges priority queue and initialize it
        queue = variablepriorityqueue.VariablePriorityQueue()
        for line in vertexes :
            for column in vertexes :
                if line != column :
                    weight = self.matrix.item(line, column)
                    edge = QueueElement(line,column,weight/2)
                    queue.add(edge)

        #define union find structure and initialize it
        C = disjointset.DisjointSet()
        f = {} #also store the values of f function (access it only through nodes representatives)
        d = {} #store nodes weights
        for vertex in vertexes :
            C.add(vertex)
            f[vertex] = 1
            d[vertex] = 0

        #store the number of trees C_r such that f(C_r) = 1
        #we need this info to check quickly for then end of the computations
        #(when it reaches 1)
        number_of_f1_trees = len(vertexes)
        
        #this is the function that computes priority of an edge
        def epsilon(i, j, d, f) :
            sum_f = f[C.find(i)] + f[C.find(j)]
            if sum_f == 0 :
               sum_f = 0.00001 # one divided by one hundred thousand
            return (self.matrix.item(i,j) -d[i] -d[j])/(sum_f)
        
        #we also need a matrix that holds the best edge between two trees
        #can only be accessed through nodes representatives
        best_edges = numpy.zeros(shape=(self.size,self.size), dtype=numpy.int)
        for line in vertexes :
            for column in vertexes :
                if line != column :
                    #we encode the edge to be able to store it in an integer
                    best_edges[line][column] = line*self.size + column

        #initialize lower bound
        LB = 0

        #initialize time
        #in priority queue the real priority we want
        #is epsilon(i,j)
        #however, it is too costly to store since any event
        #would require modifying all priorities
        #we instead use a global time T
        #and define priorities as epsilon(i,j) + T
        #see paper for more details
        T = 0

        #main loop
        while number_of_f1_trees > 1 :
            #take best edge
            best_edge = queue.get()

            if C.find(best_edge.i) == C.find(best_edge.j) :
                continue #skip edges in same tree

            if best_edges[C.find(best_edge.i)][C.find(best_edge.j)] != best_edge.i*self.size+best_edge.j :
                continue #skip edges which we do not update because not best ones

            #add it to solution
            F.append([best_edge.i, best_edge.j])
            e = epsilon(best_edge.i, best_edge.j, d, f)
            #update time
            T += e

            #update d
            for v in vertexes :
                d[v] += e * f[C.find(v)]

            #update lower bound
            for r in C.retrieve() :
                LB += e * f[r]

            #now last and complex part of the algorithm, do the fusion and weights update
            #start by updating number of trees with f=1
            r_i = C.find(best_edge.i)
            r_j = C.find(best_edge.j)
            v_i = f[r_i]
            v_j = f[r_j]

            new_f_value = cf(v_i, v_j)

            #compute new amount of trees with f value equel to 1
            if v_i + v_j > new_f_value :
                number_of_f1_trees -= 2 - cf(1, 1) #one or two less
            #now, fuse
            C.union(r_i, r_j)
            #update f
            r_union = C.find(r_i)

            f[r_union] = new_f_value

            #modify best edges going out of r_union tree
            #start with line
            for column in vertexes :
                coded_edge1 = best_edges[r_i][column]
                coded_edge2 = best_edges[r_j][column]
                edge1 = [coded_edge1 // self.size, coded_edge1 % self.size]
                edge2 = [coded_edge2 // self.size, coded_edge2 % self.size]
                e1 = epsilon(edge1[0], edge1[1], d, f)
                e2 = epsilon(edge2[0], edge2[1], d, f)
                if e1 < e2 :
                    best_edges[r_union][column] = coded_edge1
                else :
                    best_edges[r_union][column] = coded_edge2

            #continue with column
            for line in vertexes :
                coded_edge1 = best_edges[line][r_i]
                coded_edge2 = best_edges[line][r_j]
                edge1 = [coded_edge1 // self.size, coded_edge1 % self.size]
                edge2 = [coded_edge2 // self.size, coded_edge2 % self.size]
                e1 = epsilon(edge1[0], edge1[1], d, f)
                e2 = epsilon(edge2[0], edge2[1], d, f)
                best_edges[line][r_union] = coded_edge1 if e1 < e2 else coded_edge2

            #final step : modify priority queue
            #loop on all best edges leaving r_union tree
            #all edges leaving r_union tree which are not best
            #are therefore NOT UPDATED and will have ERRONEOUS priorities
            #this is why we filter edges at beginning of the loop to keep only best ones
            for column in vertexes :
                if (column != r_i) and (column != r_j) and C.find(column) != r_union :
                    coded_edge = best_edges[r_union][column]
                    edge = [coded_edge // self.size, coded_edge % self.size]
                    priority = epsilon(r_union, column, d, f) + T
                    changing_edge = QueueElement(r_union,column,priority)
                    queue.update(changing_edge)
        return F

    def solve_spanning_tree(self, vertexes) :
        return self.goemans(lambda f1,f2 : 1, vertexes)

    def solve_perfect_matching(self, vertexes) :
        F = self.goemans(lambda f1,f2 : (f1+f2)%2, vertexes)

        # create 
        edges_by_vertex = {}
        for edge in F :
            for vertex in edge :
                if vertex not in edges_by_vertex :
                    edges_by_vertex[vertex] = [edge]
                else : edges_by_vertex[vertex].append(edge)

        #returns from a graph the number of components
        #containing an even number of vertices
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
        #now try to remove one by one each edge
        #for each edge we simply remove it
        #and see on remaining graph if there exist
        #any connected component of odd size
        #if yes, keep this edge else delete it
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
            while len(edges_by_vertex[vertex]) >= 3 :
                edge1 = edges_by_vertex[vertex][-1]
                edge2 = edges_by_vertex[vertex][-2]
                if edge1[0] == edge2[0] : shortcut = [edge1[1], edge2[1]]
                elif edge1[0] == edge2[1] : shortcut = [edge1[1], edge2[0]]
                elif edge1[1] == edge2[0] : shortcut = [edge1[0], edge2[1]]
                elif edge1[1] == edge2[1] : shortcut = [edge1[0], edge2[0]]
                edges_by_vertex[vertex].pop()
                edges_by_vertex[vertex].pop()
                F.append(shortcut)
                F.remove(edge1)
                F.remove(edge2)
        return F

    #print dot file
    def display_selected_edges(self, edges, name="goemans") :
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
                                edges_file.write("\t%s -- %s [label=%s, color=red];\n" % (line, column, self.matrix[line][column]))
                                sol_file.write("\t%s -- %s;\n" % (line, column))
                            elif column > line :
                                edges_file.write("\t%s -- %s [label=%s];\n" % (line, column, self.matrix[line][column]))
                edges_file.write("}\n")
                sol_file.write("}\n")
