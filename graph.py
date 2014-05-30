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
    def display(self) :
        print self.matrix
    #see 'a general approximation technique for constrained forest problems' by goemans et williamson
    #we use it both for min spanning tree and min perfect matching
    #we need a cf function such that f(C1 U C2) = cf(f(C1),f(C2)) #see paper for what is f
    def goemans(self, cf) :
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
        for line in xrange(self.size) :
            for column in xrange(self.size) :
                if line != column :
                    weight = self.matrix.item(line, column)
                    edge = QueueElement(line,column,weight/2)
                    queue.add(edge)
        queue.disp()
        #define union find structure and initialize it
        C = disjointset.DisjointSet()
        f = [] #also store the values of f function (access it only through nodes representatives)
        d = [] #store nodes weights
        for vertex in range (self.size) :
            C.add(vertex)
            f.append(1)
            d.append(0)

        #store the number of trees C_r such that f(C_r) = 1
        #we need this info to check quickly for then end of the computations
        #(when it reaches 1)
        number_of_f1_trees = self.size
        
        #this is the function that computes priority of an edge
        def epsilon(i, j, d, f) :
            sum_f = f[C.find(i)] + f[C.find(j)]
            if sum_f == 0 :
               sum_f = 0.000001
            return (self.matrix.item(i,j) -d[i] -d[j])/(sum_f)
        
        #we also need a matrix that holds the best edge between two trees
        #can only be accessed through nodes representatives
        best_edges = numpy.zeros(shape=(self.size,self.size), dtype=numpy.int)
        for line in xrange(self.size) :
        	for column in xrange(self.size) :
        		if line != column :
                    #we encode the edge to be able to store it in an integer
        			best_edges[line][column] = line*self.size + column

        #initialize lower bound
        LB = 0

        #main loop
        while number_of_f1_trees > 1 :
            #take best edge
            best_edge = queue.get()
            print "Best edge found :", best_edge
            if C.find(best_edge.i) == C.find(best_edge.j) :
                print "But edge inside a tree"
                continue #skip edges in same tree
            if best_edges[C.find(best_edge.i)][C.find(best_edge.j)] != best_edge.i*self.size+best_edge.j :
                print "But not really best_edge indeed"
                continue #skip edges which we do not update because not best ones
            #add it to solution
            F.append([best_edge.i, best_edge.j])
            e = epsilon(best_edge.i, best_edge.j, d, f)
            print "Epsilon computed :", e
            #update d
            old_d = list(d) #backup because we will need the old info at the end
            print "Adding weight"
            for v in range(self.size) :
                d[v] += e * f[C.find(v)]
            print "Weigths before : ", old_d
            print "Weights after : ", d
            #update lower bound
            for r in C.retrieve() :
                LB += e * f[r]
            print "Lower Bound updated :", LB
            #now last and complex part of the algorithm, do the fusion and weights update
            #start by updating number of trees with f=1
            r_i = C.find(best_edge.i)
            r_j = C.find(best_edge.j)
            v_i = f[r_i]
            v_j = f[r_j]
            print "Representative of i :", r_i
            print "Value :", v_i
            print "Representative of j :", r_j
            print "Value :", v_j
            new_f_value = cf(v_i, v_j)
            print "Composition :", new_f_value
            if v_i + v_j > new_f_value :
                number_of_f1_trees -= 2 - cf(1, 1) #one or two less
                print "f1 trees less :", 2-cf(1,1)
            #now, fuse
            C.union(r_i, r_j)
            #update f
            r_union = C.find(r_i)
            print "Representative of ij :", r_union
            old_f = list(f) #backup because we will need the old info at the end
            f[r_union] = new_f_value
            print "Trees :", f
            #modify best edges going out of r_union tree
            #start with line
            for column in range(self.size) :
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
            for line in range(self.size) :
                coded_edge1 = best_edges[line][r_i]
                coded_edge2 = best_edges[line][r_j]
                edge1 = [coded_edge1 // self.size, coded_edge1 % self.size]
                edge2 = [coded_edge2 // self.size, coded_edge2 % self.size]
                e1 = epsilon(edge1[0], edge1[1], d, f)
                e2 = epsilon(edge2[0], edge2[1], d, f)
                if e1 < e2 :
                    best_edges[line][r_union] = coded_edge1
                else :
                    best_edges[line][r_union] = coded_edge2

            #final step : modify priority queue
            #loop on all best edges leaving r_union tree
            #all edges leaving r_union tree which are not best
            #are therefore NOT UPDATED and will have ERRONEOUS priorities
            #this is why we filter edges at beginning of the loop to keep only best ones
            print r_i, r_j
            for column in range(self.size) :
                if (column != r_i) and (column != r_j) and C.find(column) != r_union :
                    coded_edge = best_edges[r_union][column]
                    edge = [coded_edge // self.size, coded_edge % self.size]
                    #compute by how much the corresponding epsilon will change
                    old_e = epsilon(edge[0], edge[1], old_d, old_f)
                    new_e = epsilon(edge[0], edge[1], d, f)
                    #normaly we should have a new priority of new_e
                    #but we do not store epsilons as priorities
                    #the reason is that at each step d increases for many edges
                    #so all 'interesting' edges not concerned by the fusion
                    #should see their epsilon decrease by e
                    #we do not want to update so many priorities but still
                    #want the right edges order
                    #so instead we add e to everyone's priority
                    #this way, we do not need to update all edges with a priority change
                    # of exactly -e
                    priority_delta = new_e - old_e + e 
                    changing_edge = QueueElement(r_union,column,priority_delta)
                    queue.update(changing_edge)
        return F
    def solve_spanning_tree(cf) :
        return self.goemans(cf)
    def solve_perfect_matching(cf) :
        F = self.goemans(cf)
        # now we must remove any edge that's not needed
        edges_by_vertices = {}
        for edge in F :
            if edge[0] not in edges_by_vertices :
                edges_by_vertices[edge[0]] = [edge]
            else :
                edges_by_vertices[edge[0]].append(edge)
            if edge[1] not in edges_by_vertices :
                edges_by_vertices[edge[1]] = [edge]
            else :
                edges_by_vertices[edge[1]].append(edge)
    #print dot file
    def display_selected_edges(self, edges, name="goemans.dot") :
        with open(name, 'w') as edges_file :
            h = {}
            for edge in edges :
                h[edge[0]*self.size+edge[1]] = 1
            edges_file.write("graph {\n")
            for line in range(self.size) :
                for column in range(line, self.size) :
                    if line != column :
                        if line*self.size+column in h :
                            edges_file.write("\t%s -- %s [label=%s, color=red];\n" % (line, column, self.matrix[line][column]))
                        else :
                            edges_file.write("\t%s -- %s [label=%s];\n" % (line, column, self.matrix[line][column]))
            edges_file.write("}\n")
