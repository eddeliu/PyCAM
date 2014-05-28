import numpy
import variablepriorityqueue
import disjointset


class Graph(object) :
    def read_file(self, filename) :
        file = open(filename, 'r')
        self.size = int(file.readline().rstrip())
        self.matrix = numpy.zeros(shape=(self.size,self.size))
        current_line_number = 0
        for line in file :
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
        class QueueElement(object) :
            def __init__(self, i, j, priority) :
                self.i = i
                self.j = j
                self.priority = priority
            def __lt__(self, other) :
                return self.priority < other.priority
        #declare edges priority queue and initialize it
        queue = VariablePriorityQueue()
        for line in range (self.size) :
            for column in range (self.size) :
                weight = self.matrix.item(line, column)
                queue.add(QueueElement(line,column,weight/2))
        #define union find structure and initialize it
        C = DisjointSet()
        f = [] #also store the values of f function
        for vertex in range (self.size) :
            C.add(vertex)
            f.append(1)
        #store the number of trees C_r such that f(C_r) = 1
        #we need this info to check quickly for then end of the computations
        #(when it reaches 0)
        number_of_f1_trees = self.size

        #initialize lower bound
        LB = 0

        #main loop
        while number_of_f1_trees :
            #take best edge between two different trees
            best_edge = queue.get()
            if C.find(best_edge.i) == C.find(best_edge.j) :
                continue
            #add it to solution
            F.append([best_edge.i, best_edge.j])

