# -*- coding: utf-8 -*-
"""
stuff
"""

from pycam.PathProcessors import BasePathProcessor
# use Grid and Box objects from Contour2dCutter
from pycam.Geometry.Path import Path
import heapq # for Dijkstra's algorithm optimisation (using priority queue)
from solvinggraph import SolvingGraph
from numpy import zeros as numpy_zeros
from sys import getsizeof

class Pocket(object) :
    id = 1

    def __init__(self) :
        self.id = Pocket.id
        Pocket.id += 1
        self.boxes = []
        self.pocket_above = None
        self.pockets_underneath = []
        self.matrix_distance = None
        self.matrix_predecessor = None
        self.solved_path = None

    @property
    def bottom(self) : return not self.pockets_underneath

    @property
    def size(self) : return len(self.boxes)

    def free(self) :
        #print "size of", self, getsizeof(self)
        #print "\tsize of boxes", getsizeof(self.boxes)
        #print "\tsize of above", getsizeof(self.pocket_above)
        #print "\tsize of underneath", getsizeof(self.pockets_underneath)
        #print "\tsize of distance", getsizeof(self.matrix_distance)
        #print "\tsize of predecessor", getsizeof(self.matrix_predecessor)
        #print "\tsize of path", getsizeof(self.solved_path)
        #print "\ttotal size :", sum(map(getsizeof, (self, self.boxes, self.pocket_above, self.pockets_underneath, self.matrix_distance, self.matrix_predecessor, self.solved_path))) 
        del self.boxes[:]
        del self.boxes
        del self.pocket_above
        del self.pockets_underneath[:]
        del self.pockets_underneath
        del self.matrix_distance
        for line in self.matrix_predecessor :
            del line[:]
        del self.matrix_predecessor
        del self.solved_path[:]
        del self.solved_path

    def __str__(self) :
        return("Pocket %s (%s / %s) [%s]" % (self.id, self.pocket_above.id if self.pocket_above else "0", \
            ','.join(str(pocket.id) for pocket in self.pockets_underneath), len(self.boxes)))

    def compute_dijkstra(self) :
        # build matrices to hold results
        length = len(self.boxes)
        matrix_distance = numpy_zeros(shape=(length,length))
        matrix_predecessor = [[None]*length for _ in xrange(length)]

        # assign to each box his index in matrices
        i = 0
        while i < len(self.boxes) :
            self.boxes[i].index = i
            i += 1

        # some computation time optimizations
        pop = heapq.heappop
        push = heapq.heappush

        for box in self.boxes :
            # then here come Dijkstra :
            box_index = box.index
            matrix_distance[box_index][box_index] = -1
            matrix_predecessor[box_index][box_index] = box
            neighbours_left = [box] # is a heap
            box._distance = 0
            while neighbours_left :
                nearest = pop(neighbours_left)
                new_distance = nearest._distance + 1
                for neighbour in nearest.free_neighbours :
                    previous_distance = matrix_distance.item(box_index, neighbour.index)
                    if previous_distance == 0 :
                        neighbour._distance = new_distance
                        matrix_distance[box_index][neighbour.index] = new_distance
                        matrix_predecessor[box_index][neighbour.index] = nearest
                        push(neighbours_left, neighbour)
                    # else the path would be anyway equally long or longer

        # finally save results
        self.matrix_distance = matrix_distance
        self.matrix_predecessor = matrix_predecessor

    def solve_path(self) :
        print "dijkstra computation begin ..."
        self.compute_dijkstra()
        print "dijkstra computation ended"
        graph = SolvingGraph(self.matrix_distance)
        index_path = graph.christofides()
        box_path = [self.boxes[index[0]] for index in index_path] + [self.boxes[index_path[0][0]]]
        real_path = []
        i = 0
        while i < len(box_path) :
            current_box = box_path[i]
            real_path.append(current_box)
            next_box = box_path[(i+1)%len(box_path)]
            partial_path = []
            step_box = self.matrix_predecessor[current_box.index][next_box.index]
            while step_box != current_box :
                partial_path.insert(0, step_box)
                step_box = self.matrix_predecessor[current_box.index][step_box.index]
            real_path.extend(partial_path)
            i += 1
        self.solved_path = real_path

    def solve_path_from(self, box) :
        self.solve_path()
        cut = self.solved_path.index(box)
        return self.solved_path[cut:] + self.solved_path[:cut]

class Contour2dFlooder(BasePathProcessor) :

    def __init__(self) :
        super(Contour2dFlooder, self).__init__()
        self.maxx = None
        self.maxy = None
        self.maxz = None
        self.grid = None
        self._layer = None

    def initialise(self, grid) :
        self.maxx = grid.nb_lines - 1
        self.maxy = grid.nb_columns - 1
        self.maxz = grid.nb_layers - 1
        self.grid = grid

    def do_path(self) :
        for layer in xrange(self.maxz+1) :
            self.do_layer(self.maxz - layer)
        #self.draw_pocket_PBM()
        box_path = self.find_path_from(self.grid.get_box(0, 0, self.grid.nb_layers-1)) # TODO: skypocket
        self.paths.append(Path())
        map(self.paths[0].append, (box.get_center() for box in box_path))

    def find_path_from(self, start_box) :
        #print "find path from", start_box
        start_pocket = start_box.pocket
        print "pocket size :", len(start_pocket.boxes)
        pocket_path = start_pocket.solve_path_from(start_box)
        if start_pocket.bottom :
            path = pocket_path
        else :
            underneath_pockets = {}
            encountered_boxes = {}
            for underneath_pocket in start_pocket.pockets_underneath :
                underneath_pockets[underneath_pocket] = underneath_pocket.size
            path = []
            for box in pocket_path :
                path.append(box)
                if box not in encountered_boxes :
                    encountered_boxes[box] = 1
                    underneath_box = box.get_neighbour(0, 0, -1)
                    if underneath_box.pocket is not None :
                        underneath_pockets[underneath_box.pocket] -= 1
                        if underneath_pockets[underneath_box.pocket] == 0 :
                            path.extend(self.find_path_from(underneath_box) + [box])
                            path.append(box)
                        # else the pocket cannot be milled yet
                    # else the cutter is not above pocket
                # else the box has already been milled
        start_pocket.free()
        del start_pocket
        return path

    def do_layer(self, layer) :
        self._layer = self.grid.layers[layer]
        outside = self.do_outside()
        inside = self.do_inside()
        self.link_pockets(outside+inside)
        self.do_overhang()

    def display_layer(self) :
        print '\n'*20
        for line in self._layer :
            for column in line :
                print "%s"%column.sq(),
            print ''

    def do_overhang(self) :
        for x in xrange(self.maxx) :
            for y in xrange(self.maxy) :
                box = self._layer[x][y]
                if not box.free :
                    while box.z != 0 :
                        box = box.get_neighbour(0, 0, -1)
                        box.scan = True

    def do_outside(self) :
        outside = []
        box = self._layer[0][0]
        for x in xrange(self.maxx+1-1) :
            if box.free and not box.scan :
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(1, 0, 0)
        for y in xrange(self.maxy+1-1) :
            if box.free and not box.scan :
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(0, 1, 0)
        for x in xrange(self.maxx+1-1) :
            if box.free and not box.scan :
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(-1, 0, 0)
        for y in xrange(self.maxy+1-2) : # -2 to not iterate again on [0][0]
            if box.free and not box.scan :
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(0, -1, 0)
        return outside

    def do_inside(self) :
        inside = []
        for x in xrange(1, self.maxx) :
            for y in xrange(1, self.maxy) :
                box = self._layer[x][y]
                if not box.scan and box.free :
                    if not box.inside :
                        inside.append(Pocket())
                        self.flood(box, inside[-1].boxes)
                    else :
                        self.flood(box, [])
        return inside

    def flood(self, box, pocket) :
        # see en.wikipedia.org/wiki/Flood_fill
        # alternative implementation of the algorithm
        # using a stack and 2 loops (east and west)
        if not box.free : return
        stack = []
        stack.append(box)
        while stack :
            center = stack.pop()
            if (not center.scan) and center.free :
                west = east = center
                while True :
                    if west.y == 0 : break
                    _west = west.get_neighbour(0, -1, 0)
                    if (not _west.free) or _west.scan : break
                    else : west = _west
                while True :
                    if east.y == self.maxy : break
                    _east = east.get_neighbour(0, 1, 0)
                    if (not _east.free) or _east.scan : break
                    else : east = _east
                for y in range(west.y, east.y+1) :
                    b = box.get_box(center.x, y, center.z)
                    b.scan = True
                    pocket.append(b)
                    if b.x != self.maxx :
                        north = b.get_neighbour(1, 0, 0)
                        if north.free : stack.append(north)
                    if b.x != 0 :
                        south = b.get_neighbour(-1, 0, 0)
                        if south.free : stack.append(south)

    def link_pockets(self, pockets) :
        for pocket in pockets :
            if pocket.boxes[0].z != self.maxz :
                pocket.pocket_above = pocket.boxes[0].get_neighbour(0, 0, 1).pocket
                pocket.pocket_above.pockets_underneath.append(pocket)
                self.link_pocket(pocket)
            else :
                pocket.pocket_above = None
                self.link_pocket(pocket)

    def link_pocket(self, pocket) :
        for box in pocket.boxes :
            box.pocket = pocket
            box.know_your_free_neighbours()

    def draw_pocket_PBM(self) :
        for nb_layer in xrange(self.grid.nb_layers) :
            with open("pocket_%d.pbm"%(nb_layer+1), "w") as pocketfile :
                pocketfile.write("P1\n")
                pocketfile.write("# %d/%d\n"%(nb_layer+1, self.grid.nb_layers))
                pocketfile.write("%d %d\n"%(self.grid.nb_columns, self.grid.nb_lines))
                for line in range(self.grid.nb_lines) :
                    pocketfile.write(' '.join(['0' if box.pocket else '1'
                            for box in self.grid.layers[nb_layer][line]])+'\n')

