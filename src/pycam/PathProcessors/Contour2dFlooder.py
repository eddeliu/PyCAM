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

from numpy import zeros as numpy_zeros
import heapq # for Dijkstra's algorithm optimisation (using priority queue)

from pycam.PathProcessors import BasePathProcessor
from pycam.Geometry.Path import Path
from pycam.Utils.SolvingGraph import SolvingGraph
from pycam.Utils import ProgressCounter

import pycam.Utils.log
log = pycam.Utils.log.get_logger()

class Pocket(object) :
    """A block of boxes that all have to be milled

    That class holds its boxes, above and underneath pockets to
    simplify toolpathing, and realizes it by solving eulerian path
    within the graph formed by its boxes."""

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
    def bottom(self) :
        """If the pocket has not underlying pockets"""
        return not self.pockets_underneath

    @property
    def size(self) :
        """The size of the pocket"""
        return len(self.boxes)

    def free_memory(self) :
        """Attempt to free RAM"""
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
        """Return a not-so-awful visualisation of a pocket"""
        return("Pocket %s (%s / %s) [%s]" % (self.id, self.pocket_above.id if self.pocket_above else "0", \
            ','.join(str(pocket.id) for pocket in self.pockets_underneath), len(self.boxes)))

    def compute_dijkstra(self) :
        """Computes Dijkstra's algorithme for each box in it

        It permits to fill the two matrices that respectively hold
        the minimal distance between two boxes and the box to do so.
        It is used after by SolvingGraph to solve the path."""
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
        """Creates the path for the upcoming find_path_from

        It computes Disjkstra, feed Christofides with it,
        but what it gets is the solution for a complete graph.
        So we have now to process the edges found into actual
        box path. For this, we use the predecessor matrix
        we have filled formerly."""
        log.info("dijkstra computation begin ...")
        self.compute_dijkstra()
        log.info("dijkstra computation ended")
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
        """Solve the path, and start it from a given box"""
        self.solve_path()
        cut = self.solved_path.index(box)
        return self.solved_path[cut:] + self.solved_path[:cut]


class Contour2dFlooder(BasePathProcessor) :
    """PathProcessor for Contour2dCutter : pocketing, graph solving

    It works on the grid Contour2dCutter created and filled with free
    or empty boxes. Layer by layer, it detects the pockets, solve path
    for each and then find the global path accross all pockets."""

    def __init__(self) :
        """Creates the path processor

        Indeed, it doesn't initialise it because it is created way
        before it is used, so there is a real method "initialise" that
        do the real job for when the path processor is actually needed"""
        super(Contour2dFlooder, self).__init__()
        self.maxx = None
        self.maxy = None
        self.maxz = None
        self.grid = None
        self._layer = None
        self._callback = None
        self._nb_pockets = None
        self._current_nb_pocket = None
        self._quit_requested = None
        self._progress_counter = None

    def initialise(self, grid, callback) :
        """Initialise the path processor from grid"""
        self.maxx = grid.nb_lines - 1
        self.maxy = grid.nb_columns - 1
        self.maxz = grid.nb_layers - 1
        self.grid = grid
        self._callback = callback
        self._nb_pockets = 0
        self._current_nb_pocket = 1
        self._quit_requested = False

    def do_path(self) :
        """Computes the toolpath, starting from top [0][0] of grid"""
        # TODO: take care of the top pocket which could be filled,
        # particulary top[0][0] or can be several pockets
        # => don't forget to remove the "+1" for actual Grid creation
        progress_counter = ProgressCounter(self.maxz, self._callback)
        for layer in xrange(self.maxz+1) :
            if self._callback and self._callback(text="Contour2dCutter: pocketing " \
                        + "layer %d/%d" % (layer+1, self.maxz)):
                # cancel requested
                return
            self.do_layer(self.maxz - layer)
            progress_counter.increment()
        #self.draw_pocket_PBM()
        self._progress_counter = ProgressCounter(self._nb_pockets, self._callback)
        box_path = self.find_path_from(self.grid.get_box(0, 0, self.grid.nb_layers-1))
            # there is the hack part --^
        self.paths.append(Path())
        map(self.paths[0].append, (box.get_center() for box in box_path))

    def find_path_from(self, start_box) :
        """Return the path within the pocket, starting from this box"""
        if self._callback and self._callback(text="Contour2dCutter: toolpathing " \
                    + "pocket %d/%d" % (self._current_nb_pocket, self._nb_pockets)):
            # cancel requested
            self._quit_requested = True
            return []
        self._current_nb_pocket += 1
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
                            if self._quit_requested : return []
                        # else the pocket cannot be milled yet
                    # else the cutter is not above pocket
                # else the box has already been milled
        start_pocket.free_memory()
        del start_pocket
        self._progress_counter.increment()
        return path

    def do_layer(self, layer) :
        """Proceeds pockets of one layer"""
        self._layer = self.grid.layers[layer]
        outside = self.do_outside()
        inside = self.do_inside()
        self.link_pockets(outside+inside)
        self.do_overhang()

    def display_layer(self) :
        """Display one layer, for debugging"""
        print '\n'*20
        for line in self._layer :
            for column in line :
                print "%s"%column.sq(),
            print ''

    def do_overhang(self) :
        """Ensure that there is no filled box above free ones

        This is an optionnal safety guard. But if it is not
        applied and the model have overhangs, you will get
        toolpath such the tool handle will collide with objet,
        resulting in very costly and annoying results"""
        for x in xrange(self.maxx) :
            for y in xrange(self.maxy) :
                box = self._layer[x][y]
                if not box.free :
                    while box.z != 0 :
                        box = box.get_neighbour(0, 0, -1)
                        box.scan = True

    def do_outside(self) :
        """Proceed "outside" pockets

        An outside pocket is a pocket whose at least one box
        is on the edge of the grid. Empirically there is a ton
        of this pocket, so we use a very simple method to check
        them out, so it have a low cost. It is indeed just four
        for loops with the exact same body."""
        outside = []
        box = self._layer[0][0]
        for x in xrange(self.maxx+1-1) :
            if box.free and not box.scan :
                self._nb_pockets += 1
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(1, 0, 0)
        for y in xrange(self.maxy+1-1) :
            if box.free and not box.scan :
                self._nb_pockets += 1
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(0, 1, 0)
        for x in xrange(self.maxx+1-1) :
            if box.free and not box.scan :
                self._nb_pockets += 1
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(-1, 0, 0)
        for y in xrange(self.maxy+1-2) : # -2 to not iterate again on [0][0]
            if box.free and not box.scan :
                self._nb_pockets += 1
                outside.append(Pocket())
                self.flood(box, outside[-1].boxes)
            else :
                box.scan = True
            box = box.get_neighbour(0, -1, 0)
        return outside

    def do_inside(self) :
        """Finds pockets inside the object

        This method is a bit more costly because it have to iterate
        over every box to check if it has be scanned, and if no if
        it belongs to a box."""
        inside = []
        for x in xrange(1, self.maxx) :
            for y in xrange(1, self.maxy) :
                box = self._layer[x][y]
                if not box.scan and box.free :
                    if not box.inside :
                        self._nb_pockets += 1
                        inside.append(Pocket())
                        self.flood(box, inside[-1].boxes)
                    else :
                        self.flood(box, [])
        return inside

    def flood(self, box, pocket) :
        """Floods all boxes of a pocket

        This is the algorithm FloodFill (see
        en.wikipedia.org/wiki/Flood_fill)
        But it has an alternative implementation :
        it uses a stack and two loops (east and west)"""
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
        """Links the pockets with each other, for upcoming path finding"""
        for pocket in pockets :
            if pocket.boxes[0].z != self.maxz :
                pocket.pocket_above = pocket.boxes[0].get_neighbour(0, 0, 1).pocket
                pocket.pocket_above.pockets_underneath.append(pocket)
                self.link_pocket(pocket)
            else :
                pocket.pocket_above = None
                self.link_pocket(pocket)

    def link_pocket(self, pocket) :
        """Links the boxes to their pocket, for upcoming path finding"""
        for box in pocket.boxes :
            box.pocket = pocket
            box.know_your_free_neighbours()

    def draw_pocket_PBM(self) :
        """It creates representation of the pocket for each layer

        The output is a PBM file, you can look it graphically""" 
        for nb_layer in xrange(self.grid.nb_layers) :
            with open("pocket_%d.pbm"%(nb_layer+1), "w") as pocketfile :
                pocketfile.write("P1\n")
                pocketfile.write("# %d/%d\n"%(nb_layer+1, self.grid.nb_layers))
                pocketfile.write("%d %d\n"%(self.grid.nb_columns, self.grid.nb_lines))
                for line in range(self.grid.nb_lines) :
                    pocketfile.write(' '.join(['0' if box.pocket else '1'
                            for box in self.grid.layers[nb_layer][line]])+'\n')

