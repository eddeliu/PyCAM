# -*- coding: utf-8 -*-
"""
stuff
"""

from pycam.PathProcessors import BasePathProcessor
# use Grid and Box objects from Contour2dCutter


class Pocket(object) :
    id = 1

    def __init__(self) :
        self.id = Pocket.id
        Pocket.id += 1
        self.pocket_above = None
        self.boxes = []
        self.pockets_underneath = []

    @property
    def bottom(self) : return not self.pockets_underneath

    def __str__(self) :
        return("Pocket %s (%s / %s)\n" % (self.id, self.pocket_above.id if self.pocket_above else "0", \
            ','.join(str(pocket.id) for pocket in self.pockets_underneath)) \
            + '\n\t'.join(str(box) for box in self.boxes))

class DisjointSet(object) :
    def find(self, value) :
        if value != value._repr :
            value._repr = self.find(value._repr)
        return value._repr
    def union(self, value, other_value) : # 2nd value's repr is 1st value one's
        self.find(other_value)._repr = self.find(value)._repr
    def add(self, value) :
        value._repr = value

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
        self.draw_pocket_PBM()

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
                self.link_pocket(pocket)
            else :
                pocket.pocket_above = 0
                self.link_pocket(pocket)
            #pocket.boxes[0].export_your_graph_of_shortest_paths()

    def link_pocket(self, pocket) :
        for box in pocket.boxes :
            box.pocket = pocket
            box.know_your_free_neighbours()
        for box in pocket.boxes :
            box.know_your_pocket_neighbourhood()

    def draw_pocket_PBM(self) :
        for nb_layer in xrange(self.grid.nb_layers) :
            with open("pocket_%d.pbm"%(nb_layer+1), "w") as pocketfile :
                pocketfile.write("P1\n")
                pocketfile.write("# %d/%d\n"%(nb_layer+1, self.grid.nb_layers))
                pocketfile.write("%d %d\n"%(self.grid.nb_columns, self.grid.nb_lines))
                for line in range(self.grid.nb_lines) :
                    pocketfile.write(' '.join(['0' if box.pocket else '1'
                            for box in self.grid.layers[nb_layer][line]])+'\n')

