# -*- coding: utf-8 -*-
"""
stuff
"""

from pycam.Geometry.Plane import Plane
from pycam.Geometry.Point import Point, Vector
from pycam.Geometry.Line import Line
from pycam.Geometry.utils import INFINITE, epsilon, ceil, sqrt
def _is_near(x, y):
    return abs(x - y) < epsilon
import heapq # for Dijkstra's algorithm optimisation (using priority queue)

from pycam.Utils import ProgressCounter
import pycam.Utils.log

log = pycam.Utils.log.get_logger()

class Box(object) :
    # "Abstract" methods :
    #   - get_center
    #   - get_location (static)
    #   - know_your_free_neighbours
    def __init__(self, x, y, z) :
        self.x = x
        self.y = y
        self.z = z
        self.index = 0
        self._distance = 0 # used by Dijkstra's algorithm
        self.lines = []
        self.scan = False
        self.pocket = None
        self.free_neighbours = []
    def __lt__(self, other) :
        # used in Dijkstra's algorithm implementation by heapq.heappush
        return self._distance < other._distance
    def free_memory(self) :
        del self.pocket
        del self.free_neighbours
    @property
    def free(self) : return len(self.lines) == 0
    @property
    def filled(self) : return self.poche == None
    @property
    def inside(self) :
        p1 = self.get_center()
        p2 = Point(Box.grid.minx-1, p1.y, p1.z) # TODO: implement algo to find nearest border
        escaping_line = Line(p1, p2)
        cuts = 0
        box = self
        lines = set()
        while box.x != 0 :
            box = box.get_neighbour(-1, 0, 0)
            for line in box.lines :
                if line.id not in lines :
                    cut, d = escaping_line.get_intersection(line)
                    if cut is not None :
                        cuts += 1
                    lines.add(line.id)
        return cuts%2
    def get_neighbour(self, line, column, layer) :
        return Box.grid.get_box(self.x+line, self.y+column, self.z+layer)
    def know_your_free_neighbours(self) :
        raise NotImplementedError("Abstract Grid.know_your_free_neighbours")
    def get_center(self) :
        raise NotImplementedError("Abstract Grid.get_center")
    @staticmethod
    def get_box(x, y, z) : return Box.grid.get_box(x, y, z)
    @staticmethod
    def get_location(x, y, z) :
        raise NotImplementedError("Abstract @staticmethod Grid.get_location")
    def __str__(self) :
        return("Box(%s,%s,%s)" % (self.x, self.y, self.z))
    def __repr__(self) :
        return("Box(%r, %r, %r)" % (self.x, self.y, self.z, ))
    def sq(self) :
        return("⬜" if self.free else "⬛")
    def viz(self) :
        return('"' + ' '.join((str(self.x), str(self.y), str(self.z))) + '"')

class SquareBox(Box) :
    def get_center(self) :
        return Point((self.x+0.5)*Box.grid.padx+Box.grid.minx, \
                    (self.y+0.5)*Box.grid.pady+Box.grid.miny, \
                    (self.z)*Box.grid.padz+Box.grid.minz) # get the "floor" of the box
    @staticmethod
    def get_location(x, y, z) :
        return (int((x-Box.grid.minx)//Box.grid.padx), \
                int((y-Box.grid.miny)//Box.grid.pady), \
                int((z-Box.grid.minz)//Box.grid.padz))
    def know_your_free_neighbours(self) :
        neighbours = []
        if self.x != 0 :
            neighbours.append(self.get_neighbour(-1, 0, 0))
        if self.x != Box.grid.nb_lines-1 :
            neighbours.append(self.get_neighbour(1, 0, 0))
        if self.y != 0 :
            neighbours.append(self.get_neighbour(0, -1, 0))
        if self.y != Box.grid.nb_columns-1 :
            neighbours.append(self.get_neighbour(0, 1, 0))
        self.free_neighbours = filter(lambda b : b.free, neighbours)

class HexaBox(Box) :
    def get_center(self) :
        return Point(self.x*Box.grid.apo2 + Box.grid.minx \
                    + (Box.grid.apo if self.y % 2 else 0), # décalage pair-impair \
                    Box.grid.rad05 + self.y*Box.grid.rad15 + Box.grid.miny, \
                    self.z*Box.grid.height + Box.grid.minz)
    @staticmethod
    def get_location(x, y, z) :
        x -= Box.grid.minx
        y -= Box.grid.miny
        z -= Box.grid.minz
        rad = Box.grid.rad
        apo = Box.grid.apo
        apo2 = Box.grid.apo2
        rad05 = Box.grid.rad05
        rad15 = Box.grid.rad15
        x_possible = int(x // apo2)
        y_possible = int(y // rad15)
        if (y % rad15) < rad : # "core" of an hexagon
            y_found = y_possible
            if y_found % 2 : # odd column (first hexagon is full)
                x_found = int(x // apo2)
            else : # even column (first hexagon is halved)
                x_found = int((x+apo) // apo2)
        else :
            # "boundaries" of hexagons
            # reduce all cases to one which solution is known
            #PATTERN_V = True
            #CASE_SLASH = True
            #SUBCASE_ABOVE = True
            pattern = bool((y // rad15) % 2)
            x_reduced = x % apo2
            case = (x_reduced < apo) != pattern
            if (pattern == case) : slope = Box.grid.slope
            else :                 slope = - Box.grid.slope
            if pattern :
                if case : y_reduced = rad05
                else :    y_reduced = - rad05
            else :
                if case : y_reduced = rad15
                else :    y_reduced = rad05
            y_reduced += x_reduced*slope + y_possible*rad15
            if y_reduced > y :
                if not pattern and not case : x_possible += 1
            else :
                y_possible += 1
                if pattern and case : x_possible += 1
            x_found = x_possible
            y_found = y_possible
        z_found = int(z // Box.grid.height)
        return (x_found, y_found, z_found)
    def know_your_free_neighbours(self) :
        neighbours = []
        if self.x != 0 :
            neighbours.append(self.get_neighbour(-1, 0, 0))
            if self.y % 2 == 0 :
                if self.y != 0 :
                    neighbours.append(self.get_neighbour(-1, -1, 0))
                if self.y != Box.grid.nb_columns-1 :
                    neighbours.append(self.get_neighbour(-1, 1, 0))
        if self.y != 0 :
            neighbours.append(self.get_neighbour(0, -1, 0))
        if self.x != Box.grid.nb_lines-1 :
            neighbours.append(self.get_neighbour(1, 0, 0))
            if self.y % 2 == 1 :
                if self.y != 0 :
                    neighbours.append(self.get_neighbour(1, -1, 0))
                if self.y != Box.grid.nb_columns-1 :
                    neighbours.append(self.get_neighbour(1, 1, 0))
        if self.y != Box.grid.nb_columns-1 :
            neighbours.append(self.get_neighbour(0, 1, 0))
        self.free_neighbours = filter(lambda b : b.free, neighbours)

class GridBox(type) :
    def __new__(self, grid) :
        Box.grid = grid
        if type(grid) is HexaGrid :
            Class = HexaBox
        else :
            Class = SquareBox
        return Class

class Grid(object) :
    # "Abstract" methods :
    #   - do_pavement
    #   - discretise_line
    def __init__(self, model, cutter) :
        self.model = model
        self.radius = cutter.radius
        self.height = self.radius * 2 # TODO: enfoncement
        self.nb_layers = ceil( (model.maxz - model.minz) / self.height) +1 # TODO: +1 skypocket
        self.layers = []
        #self.heights = Grid.floatrange(model.minz, model.maxz, self.height, False)
        self.heights = [model.minz+i*self.height for i in range(self.nb_layers)]
        self.minx, self.miny, self.minz = model.minx, model.miny, model.minz
        self.maxx, self.maxy, self.maxz = model.maxx, model.maxy, model.maxz
        self.up_vector = Vector(0, 0, 1)
        self.Box = GridBox(self)
        self.do_pavement()
    def do_pavement(self) :
        raise NotImplementedError("Abstract Grid.do_pavement")
    def free_memory(self) :
        for layer in self.layers :
            for line in layer :
                for box in line :
                    box.free_memory()
                    del box
                del line
            del layer
        del self.layers
    def iterate_on_layers(self) :
        for layer in xrange(self.nb_layers) :
            self.layers.append([[self.Box(x, y, layer) \
                                for y in xrange(self.nb_columns)] \
                                for x in xrange(self.nb_lines)])
            yield self.heights[layer]
    def get_box(self, line, column, layer) :
        try :
            return self.layers[layer][line][column]
        except IndexError :
            raise IndexError("get_box(%s, %s, %s) cannot be resolved (max : %s, %s, %s)" % \
                    (line, column, layer, self.nb_lines-1, self.nb_columns-1, self.nb_layers-1))
    def discretise_point(self, p) :
        return self.get_box(*self.Box.get_location(p.x, p.y, p.z))
    def discretise_line(self, line) :
        raise NotImplementedError("Abstract Grid.discretise_line")
    def display_complete_layer(self, z) :
        index = self.heights.index(z)
        log.info("\nLayer %d/%d : " % (index+1, self.nb_layers))
        for column in range(self.nb_columns) :
            log.info(' '.join(self.get_box(line, column, index).sq() \
                for line in range(self.nb_lines)))
    def draw_contour_PBM(self) :
        for nb_layer in range(self.nb_layers) :
            with open("cutter_%d.pbm"%(nb_layer+1), "w") as cutterfile :
                cutterfile.write("P1\n")
                cutterfile.write("# %d/%d\n"%(nb_layer+1, self.nb_layers))
                cutterfile.write("%d %d\n"%(self.nb_columns, self.nb_lines))
                for line in range(self.nb_lines) :
                    cutterfile.write(' '.join(['0' if box.free else '1'
                            for box in self.layers[nb_layer][line]])+'\n')

class SquareGrid(Grid) :
    def do_pavement(self) :
        self.diameter = self.radius * sqrt(2)
        self.rangex = [self.model.minx+i*self.diameter for i in xrange(int((self.model.maxx-self.model.minx)/self.diameter+1))]
        self.rangey = [self.model.minx+i*self.diameter for i in xrange(int((self.model.maxy-self.model.miny)/self.diameter+1))]
        self.nb_lines = len(self.rangex)
        self.nb_columns = len(self.rangey)
        self.padx, self.pady, self.padz = self.diameter, self.diameter, self.radius*2 # TODO
    def discretise_line(self, line) : # TODO: OPTIMISER !!!!! (pas de boucle, un seul cas avec 2 passages, calcul de l'intervale sans list.index)
        c1, c2 = self.discretise_point(line.p1), self.discretise_point(line.p2)
        c1.lines.append(line)
        c2.lines.append(line)
        if line.p1.x == line.p2.x and line.p1.x in self.rangex :
            x = self.rangex.index(line.p1.x)
            if c1.y > c2.y : c1, c2 = c2, c1
            for y in range(c1.y, c2.y+1) :
                self.get_box(x, y, c1.z).lines.append(line)
        elif line.p1.y == line.p2.y and line.p1.y in self.rangey :
            y = self.rangey.index(line.p1.y)
            if c1.x > c2.x : c1, c2 = c2, c1
            for x in range(c1.x, c2.x+1) :
                self.get_box(x, y, c1.z).lines.append(line)
        else :
            ordx = 1
            if c1.x > c2.x : c1, c2 = c2, c1
            for x in self.rangex[c1.x+1:c2.x] :
                sec, d = line.get_intersection(Line(Point(x, self.miny, line.p1.z), \
                                                Point(x, self.maxy, line.p1.z)))
                if sec is None : break
                dsec = self.Box.get_location(0, sec.y, 0)
                self.get_box(c1.x+ordx, dsec[1], c1.z).lines.append(line)
                ordx += 1
            ordy = 1
            if c1.y > c2.y : c1, c2 = c2, c1
            for y in self.rangey[c1.y+1:c2.y] :
                sec, d = line.get_intersection(Line(Point(self.minx, y, line.p1.z), \
                                                Point(self.maxx, y, line.p1.z)))
                if sec is None : break
                dsec = self.Box.get_location(sec.x, 0, 0)
                self.get_box(dsec[0], c1.y+ordy, c1.z).lines.append(line)
                ordy += 1

class HexaGrid(Grid) :
    def do_pavement(self) :
        self.rad = self.radius
        self.apo = sqrt(3)/2*self.rad
        self.apo2 = self.apo*2
        self.rad05 = self.rad*0.5
        self.rad15 = self.rad*1.5
        self.slope = 2./sqrt(3.)
        self.nb_lines = int((self.model.maxx - self.model.minx + self.apo) // self.apo2)+1
        self.nb_columns = int((self.model.maxy - self.model.miny + self.rad05) // self.rad15)+1
    def discretise_line(self, line) :
        if line.p1 > line.p2 :
            line.p1, line.p2 = line.p2, line.p1
        b1, b2 = self.discretise_point(line.p1), self.discretise_point(line.p2)
        b1.lines.append(line)
        b2.lines.append(line)
        delta = line.p2.sub(line.p1)
        if abs(delta.x) < epsilon : # vertical line
            slope = 10**15
            iterate_on = 'Y'
        elif abs(delta.y) < epsilon : # horizontal line
            slope = 0
            iterate_on = 'X'
        else :
            slope = delta.y/delta.x
            iterate_on = 'X' if slope < self.slope else 'Y'
        z = line.p1.z
        if iterate_on == 'X' :
            apo = self.apo
            dy = slope*apo
            next_apo = (line.p1.x//apo+1)*apo
            y = line.p1.y + (next_apo - line.p1.x)*slope
            x = next_apo
            max_x = line.p2.x
            while x < max_x :
                box = self.get_box(*self.Box.get_location(x, y, z))
                box.lines.append(line)
                x += apo
                y += dy
        else : # iterate on Y
            rad = self.rad
            dx = 1/slope/rad
            next_rad = (line.p1.y//rad+1)*rad
            x = line.p1.x + (next_rad - line.p1.y)*(1/slope)
            y = next_rad
            max_y = line.p2.y
            while y < max_y :
                box = self.get_box(*self.Box.get_location(x, y, z))
                box.lines.append(line)
                y += rad
                x += dx

class Contour2dCutter(object) :

    def __init__(self, cutter, models, path_processor) :
        self.cutter = cutter
        self.model = models[0]
        self.pa = path_processor
        #self.grid = SquareGrid(models[0], cutter)
        self.grid = HexaGrid(models[0], cutter)

    def GenerateToolPath(self, callback=None) :
        import cProfile
        from datetime import datetime
        result = [None]
        cProfile.runctx('result[0] = self.gen(callback)', globals(), locals(), 'profile_'+datetime.now().strftime('%H:%M')+'.prof')
        return result[0]

    def gen(self, callback=None) :
        triangles = self.model.triangles()
        equations = {}

        # patching Point
        Point.__hash__ = lambda self : hash(' '.join(
                                [str(self.x), str(self.y), str(self.z)]))
        Point.__eq__ = lambda self, other : _is_near(self.x, other.x) \
                                        and _is_near(self.y, other.y) \
                                        and _is_near(self.z, other.z)

        #patching Line
        def equation2D(x1, y1, x2, y2):
            if abs(x1 - x2) < epsilon : # care of floating-point imprecision
                return "x=%s" % x1
            else:
                a = (y2-y1)/(x2-x1)
                return "y=%sx+%s" % (a, y2-a*x2)
        Line.__hash__ = lambda self : hash(' '.join([
                equation2D(self.p1.x, self.p1.y, self.p2.x, self.p2.y), \
                equation2D(self.p1.x, self.p1.z, self.p2.x, self.p2.z), \
                equation2D(self.p1.y, self.p1.z, self.p2.y, self.p2.z)]))
        Line.__eq__ = lambda self, other : hash(self) == hash(other)

        num_of_layers = self.grid.nb_layers
        current_layer_num = 1
        quit_requested = False
        progress_counter = ProgressCounter(num_of_layers, callback)

        for height in self.grid.iterate_on_layers() :

            if callback and callback(text="Contour2dCutter: projecting " \
                        + "layer %d/%d" % (current_layer_num, num_of_layers)):
                # cancel requested
                quit_requested = True
                break

            equations.clear()
            planeInf = Plane(Point(0, 0, height), self.grid.up_vector)
            planeSup = Plane(Point(0, 0, INFINITE), self.grid.up_vector)
            lines = [ligne for triangle in triangles for ligne in \
                self.triangleBetweenTwoPlanes(triangle, planeInf, planeSup)]
            for line in lines :
                if line not in equations :
                    equations[line] = [line]
                else :
                    equations[line].append(line)
            to_project = []
            for equation in iter(equations) :
                lines = equations[equation]
                uniques = self.mergeRiddanceLines(lines)
                to_project.extend(uniques)
            for line in to_project :
                projection = planeInf.get_line_projection(line)
                self.grid.discretise_line(projection)
            #self.grid.display_complete_layer(height)

            progress_counter.increment()

        if quit_requested : return

        #self.grid.draw_contour_PBM()
        self.pa.initialise(self.grid, callback)
        self.pa.do_path()
        #self.grid.free_memory()
        return self.pa.paths

    @staticmethod
    def mergeRiddanceLines(lines) :
        """\
        returns not-layered part of the lines
        lines must be aligned
        """
        if len(lines) <= 1 : return lines
        points = {}
        for line in lines :
            if line.p1 > line.p2 : line.p1, line.p2 = line.p2, line.p1
            if line.p1 not in points :
                points[line.p1] = 1
            else : points[line.p1] += 1
            if line.p2 not in points :
                points[line.p2] = -1
            else : points[line.p2] -= 1
        order = sorted(points.keys())
        counter = 0
        start = None
        uniques = []
        for point in order :
            counter += points[point]
            if counter == 1 and start is None :
                start = point
            elif counter != 1 and start is not None :
                if cmp(start, point) :
                    uniques.append(Line(start, point))
                start = None
        return uniques

    @staticmethod
    def triangleBetweenTwoPlanes(triangle, planeInf, planeSup) :
        """\
        returns the polygon that is the cut of a triangle between two planes
        planes must be horizontals
        planeInf must be beneath planeSup
        belonging is exclusive for both planes
        """
        if triangle.center.z == planeInf.p.z and \
                (triangle.normal.z == 1 or triangle.normal.z == -1) :
            return []
        return Contour2dCutter.makePolygonFrom(
                Contour2dCutter.lineBetweenTwoPlanes(triangle.e1, planeInf, planeSup),
                Contour2dCutter.lineBetweenTwoPlanes(triangle.e2, planeInf, planeSup),
                Contour2dCutter.lineBetweenTwoPlanes(triangle.e3, planeInf, planeSup))

    @staticmethod
    def lineBetweenTwoPlanes(line, planeInf, planeSup) :
        """\
        returns a tuple that contains :
         - the case encountered (aa, ab, ac, bb, bc, cc)
         - the line between the two planes (cases ab-c et bb-c) else None
        """
        p1, p2 = (line.p1, line.p2) if line.p1.z <= line.p2.z else (line.p2, line.p1)
        direction = p2.sub(p1)
        if p1.z <= planeInf.p.z :
            if p2.z <= planeInf.p.z :
                return ('aa', None)
            elif p2.z < planeSup.p.z :
                sec, l = planeInf.intersect_point(direction, p1)
                return ('ab', Line(sec, p2))
            elif p2.z >= planeSup.p.z :
                sec1, l1 = planeInf.intersect_point(direction, p1)
                sec2, l2 = planeSup.intersect_point(direction, p2)
                return ('ac', Line(sec1, sec2))
        elif p1.z < planeSup.p.z :
            if p2.z < planeSup.p.z :
                return ('bb', line)
            elif p2.z >= planeSup.p.z :
                sec, l = planeSup.intersect_point(direction, p2)
                return ('bc', Line(p2, sec))
        elif p1.z >= planeSup.p.z :
            return ('cc', None)

    @staticmethod
    def makePolygonFrom(r1, r2, r3) :
        """\
        returns the polygon created from edges provided
        parameters should be results from lineBetweenTwoPlanes
        """
        case = [r1, r2, r3]
        case.sort(key=lambda r : r[0]) # lexicographic order
        if case[0][0] == 'aa' :
            if case[1][0] == 'aa' and case[2][0] == 'aa' : # Pattern 1
                return []
            if case[1][0] == 'ab' and case[2][0] == 'ab' : # Pattern 2
                return [case[1][1], case[2][1],
                        Line(case[1][1].p1, case[2][1].p1)]
            if case[1][0] == 'ac' and case[2][0] == 'ac' : # Pattern 3
                return [case[1][1], case[2][1],
                        Line(case[1][1].p1, case[2][1].p1),
                        Line(case[1][1].p2, case[2][1].p2)]
        elif case[0][0] == 'ab' :
            if case[1][0] == 'ab' and case[2][0] == 'bb' : # Pattern 4
                return [case[0][1], case[1][1], case[2][1],
                        Line(case[0][1].p1, case[1][1].p1)]
            if case[1][0] == 'ac' and case[2][0] == 'bc' : # Pattern 5
                return [case[0][1], case[1][1], case[2][1],
                        Line(case[1][1].p1, case[0][1].p1),
                        Line(case[1][1].p2, case[2][1].p2)]
        elif case[0][0] == 'ac' :
            if case[1][0] == 'ac' and case[2][0] == 'cc' : # Pattern 6
                return [case[0][1], case[1][1],
                        Line(case[0][1].p1, case[1][1].p1),
                        Line(case[0][1].p2, case[1][1].p2)]
        elif case[0][0] == 'bb' :
            if case[1][0] == 'bb' and case[2][0] == 'bb' : # Pattern 7
                return [case[0][1], case[1][1], case[2][1]]
            if case[1][0] == 'bc' and case[2][0] == 'bc' : # Pattern 8
                return [case[0][1], case[1][1], case[2][1],
                        Line(case[1][1].p2, case[2][1].p2)]
        elif case[0][0] == 'bc' :
            if case[1][0] == 'bc' and case[2][0] == 'cc' : # Pattern 9
                return [case[0][1], case[1][1],
                        Line(case[0][1].p2, case[1][1].p2)]
        elif case[0][0] == 'cc' :
            if case[1][0] == 'cc' and case[2][0] == 'cc' : # Pattern 10
                return []
