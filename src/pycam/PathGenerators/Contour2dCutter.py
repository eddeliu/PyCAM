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
import heapq # for Dijkstra's algorithm optimisation (using priority queue)

import traceback # TODO: remove

from pycam.Utils import ProgressCounter
import pycam.Utils.log

log = pycam.Utils.log.get_logger()

class Box(object) :
    # "Abstract" methods :
    #   - get_center
    #   - get_location (static)
    #   - know_your_free_neighbours
    #   - know_your_pocket_neighbourhood
    def __init__(self, x, y, z) :
        self.x = x
        self.y = y
        self.z = z
        self.index = z*Box.grid.nb_boxes_per_layer + x*Box.grid.nb_columns + y
        self._distance = 0 # used by Dijkstra's algorithm
        self.lines = []
        self.scan = False
        self.pocket = None
        self.free_neighbours = []
    def __lt__(self, other) :
        # used in Dijkstra's algorithm implementation by heapq
        return self._distance < other._distance
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
    def know_your_pocket_neighbourhood(self) :
        raise NotImplementedError("Abstract Grid.know_your_pocket_neighbourhood")
    def get_center(self) :
        raise NotImplementedError("Abstract Grid.get_center")
    @staticmethod
    def get_box(x, y, z) : return Box.grid.get_box(x, y, z)
    @staticmethod
    def get_location(x, y, z) :
        raise NotImplementedError("Abstract @staticmethod Grid.get_location")
    def __str__(self) :
        return("Box({x},{y}{z})".format(self))
    def __repr__(self) :
        return("Box(%r, %r, %r)" % (self.x, self.y, self.z, ))
    def sq(self) :
        return("⬜" if self.libre else "⬛")
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

class HexaBox(Box) :
    def get_center(self) :
        return Point(self.x*Box.grid.apo2 + Box.grid.minx \
                    + (Box.grid.apo if self.y % 2 else 0), # décalage pair-impair \
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
        self.free_neighbours = filter(lambda v : v.free, neighbours)
    def know_your_pocket_neighbourhood(self) :
        # base algorithm : Dijkstra
        # modifications :
        #  - iteration using flood
        #  - all distances are 1
        pop = heapq.heappop
        push = heapq.heappush
        m = Box.grid.matrix
        m[self.index][self.index] = (0, self)
        for neighbour in self.free_neighbours : neighbour._distance = 1
        neighbours_left = heapq.heapify(self.free_neighbours[:])
        while neighbours_left :
            nearest = pop(neighbours_left)
            new_distance = nearest._distance + 1
            for neighbour in nearest.free_neighbours :
                previous_distance = m[self.index][neighbour.index][0]
                if previous_distance is None :
                    neighbour._distance = new_distance
                    m[self.index][neighbours.index] = (new_distance, nearest)
                    push(neighbours_left, neighbour)
                # else the path would be anyway equally long or longer
    def export_your_graph_of_shortest_paths(self) :
        m = Box.grid.matrix
        this = self.viz()
        ind = self.index
        pathfile = open('dijkstra.dot', 'w')
        pathfile.write('graph G {\n')
        pathfile.write('\tnode [shape=hexagon]\n')
        if self.pocket :
            for vertex in self.pocket.boxes :
                pathfile.write('\t%s -- %s [%s];\n'%(this, vertex.viz(), m[ind][vertex.index]))
        pathfile.write('}\n')
        pathfile.close()
        gridfile = open('grid.dot', 'w')
        encountered = set()
        gridfile.write('graph G {\n')
        gridfile.write('\tnode [shape=hexagon]\n')
        if self.pocket :
            for vertex in self.pocket.boxes :
                gridfile.write('\t%s -- {%s};\n'%(vertex.viz(), ' '.join(map(Box.viz, filter(lambda n : n not in encountered, vertex.free_neighbours)))))
                encountered.add(vertex)
        gridfile.write('}\n')
        gridfile.close()

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
        self.diameter = 2*cutter.radius
        self.height = self.diameter # TODO
        self.nb_layers = ceil( (model.maxz - model.minz) / self.height) # integer division TODO: +1 ? sky pocket
        self.layers = []
        self.heights = Grid.floatrange(model.minz, model.maxz, self.height, False)
        self.minx, self.miny, self.minz = model.minx, model.miny, model.minz
        self.maxx, self.maxy, self.maxz = model.maxx, model.maxy, model.maxz
        self.up_vector = Vector(0, 0, 1)
        self.Box = GridBox(self)
        self.do_pavement()
        self.nb_boxes_per_layer = self.nb_lines*self.nb_columns
        self.nb_boxes_total = self.nb_boxes_per_layer*self.nb_layers
        self.matrix = [[(None, None)]*self.nb_boxes_total]*self.nb_boxes_total
    def do_pavement(self) :
        raise NotImplementedError("Abstract Grid.do_pavement")
    def iterate_on_layers(self) :
        for height in self.heights :
            z = self.heights.index(height) # TODO improvement
            self.layers.append([[self.Box(x, y, z) \
                                for y in range(self.nb_columns)] \
                                for x in range(self.nb_lines)])
            yield height
    def get_box(self, line, column, layer) :
        try :
            return self.layers[layer][line][column]
        except IndexError :
            raise IndexError("get_box(%s, %s, %s) cannot be resolved (max : %s, %s, %s) !" % \
                    (line, column, layer, self.nb_lines-1, self.nb_columns-1, self.nb_layers-1))
    def discretise_point(self, p) :
        return self.get_box(*self.Box.get_location(p.x, p.y, p.z))
    def discretise_line(self, line) :
        raise NotImplementedError("Abstract Grid.discretise_line")
    def display_complete_layer(self, z) :
        index = self.heights.index(z)
        print("\nLayer %d/%d : " % (index+1, self.nb_layers))
        for column in range(self.nb_columns) :
            print(' '.join(self.get_box(index, line, column).sq() \
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
    @staticmethod
    def floatrange(start, end, inc, reverse=False) :
        l = []
        n = start
        while n < end :
            l.append(n)
            n += inc
        if reverse :
            l.reverse()
        return l

class SquareGrid(Grid) :
    def calculerPavage(self) :
        self.rangex = Grille.floatrange(self.model.minx, self.model.maxx, self.diametre)
        self.rangey = Grille.floatrange(self.model.miny, self.model.maxy, self.diametre)
        self.nb_lines = len(self.rangex)
        self.nb_columns = len(self.rangey)
        self.padx, self.pady, self.padz = self.diametre, self.diametre, self.diametre # TODO
    def discretise_line(self, line) : # TODO: OPTIMISER !!!!! (pas de boucle, un seul cas avec 2 passages, calcul de l'intervale sans list.index)
        # on discrétise les deux extrémités du segment
        c1, c2 = self.discretiser_point(ligne.p1), self.discretiser_point(ligne.p2)
        c1.lignes.append(ligne)
        c2.lignes.append(ligne)
        if ligne.p1.x == ligne.p2.x and ligne.p1.x in self.rangex : # la ligne appartient au pavage vertical
            # alors on va ajouter toutes les paires de cases qu'elle rencontre
            x = self.rangex.index(ligne.p1.x)
            if c1.y > c2.y : c1, c2 = c2, c1
            for y in range(c1.y, c2.y+1) :
                self.get_case(c1.z, x, y).lignes.append(ligne)
        elif ligne.p1.y == ligne.p2.y and ligne.p1.y in self.rangey : # la ligne appartient au pavage horizontal
            y = self.rangey.index(ligne.p1.y)
            if c1.x > c2.x : c1, c2 = c2, c1
            for x in range(c1.x, c2.x+1) :
                self.get_case(c1.z, x, y).lignes.append(ligne)
        else : # sinon la ligne est simplement verticale, horizontale ou oblique
            # alors on va discrétiser toutes les cases traversées
            # d'abord les intersections verticales
            ordx = 1
            if c1.x > c2.x : c1, c2 = c2, c1 # pour itérer dans le bon sens
            for x in self.rangex[c1.x+1:c2.x] :
                # pour cela on calcule l'intersection de la ligne avec la grille des X
                sec, d = ligne.get_intersection(Line(Point(x, self.miny, ligne.p1.z), \
                                                Point(x, self.maxy, ligne.p1.z)))
                if sec is None : # il n'y a pas intersection : les deux lignes sont parallèles et distinctes
                    break # la ligne sera traitée par l'itération verticale
                # puis on trouve à quelle colonne l'intersection appartient
                dsec = self.Case.get_emplacement(0, sec.y, 0)
                # et on discrétise alors la case correspondante
                self.get_case(c1.z, c1.x+ordx, dsec[1]).lignes.append(ligne)
                # le [1] correspond à la composante Y de la coordonnée
                ordx += 1
            # puis celles horizontales de la même manière
            ordy = 1
            if c1.y > c2.y : c1, c2 = c2, c1
            for y in self.rangey[c1.y+1:c2.y] :
                sec, d = ligne.get_intersection(Line(Point(self.minx, y, ligne.p1.z), \
                                                Point(self.maxx, y, ligne.p1.z)))
                if sec is None :
                    break # la ligne est horizontale, et donc déjà traitée
                dsec = self.Case.get_emplacement(sec.x, 0, 0)
                self.get_case(c1.z, dsec[0], c1.y+ordy).lignes.append(ligne)
                ordy += 1

class HexaGrid(Grid) :
    def do_pavement(self) :
        self.rad = self.diameter/2
        self.apo = sqrt(3)/2*self.rad
        self.apo2 = self.apo*2
        self.rad05 = self.rad*0.5
        self.rad15 = self.rad*1.5
        self.slope = 2./sqrt(3.)
        self.nb_lines = int((self.model.maxx - self.model.minx + self.apo) // self.apo2)+1
        self.nb_columns = int((self.model.maxy - self.model.miny + self.rad05) // self.rad15)+1
        self.height = self.diameter # TODO: take account of overlapping, ...
    def discretise_line(self, line) :
        if line.p1 > line.p2 : line.p1, line.p2 = line.p2, line.p1
        b1, b2 = self.discretise_point(line.p1), self.discretise_point(line.p2)
        b1.lines.append(line)
        b2.lines.append(line)
        delta = line.p2.sub(line.p1)
        if abs(delta.x) < epsilon : # vertical line
            slope = 0
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
            dx = slope*rad
            next_rad = (line.p1.y//rad+1)*rad
            x = line.p1.x + (next_rad - line.p1.y)*slope
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
        Line.__hash__ = lambda self : hash(' '.join(
                [equation2D(self.p1.x, self.p1.y, self.p2.x, self.p2.y), \
                equation2D(self.p1.x, self.p1.z, self.p2.x, self.p2.z), \
                equation2D(self.p1.y, self.p1.z, self.p2.y, self.p2.z)]))
        Line.__eq__ = lambda self, other : hash(self) == hash(other)

        for height in self.grid.iterate_on_layers() :
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
        self.grid.draw_contour_PBM()
        self.pa.initialise(self.grid)
        self.pa.do_path()

        return None

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
