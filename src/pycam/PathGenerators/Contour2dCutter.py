# -*- coding: utf-8 -*-
"""
trucs
"""

from pycam.Geometry.Plane import Plane
from pycam.Geometry.Point import Point, Vector
from pycam.Geometry.Line import Line
from pycam.Geometry.utils import INFINITE, epsilon, ceil, sqrt
def _is_near(x, y):
    return abs(x - y) < epsilon

from pycam.Utils import ProgressCounter
import pycam.Utils.log

log = pycam.Utils.log.get_logger()

class Case(object) :
    # Méthodes "abstraites" :
    #   - interieure (property)
    #   - get_voisin
    #   - get_centre
    #   - get_emplacement (static)
    def __init__(self, x, y, z) :
        self.x = x
        self.y = y
        self.z = z
        self.lignes = []
        self.parcourue = False
        self.poche = None
    @property
    def libre(self) : return len(self.lignes)==0
    @property
    def pleine(self) : return self.poche == None
    @property
    def interieure(self) :
        pass
    @staticmethod
    def ajouter_ligne(self, newline) :
        self.lignes.append(newline)
    def get_voisin(self, tranche, ligne, colonne) :
        pass
    def get_centre(self) :
        pass
    @staticmethod
    def get_case(x, y, z) : return Case.grille.get_case(z, x, y)
    @staticmethod
    def get_emplacement(x, y, z) :
        pass
    def __repr__(self) :
        return("Case(%r, %r, %r, %r, %r, %r)" % \
            (self.x, self.y, self.z, len(self.lignes), self.parcourue, "[P]" if self.poche else "[]"))
    def __str__(self) :
        return("⬜" if self.libre else "⬛")

class QuadraCase(Case) :
    def get_voisin(self, tranche, ligne, colonne) :
        return Case.grille.get_case(self.z+tranche, self.x+ligne, self.y+colonne)
    def get_centre(self) :
        return Point((self.x+0.5)*Case.grille.padx+Case.grille.minx, \
                    (self.y+0.5)*Case.grille.pady+Case.grille.miny, \
                    (self.z)*Case.grille.padz+Case.grille.minz) # obtention du centre du "plancher"
    @property
    def interieure(self) :
        #print("Case : %r" % self)
        p1 = self.get_centre()
        #print("P1 : %r" % p1)
        p2 = Point(Case.grille.minx-1, p1.y, p1.z)
        #print("P2 : %r" % p2)
        fuyante = Line(p1, p2)
        sections = 0
        case = self
        lignes = set() # gère la duplicité des lignes fondues au pavage
        while case.x != 0 :
            case = case.get_voisin(0, -1, 0)
            #print("Etape : %r" % case)
            for ligne in case.lignes :
                if ligne.id not in lignes :
                    #print("\tLigne : %s" % ligne)
                    sec, d = fuyante.get_intersection(ligne)
                    if sec is not None :
                        #print("\t\tIntersection : %s" % sec)
                        sections += 1
                    lignes.add(ligne.id)
                #else :
                    #print("\tLigne deja rencontre : %s" % ligne)
        return sections%2 == 1
    @staticmethod
    def get_emplacement(x, y, z) :
        # on va discrétiser par l'usage de pavé semi-ouverts
        # pour renvoyer les coordonnées entières dans le pavage
        return (int((x-Case.grille.minx)//Case.grille.padx), \
                int((y-Case.grille.miny)//Case.grille.pady), \
                int((z-Case.grille.minz)//Case.grille.padz))

class HexaCase(Case) :
    def get_voisin(self, tranche, ligne, colonne) :
        return Case.grille.get_case(self.x + ligne if self.y%2 == 0 or colonne == 0 else 0, \
                                    self.y+colonne, self.z+tranche)
    def get_centre(self) :
        pass
    @property
    def interieure(self) :
        pass
    @staticmethod
    def get_emplacement(x, y, z) :
        pass

class CaseDeGrille(type) :
    def __new__(self, grille) :
        Case.grille = grille
        if type(grille) is HexaGrille :
            Classe = HexaCase
        else :
            Classe = QuadraCase
        return Classe

class Grille(object) :
    # Méthodes "abstraites" :
    #   - calculerPavage
    def __init__(self, model, cutter) :
        self.model = model
        # on suppose un cutter cubique pour le moment
        self.diametre = 2*cutter.radius # TODO: gérer l'overlapping et le grain
        # découpage du modèle en tranche de hauteur le diamètre du cutter
        self.nb_tranches = ceil( (model.maxz - model.minz) / self.diametre) # division entière
        self.tranches = []
        # génération de la hauteur basse pour chacune des tranches, en vue de la découpe
        # par ordre descendant
        self.hauteurs = Grille.floatrange(model.minz, model.maxz, self.diametre, False)
        self.minx, self.miny, self.minz = model.minx, model.miny, model.minz
        self.maxx, self.maxy, self.maxz = model.maxx, model.maxy, model.maxz
        self.vecteur_haut = Vector(0, 0, 1)
        self.Case = CaseDeGrille(self)
        self.calculerPavage()
    def calculerPavage(self) :
        pass
    def iterer_tranches(self) :
        for hauteur in self.hauteurs :
            z = self.hauteurs.index(hauteur)
            self.tranches.append([[self.Case(x, y, z) \
                                for y in range(self.nb_colonnes)] \
                                for x in range(self.nb_lignes)])
            yield hauteur
    def get_case(self, tranche, ligne, colonne) :
        try :
            return self.tranches[tranche][ligne][colonne]
        except IndexError :
            modif = False
            if   ligne == self.nb_lignes : ligne -= 1; modif = True
            elif ligne == -1 :             ligne = 0;  modif = True
            if   colonne == self.nb_colonnes : colonne -= 1; modif = True
            elif colonne == -1 :               colonne = 0;  modif = True
            if not modif : print("Erreur : get_case(%r, %r, %r, %r) ne peut être résolu" % \
                                    (self, tranche, ligne, colonne)); return
            return self.get_case(tranche, ligne, colonne)
    def discretiser_point(self, p) :
        x, y, z = self.Case.get_emplacement(p.x, p.y, p.z)
        return self.get_case(z, x, y)
    def afficherTrancheComplete(self, z) :
        index = self.hauteurs.index(z)
        print("\nCouche %d/%d : " % (index+1, self.nb_tranches))
        for ligne in range(self.nb_lignes) :
            print ' '.join(str(self.get_case(index, ligne, colonne)) \
                for colonne in range(self.nb_colonnes))
    def dessinerContourPBM(self) :
        for nb_tranche in range(self.nb_tranches) :
            with open("cutter_%d.pbm"%(nb_tranche+1), "w") as fichier :
                fichier.write("P1\n")
                fichier.write("# %d/%d\n"%(nb_tranche+1, self.nb_tranches))
                fichier.write("%d %d\n"%(self.nb_colonnes, self.nb_lignes))
                for ligne in range(self.nb_lignes) :
                    fichier.write(' '.join(['0' if case.libre else '1'
                            for case in self.tranches[nb_tranche][ligne]])+'\n')
    def discretiser_ligne(self, l1) :
        # on discrétise les deux extrémités du segment
        c1, c2 = self.discretiser_point(l1.p1), self.discretiser_point(l1.p2)
        c1.lignes.append(l1)
        c2.lignes.append(l1)
        if l1.p1.x == l1.p2.x and l1.p1.x in self.rangex : # la ligne appartient au pavage vertical
            # alors on va ajouter toutes les paires de cases qu'elle rencontre
            x = self.rangex.index(l1.p1.x)
            if c1.y > c2.y : c1, c2 = c2, c1
            for y in range(c1.y, c2.y+1) :
                self.get_case(c1.z, x, y).lignes.append(l1)
        elif l1.p1.y == l1.p2.y and l1.p1.y in self.rangey : # la ligne appartient au pavage horizontal
            y = self.rangey.index(l1.p1.y)
            if c1.x > c2.x : c1, c2 = c2, c1
            for x in range(c1.x, c2.x+1) :
                self.get_case(c1.z, x, y).lignes.append(l1)
        else : # sinon la ligne est simplement verticale, horizontale ou oblique
            # alors on va discrétiser toutes les cases traversées
            # d'abord les intersections verticales
            ordx = 1
            if c1.x > c2.x : c1, c2 = c2, c1 # pour itérer dans le bon sens
            for x in self.rangex[c1.x+1:c2.x] :
                # pour cela on calcule l'intersection de la ligne avec la grille des X
                sec, d = l1.get_intersection(Line(Point(x, self.miny, l1.p1.z), \
                                                Point(x, self.maxy, l1.p1.z)))
                if sec is None : # il n'y a pas intersection : les deux lignes sont parallèles et distinctes
                    break # la ligne sera traitée par l'itération verticale
                # puis on trouve à quelle colonne l'intersection appartient
                dsec = self.Case.get_emplacement(0, sec.y, 0)
                # et on discrétise alors la case correspondante
                self.get_case(c1.z, c1.x+ordx, dsec[1]).lignes.append(l1)
                # le [1] correspond à la composante Y de la coordonnée
                ordx += 1
            # puis celles horizontales de la même manière
            ordy = 1
            if c1.y > c2.y : c1, c2 = c2, c1
            for y in self.rangey[c1.y+1:c2.y] :
                sec, d = l1.get_intersection(Line(Point(self.minx, y, l1.p1.z), \
                                                Point(self.maxx, y, l1.p1.z)))
                if sec is None :
                    break # la ligne est horizontale, et donc déjà traitée
                dsec = self.Case.get_emplacement(sec.x, 0, 0)
                self.get_case(c1.z, dsec[0], c1.y+ordy).lignes.append(l1)
                ordy += 1
    @staticmethod
    def floatrange(debut, fin, increment, inverse=False) :
        l = []
        n = debut
        while n < fin :
            l.append(n)
            n += increment
        if inverse :
            l.reverse()
        return l

class QuadraGrille(Grille) :
    def calculerPavage(self) :
        self.rangex = Grille.floatrange(self.model.minx, self.model.maxx, self.diametre)
        self.rangey = Grille.floatrange(self.model.miny, self.model.maxy, self.diametre)
        self.nb_lignes = len(self.rangex)
        self.nb_colonnes = len(self.rangey)
        self.padx, self.pady, self.padz = self.diametre, self.diametre, self.diametre # TODO

class HexaGrille(Grille) :
    def calculerPavage(self) :
        # ajouter les diagonales ?
        rayon = self.diametre/2
        apotheme = sqrt(3)/2*rayon
        tronque = rayon - sqrt(rayon**2 - apotheme)
        dia2tronque = self.diametre - 2*tronque
        decalage_y = - tronque
        decalage_x = - rayon
        rangex = []
        rangey = []
        pos_y = model.miny + decalage_y
        self.nb_lignes = 0
        while pos_y < model.maxx :
            pos_y += dia2tronque
            rangey.append(pos_y)
            pos_y += tronque
            rangey.append(pos_y)
            self.nb_colonnes += 1
        pos_x = model.minx + decalage_x
        while pos_x < model.maxx :
            pos_x += rayon
            rangex.append(pos_x)
            pos_x += rayon
            rangex.append(pos_x)
            self.nb_lignes += 1

class Contour2dCutter(object) :

    def __init__(self, cutter, models, path_processor) :
        self.cutter = cutter
        self.model = models[0]
        self.pa = path_processor
        # création de la grille corespondante au modèle en fonction du cutter
        self.grille = QuadraGrille(models[0], cutter)

    def GenerateToolPath(self, callback=None) :
        # le modèle ne permet pas de sélectionner des triangles par Z
        # donc récupération de tous les triangles
        triangles = self.model.triangles()
        #for t in triangles : print "%s\n\t%s\n\t%s\n\t%s" % (t, t.e1, t.e2, t.e3)
        equations = {}

        Point.__hash__ = lambda self : hash(" ".join([str(self.x), str(self.y), str(self.z)])) # TODO
        Point.__eq__ = lambda self, other : _is_near(self.x, other.x) \
                                        and _is_near(self.y, other.y) \
                                        and _is_near(self.z, other.z)

        #i = 1
        for hauteur in self.grille.iterer_tranches() :
            #print ("\n"+"*"*30)*5
            #print i; i+=1
            equations.clear()
            # création du plan qui servira d'espace de projection
            planInf = Plane(Point(0, 0, hauteur), self.grille.vecteur_haut)
            # et du second plan qui sert de borne supérieure
            #planSup = Plane(Point(0, 0, hauteur+self.grille.padz), self.grille.vecteur_haut)
            planSup = Plane(Point(0, 0, INFINITE), self.grille.vecteur_haut)

            # puis section de chaque triangle entre les deux plans
            lignes = [ligne for triangle in triangles for ligne in self.sectionDeTriangleEntreDeuxPlans(triangle, planInf, planSup)]
            # hashage des équations de droite de toutes les lignes pour les regrouper
            for ligne in lignes :
                if ligne not in equations :
                    equations[ligne] = [ligne]
                else :
                    equations[ligne].append(ligne)
            # puis simplification des lignes en vue de leur projection
            a_projeter = []
            
            for equation in iter(equations) :
                compliquees = equations[equation]
                #print "Equation : ", equation
                #print "\tCompliquees : ", compliquees
                simples = self.simplifierLignes(compliquees)
                #print "\tSimples : ", simples
                a_projeter.extend(simples)
            # finalement projection des segments
            for ligne in a_projeter :
                projete = planInf.get_line_projection(ligne)
                #print "Ligne %s a donne %s" % (ligne, projete)
                self.grille.discretiser_ligne(projete)
        self.grille.dessinerContourPBM()
        self.pa.initialiser(self.grille)
        self.pa.traiterGrille()

        return None # TODO!!!

        

    @staticmethod
    def simplifierLignes(lignes) :
        if len(lignes) <= 1 : return lignes
        points = {}
        for ligne in lignes :
            if ligne.p1 > ligne.p2 : ligne.p1, ligne.p2 = ligne.p2, ligne.p1
            if ligne.p1 not in points :
                points[ligne.p1] = 1
            else : points[ligne.p1] += 1
            if ligne.p2 not in points :
                points[ligne.p2] = -1
            else : points[ligne.p2] -= 1
        #print "\t\t%s" % points
        #for p in points : print "\t\t\t%s" % p
        ordre = sorted(points.keys())
        compteur = 0
        segments = []
        #print "\t\t%s" % ordre
        for point in ordre :
            tmpcomp = compteur
            compteur += points[point]
            if tmpcomp == 0 and compteur == 1 \
                or tmpcomp > 1 and compteur == 1 :
                debut = point
            elif tmpcomp == 1 and compteur == 0 :
                segments.append(Line(debut, point))
        return segments

    @staticmethod
    def ligneEntreDeuxPlans(ligne, planInf, planSup) :
        # les plans DOIVENT être horizontaux
        # le planInférieur DOIT être en-dessous du planSupérieur
        # l'appartenance est inclusive
        # la fonction return un tuple qui contient :
        #   - le cas rencontré (aa-b-c, bb-c, cc)
        #   - la ligne comprise entre les deux plans
        #     (cas ab-c et bb-c) sinon None
        #     avec ligne.p1.z <= ligne.p2.z

        direction = ligne.p2.sub(ligne.p1)
        # pour simplifier on va mettre p1 en-dessous de p2
        if ligne.p1.z > ligne.p2.z : ligne.p1, ligne.p2 = ligne.p2, ligne.p1
        # puis tester tous les cas de figure
        if ligne.p1.z < planInf.p.z : # le premier point est sous le plan 1
            if ligne.p2.z < planInf.p.z : # le second point aussi (cas aa)
                return ('aa', None) # donc la ligne est ignorée
            elif ligne.p2.z <= planSup.p.z : # le second point est entre les deux plans (cas ab)
                # alors on va récupérer la moitié qui nous intéresse
                sec, l = planInf.intersect_point(direction, ligne.p1)
                return ('ab', Line(sec, ligne.p2))
            elif ligne.p2.z > planSup.p.z : # le second point est au-dessus du plan 2 (cas ac)
                # alors on va récupérer le tiers qui nous intéresse
                sec1, l1 = planInf.intersect_point(direction, ligne.p1)
                sec2, l2 = planSup.intersect_point(direction, ligne.p2)
                return ('ac', Line(sec1, sec2))
        elif ligne.p1.z <= planSup.p.z : # le premier point est entre les deux plans
            if ligne.p2.z <= planSup.p.z : # le second point aussi (cas bb)
                # alors toute cette ligne est entre les deux plans
                return ('bb', ligne)
            elif ligne.p2.z > planSup.p.z : # le second point est au-dessus du plan 2 (cas bc)
                # alors on va récupérer la moitié qui nous intéresse
                sec, l = planSup.intersect_point(direction, ligne.p2)
                return ('bc', Line(ligne.p2, sec))
        elif ligne.p1.z > planSup.p.z : # le premier point est au-dessus du plan 2 (cas cc)
            return ('cc', None) # donc le second aussi, la ligne est ignorée

    @staticmethod
    def sectionnerTriangle(r1, r2, r3) :
        # r est un tuple qui contient le cas et sa ligne
        cas = [r1, r2, r3]
        cas.sort(key=lambda r : r[0]) # ordre lexicographique : aa,ab,ac,bb,bc,cc
        if cas[0][0] == 'aa' :
            if cas[1][0] == 'aa' and cas[2][0] == 'aa' : # Configuration 1
                return []
            if cas[1][0] == 'ab' and cas[2][0] == 'ab' : # Configuration 2
                return [cas[1][1], cas[2][1],
                        Line(cas[1][1].p1, cas[2][1].p1)]
            if cas[1][0] == 'ac' and cas[2][0] == 'ac' : # Configuration 3
                return [cas[1][1], cas[2][1],
                        Line(cas[1][1].p1, cas[2][1].p1),
                        Line(cas[1][1].p2, cas[2][1].p2)]
        elif cas[0][0] == 'ab' :
            if cas[1][0] == 'ab' and cas[2][0] == 'bb' : # Configuration 4
                return [cas[0][1], cas[1][1], cas[2][1],
                        Line(cas[0][1].p1, cas[1][1].p1)]
            if cas[1][0] == 'ac' and cas[2][0] == 'bc' : # Configuration 5
                return [cas[0][1], cas[1][1], cas[2][1],
                        Line(cas[1][1].p1, cas[0][1].p1),
                        Line(cas[1][1].p2, cas[2][1].p2)]
        elif cas[0][0] == 'ac' :
            if cas[1][0] == 'ac' and cas[2][0] == 'cc' : # Configuration 6
                return [cas[0][1], cas[1][1],
                        Line(cas[0][1].p1, cas[1][1].p1),
                        Line(cas[0][1].p2, cas[1][1].p2)]
        elif cas[0][0] == 'bb' :
            if cas[1][0] == 'bb' and cas[2][0] == 'bb' : # Configuration 7
                return [cas[0][1], cas[1][1], cas[2][1]]
            if cas[1][0] == 'bc' and cas[2][0] == 'bc' : # Configuration 8
                return [cas[0][1], cas[1][1], cas[2][1],
                        Line(cas[1][1].p2, cas[2][1].p2)]
        elif cas[0][0] == 'bc' :
            if cas[1][0] == 'bc' and cas[2][0] == 'cc' : # Configuration 9
                return [cas[0][1], cas[1][1],
                        Line(cas[0][1].p2, cas[1][1].p2)]
        elif cas[0][0] == 'cc' :
            if cas[1][0] == 'cc' and cas[2][0] == 'cc' : # Configuration 10
                return []
        raise NotImplementedError("Un cas de section de triangle n'a pas " +
                                    "été pris en compte : %s-%s-%s\n%s\n%s\n%s"
                                    % (cas[0][0], cas[1][0], cas[2][0], \
                                        cas[0][1], cas[1][1], cas[2][1]))

    @staticmethod
    def sectionDeTriangleEntreDeuxPlans(triangle, planInf, planSup) :
        # les plans DOIVENT être horizontaux
        # le planInférieur DOIT être en-dessous du planSupérieur
        # l'appartenance est inclusive pour planSup,
        # et exclusive pour planInf
        #print triangle
        #print triangle.e1
        #print triangle.e2
        #print triangle.e3
        if triangle.center.z == planInf.p.z and (triangle.normal.z == 1 or triangle.normal.z == -1) :
            #print "Triangle ignoré pour cause d'horizontalité : %s" % triangle
            return []
        return Contour2dCutter.sectionnerTriangle(
                Contour2dCutter.ligneEntreDeuxPlans(triangle.e1, planInf, planSup),
                Contour2dCutter.ligneEntreDeuxPlans(triangle.e2, planInf, planSup),
                Contour2dCutter.ligneEntreDeuxPlans(triangle.e3, planInf, planSup))







