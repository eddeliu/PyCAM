# -*- coding: utf-8 -*-
"""
trucs
"""

from pycam.PathProcessors import BasePathProcessor
# utilise GrilleDiscrete et CaseDeGrille//Case de Contour2dCutter


class Poche(object) :

    id = 1

    def __init__(self) :
        self.id = Poche.id
        Poche.id += 1
        self.sur_poche = None
        self.cases = []
        self.sous_poches = []

    @property
    def fond(self) : return not self.sous_poches

    def __str__(self) :
        return("Poche%s (%s / %s)\n" % (self.id, self.sur_poche.id if self.sur_poche else "None", \
            ','.join(str(sous.id) for sous in self.sous_poches)) \
            + '\n\t'.join(repr(case) for case in self.cases))

class Contour2dFlooder(BasePathProcessor) :

    def __init__(self) :
        super(Contour2dFlooder, self).__init__()
        self.maxx = None
        self.maxy = None
        self.maxz = None
        self.grille = None
        self._tranche = None

    def initialiser(self, grille) :
        self.maxx = grille.nb_lignes - 1
        self.maxy = grille.nb_colonnes - 1
        self.maxz = grille.nb_tranches - 1
        self.grille = grille

    def traiterGrille(self) :
        for tranche in range(self.maxz+1) :
            self.traiterTranche(self.maxz-tranche)
            #break
        self.dessinerPochesPBM() # TODO

    def traiterTranche(self, tranche) :
        self._tranche = self.grille.tranches[tranche]
        #self.afficherTranche()
        #print("Surplomb :");
        #self.traiterSurplombs()
        #self.afficherTranche()
        #print("Exterieurs :");
        exterieures = self.traiterExterieures()
        #print("Poches interieures de la tranche :", tranche)
        interieures = self.traiterInterieures()
        self.lierPoches(exterieures+interieures)
        #print("*" * 50 + "\n")
        #print("POCHES EX TERIEURES")
        #for poche in exterieures :
            #print(poche)
        #print("\nPOCHES IN TERIEURES")
        #for poche in poches :
            #print(poche)
        #print("\n" + "*" * 50)
        #self.afficherTranche()

    def afficherTranche(self) :
        print '\n'*20
        for ligne in self._tranche :
            for colonne in ligne :
                print "%s"%colonne,
            print ''

    def traiterSurplombs(self) : # TODO: utile ?
        for x in range(self.maxx) :
            for y in range(self.maxy) :
                case = self._tranche[x][y]
                if not case.libre :
                    #print("%r est un surplomb" % case)
                    case.parcourue = True
                    while case.z != 0 :
                        case = self.grille.get_case(case.z-1, x, y)
                        case.parcourue = True
                        #print("\t%r est surplombe" % case)

    def traiterExterieures(self) :
        exterieures = []
        case = self._tranche[0][0]
        direction = 0
        # direction 0 : incrémenter x
        # direction 1 : incrémenter y
        # direction 2 : décrémenter x
        # direction 3 : décrémenter y
        while direction != 4 : # TODO: remplacer par 4 for directionnelles ?
            #print("%r case exterieure ?" % case)
            if case.libre and not case.parcourue :
                #print("\t%r est une nouvelle poche" % case)
                exterieures.append(Poche())
                self.flooder(case, exterieures[-1].cases)
            else :
                case.parcourue = True
                case = case.get_voisin(*((0, 1, 0),
                                        (0, 0, 1),
                                        (0, -1, 0),
                                        (0, 0, -1))[direction])
                if direction == 0 and case.x == self.maxx or \
                        direction == 1 and case.y == self.maxy or \
                        direction == 2 and case.x == 0 or \
                        direction == 3 and case.y == 0 :
                    direction += 1
        return exterieures

    def traiterInterieures(self) :
        interieures = []
        for x in range(self.maxx) : # TODO: modifier les range pour éviter de parcourir à nouveau les bordures
            for y in range(self.maxy) :
                case = self._tranche[x][y]
                if not case.parcourue and case.libre :
                    print "Case %r appartient à une poche "%case
                    if not case.interieure :
                        print "\nEX TERIEURE\n\n\n"
                        interieures.append(Poche())
                        self.flooder(case, interieures[-1].cases)
                    else :
                        print "\nIN TERIEURE\n\n\n"
                        self.flooder(case, [])
        print("Fin de traiter poches : %s poches trouvées" % len(interieures))
        return interieures

    def flooder(self, case, poche) :
        # via http://en.wikipedia.org/wiki/Flood_fill#Alternative_implementations
        # implémentation de l'algorithme de remplissage par diffusion
        # avec utilisation d'une pile et emploi de deux boucles (Est et Ouest)
        pile = []
        if not case.libre : return
        pile.append(case)
        while pile :
            centre = pile.pop()
            if centre.libre and not centre.parcourue :
                ouest = est = centre
                while True :
                    if ouest.y == 0 : break
                    _ouest = ouest.get_voisin(0, 0, -1)
                    if not _ouest.libre or _ouest.parcourue : break # TODO: enlever le test parcourue -> rapidité ? exactitude ?
                    else : ouest = _ouest
                while True :
                    if est.y == self.maxy : break
                    _est = est.get_voisin(0, 0, 1)
                    if not _est.libre or _est.parcourue : break # TODO: idem
                    else : est = _est
                for y in range(ouest.y, est.y+1) :
                    c = case.get_case(centre.x, y, centre.z)
                    c.parcourue = True
                    poche.append(c)
                    nord = c.get_voisin(0, 1, 0)
                    if nord.libre : pile.append(nord)
                    sud = c.get_voisin(0, -1, 0)
                    if sud.libre : pile.append(sud)

    def lierPoches(self, poches) :
        for poche in poches :
            if poche.cases[0].z != self.maxz :
                poche.sur_poche = poche.cases[0].get_voisin(1, 0, 0).poche
            for case in poche.cases :
                case.poche = poche

    def dessinerPochesPBM(self) :
        for nb_tranche in range(self.grille.nb_tranches) :
            with open("poche_%d.pbm"%(nb_tranche+1), "w") as fichier :
                fichier.write("P1\n")
                fichier.write("# %d/%d\n"%(nb_tranche+1, self.grille.nb_tranches))
                fichier.write("%d %d\n"%(self.grille.nb_colonnes, self.grille.nb_lignes))
                for ligne in range(self.grille.nb_lignes) :
                    fichier.write(' '.join(['0' if case.poche else '1'
                            for case in self.grille.tranches[nb_tranche][ligne]])+'\n')

