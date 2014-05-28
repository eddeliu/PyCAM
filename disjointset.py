class DisjointSet(object) :
    def __init__(self) :
        self.values = set()
        self.representatives = {}
    def add(self, value) :
        self.values.add(value)
        self.representatives[value] = value
    def union(self, value, other_value) :
        # repr of 1st value become repr of 2nd value
        self.representatives[self.find(other_value)] = \
            self.representatives[self.find(value)]
    def find(self, value) :
        if value != self.representatives[value] :
            self.representatives[value] = self.find(self.representatives[value])
        return self.representatives[value]
    def select(self, value) :
        rep = self.find(value)
        return filter(lambda val : self.find(val) == rep, self.values)
    #return all representatives (one per distinct set)
    def get_representatives(self) :
        return filter(lambda val : self.find(val) == val, self.values)
    def card(self, value) :
        return len(self.select(value))
    def disp(self) :
        print('Disjoint set')
        print(' '.join('{0}'.format(self.representatives[value], 3) for value in self.values))
        print(' '.join('{0}'.format(value, 3) for value in self.values))

class Value(int) : pass

values = {}
ds = DisjointSet()
choice = None
while choice != "quit" :
	print ''
	ds.disp()
	print ''
	choice = raw_input("Actions : add find union select card disp quit > ")
	if choice == "add" :
		val = input("Value > ")
		if val not in values :
			value = Value(val)
			values[val] = value
			ds.add(value)
			print 'Added !'
		else :
			print 'Already added.'
	elif choice == "find" :
		val = input("Value > ")
		if val in values :
			print 'Representative :', ds.find(values[val])
		else :
			print 'Add first'
	elif choice == "union" :
		val1 = input("Value > ")
		if val1 in values :
			val2 = input("Value > ")
			if val2 in values :
				ds.union(values[val1], values[val2])
				print 'Union done !'
			else :
				print 'Add first'
		else :
			print 'Add first.'
	elif choice == "select" :
		val = input("Value > ")
		if val in values :
			print 'Set :', ds.select(values[val])
		else :
			print 'Add first'
	elif choice == "card" :
		val = input("Value > ")
		if val in values :
			print 'Card: ', ds.card(values[val])
		else :
			print 'Add first'
	elif choice == "disp" :
		pass
	elif choice == "quit" :
		pass
	else :
		print 'Commande inconnue !'

#one = Value(1)
#two = Value(2)
#three = Value(3)
#four = Value(4)
#five = Value(5)
#six = Value(6)
#seven = Value(7)
#eight = Value(8)

#ds = DisjointSet()
#for n in (one,three,seven,five,two,four) :
	#ds.add(n)
#ds.disp()

#for v,ov in ((two,three), (four,five), (three,five)) :
	#ds.union(v,ov)
	#ds.disp()

#assert ds.find(four) == two
#assert five._representative == four
#assert ds.find(five) == two
#assert five._representative == two

#assert sorted(ds.select(five)) == sorted([two,three,four,five])
#assert sorted(ds.select(one)) == sorted([one])

