class DisjointSet(object) :
	def __init__(self) :
		self.values = set()
	def find(self, value) :
		if value != value._representative :
			value._representative = self.find(value._representative)
		return value._representative
	def union(self, value, other_value) : # repr of 1st value become repr of 2nd value
		self.find(other_value)._representative = self.find(value)._representative
	def add(self, value) :
		self.values.add(value)
		value._representative = value
	def disp(self) :
		print('Disjoint set')
		print(' '.join('{0}'.format(value._representative, 3) for value in self.values))
		print(' '.join('{0}'.format(value, 3) for value in self.values))

class Value(int) : pass

one = Value(1)
two = Value(2)
three = Value(3)
four = Value(4)
five = Value(5)
six = Value(6)
seven = Value(7)
eight = Value(8)

ds = DisjointSet()
for n in (one,three,seven,five,two,four) :
	ds.add(n)
ds.disp()

for v,ov in ((two,three), (four,five), (three,five)) :
	ds.union(v,ov)
	ds.disp()

assert ds.find(four) == two
assert five._representative == four
assert ds.find(five) == two
assert five._representative == two

