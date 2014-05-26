import heapq

class Box(object) :
	def __init__(self, index) :
		self._distance = 0
		self.free_neighbours = []
		self.index = index
	def know_your_pocket_neighbourhood(self) :
		# base algorithm : Dijkstra
		# modifications :
		#  - iteration using flood
		#  - all distances are 1
		pop = heapq.heappop
		push = heapq.heappush
		m = Box.grid.matrix
		m[self.index][self.index] = (0, self)
		neighbours_left = [self]
		self._distance = 0
		while neighbours_left :
			nearest = pop(neighbours_left)
			new_distance = nearest._distance + 1
			for neighbour in nearest.free_neighbours :
				previous_distance = m[self.index][neighbour.index][0]
				if previous_distance is None :
					neighbour._distance = new_distance
					m[self.index][neighbour.index] = (new_distance, nearest)
					push(neighbours_left, neighbour)
				# else the path would be anyway equally long or longer
	def __str__(self) :
		return("Box%s(%s)" % (self.index, self._distance))
	def __repr__(self) :
		return("Box%s(%s)" % (self.index, self._distance))

class Grid(object) :
	def __init__(self) :
		self.matrix = list([[(None, None) for line in xrange(25)] for column in xrange(25)])
	def disp(self) :
		print '\n'
		print '*'*30
		print '\n'*3
		print '*'*30
		print '\n'
		for line in self.matrix :
			for column in line :
				if column != (None, None) :
					print column[0],# column[1],
				else : print '.',
			print ''



grid = Grid()
Box.grid = grid

box3 = Box(3)
box5 = Box(5)
box6 = Box(6)
box8 = Box(8)
box10 = Box(10)
box11 = Box(11)
box12 = Box(12)
box13 = Box(13)
box14 = Box(14)
box16 = Box(16)
box18 = Box(18)
box19 = Box(19)
box20 = Box(20)

box3.free_neighbours = [box6, box8]
box5.free_neighbours = [box6, box10]
box6.free_neighbours = [box3, box5, box10]
box8.free_neighbours = [box3, box12]
box10.free_neighbours = [box5, box6, box11, box14]
box11.free_neighbours = [box10, box12, box14, box16]
box12.free_neighbours = [box8, box11, box16]
box13.free_neighbours = [box14, box18]
box14.free_neighbours = [box10, box11, box13, box18]
box16.free_neighbours = [box11, box12, box20]
box18.free_neighbours = [box13, box14, box19]
box19.free_neighbours = [box18, box20]
box20.free_neighbours = [box16, box19]

boxes = [box3, box5, box6, box8, box10, box11, box12, box13, box14, box16, box18, box19, box20]

for box in boxes :
	#grid.disp()
	box.know_your_pocket_neighbourhood()
#grid.disp()

print 'Shortest path'
start_index = locals()['box'+raw_input('from : ')].index
end_index = locals()['box'+raw_input('to : ')].index
distance, step = grid.matrix[start_index][end_index]
path = [locals()['box'+str(end_index)]]
while step.index != start_index :
	path.append(step)
	step = grid.matrix[start_index][step.index][1]
path.append(step)
print 'Distance :', distance
print 'Path : ' + '-'.join(str(box.index) for box in path[::-1])

