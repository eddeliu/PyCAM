
class VariablePriorityQueue(object) :
	def __init__(self) :
		self.queue = []
		self.position = {}
	def make(self, l) :
		for elem in l :
			self.disp()
			self.add(elem)
		self.disp()
	def add(self, elem) :
		self.position[elem] = len(self.queue)
		self.queue.append(elem)
		self.up(elem)
	def get(self) :
		self.swap(0, len(self.queue)-1)
		elem = self.queue.pop()
		self.down(self.queue[0])
		return elem
	def update(self, elem, value) :
		if value > 0 : self.down(elem)
		elif value < 0 : self.up(elem)
	def up(self, elem) : # TODO FIXME: choose the highest of the two childrens !!!!
		index_elem = self.position[elem]
		while True :
			index_parent = ((index_elem-2 if not index_elem%2 else index_elem)//2) \
														if index_elem>0 else 0
			print '\tindex_parent =', index_parent
			parent = self.queue[index_parent]
			print '\tparent =', parent
			print '\telem =', elem
			if parent > elem :
				self.swap(index_parent, index_elem)
				index_elem = index_parent
			else : break
	def down(self, elem) :
		length = len(self.queue)-1
		index_elem = self.position[elem]
		while True :
			index_child1 = index_parent*2+1
			index_child2 = index_child1+1
			if index_child1 > length : # no child
				return
			elif index_child2 > length : # only left child
				index_child, child = index_child1, child1
			else : # both child
				index_child, child = (index_child1, child1) if child1 > child2 \
								else (index_child2, child2)
			if child < parent :
				self.swap(index_elem, index_child)
				index_elem = index_child
			else : break
	def swap(self, index1, index2) :
		self.queue[index1], self.queue[index2] = self.queue[index2], self.queue[index1]
		self.position[index1], self.position[index2] = index2, index1
	def disp(self) :
		print ' '.join(str(task) for task in self.queue)

class Task(object) :
	def __init__(self, priority) :
		self.priority = priority
	def __lt__(self, other) :
		return self.priority < other.priority
	def __str__(self) :
		return str(self.priority)


queue = VariablePriorityQueue()
#queue.make([5,9,10,5,6,7,5,9,6,4,1,3,4,8,6,4,6,45,14])# [1, 4, 3, 4, 5, 4, 5, 5, 6, 9, 6, 7, 10, 8, 6, 9, 6, 45, 14]
queue.make([8,12,9,7,22,3,26,14,11,15,22]) # [3, 7, 8, 11, 15, 9, 26, 14, 12, 22, 22]
queue.disp()

	
