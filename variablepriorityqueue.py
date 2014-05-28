
class VariablePriorityQueue(object) :
	def __init__(self) :
		self.queue = []
		self.position = {}
		self.n_export = 0
	def make(self, l) :
		for elem in l :
			self.add(elem)
	def add(self, elem) :
        position = len(self.queue)
		self.position[elem] = position
		self.queue.append(elem)
		self.up(position)
	def get(self) :
		self.swap(0, len(self.queue)-1)
		elem = self.queue.pop()
		if len(self.queue) :
			self.down(0)
		return elem
	def update(self, elem) :
        position = self.position[elem]
		self.queue[position].priority += elem.priority
		if elem.priority > 0 : self.down(position)
		elif elem.priority < 0 : self.up(position)
	def up(self, index_elem) :
		while True :
			index_parent = ((index_elem-2 if not index_elem%2 \
											else index_elem)//2) \
							if index_elem>0 else 0
			parent = self.queue[index_parent]
			if parent > elem :
				self.swap(index_parent, index_elem)
				index_elem = index_parent
			else : break
	def down(self, index_elem) :
		length = len(self.queue)-1
		while True :
			index_child1 = index_elem*2+1
			index_child2 = index_child1+1
			if index_child1 > length : # no child
				return
			elif index_child2 > length : # only left child
				index_child, child = index_child1, child1
			else : # two child
				child1 = self.queue[index_child1]
				child2 = self.queue[index_child2]
				index_child, child = (index_child1, child1) \
								if child1 < child2 \
								else (index_child2, child2)
			if child < elem :
				self.swap(index_elem, index_child)
				index_elem = index_child
			else : break
	def swap(self, index1, index2) :
		self.queue[index1], self.queue[index2] = \
			self.queue[index2], self.queue[index1]
		self.position[self.queue[index1]], \
			self.position[self.queue[index2]] = index1, index2
	def disp(self) :
		print ' '.join(str(task.priority) for task in self.queue)
	def exportPNG(self) :
		from os import system
		with open('queue'+str(self.n_export)+'.dot', 'w') as dot_file :
			length = len(self.queue)
			dot_file.write('graph G {\n')
			for pos in range(length) :
				child = ''
				left = pos*2+1
				if left < length :
					child += str(self.queue[left]) # interface
				if left+1 < length :
					child += ' '+str(self.queue[left+1]) # interface
				dot_file.write('\t%s -- {%s};\n' % (self.queue[pos], child))
			dot_file.write('}\n')
		system("dot -Tpng queue"+str(self.n_export)+".dot -o queue"+str(self.n_export)+".png")
		self.n_export += 1

class Task(object) :
	num = 0
	def __init__(self, priority, num=None) :
		if num is None :
			self.num = Task.num
			Task.num += 1
		else : self.num = num
		self.priority = priority
	def __lt__(self, other) :
		return self.priority < other.priority
	def __repr__(self) :
		return 'Task('+str(self.priority)+','+str(self.num)+')'
	def __str__(self) :
		return 'Task_'+str(self.priority)+'_'+str(self.num)+''


queue = VariablePriorityQueue()
tasks = []
#                0 1 2  3 4 5 6 7 8 9 10  12  14  16   18
for priority in (5,9,10,5,6,7,5,9,6,4,1,3,4,8,6,4,6,45,14) :
	tasks.append(Task(priority))
queue.make(tasks)
queue.disp()
queue.exportPNG()
#raw_input("...")
queue.update(tasks[9], 10)
queue.disp()
queue.exportPNG()
queue.update(tasks[9], -10)
queue.disp()
queue.exportPNG()

	
