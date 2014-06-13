# -*- coding: utf-8 -*-
"""
This file is part of PyCAM.

PyCAM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PyCAM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PyCAM.  If not, see <http://www.gnu.org/licenses/>.
"""

class VariablePriorityQueueElement(object) :
    """Utility class for elements of VariablePriorityQueue

    It is aimed to be inherited, and defines minimal behaviour
    for the VariablePriorityQueue to works properly."""

    def __init__(self, priority) :
        if (not isinstance(priority, float)) or \
            (not isinstance(priority, int)) :
            raise TypeError("You tried to instanciate an VariablePriority" + \
                "QueueElement with wrong priority : " + str(type(priority)))
        else :
            self.priority = priority

    def __eq__(self, other) :
        raise NotImplementedError("You tried to compare two Variable" + \
            "PriorityQueueElement withoud defining such method : __eq__")

    def __hash__(self) :
        raise NotImplementedError("You tried to hash an VariablePriority" + \
            "QueueElement without defining such function : __hash__")


class VariablePriorityQueue(object) :
    """Implementation of a priority queue whose elements can change priority

    You can use it as a normal queue, but also change priority
    of elements within, and changes will be reflected back.
    But the elemnts have to be a bit tricky, because for sack of
    speed we use a dictionnary, which will use both __eq__(self, other)
    and __hash__(self) from the element class, so both must be
    correctly implemented. Be careful ! Else there will be really
    any warranty of the outcome.
    Internally, this is a minimal heap, and index of elements inside
    the heap (seen as a list) are held in the dictionnary to retrieve
    them really fast."""

    def __init__(self) :
        """Do nothing special"""
        self.queue = []
        self.position = {}

    def make(self, l) :
        """Add a list of elements to the queue"""
        for elem in l :
            self.add(elem)

    def add(self, elem) :
        """Add a single element to the queue"""
        position = len(self.queue)
        self.position[elem] = position
        self.queue.append(elem)
        self._up(position)

    def get(self) :
        """Return min priority element"""
        self._swap(0, len(self.queue)-1)
        elem = self.queue.pop()
        if len(self.queue) :
            self._down(0)
        del self.position[elem]
        return elem

    def update(self, elem) :
        """Update the priority of an element

        This is were __eq__ is specially needed because
        for the sack of speed we cannot pop the old element
        and insert the new one, so the parameter elem has
        to be some king of signal which holds the priority bonus/malus
        and will hash on the real element. This is why __hash__ and
        __eq__ have to be well defined for elements."""
        position = self.position[elem]
        priority_delta = elem.priority - self.queue[position].priority
        self.queue[position].priority = elem.priority
        if priority_delta > 0 : self._down(position)
        elif priority_delta < 0 : self._up(position)

    def _up(self, index_elem) :
        """Private method that move up an element in the heap"""
        elem = self.queue[index_elem]
        while True :
            index_parent = ((index_elem-2 if not index_elem%2 \
                                            else index_elem)//2) \
                            if index_elem>0 else 0
            parent = self.queue[index_parent]
            if parent > elem :
                self._swap(index_parent, index_elem)
                index_elem = index_parent
            else : break

    def _down(self, index_elem) :
        """Private method that move down an element in the heap"""
        elem = self.queue[index_elem]
        length = len(self.queue)-1
        while True :
            index_child1 = index_elem*2+1
            index_child2 = index_child1+1
            if index_child1 > length : # no child
                return
            elif index_child2 > length : # only left child
                index_child, child = index_child1, self.queue[index_child1]
            else : # two child
                child1 = self.queue[index_child1]
                child2 = self.queue[index_child2]
                index_child, child = (index_child1, child1) \
                                if child1 < child2 \
                                else (index_child2, child2)
            if child < elem :
                self._swap(index_elem, index_child)
                index_elem = index_child
            else : break

    def _swap(self, index1, index2) :
        """Private method which switch two elements inside the heap

        This does not garanty that the queue at the end is a heap
        because all it does is switching.
        It is to the responsability of calling methods to check
        for the queue remains a heap.
        This method also update the disctionnary, to not forgive it."""
        self.queue[index1], self.queue[index2] = \
            self.queue[index2], self.queue[index1]
        self.position[self.queue[index1]], \
            self.position[self.queue[index2]] = index1, index2

    def disp(self) :
        """Display what is in the queue, for debugging"""
        print '\n'.join(str(task) for task in self.queue)
        print self.position

    def print_heap(self) :
        """Create a DOT file that hold the queue

        The file is created, but it needs then to bo processed by :
        dot -Tpng queue.dot -o queue.png
        in order to see graphically what the heap looks like."""
        from os import system
        with open('queue.dot', 'w') as dot_file :
            length = len(self.queue)
            dot_file.write('graph G {\n')
            for pos in range(length) :
                child = ''
                left = pos*2+1
                if left < length :
                    child += str(self.queue[left])
                if left+1 < length :
                    child += ' '+str(self.queue[left+1])
                dot_file.write('\t%s -- {%s};\n' % (self.queue[pos], child))
            dot_file.write('}\n')
