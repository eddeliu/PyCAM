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

class DisjointSet(object) :
    """Hold values and permit to group them

    See Wikipedia for what it is, but basically it hold values.
    Values have each a representative, which are others' value.
    Values which have same representative belong to the same group.
    You can that way retrieve all values of a groups, just by
    knowing one value of this group."""

    def __init__(self) :
        """Just do basic initialisation"""
        self.values = set()
        self.representatives = {}

    def add(self, value) :
        """Add a value to the set, its representative is itself"""
        self.values.add(value)
        self.representatives[value] = value

    def union(self, value, other_value) :
        """Make two values have the same representative, grouping them"""
        # repr of 1st value become repr of 2nd value too
        self.representatives[self.find(other_value)] = \
            self.representatives[self.find(value)]

    def find(self, value) :
        """Get the representative of the group value belong to"""
        while self.representatives[value] != self.representatives[self.representatives[value]] :
            self.representatives[value] = self.representatives[self.representatives[value]]
        return self.representatives[value]
        # this is an other way to do it equally faster but less clearer :
        #this_repr = self.representatives[value]
        #repr_of_this_repr = self.representatives[this_repr]
        #while this_repr != repr_of_this_repr :
            #self.representatives[this_repr] = self.representatives[repr_of_this_repr]
            #this_repr = self.representatives[this_repr]
            #repr_of_this_repr = self.representatives[this_repr]
        #return this_repr

    def select(self, value) :
        """Get all the values of the group of given value"""
        rep = self.find(value)
        return filter(lambda val : self.find(val) == rep, self.values)

    def retrieve(self) :
        """Get all the groups' representatives"""
        return filter(lambda val : self.find(val) == val, self.values)

    def card(self, value) :
        """Get the size of the group value belong to"""
        return len(self.select(value))

    def disp(self) :
        print('Disjoint set')
        print(' '.join('{0}'.format(self.representatives[value], 3) for value in self.values))
        print(' '.join('{0}'.format(value, 3) for value in self.values))
