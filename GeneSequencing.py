#!/usr/bin/python3

from which_pyqt import PYQT_VER
from enum import Enum
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# Neighbor Types
UPPER = 0
LEFT = 1
DIAGONAL = 2

class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		# Trim the sequences to specified characters
		seq1 = seq1[0:align_length]
		seq2 = seq2[0:align_length]

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = self.alignHelper(seq1, seq2, banded)
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def alignHelper(self, seq1, seq2, banded):

		if seq1 == seq2:
			# The sequences are the same. This saves all the computation time along the diagonalal
			return -3 * len(seq1)

		# Add the blank character to the beginning
		seq1 = " " + seq1
		seq2 = " " + seq2

		# Else begin the algorithm for finding cost to align
		if banded:
			return self.alignBanded(seq1, seq2)

		return self.alignNoBand(seq1, seq2)


	def alignBanded(self, seq1, seq2):
		return -1

	def alignNoBand(self, seq1, seq2):

		""" Initialize the table for memoization """
		len_sequence1 = len(seq1)  # This will be the "top" word in the table
		len_sequence2 = len(seq2)  # This will be the "bottom" word in the table
		Table = [[{"cost": 0, "back_ptr": []} for i in range(len_sequence2)] for j in range(len_sequence1)]

		""" Initialize first row and first column of the table """

		# First row
		for i in range(1, len_sequence2):
			Table[0][i] = {
				"cost": 5 * i,  # Because the cost of an insert / delete is 5
				"back_ptr": []
			}

		# First column
		for j in range(1, len_sequence1):
			Table[j][0] = {
				"cost": 5 * j,  # Because the cost of an insert / delete is 5
				"back_ptr": []
			}

		"""Calculate all the costs and back pointers row-by-row"""

		UPPER = 0
		LEFT = 1
		DIAG = 2

		for i in range(1, len_sequence1):
			for j in range(1, len_sequence2):
				minCost, minCosts = self.minNeighborCosts(i, j, seq1, seq2, Table)  # Sets up an array of the costs for neighbors so min can be found
				Table[i][j]["cost"] = minCost

				# Add the minimum cost(s) as back pointers to the cell in the Table
				for k in range(len(minCosts)):
					if minCosts[k] == LEFT:
						Table[i][j]["back_ptr"].append("LEFT")
					elif minCosts[k] == UPPER:
						Table[i][j]["back_ptr"].append("UPPER")
					elif minCosts[k] == DIAG:
						Table[i][j]["back_ptr"].append("DIAG")

		return Table[len_sequence1-1][len_sequence2-1]["cost"]

	def minNeighborCosts(self, i, j, seq1, seq2, Table):
		""" Calculate the cost for left, upper, and diagonal neighbor """
		costs = [None for i in range(3)]  # There are 3 costs to compute

		costs[0] = Table[i-1][j]["cost"] + 5  # Cost of upper neighbor (cost of indel)  # TODO fixme
		costs[1] = Table[i][j-1]["cost"] + 5  # Cost of left neighbor (cost of indel)
		costs[2] = Table[i-1][j-1]["cost"] + self.diff(i, j, seq1, seq2)  # Cost of diag neighbor (-3 if char match, else 1)

		"""Find minimum cost"""
		minCost = float('inf')
		minCosts = []  # Holds the indexes for the minimum costs
						# Index 0 is left neighbor
						# Index 1 is upper neighbor
						# Index 2 is diag neighbor

		for i in range(3):
			if costs[i] < minCost:
				minCosts = []  # New minimum found, reset the array of pointers
				minCosts.append(i)  # Add the new minimum
				minCost = costs[i]
			elif costs[i] == minCost:
				minCosts.append(i)

		return costs[minCosts[0]], minCosts

	def diff(self, i, j, seq1, seq2):
		if seq2[j] == seq1[i]:
			return -3
		return 1