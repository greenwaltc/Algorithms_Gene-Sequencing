#!/usr/bin/python3

from which_pyqt import PYQT_VER
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
		score, seq1Alignment, seq2Alignment = self.alignHelper(seq1, seq2, banded)
		alignment1 = '{}'.format(seq1Alignment[:100])
		alignment2 = '{}'.format(seq2Alignment[:100])
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def alignHelper(self, seq1, seq2, banded):

		if seq1 == seq2:
			# The sequences are the same. This saves all the computation time along the diagonalal
			return -3 * len(seq1), "", ""

		# Add the blank character to the beginning
		seq1 = " " + seq1
		seq2 = " " + seq2

		# Else begin the algorithm for finding cost to align
		if banded:
			if abs(len(seq1) - len(seq2)) > 100:  # An arbitrary cutoff for no valid solution
				return math.inf, "No Alignment Possible", "No Alignment Possible"
			cost, seq1Alignment, seq2Alignment = self.alignBanded(seq1, seq2)
			return cost, seq1Alignment, seq2Alignment

		cost, seq1Alignment, seq2Alignment = self.alignMeaty(seq1, seq2)
		return cost, seq1Alignment, seq2Alignment


	def alignBanded(self, seq1, seq2):

		"""Initialize table for memoization"""
		len_sequence1 = len(seq1)  # This will be the "top" word in the table
		len_sequence2 = len(seq2)  # This will be the "bottom" word in the table
		# In this case, d (the max distance the letters can be apart) is always 3
		# As such, k = 2d + 1 = 7
		d = 3
		k = 7
		# Our offset is -d + i, where i is the row index in the table

		self.Table = [[{"cost": None, "back_ptr": None} for j in range(k)] for i in range(len_sequence1)]

		"""Initialize first row and column of the table"""
		# First row
		for k_j in range(k):
			offset = -d + 0
			if k_j + offset < 0:
				self.Table[0][k_j] = {
					"cost": None,
					"back_ptr": None
				}

			else:
				self.Table[0][k_j] = {
					"cost": 5 * (k_j + offset),
					"back_ptr": None if (k_j + offset) == 0 else "LEFT"
				}

		# First Column
		for k_i in range(1, d + 1):
			index = 3 - k_i
			self.Table[k_i][index] = {
				"cost": 5 * k_i,
				"back_ptr": "UPPER"
			}

		"""Calculate all min neighbor costs"""
		# For banded, left neighbor is (i, j - 1), upper is (i - 1, j + 1), diagonal is (i - 1, j)
		for i in range(1, len_sequence1):
			offset = -d + i
			for j in range(k):
				if (j + offset) <= 0 or (j + offset) > len_sequence2 - 1:
					continue
				else:
					self.minNeighborCosts_banded(i, j, offset, seq1, seq2)

		"""Get the minimum cost from the last row"""
		j = 6
		cost = float('inf')
		while self.Table[len_sequence1 - 1][j]["cost"] == None:
			j -= 1

		cost = self.Table[len_sequence1 - 1][j]["cost"]
		alignmentExits, seq1Alignment, seq2Alignment = self.makeAlignments_banded(seq1, seq2, i, j, d)
		return cost, seq1Alignment, seq2Alignment

	def makeAlignments_banded(self, seq1, seq2, i, j, d):
		seq1Alignment = seq1
		seq2Alignment = seq2
		alignmentExists = True

		while True:
			offset = -d + i
			if i > 0 and j > 2 and self.Table[i][j]["back_ptr"] is None:
				alignmentExists = False
				break
			elif i == 0 and j == 3:
				break
			else:
				if self.Table[i][j]["back_ptr"] == "DIAG":
					i -= 1
					continue
				elif self.Table[i][j]["back_ptr"] == "LEFT":
					seq1Alignment = seq1Alignment[:i+1] + "-" + seq1Alignment[i+1:]
					j -= 1
					continue
				elif self.Table[i][j]["back_ptr"] == "UPPER":
					seq2Alignment = seq2Alignment[:j+1+offset] + "-" + seq2Alignment[j+1+offset:]
					i -= 1
					j += 1
		return alignmentExists, seq1Alignment[1:], seq2Alignment[1:]


	def minNeighborCosts_banded(self, i, j, offset, seq1, seq2):
		cost = float('inf')
		back_ptr = None

		# Left
		if i >= 0 and (j - 1) < 7 \
			and self.Table[i][j - 1]["cost"] != None:
			cost = self.Table[i][j - 1]["cost"] + 5
			back_ptr = "LEFT"

		# Upper
		if (i - 1) >= 0 and (j + 1) < 7 \
				and self.Table[i - 1][j + 1]["cost"] != None \
				and self.Table[i - 1][j + 1]["cost"] + 5 < cost:
			cost = self.Table[i - 1][j + 1]["cost"] + 5
			back_ptr = "UPPER"

		# Diagonal
		if (i - 1) >= 0 and j < 7 \
			and self.Table[i - 1][j]["cost"] != None \
			and (self.Table[i - 1][j]["cost"] + self.diff_banded(i, j, offset, seq1, seq2)) < cost:
			cost = self.Table[i - 1][j]["cost"] + self.diff_banded(i, j, offset, seq1, seq2)
			back_ptr = "DIAG"

		self.Table[i][j]["cost"] = cost
		self.Table[i][j]["back_ptr"] = back_ptr

	def diff_banded(self, i, j, offset, seq1, seq2):
		if seq2[j + offset] == seq1[i]:
			return -3
		return 1

	def alignMeaty(self, seq1, seq2):

		""" Initialize the table for memoization """
		len_sequence1 = len(seq1)  # This will be the "top" word in the table
		len_sequence2 = len(seq2)  # This will be the "bottom" word in the table
		self.Table = [[{"cost": None, "back_ptr": None} for i in range(len_sequence2)] for j in range(len_sequence1)]

		""" Initialize first row and first column of the table """

		# First row
		for i in range(1, len_sequence2):
			self.Table[0][i] = {
				"cost": 5 * i,  # Because the cost of an insert / delete is 5
				"back_ptr": "LEFT"
			}

		# First column
		for j in range(1, len_sequence1):
			self.Table[j][0] = {
				"cost": 5 * j,  # Because the cost of an insert / delete is 5
				"back_ptr": "UPPER"
			}

		# First Cell
		self.Table[0][0] = {
			"cost": 0,  # Because the cost of an insert / delete is 5
			"back_ptr": None
		}

		"""Calculate all the costs and back pointers row-by-row"""
		for i in range(1, len_sequence1):
			for j in range(1, len_sequence2):
				minCost, backPtr = self.minNeighborCosts(i, j, seq1, seq2)  # Sets up an array of the costs for neighbors so min can be found
				self.Table[i][j]["cost"] = minCost
				self.Table[i][j]["back_ptr"] = backPtr

		"""Align the sequences"""
		alignmentExits, seq1Alignment, seq2Alignment = self.makeAlignments(seq1, seq2)

		return self.Table[len_sequence1-1][len_sequence2-1]["cost"], seq1Alignment[1:], seq2Alignment[1:]  # Remove leading space in alignments

	def minNeighborCosts(self, i, j, seq1, seq2):
		""" Calculate the cost for left, upper, and diagonal neighbor """

		# Since the tie breaker rule is always true, there only ever needs to be one back pointer: the minimum
		# of the first LEFT, UPPER, or DIAG
		cost = self.Table[i][j - 1]["cost"] + 5  # Cost of left neighbor (cost of indel)
		back_ptr = "LEFT"

		if self.Table[i-1][j]["cost"] + 5 < cost:
			cost = self.Table[i-1][j]["cost"] + 5
			back_ptr = "UPPER"

		if self.Table[i-1][j-1]["cost"] + self.diff(i, j, seq1, seq2) < cost:
			cost = self.Table[i-1][j-1]["cost"] + self.diff(i, j, seq1, seq2)
			back_ptr = "DIAG"

		return cost, back_ptr

	def diff(self, i, j, seq1, seq2):
		if seq2[j] == seq1[i]:
			return -3
		return 1

	def makeAlignments(self, seq1, seq2):
		i = len(seq1) - 1
		j = len(seq2) - 1
		seq1Alignment = seq1
		seq2Alignment = seq2
		alignmentExists = True

		while True:
			if i != 0 and j != 0 and self.Table[i][j]["back_ptr"] is None:
				alignmentExists = False
				break
			elif i == 0 and j == 0:
				break
			else:
				if self.Table[i][j]["back_ptr"] == "DIAG":
					i = i - 1
					j = j - 1
					continue
				elif self.Table[i][j]["back_ptr"] == "LEFT":
					seq1Alignment = seq1Alignment[:i+1] + "-" + seq1Alignment[i+1:]
					j = j - 1
					continue
				elif self.Table[i][j]["back_ptr"] == "UPPER":
					seq2Alignment = seq2Alignment[:j+1] + "-" + seq2Alignment[j+1:]
					i = i - 1
		return alignmentExists, seq1Alignment, seq2Alignment