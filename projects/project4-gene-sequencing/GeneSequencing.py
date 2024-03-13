#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

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

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		#score = random.random()*100
		#alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			#len(seq1), align_length, ',BANDED' if banded else '')
		#alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			#len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################
		
		if banded:
			alignment1 = "abc"
			alignment2 = "def"
			a1 = seq1[:align_length]
			a2 = seq2[:align_length]
			table = [[0 for i in range(len(a1) + 1)] for j in range(7)]
			for i in range(len(a1) + 1):
				table[0][i] = i * INDEL
			for i in range(1, 7):
				table[i][0] = i * INDEL
			

			for i in range(1, len(a1) + 1):
				for j in range(1, 7):
					if a1[i - 1] == a2[j - 1]:
						table[i][j] = min(table[i-1][j-1] + MATCH, table[i-1][j] + INDEL, table[i][j-1] + INDEL)
					else:
						table[i][j] = min(table[i-1][j-1] + SUB, table[i-1][j] + INDEL, table[i][j-1] + INDEL)
			score = table[len(a1)][6]


		else:
			a1 = seq1[:align_length]
			a2 = seq2[:align_length]

			#build table with infinity for all values
			table = [[0 for i in range(len(a2)+1)] for j in range(len(a1)+1)]

			#main algorithm
			for i in range(len(a1) + 1):
				table[i][0] = i * INDEL
			for j in range(len(a2) + 1):
				table[0][j] = j * INDEL
			for i in range(1, len(a1) + 1):
				for j in range(1, len(a2) + 1):
					if a1[i - 1] == a2[j - 1]:
						table[i][j] = min(table[i-1][j-1] + MATCH, table[i-1][j] + INDEL, table[i][j-1] + INDEL)
					else:
						table[i][j] = min(table[i-1][j-1] + SUB, table[i-1][j] + INDEL, table[i][j-1] + INDEL)
			score = table[len(a1)][len(a2)]

			#backtrack to find the alignment
			i = len(a1)
			j = len(a2)
			alignment1 = ""
			alignment2 = ""

			while i > 0 and j > 0:
				if a1[i-1] == a2[j-1]:
					alignment1 = a1[i-1] + alignment1
					alignment2 = a2[j-1] + alignment2
					i -= 1
					j -= 1
				else:
					if table[i][j] == table[i-1][j] + INDEL:
						alignment1 = a1[i-1] + alignment1
						alignment2 = "-" + alignment2
						i -= 1
					elif table[i][j] == table[i][j-1] + INDEL:
						alignment1 = "-" + alignment1
						alignment2 = a2[j-1] + alignment2
						j -= 1
					elif table[i][j] == table[i-1][j-1] + SUB:
						alignment1 = a1[i-1] + alignment1
						alignment2 = a2[j-1] + alignment2
						i -= 1
						j -= 1
			while i > 0:
				alignment1 = a1[i-1] + alignment1
				alignment2 = "-" + alignment2
				i -= 1
			while j > 0:
				alignment1 = "-" + alignment1
				alignment2 = a2[j-1] + alignment2
				j -= 1


			alignment1 = alignment1[:100]
			alignment2 = alignment2[:100]

		#return the results

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
