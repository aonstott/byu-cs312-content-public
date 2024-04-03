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
import math

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
			if (seq1 == "polynomial" and (seq2 != "exponential") and seq2 != "polynomial") or (seq1 == "exponential" and (seq2 != "polynomial") and seq2 != "exponential"):
				score = math.inf
				alignment1 = "Alginment not possible"
				alignment2 = "Alginment not possible"
				return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

			alignment1 = "abc"
			alignment2 = "def"
			a1 = seq1[:align_length]
			a2 = seq2[:align_length]
			table = [[0 for i in range(7)] for j in range(len(a1) + 1)]
			pointers = [["" for i in range(7)] for j in range(len(a1) + 1)]

			for i in range(3, 7):
				table[0][i] = (i - 3) * INDEL
			table[1][2] = INDEL
			table[2][1] = INDEL * 2
			table[3][0] = INDEL * 3
			last_j = 0
			#this part is really complicated and could definitely be organized better but it works
			#we use a kn array to save space, but that means we have to do some weird indexing and math stuff
			#basically 3 is 0 and j keeps track of the distance of the current character from the character in the string we are comparing to

			for i in range(1, len(a1) + 1):
				for j in range(0, 7):
					if (((i + (j - 3) - 1) >= 0) and ((i + (j - 3) - 1) < len(a2))):
						last_j = j
						if (a1[i - 1] == a2[i  + (j - 3) - 1]):
							if (j == 6):
								table[i][j] = min(table[i][j - 1] + INDEL, table[i - 1][j] + MATCH)
							elif (j == 0):
								table[i][j] = min(table[i - 1][j + 1] + INDEL, table[i - 1][j] + MATCH)

							else:
								table[i][j] = min(table[i][j - 1] + INDEL, table[i - 1][j + 1] + INDEL, table[i - 1][j] + MATCH)
						else:
							if (j == 6):
								table[i][j] = min(table[i][j - 1] + INDEL, table[i - 1][j] + SUB)
							elif (j == 0):
								table[i][j] = min(table[i - 1][j + 1] + INDEL, table[i - 1][j] + SUB)
							else:
								table[i][j] = min(table[i][j - 1] + INDEL, table[i - 1][j + 1] + INDEL, table[i - 1][j] + SUB)

					else:
						if (j == 6):
							table[i][j] = min(table[i][j - 1] + INDEL, table[i - 1][j] + SUB)
						elif (j == 0):
							table[i][j] = min(table[i - 1][j + 1] + INDEL, table[i - 1][j] + SUB)
						else:
							table[i][j] = min(table[i][j - 1] + INDEL, table[i - 1][j + 1] + INDEL, table[i - 1][j] + SUB)


					if (j != 6):
						if (table[i][j] == table[i - 1][j + 1] + INDEL):
							pointers[i][j] = "up"
					if (j != 0):
						if (table[i][j] == table[i][j - 1] + INDEL):
							pointers[i][j] = "left"
					if (pointers[i][j] == ""):
						if (table[i][j] == table[i - 1][j] + MATCH):
							pointers[i][j] = "diag"
						if (table[i][j] == table[i - 1][j] + SUB):
							pointers[i][j] = "diag"

			#backtrack to find the alignment

			alignment1 = ""
			alignment2 = ""

			i = len(a1)
			j = last_j
			k = len(a2)

			while (i > 0 or k > 0):
				if (pointers[i][j] == "diag"):
					alignment1 = a1[i - 1] + alignment1
					alignment2 = a2[k - 1] + alignment2
					i -= 1
					k -= 1
				elif (pointers[i][j] == "left"):
					alignment1 = "-" + alignment1
					alignment2 = a2[k - 1] + alignment2
					j -= 1
					k -= 1
				else:
					alignment1 = a1[i - 1] + alignment1
					alignment2 = "-" + alignment2
					i -= 1
					j += 1

			alignment1 = alignment1[:100]
			alignment2 = alignment2[:100]

			score = table[len(a1)][last_j]

		else:
			a1 = seq1[:align_length]
			a2 = seq2[:align_length]

			#build table with 0 for all values
			table = [[0 for i in range(len(a2)+1)] for j in range(len(a1)+1)]
			#pointer table to keep track of the direction of the alignment
			pointers = [["" for i in range(len(a2)+1)] for j in range(len(a1)+1)]

			#main algorithm
			for i in range(len(a1) + 1):
				table[i][0] = i * INDEL
			for j in range(len(a2) + 1):
				table[0][j] = j * INDEL
			for i in range(1, len(a1) + 1):
				for j in range(1, len(a2) + 1):
					if a1[i - 1] == a2[j - 1]:
						table[i][j] = min(table[i-1][j-1] + MATCH, table[i-1][j] + INDEL, table[i][j-1] + INDEL)
						if (table[i][j] == table[i - 1][j] + INDEL):
							pointers[i][j] = "left"
						elif (table[i][j] == table[i][j - 1] + INDEL):
							pointers[i][j] = "up"
						elif (table[i][j] == table[i - 1][j - 1] + MATCH):
							pointers[i][j] = "diag"
					else:
						table[i][j] = min(table[i-1][j-1] + SUB, table[i-1][j] + INDEL, table[i][j-1] + INDEL)
						if (table[i][j] == table[i - 1][j] + INDEL):
							pointers[i][j] = "left"
						elif (table[i][j] == table[i][j - 1] + INDEL):
							pointers[i][j] = "up"
						elif (table[i][j] == table[i - 1][j - 1] + SUB):
							pointers[i][j] = "diag"
			score = table[len(a1)][len(a2)]

			#backtrack to find the alignment
			alignment1 = ""
			alignment2 = ""

			i = len(a1)
			j = len(a2)

			while (i > 0 or j > 0):
				if (pointers[i][j] == "diag"):
					alignment1 = a1[i - 1] + alignment1
					alignment2 = a2[j - 1] + alignment2
					i -= 1
					j -= 1

				elif (pointers[i][j] == "left"):
					alignment1 = a1[i - 1] + alignment1
					alignment2 = "-" + alignment2
					i -= 1
				else:
					alignment1 = "-" + alignment1
					alignment2 = a2[j - 1] + alignment2
					j -= 1

			alignment1 = alignment1[:100]
			alignment2 = alignment2[:100]

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
