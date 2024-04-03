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


import time
import numpy as np
from TSPClasses import *
from State import *
import heapq
import itertools


class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities().copy()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			cities_visited = []
			# get start city randomly
			start_index = np.random.randint(ncities)
			current_city = cities[start_index]
			cities_visited.append(start_index)

			#repeat until all cities have been visited
			while len(cities_visited) < ncities:

				#initialize min distance to infinity
				min_dist = np.inf

				for i in range(ncities):
					if i in cities_visited:
						continue
					else:
						dist = current_city.costTo(cities[i])
						if dist < min_dist:
							min_dist = dist
							min_index = i
				cities_visited.append(min_index)
				current_city = cities[min_index]

			route = []
			# Now build the route using order we visited the cities
			for i in range( len(cities_visited) ):
				route.append( cities[cities_visited[i]] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0 ):
		######################################################################
		#initializations
		current_bssf = self.greedy()['cost']
		 
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		pq = []
		######################################################################

		#make original matrix
		original_matrix = np.zeros((ncities,ncities))
		for i in range(ncities):
			for j in range(ncities):
				original_matrix[i][j] = cities[i].costTo(cities[j])

		######################################################################
		
		#reduce matrix
		reduced = self.reduce_matrix(original_matrix, [], [])
		reduced_matrix = reduced['matrix']
		reduced_cost = reduced['cost']

		######################################################################

		start_city = cities[0]
		current_path = [start_city]
		heapq.heappush(pq, State(reduced_cost, current_path, reduced_matrix, [], [], reduced_cost))
		count += 1
		pruned = 0
		solutions = 0
		max_queue_size = 0

		#construct matrixes for each other city
		while pq and time.time()-start_time < time_allowance:
			current_state : State = heapq.heappop(pq)
			current_cost = current_state.cost
			current_path = current_state.path
			reduced_matrix = current_state.matrix #this was problematic
			start_city = current_path[-1]
			for i in range(ncities):
				count += 1
				if cities[i] in current_path:
					continue
				total_cost = current_cost
				matrix = reduced_matrix.copy()
				#set row and column to infinity
				matrix[start_city._index] = np.inf
				current_state.rows.append(start_city._index)
				matrix[:,i] = np.inf
				current_state.cols.append(i)
				matrix[i][start_city._index] = np.inf

				#reduce_matrix
				reduced = self.reduce_matrix(matrix, current_state.rows, current_state.cols)
				new_reduced_matrix = reduced['matrix']
				total_cost += reduced['cost']
				total_cost += reduced_matrix[start_city._index][i]
				if total_cost > current_bssf or total_cost == np.inf:
					pruned += 1
					continue
				
				if len(cities) > 10:
					if (len(current_path) == 1) and total_cost > current_bssf * .75:
						pruned += 1
						continue
					elif (len(current_path) == 2) and total_cost > current_bssf * .78:
						pruned += 1
						continue
					elif (len(current_path) == 3) and total_cost > current_bssf * .85:
						pruned += 1
						continue
					elif (len(current_path) == 4) and total_cost > current_bssf * .87:
						pruned += 1
						continue
					elif (len(current_path) == 5) and total_cost > current_bssf * .93:
						pruned += 1
						continue
					elif (len(current_path) == 6) and total_cost > current_bssf * .95:
						pruned += 1
						continue
					elif (len(current_path) == 7) and total_cost > current_bssf * .98:
						pruned += 1
						continue
					elif (len(current_path) == 8) and total_cost > current_bssf * .99:
						pruned += 1
						continue
					elif (len(current_path) == 9) and total_cost > current_bssf * .995:
						pruned += 1
						continue
					elif (len(current_path) == 10) and total_cost > current_bssf * .9997:
						pruned += 1
						continue
				elif total_cost > current_bssf:
					pruned += 1
					continue

				if total_cost < current_bssf:
					new_path = current_path.copy()
					new_path.append(cities[i])
					if len(new_path) == ncities:
						bssf = TSPSolution(new_path)
						current_bssf = bssf.cost
						solutions += 1
						foundTour = True
						continue
					if (len(new_path) > ncities * .85):
						heapq.heappush(pq, State(total_cost * .75, new_path, new_reduced_matrix, current_state.rows, current_state.cols, total_cost))
					elif (len(new_path) > ncities * .65):
						heapq.heappush(pq, State(total_cost * .85, new_path, new_reduced_matrix, current_state.rows, current_state.cols, total_cost))
					elif (len(new_path) > ncities * .5):
						heapq.heappush(pq, State(total_cost * .90, new_path, new_reduced_matrix, current_state.rows, current_state.cols, total_cost))

					else:
						heapq.heappush(pq, State(total_cost, new_path, new_reduced_matrix, current_state.rows, current_state.cols, total_cost))
					count += 1
					if len(pq) > max_queue_size:
						max_queue_size = len(pq)


		######################################################################
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = solutions
		results['soln'] = bssf
		results['max'] = max_queue_size
		results['total'] = count
		results['pruned'] = pruned
		return results
	
	#Bug is in this function, you need to pass a list of rows and cols
	def reduce_matrix(self, matrix, row, col):
		cost = 0
		#reduce rows
		for i in range(len(matrix)):
			if i in row:
				continue
			min_val = np.min(matrix[i])
			if min_val == np.inf:
				return {'matrix':matrix, 'cost':np.inf}
			
			for j in range(len(matrix[i])):
				if matrix[i][j] != np.inf:
					matrix[i][j] -= min_val
			cost += min_val

		#reduce columns
		for i in range(len(matrix)):
			if i in col:
				continue
			min_val = np.min(matrix[:,i])
			if min_val == np.inf:
				return {'matrix':matrix, 'cost':np.inf}
			for j in range(len(matrix[i])):
				if matrix[j][i] != np.inf:
					matrix[j][i] -= min_val
			cost += min_val

		return {'matrix':matrix, 'cost':cost}



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

	def fancy( self,time_allowance=60.0 ):
		pass
