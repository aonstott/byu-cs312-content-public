class State:
    def __init__(self, pq_val, path, matrix, rows, cols, cost):
        self.pq_val = pq_val
        self.cost = cost
        self.path = path
        self.matrix = matrix
        self.rows = rows
        self.cols = cols

    def __lt__(self, other):
        return self.pq_val < other.pq_val