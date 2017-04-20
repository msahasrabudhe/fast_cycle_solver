#! /usr/bin/python

# Python implementation of a fast cycle solver [1]. 
#
#
# References:
# [1]: A Fast and Exact Energy Minimization Algorithm for Cycle MRFs, H. Wang, D. Koller, ICML 2013.

class Cycle:
	# A class for a cycle MRFs: MRFs defined on graphs that are simple cycles: 
	#   X1 --- X2 --- ... --- XN --- X1

	def __init__(self, n_nodes, n_labels):
		# We need at least three nodes in the graph for this algorithm. 
		if n_nodes < 3:
			print 'Cycle.__init__(): A cycle must have at least three nodes.'
			raise ValueError

		# The number of nodes in the graph. 
		self.n_nodes = n_nodes

		# The number of labels. If n_labels is an int, all nodes get the same
		#    number of labels. Else, each node i gets the number of labels indicated
		#    by n_labels[i]
		if type(n_labels) == np.int:
			self.n_labels = n_labels + np.zeros(self.n_nodes)
		elif n_labels.size == self.n_nodes:
			self.n_labels = n_labels

		# The maximum label a node can take.
		self.max_n_labels = np.max(self.n_labels)

		# Initialise messages.

		# Clique tree messages, \Delta^*_i in the paper. 
		self.Delta     = np.zeros((self.n_nodes - 1, self.n_labels[0]*self.max_n_labels))
		# Original energy terms. \theta_i in the paper. 
		self.theta     = np.zeros((self.n_nodes, self.max_n_labels*self.max_n_labels))
		# Indicator variables to denote whether a message is active or not. These variables
		#    respectively determine \Delta^+_i and \theta^+_i in the paper.
		self.D_active  = np.zeros((self.n_nodes - 1, self.n_labels[0]*self.max_n_labels), dtype=np.bool)  
		self.t_active  = np.zeros((self.n_nodes, self.n_labels[0]*self.max_n_labels), dtype=np.bool)

		# Messages on cycle edges
		self.delta     = np.zeros((self.n_nodes - 1, self.max_n_labels))

		# The optimal labelling.
		self.l_optimal = np.zeros(n_nodes)
		# The optimal energy. 
		self.Theta_hat = np.inf

	
	def compute_deltas(self):
		# Compute the messages \delta_i for i \in [1, self.n_nodes)
		
		# We first compute \delta_{self.n_nodes-1} (\delta_{N} in the paper).
		# However, we actually compute self.delta[self.n_nodes-2, 0:self.n_labels[self.n_nodes-1]] 
		#    as the \delta-s skip the first index. 
		n_lbl_N, n_lbl_1 = self.n_labels[self.n_nodes - 1], self.n_labels[0]

		_theta    = np.reshape(self.theta[self.n_nodes - 1, 0:n_lbl_N*n_lbl_1], [n_lbl_N, n_lbl_1])
		self.delta[self.n_nodes-2, 0:self.n_labels[self.n_nodes-1]] = np.min(self.theta[self.n_nodes - 1, 

	def minsum(self, i):
		# Calculates the minsum for \Delta^*_{i-1} and \theta_{i-1}. 
		# Returns the minimum as well as the index of the minimum. 
		
		# The input index, i, cannot be less than 1 as the messages \Delta^* start from
		#    index 1. 
		if i < 1:
			print 'Cycle.minsum(): Invalid index: %d. The input index must be greater than 0.' %(i)
			raise ValueError
		
		# We decrease i by 1, as the class does not store anything for \Delta^*_1. 
		i = i - 1
		# The arrays for which we perform a minsum. 
		_d       = self.Delta[i, :]
		_t       = self.theta[i, :]
		# The minsum is just the sum of the minimum values in these two arrays. 
		min_l_d  = np.argmin(_d)
		min_l_t  = np.argmin(_t)
		m_sum    = _d[min_l_d] + _t[min_l_t]

		return m_sum, min_l_d, min_l_t
		
	def partial_minsum(self, i):
		# Calculates the partial minsum for \Delta^+{i-1} and \theta^+_{i-1}
		# Returns the minimum as well as the index of the minimum. However, this
		#    is in agreement with self.Delta and self.theta rather than 
		#    self.Delta[:, _d_plus] and self.theta[:, _t_plus], where _d_plus 
		#    and _t_plus are indices where self.D_active and self.t_active are true
		#    for a particular index. That is to say, the indices returned 
		#    can be used directly in self.Delta[:,] and self.theta[:,].
		
		# The input index, i, cannot be less than 1 as the messages \Delta^* start
		#    from index 1. 
		if i < 1:
			print 'Cycle.partial_minsum(): Invalid index: %d. The input index must be greater than 0.' %(i)
			raise ValueError
	
		# We decrease i by 1, as the class does not store anything for \Delta^*_1.
		i = i - 1

		# The arrays for we perform minsum. 
		_d_plus  = np.where(self.D_active[i,:] == True)[0]
		_t_plus  = np.where(self.t_active[i,:] == True)[0]
		_d       = self.Delta[i, _d_plus]
		_t       = self.theta[i, _t_plus]

		min_l_d  = np.argmin(_d)
		min_l_t  = np.argmin(_t)
		m_sum    = _d[min_l_d] + _t[min_l_t]

		return m_sum, _d_plus[min_l_d], _t_plus[min_l_t]
