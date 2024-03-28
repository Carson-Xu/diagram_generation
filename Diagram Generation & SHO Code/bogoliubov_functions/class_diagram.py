"""
-=-=-=-

Functions to generate single oscillator diagrams.


Changelog:

Version 1.0
Notes:
- First working version.
- Updated comments.

Version 1.1
Notes:
- Changed format for final equations in code_(Vn,Tm).py files.
- Cancelled large redundant factorials.

	Original form:
	invf(#creation) * invf(#annihilation) * C * (#creation)!/(...!...!...!) * tn1 * tn2 * tn3

	Updated form:
					  invf(#annihilation) * C * 1/(...!...!...!) * tn1 * tn2 * tn3

- "symm_coeff" has been replaced with "reduced_symm_coeff" to account for this change.


Version 2.0
Notes:
- Complete freaking overhaul of the entire diagrammatic system.
- Cut 1h worth of processing time into roughly 10s by using braincells.

- Decided to use classes.  Should've done that originally.
- Fixed diagrammatic connections to work more efficiently.
- Fixed testing and filtering to work more efficiently.
- Pretty much everything's changed and is more efficient as a result.
- Why did I not optimize this code before having to write a paper on it!?

- Updated comments.

-=-=-=-
"""



### Import Packages.

import os
import numpy as np
from scipy.special import factorial as fac


class Diagram:
	"""
	This is a class for building and manipulating harmonic oscillator diagrams.

	Get Attributes
		> .get_perturbation()
		> .get_amplitudes()
		> .get_connections()
		> .get_remaining()
		> .get_components()
		> .get_unique()

	Calculate Values
		> .calc_symmetrizing()
		> .calc_factorials()
		> .calc_coefficient()
		> .calc_coefficient_for_comparison()
	
	Representations
		> .convert_dict()
		> .convert_str()
	
	Static Methods
		> connect(diagram, new_amplitudes, new_connections)
		> populate_connections(annihilation_count, amplitudes)
	"""

	def __init__(self, perturbation, amplitudes, connections):
		"""
		The constructor for the Diagram class.

		*-*-*-*-*

		Inputs are as follows.

		perturbation:	3-tuple		|	(Y, Z, "D")
		- Gives the number of creation and annihilation operators in the normal ordered perturbation.	D (a†)^Y (a)^Z
		- Takes in two ints and a string - the final value being the coefficient.

		amplitudes:		tuple		|	(i, j, k, ...)
		- Each index represents a cluster amplitude.		t_i, t_j, t_k, ...
		- Indices are forced to be in ascending order.		i < j < k < ...

		connections:	tuple		|	(p, q, r, ...)
		- Shows connections between the perturbation and the respective cluster amplitudes.				D(i:p, j:q, k:r, ...)
		- The total number of connections must not exceed the number of available operators.			p + q + r + ... ≤ Z
		- Individual connections must not exceed their respective cluster amplitude.					p ≤ i, q ≤ j, r ≤ k
		"""

		valid_diagram = True

		if np.any([(isinstance(k, int) == False) for k in (perturbation[0], perturbation[1], *amplitudes, *connections)]):
			valid_diagram = False
			raise ValueError("Non-integers found in input!")

		if np.any([(k < 0) for k in (perturbation[0], perturbation[1])]):
			valid_diagram = False
			raise ValueError("Negative exponents found in input!")

		if np.any([(k <= 0) for k in (*amplitudes, *connections)]):
			valid_diagram = False
			raise ValueError("Negative or zero-value indices found in input!")

		if len(amplitudes) > perturbation[1]:
			valid_diagram = False
			raise ValueError("Too many amplitudes given for the perturbation!")

		if sum(connections) > perturbation[1]:
			valid_diagram = False
			raise ValueError("Too many connections given for the perturbation!")

		if len(connections) != len(amplitudes):
			valid_diagram = False
			raise ValueError("Each amplitude requires a singular associated connection!")

		if np.any([(amplitudes[i] < connections[i]) for i in range(0, len(amplitudes))]):
			valid_diagram = False
			raise ValueError("All connections must be applicable!")

		if valid_diagram == True:
			self.perturbation = tuple(perturbation)
			self.amplitudes = tuple(amplitudes)
			self.connections = tuple(connections)

			self.remaining = (perturbation[0] + sum(amplitudes) - sum(connections), perturbation[1] - sum(connections))

			unsorted_components = list(zip(amplitudes, connections))
			components = sorted(unsorted_components, key=lambda ampl_conn_pair:(ampl_conn_pair[0], ampl_conn_pair[1]))

			self.components = tuple(components)

			collapsed_components = sorted(set(components), key=lambda ampl_conn_pair:(ampl_conn_pair[0], ampl_conn_pair[1]))
			remaining_amplitudes = [item[0] for item in collapsed_components]

			if len(amplitudes) == len(set(amplitudes)):
				self.uniqueness = "Pure"

			elif len(remaining_amplitudes) == len(set(remaining_amplitudes)):
				self.uniqueness = "Pseudo"

			else:
				self.uniqueness = "None"


	def get_perturbation(self):
		"""
		Obtains the perturbation term from the diagram.
		"""
		return self.perturbation

	def get_amplitudes(self):
		"""
		Obtains the cluster amplitudes from the diagram.
		"""
		return self.amplitudes

	def get_connections(self):
		"""
		Obtains the connections from the diagram.
		"""
		return self.connections

	def get_remaining(self):
		"""
		Obtains the numbers of remaining creation and annihilation operators in the diagram.
		"""
		return self.remaining

	def get_components(self):
		"""
		Obtains the uniquely sorted 'components' of the diagram.
		"""
		return self.components

	def get_uniqueness(self):
		"""
		Obtains the type of uniquess determined by the amplitudes and connections at different depths.
		Used to optimize the process of filtering out degenerate diagrams.
		Possible values include:

		'Pure':		t_i, t_j, t_k, ...		where	i ≠ j ≠ k ≠ ...
		- Returned when a diagram consists exclusively of unique orders of cluster amplitudes.
		- In other words, there are no t_i^2, t_i^3, t_i^4, and so on involved.
		- These terms are not checked for degeneracy, as there should be no confusion in the order of amplitudes.

		'Pseudo':	D(i:p, i:q, ..., j:r, j:s, ...)		where	i:p = i:q = ..., j:r, j:s = ..., etc.
		- Returned when a diagram may not consist of unique amplitudes, but every duplicate is connected equivalently.
		- For example, D(t3:2, t3:2, t4:1) consists of duplicate t3's, but they are both doubly connected.
		- These terms are not checked for degeneracy, as swapping around either of the t3:2 terms makes no change in the diagram.

		'None':		
		- The default result when none of the above requirements for uniqueness are met.
		- For example, D(t3:2, t3:2, t4:1, t4:2) may meet the conditions for t3, but t4 connections could be permuted to form a duplicate diagram.
		- These terms are the only ones checked for degeneracy, as an equivalent diagram with D(t3:2, t3:2, t4:2, t4:1) may exist.
		"""
		return self.uniqueness


	def calc_symmetrizing(self):

		if len(self.components) == 0:
			symmetrizing = None

		else:
			open_amplitudes = [(pair[0] - pair[1]) for pair in self.components]

			symmetrizing = (fac(self.remaining[0]))/(fac(self.perturbation[0]))
			for open_ampl_i in open_amplitudes:
				symmetrizing /= fac(open_ampl_i)

			if symmetrizing.is_integer():
				symmetrizing = int(symmetrizing)

		return symmetrizing

	def calc_factorials(self):

		duplicate_amplitudes = list(np.unique(np.array(self.amplitudes), return_counts=True)[1])
			
		factorials = tuple([i for i in (*self.connections, *duplicate_amplitudes) if i>1])

		return factorials

	def calc_coefficient(self):

		if len(self.components) == 0:
			coefficient = None

		else:
			open_amplitudes = [(pair[0] - pair[1]) for pair in self.components]

			coefficient = 1/(fac(self.perturbation[0]))
			for open_ampl_i in open_amplitudes:
				coefficient /= fac(open_ampl_i)

			if coefficient.is_integer():
				coefficient = int(coefficient)

		return coefficient

	def calc_coefficient_for_comparison(self):

		coefficient = Diagram.calc_coefficient(self)

		if coefficient == None:
			coefficient = 1
			
		return coefficient


	def convert_dict(self):

		display = dict({"perturbation": self.perturbation,
				"amplitudes": self.amplitudes,
				"connections": self.connections, 
				"remaining": self.remaining, 
				"factorials": Diagram.calc_factorials(self), 
				"symmetrizing": Diagram.calc_symmetrizing(self),
				"coefficient": Diagram.calc_coefficient(self),
				"uniqueness": self.uniqueness})

		return display

	def convert_str(self):
		# Factorials for remaining operators.
		diag_str = []

		# Only show for unconnected terms (i.e., directly from Hamiltonian).		
		if len(self.components) == 0 and self.remaining[0] > 1:
			diag_str.append(f"invf{self.remaining[0]}")

		#Any leftover annihilation operators.
		if self.remaining[1] > 1:
			diag_str.append(f"invf{self.remaining[1]}")
		
		# Desired perturbation coefficient.
		diag_str.append(self.perturbation[2])

		# Other associated factorials.
		diag_str.extend([f"invf{i}" for i in self.calc_factorials() if i>1])
			
		# Symmetrizing coefficient.
		if self.calc_coefficient() != None:
			diag_str.append(str(self.calc_coefficient()))
		
		# Cluster amplitudes.
		diag_str.extend([f"t{n}" for n in self.amplitudes])
		
		result_str = "*".join(diag_str)

		return result_str

	def __str__(self):
		return Diagram.convert_str(self)

	def equality_test(self, other):

		equals = np.all([self.perturbation == other.perturbation,
						self.components == other.components])

		return equals

	def __eq__(self, other):
		return Diagram.equality_test(self, other)


	@staticmethod
	def connect(diagram, new_amplitudes, new_connections):
		"""
		Connect further amplitudes to an existing diagram, resulting in a new diagrammatic object.

		*-*-*-*-*

		diagram:			Diagram	|	<...>
		- An existing diagram object, to build on top of.

		new_amplitudes:		tuple	|	(l, m, ...)
		- Amplitudes to connect to the existing diagram.
		
		new_connections:	tuple	|	(u, v, ...)
		- Connections to add, associated with the below amplitudes respectively.

		*-*-*-*-*
		
		Returns a new diagrammatic object with the given connections added on.
		"""
		new_diag = Diagram(perturbation =   diagram.perturbation,
							amplitudes  = (*diagram.amplitudes, *tuple(new_amplitudes)),
							connections = (*diagram.connections, *tuple(new_connections)))
		return new_diag

	@staticmethod
	def populate_connections(annihilation_count, amplitudes):
		"""
		Calculates all possible allowed connections for a set of amplitudes.

		*-*-*-*-*
		
		annihilation_count:		int		|	Z
		- Remaining annihilation operators in the diagram.

		amplitudes:				tuple	|	(i, j, k, ...)
		- Amplitudes to connect to the diagram.
		
		*-*-*-*-*

		Generates an full list of allowed connections by checking for the highest limit of each 
		amplitude, defined as the minimum between:

		- Ensuring generated connections are valid for the associated cluster operator.

			p ≤ i

		- Ensuring all other amplitudes connect at least once.

			p + q + r + ... ≤ Z
			p ≤ Z - (q + r + ...) = Z - (1 + 1 + ...)
		"""

		higher_limit = min(amplitudes[0], annihilation_count - (len(amplitudes) - 1))

		if len(amplitudes) > 1:
			for nth_connection in range(1, higher_limit+1):
				for next_connection in Diagram.populate_connections(annihilation_count - nth_connection, amplitudes[1:]):
					yield (nth_connection,) + next_connection

		if len(amplitudes) == 1:
			for nth_connection in range(1, higher_limit+1):
				yield (nth_connection,)

	@staticmethod
	def filter_degeneracy(diagram_group):
		"""
		The constructor for the Diagram class.

		*-*-*-*-*

		Inputs are as follows.

		perturbation:	3-tuple		|	(Y, Z, "D")
		- Gives the number of creation and annihilation operators in the normal ordered perturbation.	D (a†)^Y (a)^Z
		- Takes in two ints and a string - the final value being the coefficient.

		amplitudes:		tuple		|	(i, j, k, ...)
		- Each index represents a cluster amplitude.		t_i, t_j, t_k, ...
		- Indices are forced to be in ascending order.		i < j < k < ...

		connections:	tuple		|	(p, q, r, ...)
		- Shows connections between the perturbation and the respective cluster amplitudes.				D(i:p, j:q, k:r, ...)
		- The total number of connections must not exceed the number of available operators.			p + q + r + ... ≤ Z
		- Individual connections must not exceed their respective cluster amplitude.					p ≤ i, q ≤ j, r ≤ k
		"""

		unique_diagrams = []
		pseudo_diagrams = []

		degeneracy_save = []
		degeneracy_kill = []

		mono_degenerate = []

		for diagram in diagram_group:

			if diagram.get_uniqueness() == "Pure":
				unique_diagrams.append(diagram)

			elif diagram.get_uniqueness() == "Pseudo":
				pseudo_diagrams.append(diagram)

			elif diagram.get_uniqueness() == "None":

				max_ord_amplitude = max(diagram.get_amplitudes())
				max_ord_components = [ampl_conn_pair for ampl_conn_pair in diagram.get_components() if (ampl_conn_pair[0] == max_ord_amplitude)]

				if len(set(max_ord_components)) == 1:
					mono_degenerate.append(diagram)


				if diagram in degeneracy_save:
					degeneracy_kill.append(diagram)
					# print(f"Killed diagram with {diagram.get_perturbation()} and {diagram.get_components()}: {diagram.get_connections()}")

				if diagram not in degeneracy_save:
					degeneracy_save.append(diagram)
					# print(f"Saved diagram with {diagram.get_perturbation()} and {diagram.get_components()}: {diagram.get_connections()}")

		filtered_dict = dict({
							"unique_diagrams": unique_diagrams,
							"pseudo_diagrams": pseudo_diagrams,
							"degeneracy_save": degeneracy_save,
							"degeneracy_kill": degeneracy_kill,
							"mono_degenerate": mono_degenerate
							})

		return filtered_dict


"""
-=-=-=-

End of function file.

-=-=-=-
"""
