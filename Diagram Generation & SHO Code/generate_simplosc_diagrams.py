"""
-=-=-=-

Generate single oscillator diagrams.


Changelog:

Version 1.0
Notes:

- Written using the diagrammatic class file at bogoliubov_functions/class_diagram.py.
- Code generates a python file with all single oscillator equations up to an arbitrary Vm and Tn for 
	later use in any coupled cluster calculations (i.e., writing summary files or graphing.)
- This file is located at bogoliubov_functions/diagram_equations/loadeqns_sho_V{m}_T{n}.py.


-=-=-=-


Inputs are as follows.

Vmax:		int		|	m
- Gives the maximum order perturbation to use in generating diagrams.		V = V1 + V2 + V3 + ... + Vm

Tmax:		int		|	n
- Gives the maximum order amplitude to use in generating diagrams.			T = T1 + T2 + T3 + ... + Tn
- Exponential cluster operator e^T is used in this code, to give all combinations of amplitudes below the given maximum.

minimum_coefficient:	float  or  None		|	q  or  None
- Removes diagrams if their associated factorial coefficient is below a certain value.
- To disable this feature, 'None' is considered a valid input.

	
-=-=-=-
"""

# Input Parameters.

Vmax = 4
Tmax = 20
minimum_coefficient = None

"""
-=-=-=-

Inputs Above.

#
#
#

Code Below.

-=-=-=-
"""

import bogoliubov_functions.loadfuncs_generate_simplosc

bogoliubov_functions.loadfuncs_generate_simplosc.generate_equations(Vmax = Vmax, 
																	Tmax = Tmax, 
																	minimum_coefficient = minimum_coefficient)


