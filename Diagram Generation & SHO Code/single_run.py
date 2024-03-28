"""
-=-=-=-

User-input for calculating (ground-state and excited-state) eigenvalues of single oscillator systems.


Changelog:

Version 1.0
Notes:
- First working version.
- Updated comments.

Version 1.1
Notes:
- Fixed final output.
- Added starting amplitude input for t1 and t2.


-=-=-=-


Inputs are as follows.

a_lin:		float	|	a
- Linear perturbation coefficient.

b_qdr:		float	|	b
- Quadratic perturbation coefficient.

c_cub:		float	|	c
- Cubic perturbation coefficient.

d_qrt:		float	|	d
- Quartic perturbation coefficient.


angu_freq:		float	|	w
- Angular frequency of the unperturbed oscillator.  Default is 1.

dim_range:		list	|	[N1, N2, N3, ...]
- Values of matrix dimensions to iterate over, if making a number of files.


bogo_first:		bool	|	True or False
- If True, begins the iterative cycle by first running a Bogoliubov transform with the below starting amplitudes.  
- If False, skips this adjustment and simply begins the cycle with CC convergence.

starting_amplitude_t1:				float		|	t1
- Starting amplitude for t1, if bogo_first = True.

starting_amplitude_t2_range:		2-tuple		|	(t2_min, t2_max)
- Range of values for the t2 starting amplitude, used for the lower and upper limits of the x-axis.

steps_t2:		int
- Number of steps, or x-ticks, to use on the graph.


gs_window:		2-tuple		|	(gs_min, gs_max)
- Range of values for the ground state energy, used for the lower and upper limits of the y-axis.


config_params:		dict	|	{...}
- An extra dictionary of configurative flags to control how the code runs.

	Tmax: 		int 	|	n
	- The maximum order of amplitudes Tn expected.
	- Ensure the top line in bogoliubov_functions/loadfuncs_bogoliubov.py reflects this!

	omit_amplitudes: 	list 	|	[i, j, k, ...]
	- A list of amplitudes ti, tj, tk, ... to ignore in calculations.

	cc_conv: 			float
	- Convergence criteria for the coupled cluster iterations.

	cc_denominator: 	bool 	|	True or False
	- If True, uses adapted equations with extra denominator diagrams for faster convergence.
	- If False, uses simple equations with no extra denominator diagrams.
	
	dampen_oscill: 		bool 	|	True or False
	- If True, reduces iterative oscillations by averaging out the previous two values.				t(n) from (t(n-1) + t(n-2))/2
	- If False, uses simple iterative behavior to calculate directly from the previous value.		t(n) from t(n-1)

	show_all_cc: 		bool 	|	True or False
	- If True, writes cluster amplitudes at all iterations in the final output file.

	max_loops_cc: 		int
	- Cuts coupled cluster process after a certain number of iterations if necessary.

	max_loops_bogo: 	int
	- Cuts Bogoliubov process after a certain number of transformations if necessary.

	bogo_conv: 			float
	- Convergence criteria for the Bogoliubov transformations.

	print_debug: 		bool 	|	True or False
	- If True, prints relevant values for each iteration.  Should only be used for debugging purposes.


-=-=-=-
"""


# Input Parameters.

# a_lin = 0.
# b_qdr = 0.125
# c_cub = 0.
# d_qrt = 0.

# a_lin = 0.
# b_qdr = 0.1
# c_cub = 0.
# d_qrt = 0.12

#a_lin = 0.
#b_qdr = 0.1
#c_cub = 0.
#d_qrt = 0.

a_lin = 0.1
b_qdr = 0.075
c_cub = -0.03
d_qrt = 0.09

angu_freq = 1
dim_range = [6, 8, 10, 30]

Tmax = 2
minimum_coefficient = None

bogo_first = True
starting_amplitude_t1 = 0
starting_amplitude_t2 = 0

gs_window = (0.58, 0.68)

config_params = dict({
	"omit_amplitudes": [],
	"cc_conv": 1e-10,
	"cc_denominator": True,
	"dampen_oscill": True,
	"show_all_cc": True,
	"max_loops_cc": 10000,
	"max_loops_bogo": 10,
	"bogo_conv": 1e-10,
	"print_debug": False
	})


"""
-=-=-=-

Inputs Above.

#
#
#

Code Below.

-=-=-=-
"""


### Import Packages.

import os
import numpy as np
import math as math
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import figure
from scipy.special import factorial as fac

import importlib
import bogoliubov_functions.loadfuncs_general as gnrl_func
import bogoliubov_functions.loadfuncs_generate_simplosc

from configparser import ConfigParser
config_object = ConfigParser()


Vmax = 4
config_params["Vmax"] = Vmax
config_params["Tmax"] = Tmax

config_params["angu_freq"] = angu_freq

run_subtitle = f"generatedt{Tmax}"
plot_title_suffix = f" | Maximum Order $t_{{{Tmax}}}$"

equation_path = f"bogoliubov_functions/diagram_equations/loadeqns_sho_V{Vmax}_T{Tmax}.py"

if os.path.exists(equation_path):

	print(f"Equation file loadeqns_sho_V{Vmax}_T{Tmax}.py successfully found!\n")

else:

	print(f"Equation file loadeqns_sho_V{Vmax}_T{Tmax}.py not found.  Generating equations...\n")

	bogoliubov_functions.loadfuncs_generate_simplosc.generate_equations(Vmax = Vmax, 
																	Tmax = Tmax, 
																	minimum_coefficient = minimum_coefficient)

	print(f"")



#Run iterations until values of T are within...
#convergence_bogo = convergences_cc[0]

for dim in dim_range:

	collect_amplit = []
	collect_energy = []

	config_params["starting_amplitude_t1"] = starting_amplitude_t1
	config_params["starting_amplitude_t2"] = starting_amplitude_t2

	config_object["Parameters"] = config_params

	#Write the above sections to config.ini file
	with open('config.ini', 'w') as conf:
		config_object.write(conf)

	#Forces function file to reload itself with the newly updated config.
	import bogoliubov_functions.loadfuncs_filewriting as main
	importlib.reload(main)
	
	main.write_full_file(perturbations = [a_lin, b_qdr, c_cub, d_qrt], 
							dim = dim, 
							angu_freq = angu_freq, 
							bogo_first = bogo_first, 
							starting_amplitudes = [starting_amplitude_t1, starting_amplitude_t2], 
							run_subtitle = run_subtitle)

	output_file_directory = f"Data Files/Saved Bogoliubov Data (Summary)/Bogoliubov Summary ({run_subtitle}_dim{dim})/"
	output_file_name = f"bogoliubov_data_summary (a={a_lin}, b={b_qdr}, c={c_cub}, d={d_qrt}, size={dim}).txt"
	file_path = f"{output_file_directory}{output_file_name}"

	#pull_amplit = pull_nth_amplitudes(file_path, 0)
	#pull_energy = pull_nth_gs_energy(file_path, 0)

	pull_t1_t2 = gnrl_func.pull_from_file(file_path=file_path, search_term="Amplitudes: ", instance=1, show_line=False)
	pull_energy_clos = gnrl_func.pull_from_file(file_path=file_path, search_term="GS Energy (No Open Diagrams): ", instance=0, show_line=False)
	pull_energy_zero = gnrl_func.pull_from_file(file_path=file_path, search_term="Eig. #0: ", instance=0, show_line=False)
	pull_energy_conv = gnrl_func.pull_from_file(file_path=file_path, search_term="Eig. #0: ", instance=1, show_line=False)

	pull_amplit = dict({"t1_amplitude": pull_t1_t2[0], 
						"t2_amplitude": pull_t1_t2[1]})

	pull_energy = dict({"closed_only": pull_energy_clos[0], 
						"zero_ampl": pull_energy_zero[2], 
						"conv_ampl": pull_energy_conv[2]})

	print(f"t1_amplitude:\t{pull_amplit['t1_amplitude']}\t\t|\tt2_amplitude:\t{pull_amplit['t2_amplitude']}")
	print(f"closed_only: \t{pull_energy['closed_only']}\t\t|\tzero_ampl:   \t{pull_energy['zero_ampl']}\t\t|\tconv_ampl:\t{pull_energy['conv_ampl']}")
	print("")


