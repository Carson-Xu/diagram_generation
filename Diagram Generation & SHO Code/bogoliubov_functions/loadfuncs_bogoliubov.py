"""
-=-=-=-

Functions for calculating (ground-state and excited-state) eigenvalues of single oscillator systems.


Changelog:

Version 1.0
Notes:
- First working version.
- Updated comments.


-=-=-=-


No inputs used.


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
from scipy.special import factorial as fac
from configparser import ConfigParser

import importlib
import bogoliubov_functions.loadfuncs_general as genl_func



### Load in config.
config = ConfigParser()		
config.read("config.ini")
params = config["Parameters"]



angu_freq = float(params["angu_freq"])
starting_amplitude_t1 = float(params["starting_amplitude_t1"])
starting_amplitude_t2 = float(params["starting_amplitude_t2"])

omit_amplitudes = genl_func.tuple_from_string(params["omit_amplitudes"], int)

cc_conv = float(params["cc_conv"])
cc_denominator = genl_func.convert_bool(params["cc_denominator"])
dampen_oscill = genl_func.convert_bool(params["dampen_oscill"])
show_all_cc = genl_func.convert_bool(params["show_all_cc"])
max_loops_cc = int(params["max_loops_cc"])
max_loops_bogo = int(params["max_loops_bogo"])
bogo_conv = float(params["bogo_conv"])
print_debug = genl_func.convert_bool(params["print_debug"])


Vmax = int(params["Vmax"])
Tmax = int(params["Tmax"])


equation_path = f"bogoliubov_functions.diagram_equations.loadeqns_sho_V{Vmax}_T{Tmax}"
dgrm_eqns = importlib.import_module(equation_path)


zero_array = np.zeros(shape = Tmax)
zero_array[0], zero_array[1] = starting_amplitude_t1, starting_amplitude_t2

amplitude_use = [True for i in range(0, Tmax)]

for i in omit_amplitudes:
	amplitude_use[i-1] = False


#######

# In[25]:


### 1.0
#Calculate open diagrams in preparation for matrix element calculation.

invf2 = 1/fac(2)
invf3 = 1/fac(3)
invf4 = 1/fac(4)



#max_loops_conv#' (int) - Maximum number of iterations to run through before ending the cycle.
#dampen_oscill#' (bool) - Dampen oscillations.
#convergences#' (float list) - List of desired convergences for each [t1, t2, t3, t4].
#use_ampl#' (bool list) - Effectively turns on or off use of [t1, t2, t3, t4] while iterating.
#angu_freq#' - Oscillator frequency, ω.
#max_loops_bogo#

def file_break(file_path, filephase):

	if (filephase == "intialize"):
		with open(file_path, "w") as file:
			file.write("")

	elif (filephase == "energy_calc"):
		with open(file_path, "a") as file:
			file.write("=-= ################################### =-=\n")
			file.write("=-= \n")
			file.write("=-= Energy Calculation Input \n")
			file.write("=-= \n")
			file.write("=-= ################################### =-=\n")
			file.write("\n\n")

	elif (filephase == "config_params"):
		with open(file_path, "a") as file:
			file.write(f"Angular frequency: {angu_freq}\n\n")

			file.write(f"Amplitudes used: [{genl_func.join_strings(string_array = amplitude_use, stagger = 5, separator = ', ')}]\n")

			file.write(f"CC convergence criteria: {cc_conv}\n")

			file.write(f"Starting amplitudes: {starting_amplitude_t1}, {starting_amplitude_t2}\n")
			
			file.write(f"Bogo convergence criteria: {bogo_conv}\n")
			file.write(f"CC denominator diagrams: {cc_denominator}\n")
			file.write(f"Maximum CC iterations allowed: {int(max_loops_cc)}\n")
			file.write(f"Maximum Bogo iterations allowed: {int(max_loops_bogo)}\n")
			file.write(f"Dampened oscillations for convergence: {dampen_oscill}\n")

	elif (filephase == "bogo_transform"):
		with open(file_path, "a") as file:
			file.write("\n\n\n\n=-= ################################### =-=\n")
			file.write("=-= \n")
			file.write("=-= Bogoliubov Transformations \n")
			file.write("=-= Unperturbed: [Delta, Lambda, Gamma] \n")
			file.write("=-= Perturbed: [a, b, c, d] \n")
			if show_all_cc == True:
				file.write("\n=-= Individual CC Iterations:\n")
				file.write("=-= #N [t1(N), t2(N), t3(N), ...] \n")
				file.write("=-= dN [t1(N) - t1(N-1), t2(N) - t2(N-1), t3(N) - t3(N-1), ...] \n\n")
			file.write("=-= Amplitudes: [t1, t2, t3, ...] \n")
			file.write("=-= Bogo. Coeff.: [F, G, D] \n")
			file.write("=-= GS Energy: Using zero & converged amplitude calculations. \n")
			file.write("=-= \n")
			file.write("=-= ################################### =-=\n\n")

	elif (filephase == "final_output"):
		with open(file_path, "a") as file:
			file.write("\n\n=-= ################################### =-=\n")
			file.write("=-= \n")
			file.write("=-= Final Output \n")
			file.write("=-= \n")
			file.write("=-= ################################### =-=\n\n")


def normal_order_coeff(pert):
	"""
	Normal orders coefficients given a quartic perturbation of the following form.
	
	V = αx + βx^2 + γx^3 + δx^4
	
	The normal ordered potential can then be split as follows.
	
	V0 = constant
	V1 = (A* a†) + (A a)
	V2 = 1/2!(B* a†a†) + (J a†a) + 1/2!(B aa)
	V3 = 1/3!(C* a†a†a†) + 1/2!(K* a†a†a) + 1/2!(K a†aa) + 1/3!(C aaa)
	V4 = 1/4!(D* a†a†a†a†) + 1/3!(L* a†a†a†a) + 1/2! 1/2!(M a†a†aa) 
		 + 1/3!(L a†aaa) + 1/4!(D aaaa)
	
	Input is as follows.
	'pert' - An array of the original perturbation, [α, β, γ, δ].
	
	Returns a tuple of normal ordered coefficients in the following order.
	(V0, A, B, J, C, K, D, L, M)
	"""
	
	α, β, γ, δ = pert
	
	## Update normal-ordered coefficients from [α, β, γ, δ].
	V0 = β/2 + 3*δ/4
	
	#V1 = A* a† + A a
	A = α/np.sqrt(2) + 3*γ/(2*np.sqrt(2))
	
	#x^2 = 1/2(a†a† + 2a†a + aa + 1)
	B = fac(2)*(β/2 + 6*δ/4)
	J = 2*β/2 + 12*δ/4
	
	#V3 = 1/3!(C* a†a†a†) + 1/2!(K* a†a†a) + 1/2!(K a†aa) + 1/3!(C aaa)
	C = fac(3)*γ/(2*np.sqrt(2))
	K = fac(2)*3*γ/(2*np.sqrt(2))

	#x^4 = 1/4(a†a†a†a† + 4a†a†a†a + 6a†a†aa + 4a†aaa + aaaa + 6a†a† + 12a†a + 6aa + 3)
	D = fac(4)*δ/4
	L = fac(3)*4*δ/4
	M = fac(2)*fac(2)*6*δ/4

	normal_ordered_coeffs = np.array([[V0, A,  B,  C,  D],
									 [A,  J,  K,  L,  0],
									 [B,  K,  M,  0,  0],
									 [C,  L,  0,  0,  0],
									 [D,  0,  0,  0,  0]])
	
	return normal_ordered_coeffs



def open_diagrams(unpert, pert, amplitudes):
	"""
	Calculate open diagrams in preparation for matrix element calculation.
	
	Inputs are as follows.
	'unpert' - An array of the unperturbed coefficients, [Δ, Λ, Γ].
	'pert' - An array of the original perturbation, [α, β, γ, δ].
	'amplitudes' - An array of the amplitudes, [t1, t2, t3, t4].
	
	Starting amplitudes are typically all set to zero, though the process 
	can be done far more quickly if the amplitudes are well-known for an 
	approximately similar perturbation.
	
	References the following functions.
	> normal_order_coeff
	
	Returns totaled diagrams of the form (W (a†)^k (a)^l) for later use, 
	represented as a list of tuples (W, k, l).
	"""

	## Prepare inputs to use in calculation.
	contrib = dgrm_eqns.get_diagrams(unpert = unpert, 
								normal_ordered_coefficients = normal_order_coeff(pert), 
								amplitudes = amplitudes, 
								angu_freq = angu_freq, 
								mode = "open")
	
	## ...where all open diagrams are sorted in one array.
	diagram_list = [tuple([contrib[k], int(k.split(',')[0]), int(k.split(',')[1])]) for k in contrib]
	
	if print_debug == True:
		b = [tuple([k[0], type(k[0]), k[1], k[2]]) for k in diagram_list]
		print('=-=')
		for i in b:
			print(i)

	#print(diagram_list)
	return diagram_list


#Find matrix elements <m|W|n> in preparation to build matrix.
def element(m, diagram_list, n):
	"""
	Calculate matrix element for <m|H|n> from a list of diagrams in preparation 
	to build the larger Hamiltonian matrix.
	
	Inputs are as follows.
	'm' - Row index of the desired Hamiltonian element.
	'diagram_list' - List of all diagrams to consider.
	'n' - Column index of the desired Hamiltonian element.
	
	Returns a total of all non-zero contributions to the given element (m,n).
	"""

	#tot_array = []
	tot_result = 0.
	result = 0.
	
	#For each total entry in the diagram_list (of the form W (a†)^k (a)^l),...
	for diagram in diagram_list:
		
		W = diagram[0]
		k = diagram[1]
		l = diagram[2]
		
		#These are all open diagrams, so <0|W|0> = 0.
		if m == 0 and n == 0:
			result = 0.
			
		elif m - k == n - l:
			if l <= n:
				result = W * np.sqrt(fac(m) * fac(n))/(fac(n-l))
			if l > n:
				result = 0.
			
		elif m - k != n - l:
			result = 0.
		
		if np.isfinite(result) == False:
			if print_debug == True:
				print('=-=')
				print(f"Invalid value ({result}) encountered for:")
				print(f"k: {k}\nl: {l}\nm: {m}\nn: {n}\nW: {W}")
			result = 0.
		
		#... adds result to find the total matrix element resulting from all possible contributions.
		tot_result += result

		if (print_debug == True) and (result != 0.):	
			print('=-=')
			print("New contribution for:")
			print(f"k: {k}\nl: {l}\nm: {m}\nn: {n}\nW: {W}")
			print(f"result: {result}")
			print(f"tot_result: {tot_result}")
		#tot_array = np.append(tot_array, result)
		#print(k, l, tot_result, result)
	
	#print("-----")
	return tot_result


#Calculate perturbed excited state energies.
def excited_state_energies(file_path, write_output, diagram_list, dim, pert_ground_E):
	"""
	Uses list of diagrams to find excited state energies by diagonalization 
	of the effective Hamiltonian.
	
	Inputs are as follows.
	'diagram_list' - List of all diagrams to consider.
	'dim' - Dimensionality or size of the matrix to use. 
			(i.e. the 'N' of an N x N matrix.)
	'pert_ground_E' - Perturbed ground state energy.
	
	References the following functions.
	> element
	
	Returns a tuple of length 'dim' with all excited state energies, sorted in
	order from smallest to largest.  Note that this is NOT NECESSARILY the order
	one would expect if sorted by the nth excited state, especially if values
	used here are far from convergence.
	"""
	
	#Use 'dim' to define size of the Hamiltonian matrix.
	mat_elements = np.zeros( shape=(dim, dim) )
	
	#Now, create H_effective matrix elements by feeding the diagram_list into the matrix element function.
	for i in range(0, dim):
		for j in range(0, dim):
			mat_elements[i][j] = element(i, diagram_list, j)
			#print(i, j, mat_elements[i][j])
	
	#if (print_debug == True):
	#	np.savetxt(f"testmatrix.txt", mat_elements)
	#	print(mat_elements)
	#Calculate the eigenvalues/eigenvectors
	eigvalues, eigvectors = np.linalg.eig(mat_elements)
	eigvalues = sorted(eigvalues)
	
	#Add the ground state energy to find the true excited state energy.
	pert_excit_E = eigvalues + pert_ground_E

	with open(file_path, "a") as file:
		
		file.write(f"Effective Hamiltonian Matrix:\n")
		file.write(f"{genl_func.np_strings(mat_elements.round(8))}\n")

		file.write(f"Eigenvalue Number:\t(Effective Hamiltonian Eig.)\tExcited Energy\n")
		for i in range(0, len(pert_excit_E)):
			effective_eig = eigvalues[i]

			if (effective_eig.imag == 0.0):
				effective_eig = effective_eig.real
				if effective_eig >= 0.0:
					parenth = f"(+{effective_eig})"
				if effective_eig < 0.0:
					parenth = f"({effective_eig})"

			if (effective_eig.imag != 0.0):
				parenth = f"{effective_eig}"

			file.write(f"Eig. #{i}:\t\t\t{parenth}\t\t\t{pert_excit_E[i]}\n")

		file.write(f"\nExcited Energies: {tuple(pert_excit_E)}\n")

	#print(mat_elements)
	return tuple(sorted(pert_excit_E))



def calculate_tn(contrib_eqns, tcalc, torder, use_torder, cc_denominator, coeffs):
	"""
	Calculate a given tn after an iteration of Equations of Motion Coupled Cluster equations.

	*-*-*-*-*
	
	Inputs are as follows.
	
	contrib_eqns:		dict		|	{"[n,m]": Q}
	- Represents all diagrammatic contributions.
	- Sorted into keywords [n,m] with values of coefficient Q, for any total contribution of the form Q (a†)^n (a)^m.
	- Only amplitude-related contributions are required here, thus the only used keywords are those where m = 0.

	tcalc:				float		|	tn
	- The amplitude value tn to use in calculations.

	torder:				int			|	n
	- The amplitude order n of the given tn.

	use_torder:			bool		|	Config parameter.
	- Boolean parameter to toggle the inclusion of a given cluster amplitude.

	cc_denominator:		bool		|	Config parameter.
	- Boolean parameter to toggle the inclusion of extra (J x tn) and (M x tn) diagrams in the denominator.

	coeffs:				3-tuple		|	(Δ, J, M)
	- Necessary coefficients from the unperturbed Hamiltonian and perturbation.
	
	*-*-*-*-*

	Returns the appropriate value of tn, given the associated parameters and inputs.
	"""

	tcalc_diagrams = contrib_eqns[f"{torder},0"]

	hamil_delta, J, M = coeffs

	invfJ = 1/fac(torder-1)

	if use_torder == False:
		tcalc_new = 0

	if use_torder == True:
		if cc_denominator == True:
			if torder <= 1:
				tcalc_new = - (tcalc_diagrams - ((J+hamil_delta)*invfJ*tcalc))/((J+hamil_delta)*invfJ)
			if torder > 1:
				invfM = 1/fac(torder-2)
				tcalc_new = - (tcalc_diagrams - ((J+hamil_delta)*invfJ*tcalc + M*invf2*invf2*invfM*tcalc))/((J+hamil_delta)*invfJ + M*invf2*invf2*invfM)

		if cc_denominator == False:
			tcalc_new = - (tcalc_diagrams - ((hamil_delta)*invfJ*tcalc))/((hamil_delta)*invfJ)

	return tcalc_new



def converge_amplitudes(file_path, write_output, unpert, pert, amplitudes):
	"""
	Iterate multiple times through each amplitude t1 through tn.
	
	*-*-*-*-*
	
	Inputs are as follows.
	
	unpert:			3-tuple		|	(Δ, Λ, Γ)
	- A tuple of the unperturbed coefficients, otherwise called (hamil_delta, hamil_lambda, hamil_gamma).
	
	pert:			4-tuple		|	(α, β, γ, δ)
	- A tuple of the original perturbation, with V = αx + βx^2 + γx^3 + δx^4.
	
	amplitudes:		tuple		|	(t1, t2, t3, ..., tn)
	- A tuple of all starting amplitudes.
	
	Starting amplitudes are typically all set to zero, though the process can be done far 
	more quickly if the amplitudes are well-known for an approximately similar perturbation.
	
	*-*-*-*-*

	Returns an array with all iterations of t1 through tn from start to finish.
	"""
	
	# Prepare inputs to use in calculation.
	hamil_delta = unpert[0]
	
	# Update normal-ordered coefficients from [α, β, γ, δ].
	normal_ordered_coeffs = normal_order_coeff(pert)
	
	A = normal_ordered_coeffs[1,0]

	B = normal_ordered_coeffs[2,0]
	J = normal_ordered_coeffs[1,1]

	C = normal_ordered_coeffs[3,0]
	K = normal_ordered_coeffs[2,1]

	D = normal_ordered_coeffs[4,0]
	L = normal_ordered_coeffs[3,1]
	M = normal_ordered_coeffs[2,2]

	# Generate the dictionary and iterable to keep up with different tn_tracker 'variable names'.
	amplitude_tracker = dict({})
	tn_order = len(amplitudes)
	all_tn_orders = range(1, tn_order+1)
	tn_calc_array = amplitudes

	for order in all_tn_orders:
		amplitude_tracker[f"t{order}_array"] = np.array(amplitudes[order-1])
		amplitude_tracker[f"t{order}_diff"] = cc_conv + 1
	
	# Looping for each coupled cluster calculation.
	for loop_num in range(0, max_loops_cc):

		# Toggle oscillator dampening for different behavior used.
		if dampen_oscill == True:
			if loop_num > 0:

				# If True, perform calculations with the average of the two previous results.
				tn_calc_array = []

				for order in all_tn_orders:
					tn_array = amplitude_tracker[f"t{order}_array"]
					toggle_dampened_tn = (tn_array[-1] + tn_array[-2])/2
					tn_calc_array.append(toggle_dampened_tn)

		if dampen_oscill == False:
			if loop_num > 0:
				
				# If False, perform calculations with the previous result.
				tn_calc_array = []
				
				for order in all_tn_orders:
					tn_array = amplitude_tracker[f"t{order}_array"]
					toggle_dampened_tn = tn_array[-1]
					tn_calc_array.append(toggle_dampened_tn)

		# Calculate diagrammatic contributions from previous iteration cycle.
		contrib = dgrm_eqns.get_diagrams(unpert = unpert, 
								normal_ordered_coefficients = normal_order_coeff(pert), 
								amplitudes = tn_calc_array, 
								angu_freq = angu_freq, 
								mode = "amplitude")


		# Calculate next iteration of all amplitudes here.
		tn_new_array = []
		for order in all_tn_orders:
			tn_new = calculate_tn(contrib_eqns = contrib, 
								tcalc = tn_calc_array[order-1], 
								use_torder = amplitude_use[order-1], 
								torder = order, 
								cc_denominator = cc_denominator, 
								coeffs = (hamil_delta, J, M))
			tn_new_array.append(tn_new)
		

		# Update all amplitude values in the tracker for the new cycle.
		for order in all_tn_orders:

			pull_tn_array = amplitude_tracker[f"t{order}_array"]

			tn_array = np.append(pull_tn_array, tn_new_array[order-1])
			tn_diff = np.abs(tn_array[-2] - tn_array[-1])

			amplitude_tracker[f"t{order}_array"] = tn_array
			amplitude_tracker[f"t{order}_diff"] = tn_diff


		# End iterations if any nan are found.
		if np.any(np.isnan(tn_new_array)):
			print("Found nan value while converging amplitudes!")
			break

		# End iterations if the amplitudes have all converged.
		if np.all([(amplitude_tracker[f"t{order}_diff"] < cc_conv) for order in all_tn_orders]):
			break

	# Final output and file writing.
	cc_amplitudes = np.array([amplitude_tracker[f"t{order}_array"] for order in all_tn_orders])

	if (write_output == True) and (show_all_cc == True):
		with open(file_path, "a") as file:
			file.write(f"\nIndividual CC Iterations:\n")
			file.write(f"\t#0: {cc_amplitudes[:,0]}\n")
			for j in range(1, len(cc_amplitudes[0])):
				file.write(f"\t#{j}: {cc_amplitudes[:,j]}\t\t")
				file.write(f"\td{j}: {cc_amplitudes[:,j] - cc_amplitudes[:,(j-1)]}\n")
			file.write(f"\n")
			
	return cc_amplitudes



#Create array for ground state energies.
def calc_gs_energy(file_path, write_output, unpert, pert, amplitudes):
	"""
	Calculate the ground state energy.
	
	Inputs are as follows.
	'unpert' - An array of the unperturbed coefficients, [Δ, Λ, Γ].
	'pert' - An array of the original perturbation, [α, β, γ, δ].
	'amplitudes' - An array of the converged amplitudes, [t1, t2, t3, t4].
	
	References the following functions.
	> normal_order_coeff
	
	Returns ground state energy from the given coefficients.
	"""
	
	_, hamil_lambda, hamil_gamma = unpert
	
	normal_ordered_coeffs = normal_order_coeff(pert)

	V0 = normal_ordered_coeffs[0,0]
	A = normal_ordered_coeffs[1,0]
	B = normal_ordered_coeffs[2,0]
	C = normal_ordered_coeffs[3,0]
	D = normal_ordered_coeffs[4,0]
	
	#Also have (Λa†a† + Λaa), rename to (1/2Ra†a† + 1/2Raa))
	R = 2*hamil_lambda

	#Perturbed G.S. Energy:
	contrib = dgrm_eqns.get_diagrams(unpert = unpert, 
								normal_ordered_coefficients = normal_order_coeff(pert), 
								amplitudes = amplitudes, 
								angu_freq = angu_freq, 
								mode = "closed")

	new_gs_energy = contrib["0,0"]
	
	if print_debug == True:	
		print(f"=-=")
		print(f"unpert: {unpert}")
		print(f"pert: {pert}")
		print(f"amplitudes: {amplitudes}")
		print(f"v")
		print(f"new_gs_energy: {new_gs_energy}\n")
		print(f"=-=")

	if (write_output == True):

		# Defined in case values do not exceed t4.
		tn_placeholder = (*amplitudes, 0, 0, 0, 0)
		t1, t2, t3, t4 = tn_placeholder[:4]

		with open(file_path, "a") as file:
			file.write(f"Energy Contributions:\n\n")
			file.write(f"Unpert. \t(1/2)*angu_freq: \t{(1/2)*angu_freq}\n")
			file.write(f"Unpert. \tGamma: \t{hamil_gamma}\n")
			file.write(f"Pert. \tV0: \t{V0}\n")

			file.write(f"T1   \t[A*1*t1]: \t{A*1*t1}\n")
			file.write(f"T1^2 \t[(B+R)*invf2*1*t1*t1]: \t{(B+R)*invf2*1*t1*t1}\n")
			file.write(f"T1^3 \t[C*invf3*1*t1*t1*t1]: \t{C*invf3*1*t1*t1*t1}\n")
			file.write(f"T1^4 \t[D*invf4*1*t1*t1*t1*t1]: \t{D*invf4*1*t1*t1*t1*t1}\n")

			file.write(f"T1*T2 \t[C*invf2*1*t1*t2]: \t{C*invf2*1*t1*t2}\n")
			file.write(f"T1*T3 \t[D*invf3*1*t1*t3]: \t{D*invf3*1*t1*t3}\n")

			file.write(f"T1^2*T2 \t[D*invf2*invf2*1*t1*t1*t2]: \t{D*invf2*invf2*1*t1*t1*t2}\n")

			file.write(f"T2   \t[(B+R)*invf2*1*t2]: \t{(B+R)*invf2*1*t2}\n")
			file.write(f"T2^2 \t[D*invf2*invf2*invf2*1*t2*t2]: \t{D*invf2*invf2*invf2*1*t2*t2}\n")

			file.write(f"T3   \t[C*invf3*1*t3]: \t{C*invf3*1*t3}\n")
			file.write(f"T4   \t[D*invf4*1*t4]: \t{D*invf4*1*t4}\n\n")

			file.write(f"GS Energy (No Open Diagrams): \t\t{new_gs_energy}\n")

	return new_gs_energy
	#print(ground_energy_array)



def bogoliubov_single_run(file_path, write_output, unpert, pert, mode):
	"""
	Run a number of Bogoliubov transformations.
	
	Inputs are as follows.
	'unpert' - An array of the (inital) unperturbed coefficients, [Δ, Λ, Γ].
	'pert' - An array of the (inital) perturbation, [α, β, γ, δ].
	
	References the following functions.
	> converge_amplitudes
	> calc_gs_energy
	
	Returns ground state energy from the given coefficients.
	"""
	
	hamil_delta, hamil_lambda, hamil_gamma = unpert
	α, β, γ, δ = pert

	if mode != "start":
	
		all_ampl = converge_amplitudes(file_path, True, unpert, pert, zero_array)
		next_amplitudes = all_ampl[:,-1]

		gs_energy_final = calc_gs_energy(file_path, False, unpert, pert, next_amplitudes)

		# Defined in case values do not exceed t2.
		tn_placeholder = (*next_amplitudes, 0, 0)
		t1, t2 = tn_placeholder[:2]

		ampl_next = next_amplitudes

	if mode == "start":

		t1, t2 = starting_amplitude_t1, starting_amplitude_t2

		ampl_next = None
		gs_energy_final = None

	## Update Bogoliubov coefficients from [t1, t2, t3, t4].
	F_b = (1 - t2**2)**(-1/2)
	D_b = -F_b*t1
	G_b = -F_b*t2

	#bogo_coeff = np.array([F_b, G_b, D_b])

	μ = F_b - G_b
	ρ = 2*D_b*(G_b - F_b)
	P = ρ/np.sqrt(2)

	hamil_delta_new = (F_b**2 + G_b**2)*hamil_delta - 4*F_b*G_b*hamil_lambda
	hamil_lambda_new = (F_b**2 + G_b**2)*hamil_lambda - F_b*G_b*hamil_delta
	hamil_gamma_new = hamil_gamma + (G_b**2)*hamil_delta - 2*F_b*G_b*hamil_lambda + ( α*P + (2*β + hamil_delta + 2*hamil_lambda)*(1/2)*(P**2) + γ*(P**3) + δ*(P**4) )

	α_new = (μ)*(α + (2*β + hamil_delta + 2*hamil_lambda)*P + 3*γ*(P**2) + 4*δ*(P**3))
	β_new = (μ**2)*(β + 3*γ*P + 6*δ*(P**2))
	γ_new = (μ**3)*(γ + 4*δ*P)
	δ_new = (μ**4)*δ

	#V0, A, B, J, C, K, D, L, M = normal_order_coeff(pert)
	
	unpert_next = [hamil_delta_new, hamil_lambda_new, hamil_gamma_new]
	pert_next = [α_new, β_new, γ_new, δ_new]
	bogo_coeff = [F_b, G_b, D_b]
	
	return [unpert_next, pert_next, ampl_next, bogo_coeff, gs_energy_final]#, bogo_total, bogo_gs_total, cc_iter_ampl]



def bogoliubov_converge(file_path, write_output, unpert, pert, dim):
	"""
	Run a number of Bogoliubov transformations.
	
	Inputs are as follows.
	'unpert' - An array of the (inital) unperturbed coefficients, [Δ, Λ, Γ].
	'pert' - An array of the (inital) perturbation, [α, β, γ, δ].
	
	References the following functions.
	> converge_amplitudes
	> calc_gs_energy
	
	Returns ground state energy from the given coefficients.
	"""
	
	#unpert_total = np.array([unpert])
	#pert_total = np.array([pert])
	#ampl_total = np.array([zero_array])
	#bogo_total = np.array([[0, 0, 0]])
	bogo_gs_delay = np.array([0, 2*bogo_conv])
	
	for N in range(0, max_loops_bogo+1):

		#print(unpert_total[-1])
		hamil_delta, hamil_lambda, hamil_gamma = unpert
		α, β, γ, δ = pert

		if (write_output == True):
			with open(file_path, "a") as file:
				file.write("\n=-= ################################### =-=\n") 
				file.write(f"=-= Bogo Iteration #{N} =-=\n")
				file.write("=-= ################################### =-=\n\n") 
				file.write(f"Unperturbed: [\t{hamil_delta},\t{hamil_lambda},\t{hamil_gamma}\t]\n")
				file.write(f"Perturbed:   [\t{α},\t{β},\t{γ},\t{δ}\t]\n")
		
		unpert_next, pert_next, ampl_next, bogo_coeff, gs_energy_final = bogoliubov_single_run(file_path, True, [hamil_delta, hamil_lambda, hamil_gamma], [α, β, γ, δ], mode = "standard")

		t1, t2 = ampl_next[0], ampl_next[1]

		F_b, G_b, D_b = bogo_coeff

		bogo_gs_delay[0] = bogo_gs_delay[1]
		bogo_gs_delay[1] = gs_energy_final

		#print(bogo_gs_delay)

		with open(file_path, "a") as file:
			file.write(f"Amplitudes: [\t{t1},\t{t2}\t]\n")
			file.write(f"Bogo Coeff: [\t{F_b},\t{G_b},\t{D_b}\t]\n\n")

		calc_gs_energy(file_path, True, unpert, pert, ampl_next)

		with open(file_path, "a") as file:
			file.write("\n\n=-= Zero Amplitude Calculation: =-=\n\n") 
	
		grnd_energy_zero_ampl = calc_gs_energy(file_path, False, unpert, pert, zero_array)
		exct_energy_zero_ampl = excited_state_energies(file_path, True, open_diagrams(unpert, pert, zero_array), dim, grnd_energy_zero_ampl)

		with open(file_path, "a") as file:
			file.write("\n\n=-= Converged Amplitude Calculation: =-=\n\n") 

		grnd_energy_conv_ampl = calc_gs_energy(file_path, False, unpert, pert, ampl_next)
		exct_energy_conv_ampl = excited_state_energies(file_path, True, open_diagrams(unpert, pert, ampl_next), dim, grnd_energy_conv_ampl)

		unpert = unpert_next
		pert = pert_next
		
		if (N >= 1) and (np.abs(bogo_gs_delay[1] - bogo_gs_delay[0]) <= bogo_conv):
			break
		
		with open(file_path, "a") as file:
			file.write("\n\n\n\n") 

	#return [unpert_next, pert_next, ampl_next]#, bogo_total, bogo_gs_total, cc_iter_ampl]
	return grnd_energy_zero_ampl, exct_energy_zero_ampl, grnd_energy_conv_ampl, exct_energy_conv_ampl



"""
-=-=-=-

End of function file.

-=-=-=-
"""
