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



### Import Packages.

import os
import numpy as np
import math as math
from matplotlib import pyplot as plt
from scipy.special import factorial as fac
import itertools as itt
import string
import json
import warnings
import time

import bogoliubov_functions.class_diagram as dgrm_cls
import bogoliubov_functions.loadfuncs_general as gnrl_func

global checkpoints

# Time-Keeping Stuff.

start_time = time.perf_counter()

checkpoints = [start_time]
checkpoint = start_time

def get_time_since_last():

	global checkpoints

	checkpoint = time.perf_counter()
	time_since_last = checkpoint - checkpoints[-1]
	checkpoints.append(time_since_last)

	return time_since_last



def get_time_total():

	current_time = time.perf_counter()
	time_total = current_time - start_time

	return time_total



def generate_equations(Vmax, Tmax, minimum_coefficient):

	# Create all starting perturbation indices, amplitude indices, and coefficients based on the given Vmax and Tmax.

	# Perturbations from the given V.
	V_contribs = [(k, m-k, f"Q_{k}_{m-k}") for m in range(0,Vmax+1) for k in range(0,m+1)]

	# Default terms from the generalized Hamiltonian.
	H_contribs = [(0, 0, "(1/2)*angu_freq"),
				(1, 1, "hamil_delta"),
				(0, 0, "hamil_gamma"),
				(2, 0, "R"),
				(0, 2, "R")]

	pert_hamil_contribs = (*V_contribs, *H_contribs)

	pert_hamil_coeffs = [["" for k in range(0, Vmax+1)] for l in range(0, Vmax+1)]

	for K in range(0, Vmax+1):
		for L in range(0, Vmax+1):

			relevant_contribs_K_L = tuple([i[2] for i in pert_hamil_contribs if (K == i[0] and L == i[1])])

			all_K_L = " + ".join(relevant_contribs_K_L)

			if len(relevant_contribs_K_L) > 1:
				all_K_L = f"({all_K_L})"

			pert_hamil_coeffs[K][L] = all_K_L

	pert_hamil_coeffs = np.array(pert_hamil_coeffs)

	perturbation_indices = [(creat_count, Vn - creat_count) for Vn in range(0, Vmax+1) for creat_count in range(0, Vn+1)]
	amplitude_indices = [tuple([tuple([Tn for Vi in range(0, Vn)]) for Vn in range(0, Vmax+1)]) for Tn in range(1, Tmax+1)]



	# Creates simple 'starting' diagram objects, based off the perturbation Vn.

	starting_diagrams = []

	for perturbation_term in perturbation_indices:

		creat_count = perturbation_term[0]
		annih_count = perturbation_term[1]

		next_diagram = dgrm_cls.Diagram(perturbation = (creat_count, annih_count, pert_hamil_coeffs[creat_count, annih_count]), 
									amplitudes = [], 
									connections = [])

		starting_diagrams.append(next_diagram)

	print(f"Built starting diagrams from the perturbed Hamiltonian.")



	# Connects diagrams iteratively.
	# Begins by connecting all available diagrams with powers of t1, then all resulting available diagrams with powers of t2, and so on.

	# degeneracy_save_count = 0
	degeneracy_kill_count = 0
	mono_degenerate_count = 0

	#number_generated_diagrams = []
	diagram_list = np.copy(starting_diagrams)

	for cluster_level in amplitude_indices:		# e.g. For every e^Tn.		|	e^T2 = ((), (2,), (2, 2), (2, 2, 2), (2, 2, 2, 2)).
		
		placehold_diagram_list = []

		for diagram in diagram_list:			# e.g. For every diagram. 	|	D(i:p, j:q, k:r, ...).

			annih_count = diagram.get_remaining()[1]

			applicable_amplitudes = [i for i in cluster_level if len(i) <= annih_count]

			for amplitude_term in applicable_amplitudes:  # e.g. For every individual Tn^m.		|	Tn^3 = (n, n, n)

				if len(amplitude_term) == 0:
					next_diagram = diagram
					placehold_diagram_list.append(next_diagram)

				else:
					possible_connections = dgrm_cls.Diagram.populate_connections(annihilation_count = annih_count, 
									 										amplitudes = amplitude_term)

					diagram_group = [dgrm_cls.Diagram.connect(diagram = diagram,
															new_amplitudes = amplitude_term, 
															new_connections = conn)
									for conn in possible_connections]

					filtered_dict = dgrm_cls.Diagram.filter_degeneracy(diagram_group)

					placehold_diagram_list.extend(filtered_dict["unique_diagrams"])
					placehold_diagram_list.extend(filtered_dict["pseudo_diagrams"])
					placehold_diagram_list.extend(filtered_dict["degeneracy_save"])

					# Filter out degenerate diagrams, but keep a count for the output.
					# degeneracy_save_count += len(filtered_dict["degeneracy_save"])
					degeneracy_kill_count += len(filtered_dict["degeneracy_kill"])
					mono_degenerate_count += len(filtered_dict["mono_degenerate"])

					# if len(filtered_dict["degeneracy_save"]) != 0 or len(filtered_dict["degeneracy_kill"]) != 0:
					# 	print("")



		print(f"Connected diagrams with e^(T{cluster_level[-1][0]}). ", end="")
		print(f"\t\tTotal time: {get_time_total():.3f}s", end="")
		print(f"\t\tSaved diagrams: {len(placehold_diagram_list)}", end="")
		print(f"\t\tRemoved: {degeneracy_kill_count}")
		#number_generated_diagrams.append(len(placehold_diagram_list))

		diagram_list = np.copy(placehold_diagram_list)
		#print(len(diagram_list))

	print("")



	# If minimum coefficient is defined, remove any diagrams which fail the check.

	if minimum_coefficient != None:
		print(f"Addressing minimum coefficient.  Input: {minimum_coefficient}")

		filtered_diagram_list = [dia for dia in diagram_list if dia.calc_coefficient_for_comparison() >= minimum_coefficient]

		print(f"Removed diagrams with coefficients lower than the provided value.")
		print(f"Cutting process complete.", end="")
		print(f"\t\t\t\tTotal time: {get_time_total():.3f}s", end="")
		print(f"\t\tTotal diagrams: {len(filtered_diagram_list)}", end="")
		print(f"\t\t(Removed diagrams: {len(diagram_list) - len(filtered_diagram_list)})\n")

		diagram_list = np.copy(filtered_diagram_list)



	# Rigidly filters for degenerate diagrams.
	# Accomplishes this efficiently by exclusively checking over diagrams that are eligible to be degenerate in the first place.

	#print("Addressing degeneracy.")

	unique_diagrams = []
	pseudo_diagrams = []
	degenerate_diagrams = []

	for diagram in diagram_list:
		if diagram.get_uniqueness() == "Pure":
			unique_diagrams.append(diagram)

		elif diagram.get_uniqueness() == "Pseudo":
			pseudo_diagrams.append(diagram)

		elif diagram.get_uniqueness() == "None":
			degenerate_diagrams.append(diagram)

	print(f"Found {len(unique_diagrams)} unique diagrams, {len(pseudo_diagrams)} pseudo-unique diagrams, and {len(degenerate_diagrams)} degenerate diagrams.")
	#print(f"Total: {len(diagram_list)}")

	# degenerate_save = []
	# degenerate_trash = []

	# for diagram in duplicate_diagrams:

	# 	if diagram in degenerate_save:
	# 		degenerate_trash.append(diagram)

	# 	if diagram not in degenerate_save:
	# 		degenerate_save.append(diagram)

	functionally_unique_diagrams = np.copy(diagram_list)

	# print(f"Of the duplicate diagrams, {len(degenerate_save)} instances were deemed to be unique.")
	# print(f"Filtering process complete.", end="")
	# print(f"\t\t\t\tTotal time: {get_time_total():.3f}s", end="")
	# print(f"\t\tTotal diagrams: {len(functionally_unique_diagrams)}", end="")
	
	print(f"Total removed diagrams: {degeneracy_kill_count}\n")

	number_separated_diagrams = tuple([Tmax, len(placehold_diagram_list), len(unique_diagrams), len(pseudo_diagrams), len(degenerate_diagrams), degeneracy_kill_count, mono_degenerate_count])
	print(f"Degeneracy breakdown for graphing: {number_separated_diagrams}")
	print(f"Max T, Total Generated, Unique, Semi-Unique, Kept Degenerate, Removed Degenerate, Mono Degenerate\n")
	# Sorts diagrams by order.

	def sortingfunc(diagram):

		return (sum(diagram.get_perturbation()[:2]),
		diagram.get_perturbation()[:0],
		diagram.get_perturbation()[:1],
		len(diagram.get_amplitudes()),
		(diagram.get_remaining())[0] - (diagram.get_remaining())[1], 
		(diagram.get_remaining())[0],
		(diagram.get_amplitudes()))

	functionally_unique_diagrams = sorted([dia for dia in functionally_unique_diagrams], key = sortingfunc)



	# Prepares output by sorting diagrams based on their contributions.
	# Categories include closed diagrams, amplitude-contributing diagrams (i.e., exclusively creation operators), and others.

	file_path = "bogoliubov_functions/diagram_equations/"
	gnrl_func.create_dir(file_path)
	file_name = file_path + f"loadeqns_sho_V{Vmax}_T{Tmax}.py"

	print("Beginning output.")
	print("Diagrams have been sorted into the following contributions:\n")

	contribution_set = set([dia.get_remaining() for dia in functionally_unique_diagrams])
	contributions = sorted(contribution_set, key = lambda term:(term[0]-term[1], term[1]))

	contrib_list = [(remaining_operators, [dia.convert_str() for dia in functionally_unique_diagrams if dia.get_remaining() == remaining_operators]) for remaining_operators in contributions]

	contrib_closed = [terms for terms in contrib_list if (terms[0] == (0,0))]
	contrib_amplitude = [terms for terms in contrib_list if (terms[0][0] != 0 and terms[0][0] <= Tmax and terms[0][1] == 0)]
	contrib_other = [terms for terms in contrib_list if terms not in [*contrib_closed, *contrib_amplitude]]



	# Depending on the value, defining a minimum_coefficient may cut out all contributions to a necessary Tn cluster equation.
	# This checks for any amplitudes for which we would expect a contribution, then adds a dummy contribution of 0 to ensure they get accounted for.

	existing_amplitudes = [terms[0][0] for terms in contrib_amplitude]
	expected_amplitudes = tuple(range(1, Tmax+1))

	missing_amplitudes = [Tn for Tn in expected_amplitudes if Tn not in existing_amplitudes]
	contrib_missing = [((Tn, 0), ["0"]) for Tn in missing_amplitudes]

	contrib_amplitude.extend(contrib_missing)

	#print(contrib_amplitude)



	# Print statements for output.

	print(f"GS Energy Equation (Closed):\t{sum([len(contrib_pair[1]) for contrib_pair in contrib_closed])}")
	print(f"Cluster Amplitude Equations:\t{sum([len(contrib_pair[1]) for contrib_pair in contrib_amplitude])}")
	for contrib_pair in contrib_amplitude:
		ampl_i = contrib_pair[0][0]
		if ampl_i in existing_amplitudes:
			print(f"\t- T{ampl_i} Equation:\t\t{len(contrib_pair[1])}")
		if ampl_i not in existing_amplitudes:
			print(f"\t- T{ampl_i} Equation:\t\t0")

	print(f"Found {sum([len(contrib_pair[1]) for contrib_pair in contrib_other])} other diagrams exclusively used in matrix elements.\n")



	# And finally, write the function file to use.

	with open(file_name, "w") as outfile:

		outfile.write(f"\"\"\"\n-=-=-=-\n\n")
		outfile.write(f"Equations for calculating groups of diagrammatic contributions for single oscillator systems.\n\n")
		outfile.write(f"Vmax = {Vmax}\n")
		outfile.write(f"Tmax = {Tmax}\n")
		if minimum_coefficient == None:
			outfile.write(f"Minimum Coefficient: N/A\n\n")
		else:
			outfile.write(f"Minimum Coefficient: {minimum_coefficient}\n\n")
		outfile.write(f"Generated {len(functionally_unique_diagrams)} total diagrams in {get_time_total():.3f}s.")
		outfile.write(f"\n\n-=-=-=-\n\"\"\"\n\n\n")

		outfile.write(f"from scipy.special import factorial as fac\n\n")

		for f in range(2, Vmax+1):
			outfile.write(f"invf{f} = 1/fac({f})\n")
		outfile.write(f"\n")

		outfile.write(f"def get_diagrams(unpert, normal_ordered_coefficients, amplitudes, angu_freq, mode):\n\n") 

		outfile.write(f"\tcontrib = dict({{}})\n\n")

		outfile.write(f"\thamil_delta, hamil_lambda, hamil_gamma = unpert\n")
		#Also have (Λa†a† + Λaa), rename to (1/2Ra†a† + 1/2Raa)
		outfile.write(f"\tR = 2*hamil_lambda\n")

		empty_coeffs = [["_    " for K in range(0, Vmax+1)] for L in range(0, Vmax+1)]
		
		for m in range(0, Vmax+1):
			for k in range(0, m+1):
				empty_coeffs[k][m-k] = f"Q_{k}_{m-k}"

		for K in range(0, Vmax+1):
			outfile.write(f"\t{', '.join(empty_coeffs[K])} = normal_ordered_coefficients[{K}]\n")

		outfile.write(f"\n")

		print_amplitudes = ", ".join([f"t{i}" for i in range(1, Tmax+1)])
		
		outfile.write(f"\t{print_amplitudes} = amplitudes\n\n")

		outfile.write(f"\tif mode in [\"all\", \"closed\"]:")

		for contrib in contrib_closed:
			new_categ = contrib[0]
			new_contr = gnrl_func.join_strings(string_array = contrib[1], 
										stagger = 1, 
										separator = " + ")
			new_len = len(contrib[1])
			
			outfile.write(f"\n\n")
			outfile.write(f"\t\t#Ground State | {new_len} entries.\n")

			outfile.write(f"\t\tcontrib[\"{new_categ[0]},{new_categ[1]}\"]")
			outfile.write(f" = {new_contr}\n")

		outfile.write(f"\n\n")
		outfile.write(f"\tif mode in [\"all\", \"open\", \"amplitude\"]:")

		for contrib in contrib_amplitude:
			new_categ = contrib[0]
			new_contr = gnrl_func.join_strings(string_array = contrib[1], 
										stagger = 1, 
										separator = " + ")
			new_len = len(contrib[1])

			outfile.write(f"\n\n")
			outfile.write(f"\t\t#Amplitude t{new_categ[0]} | {new_len} entries.\n")

			outfile.write(f"\t\tcontrib[\"{new_categ[0]},{new_categ[1]}\"]")
			outfile.write(f" = {new_contr}\n")

		outfile.write(f"\n\n")
		outfile.write(f"\tif mode in [\"all\", \"open\", \"other\"]:\n\n")

		offset_prev = None

		for contrib in contrib_other:
			new_categ = contrib[0]
			new_contr = gnrl_func.join_strings(string_array = contrib[1], 
										stagger = 1, 
										separator = " + ")
			new_len = len(contrib[1])

			offset = new_categ[0] - new_categ[1]
			
			if offset != offset_prev:
				outfile.write(f"\t\t#=-=\n\n")

			outfile.write(f"\t\t#Offset: {offset} | {new_len} entries.\n")
			
			outfile.write(f"\t\tcontrib[\"{new_categ[0]},{new_categ[1]}\"]")
			outfile.write(f" = {new_contr}\n\n")
		
			offset_prev = offset

		outfile.write(f"\treturn contrib\n")
			
	print(f"{file_name} successfully written in {get_time_total():.3f}s.")

	#print(number_generated_diagrams)









