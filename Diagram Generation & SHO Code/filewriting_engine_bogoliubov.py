"""
-=-=-=-

Calculating (ground-state and excited-state) eigenvalues of single oscillator systems.


Changelog:

Version 1.0
Notes:
- First working version.
- Moved main code into its own file for ease of access.

-=-=-=-
"""

from configparser import ConfigParser
config_object = ConfigParser()

import os
import numpy as np
import math as math
from matplotlib import pyplot as plt
from scipy.special import factorial as fac

#Forces function file to reload itself with the newly updated config.
import importlib
import bogoliubov_functions.loadfuncs_bogoliubov as bogo_func
import bogoliubov_functions.loadfuncs_general as genl_func
importlib.reload(bogo_func)


#Run iterations until values of T are within...
#convergence_bogo = convergences_cc[0]

def write_full_file(perturbations, dim, angu_freq, bogo_first, starting_amplitudes, run_subtitle):

	hamil_delta = angu_freq
	hamil_lambda = 0.0
	hamil_gamma = 0.0

	unperturbed = [hamil_delta, hamil_lambda, hamil_gamma]
	
	a_lin, b_qdr, c_cub, d_qrt = perturbations

	starting_amplitude_t1, starting_amplitude_t2 = starting_amplitudes

	# zero_array = np.zeros(shape = 16)
	# zero_array[0], zero_array[1] = starting_amplitudes
	# t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16 = zero_array


	output_file_directory = f"Data Files/Saved Bogoliubov Data (Summary)/Bogoliubov Summary ({run_subtitle}_dim{dim})/"
	output_file_name = f"bogoliubov_data_summary (a={a_lin}, b={b_qdr}, c={c_cub}, d={d_qrt}, size={dim}).txt"
	file_path = f"{output_file_directory}{output_file_name}"

	genl_func.create_dir(folder_path=output_file_directory)

	bogo_func.file_break(file_path, "intialize")
	bogo_func.file_break(file_path, "energy_calc")

	with open(file_path, "a") as file:
		file.write("Given Parameters: (V = ax + bx^2 + cx^3 + dx^4)\n") 
		file.write(f"a = {a_lin}\n")
		file.write(f"b = {b_qdr}\n")
		file.write(f"c = {c_cub}\n")
		file.write(f"d = {d_qrt}\n")

	bogo_func.file_break(file_path, "config_params")

	unpert_to_use = unperturbed
	pert_to_use = perturbations

	if bogo_first == True:
		output_start = bogo_func.bogoliubov_single_run(file_path, True, unperturbed, perturbations, mode = "start")
		#print(output_start[0], output_start[1])
		unpert_to_use = output_start[0]
		pert_to_use = output_start[1]

		with open(file_path, "a") as file:
			file.write(f"\n\n-=-=-\n\nAdjusted Parameters:\n\n") 

			file.write(f"starting_amplitude_t1 = {starting_amplitude_t1}\n")
			file.write(f"starting_amplitude_t2 = {starting_amplitude_t2}\n\n")

			file.write(f"Delta = {unpert_to_use[0]}\n")
			file.write(f"Lambda = {unpert_to_use[1]}\n")
			file.write(f"Gamma = {unpert_to_use[2]}\n\n")

			file.write(f"a = {pert_to_use[0]}\n")
			file.write(f"b = {pert_to_use[1]}\n")
			file.write(f"c = {pert_to_use[2]}\n")
			file.write(f"d = {pert_to_use[3]}\n\n")

			file.write(f"-=-=-\n\n")

		print("Starting parameter adjustment complete.")


	bogo_func.file_break(file_path, "bogo_transform")

	# Excited state energies are calculated here for EoM-CC and each Bogo transformation.
	output = bogo_func.bogoliubov_converge(file_path, True, unpert_to_use, pert_to_use, dim)

	grnd_energy_zero_ampl, exct_energy_zero_ampl, grnd_energy_conv_ampl, exct_energy_conv_ampl = output
		
	with open(file_path, "a") as file:
		file.write("\n\n\n\n\n") 

	bogo_func.file_break(file_path, "final_output")

	with open(file_path, "a") as file:
		file.write("-= Arrays for Graphing: =-\n\n") 

		file.write(f"- Zero Amplitude Energies:\n\n")
		file.write(f"Ground State from Purely Closed Diagrams: {grnd_energy_zero_ampl}\n\n")
		file.write(f"All States from Matrix: \n{exct_energy_zero_ampl}\n\n")
		file.write(f"- Converged Amplitude Energies:\n\n")
		file.write(f"Ground State from Purely Closed Diagrams: {grnd_energy_conv_ampl}\n\n")
		file.write(f"All States from Matrix: \n{exct_energy_conv_ampl}\n\n")

		# exct_energy_conv_ampl_string = ""
		# for i in exct_energy_conv_ampl:
		# 	exct_energy_conv_ampl_string += f"{i}\t0\n"

		# file.write(f"{exct_energy_conv_ampl_string}"

	#print(output_i)

	print(f"Completed files for perturbation {perturbations} using {dim}x{dim} matrix and starting amplitudes {starting_amplitudes}.")

