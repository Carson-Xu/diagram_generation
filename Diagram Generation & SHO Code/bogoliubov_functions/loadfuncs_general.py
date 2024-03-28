"""

List of functions for general use.

> convert_bool(parameter)
> convert_tuple(parameter)
> create_dirs(folder_path)
> all_x_in_y(x,y)
> bulk_replace(input, cut_values, sub)
> join_strings(string_array, stagger, separator)
> pull_from_file(file_path, search_term, instance, show_line)

"""

import os
import numpy as np
import math as math
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import figure
from scipy.special import factorial as fac
import re

from configparser import ConfigParser
config_object = ConfigParser()


#-=-=-=-


def convert_bool(parameter):

	#Default value.
	bool_value = True

	true_values = ["True", True, "1", 1]
	false_values = ["False", False, "0", 0]

	#Check if the parameter is any version of True or False, and assign value accordingly.
	if parameter in true_values:
		bool_value = True

	if parameter in false_values:
		bool_value = False

	#If that fails, raise error and default to True.
	if parameter not in (*true_values, *false_values):
		print(f"Encountered boolean parameter with invalid value {parameter}.  Defaulting to True.")
		bool_value = True

	return bool_value


#-=-=-=-


def convert_tuple(parameter):

	#Default value.
	tuple_value = None

	# Try converting the parameter to a tuple immediately.
	try:	tuple_value = tuple(parameter)

	# If that fails, attempt to convert [parameter] to a tuple.
	except:
		try:	tuple_value = tuple([parameter])

		# If that fails, raise error and default to None.
		except:
			print(f"Encountered tuple parameter with invalid value {parameter}.  Defaulting to None.")
			tuple_value = None

	return tuple_value


#-=-=-=-


def tuple_from_string(parameter, type_func):

	#Default value.
	tuple_value = tuple([])

	if type_func == bool:
		type_func = convert_bool

	try:
		individuals = parameter.strip("[]()").replace(" ", "").split(",")

		if individuals == [""]:
			tuple_value = tuple([])

		if individuals != [""]:
			tuple_value = tuple([type_func(i) for i in individuals])

	except:
		print(f"Encountered tuple parameter with invalid value for given type.  Defaulting to None.")
		tuple_value = tuple([])

	return tuple_value


#-=-=-=-


def create_dir(folder_path):

	#Check for existing directory, to save files.
	directory = os.path.dirname(folder_path)

	#If directory does not exist, create it.
	if not os.path.exists(directory):
		os.makedirs(directory)


#-=-=-=-


def all_x_in_y(small_set,overall_set):
	return all(i in overall_set for i in small_set)


#-=-=-=-


def bulk_replace(input, cut_values, sub):

	bulk_string = input

	for i in cut_values:	
		bulk_string = bulk_string.replace(i, sub)

	return bulk_string


#-=-=-=-


def join_strings(string_array, stagger, separator):
	"""
	Joins a list of strings together, offering a stagger customization over the usual "".join() method.
	
	*-*-*-*-*
	
	string_array:	tuple
	- A list of strings to merge.
	
	stagger:		int or None		|	m
	- Inserts a line space for every m items.
	- If None, disables this behavior
	
	separator:		string
	- Connections to add, associated with the below amplitudes respectively.
	
	*-*-*-*-*
	
	Returns a large customized string for output.
	"""
	
	full_string = ""
	
	for index in range(0,len(string_array)):
		
		if stagger != None:
			if index > 0 and index%stagger == 0:
				full_string += f"\\\n\t\t\t"

		if index < len(string_array) - 1:
			full_string += f"{string_array[index]}{separator}"

		if index == len(string_array) - 1:
			full_string += f"{string_array[index]}"
	
	return full_string


#-=-=-=-


def np_strings(numpy_array):

	full_string = ""

	for i in numpy_array:

		full_string += f"[\t"
		for j in i:
			if (j == 0.0):
				full_string += f"{0.0},\t\t\t\t"
			elif (j < 0.0):
				full_string += f"{j:.8e},\t"
			else:
				full_string += f"{j:.8e},\t\t"

		full_string += f"]\n"

	return full_string
	

#-=-=-=-


def pull_from_file(file_path, search_term, instance, show_line):

	whitespace = ["\n", "\t", "+0j"]

	search_term_list = bulk_replace(input=search_term, cut_values=whitespace, sub=" ").strip().split("<...>")

	file = open(file_path, 'r')
	lines_all = [bulk_replace(input=line, cut_values=whitespace, sub=" ").strip() for line in file.readlines()]

	lines_of_interest = [line for line in lines_all if all_x_in_y(small_set=search_term_list, overall_set=line)]
	chosen_instance = lines_of_interest[instance]

	if show_line == True:
		print(f"Selected: \"{chosen_instance}\"")
	
	numerical_instance = re.findall(pattern="(-?)(\d+)(\.?)(\d+)?(e?)(-?)(\d+)?", 
							string=chosen_instance, 
							flags=re.IGNORECASE)

	numerical_instance = tuple([float(join_strings(i, None, "")) for i in numerical_instance])

	if len(numerical_instance) == 0:
		print(f"No numerical data found for instance {instance}:\n\"{chosen_instance}\"")

	#print(numerical_instance)
	return numerical_instance


#-=-=-=-


def invfac(num):

	value = 1/fac(num)

	return value


#-=-=-=-













