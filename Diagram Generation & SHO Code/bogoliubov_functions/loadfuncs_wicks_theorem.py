
# coding: utf-8

# In[1]:


### 1. IMPORT PACKAGES

import os
import numpy as np
import math as math
import string
import itertools as itt
from matplotlib import pyplot as plt
from scipy.special import factorial as fac
from configparser import ConfigParser
config_object = ConfigParser()


# In[2]:


def sum_string(list_item):

	result = ""
	for i in list_item:
		result += str(i)

	return result



def sum_dictionaries(list_item):

	values_to_check = sorted(list(set([j["operators"] for j in normal_ordered_product])), key = lambda oper: (len(oper), oper.count("a")))
	#print(f"v:{values_to_check}")

	total_dict_list = []
	example = normal_ordered_product[0]
	for oper in values_to_check:

		if type(example["coefficient"]) == int:
			associated_value = 0
			for j in normal_ordered_product:

				if j["operators"] == oper:

					associated_value += j["coefficient"]

		if type(example["coefficient"]) == str:
			associated_value = ""
			for j in normal_ordered_product:

				if j["operators"] == oper:

					add_val = j["coefficient"]
					associated_value += f" + {add_val}"

		total_dict = {"operators": oper,
					"coefficient": associated_value}
		total_dict_list.append(total_dict)

	return total_dict_list



def normal_order(creat_count, annih_count):

	#Counts all creation and annihilation operators.
	#creat_count = operator_string.count("C")
	#annih_count = operator_string.count("a")

	#Places creation and annihilation operators in order.
	normal_ordered_string = ""

	for i in range(0, creat_count):
		normal_ordered_string += "C"

	for j in range(0, annih_count):
		normal_ordered_string += "a"

	#Returns output.
	return normal_ordered_string



def count_connections(operator_string, connect_num, debug):

	core_len = len(operator_string)

	#Find positions of all the creation and annihilation operators in the string.
	core_creat_positions = tuple([i for i in range(0, core_len) if (operator_string[i] == "C")])
	core_annih_positions = tuple([i for i in range(0, core_len) if (operator_string[i] == "a")])

	#A single connection can occur for any annihilation operator, so long as its position is
	#less than (to the left of) the creation operator it connects to.
	connections_single = [(i,j) for i in core_annih_positions for j in core_creat_positions if (i < j)]

	#Obtain all nth-order connections by forming combinations of the above single connections.
	connections_preliminary = list(itt.combinations(connections_single, r=connect_num))

	#Test to see if any of the above nth-order connections are to a duplicate operator.
	#Two connections occurring to the same operator simultaneously are not allowed.
	chained_indices = [tuple(itt.chain.from_iterable(i)) for i in connections_preliminary]
	mask = [(len(i) == len(set(i))) for i in chained_indices]
	
	#Filter out the above such that only nth-order connections which use unique positions remain. 
	connections_allowed = [connections_preliminary[i] for i in range(0, len(connections_preliminary)) if mask[i]]

	#The coefficient will be given by how many unique connections are left.
	connections_count = len(connections_allowed)

	if debug == True:
		print(f"Operator String: {operator_string}")
		print(f"Connection Number: {connect_num}")
		print(f"Creat Positions: {core_creat_positions}")
		print(f"Annih Positions: {core_annih_positions}")
		print(f"Single Connections: {connections_single}")
		print(f"Prelim Nth Connections: {connections_preliminary}")
		print(f"Chained Indices: {chained_indices}")
		print(f"Mask: {mask}")
		print(f"Allowed Nth Connections: {connections_allowed}")
		print(f"Coefficient: {connections_count}")

	return connections_count


def connect(operator_string):

	#Copy strings.
	core_to_order = operator_string

	#Removes creation operators from the left and annihilation operators from the right.
	core_to_order = core_to_order.lstrip("C").rstrip("a")

	#Finds the maximum amount of connections to be made.
	#(The lowest number between available creation and annihilation operators.)
	max_connections = min(core_to_order.count("C"), core_to_order.count("a"))

	creat_count = operator_string.count("C")
	annih_count = operator_string.count("a")

	#The '0th connection' is simple.
	operator_term = [{"operators": normal_order(creat_count, annih_count),
					"coefficient": 1,
					"order": 0}]

	#Then connect and save information for higher orders.
	for conn in range(1, max_connections+1):

		dict = {"operators": normal_order(creat_count - conn, annih_count - conn),
				"coefficient": count_connections(operator_string, conn, False),
				"order": conn}

		operator_term.append(dict)


	return operator_term


# In[3]:


Vmax = 5
coeff = list(string.ascii_lowercase)


front_values = (np.sqrt(1/2))**Vmax

operators = ['C','a']

full_normal_order = []

for Vorder in range(1, Vmax+1):
	full_product = itt.product(operators, repeat=Vorder)
	condensed_product = [sum_string(i) for i in full_product]
	print(f"V{Vorder} term with {len(condensed_product)} operators: {condensed_product}")

	normal_ordered_product = [dict_item for i in condensed_product for dict_item in connect(i)]

	term_list = sum_dictionaries(normal_ordered_product)

	for dict_item in term_list:
		dict_item["Vorder"] = coeff[Vorder-1]

		full_normal_order.append(dict_item)

print("====")
for k in full_normal_order:
	print(k)

















