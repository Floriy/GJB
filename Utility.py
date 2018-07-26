# -*- coding: utf-8 -*-

from math import sqrt
from sys import argv, stdout
import os
import subprocess as sub
import textwrap
import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class cd:
	"""Context manager for changing the current working directory"""
	def __init__(self, newPath):
		self.newPath = os.path.expanduser(newPath)

	def __enter__(self):
		self.savedPath = os.getcwd()
		os.chdir(self.newPath)

	def __exit__(self, etype, value, traceback):
		os.chdir(self.savedPath)
	
def ComputeAverage(qtty):
	dataIterator = 0
	data = np.zeros((100,3))
	data[dataIterator,0] = np.mean(qtty)
	data[dataIterator,1] = np.var(qtty)
	data[dataIterator,2] = len(qtty)
	dataIterator = 1

	#start the blocking
	length = len(qtty)
	length = np.floor_divide(length,2)
	box = np.zeros(length)
	while len(box) > 2:
		qttyIterator = 0
		boxIterator = 0
		while boxIterator < length:
			box[boxIterator] = 0.5 * (qtty[qttyIterator] + qtty[qttyIterator+1])
			qttyIterator += 2
			boxIterator += 1
		data[dataIterator,0] = np.mean(box)
		data[dataIterator,1] = np.var(box)
		data[dataIterator,2] = len(box)
		dataIterator += 1
		#qtty becomes old box
		qtty = box
		#cut box into half
		length = np.floor_divide(length,2)
		box = np.zeros(length)
	#select nonzero elements and return them
	data = data[0:dataIterator,:]
	return data

#Completer for interactive session
class Completer:
	def __init__(self, words):
		self.words = words
		self.prefix = None
	def complete(self, prefix, index):
		if str(prefix) != self.prefix:
			# we have a new prefix!
			# find all words that start with this prefix
			self.matching_words = [
				w for w in self.words if w.startswith(prefix)
				]
			self.prefix = prefix
		try:
			return self.matching_words[index]
		except IndexError:
			return None

#Removes first indentation of a string
def RemoveUnwantedIndent(indentedtext):
	indentedtext = indentedtext[indentedtext.find('\n')+1:indentedtext.rfind('\n')]
	indentedtext = textwrap.dedent(indentedtext)
	return indentedtext

#Deletes the content of a file
def deleteContent(pfile):
    pfile.seek(0)
    pfile.truncate()


def tail(filef, n):
	"""
		Function returning the n last lines of a file
	"""
	tail_cmd	= "tail -n{0} {1}".format(n, file_name)
	output		= sub.Popen(tail_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
	text = output.stdout.readlines()
	
	return text

def grep_from_file(file_name, pattern):
	"""
		Function to grep from a file lines matching a pattern.
	"""
	grep_pattern_cmd	= """cat {0} | grep {1} """.format(file_name,pattern)
	grep_pattern_proc	= sub.Popen(read_index_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
	output				= read_index_proc.stdout.read()
	
	return output

def get_box_dimensions(file_name):
	"""
		Function to get the dimensions of the box from a .gro file
	"""
	output = str(tail(file_name, 1)[0], 'utf-8').replace(repr('\n'),'')
	output = output.strip().split()
	dimensions = [float(o) for o in output]
	
	return dimensions
	
	
#Function to search text under heading
def group_by_heading(some_source, heading):
		buffer= []
		FoundHeading = False
		for line in some_source:
				if heading in line and not FoundHeading:
						if buffer: yield buffer
						buffer= [ line ]
						FoundHeading = not FoundHeading
				if FoundHeading and '#' not in line:
						buffer.append( line )
				if FoundHeading and '#' in line:
					break
		yield buffer

# Do-order script integration
#2011.11.27 - Helgi I. Ingolfsson - Fix POPC and POPE (tail order is flipped in regular Martini 2.1 itp)
def read_gro(file, atoms):
	line_counter = 0
	number_of_particles = 0
	first, second = [], []
	for line in open(file):
		if line_counter == 1:
			number_of_particles = int(line)
		elif line_counter > 1 and line_counter < number_of_particles + 2:
			if line[10:15].strip() == atoms[0]:
				first.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
			elif line[10:15].strip() == atoms[1]:
				second.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
		line_counter += 1
	return [first, second]

### REAL STUFF
def do_order(trajfile, TimeRange, trajskip, XYZ , number_of_lipids, lipid_type, GROMACS_LOC_prefixPath):
	######if len(argv) != 10:
		####### coments/usage
		######print(RemoveUnwantedIndent('''
		######Compute (second rank) order parameter, defined as:

			######P2 = 0.5*(3*<cos²(theta)> - 1)

		######where "theta" is the angle between the bond and the bilayer normal.
		######P2 = 1      perfect alignement with the bilayer normal
		######P2 = -0.5   anti-alignement
		######P2 = 0      random orientation

		######All lipids defined in the "martini_v2.0_lipids.itp" file can be analyzed
		######with this script.
		######Usage: {0} <traj file> <initial time> <final time> <skip frames> <bilayer normal - xyz> <#lipids> <lipid type>

			######> {0} traj.xtc 0 10000 5 0 1 0 64 DSPC

		######will for example read a 10ns trajectory of 64 DSPC lipids, calculating the order parameter for 
		######every 5th frame and averaging the results. P2 will be calculated relative to the y-axis.

		######WARNING script will output all frames in one go, into files called frame_dump_XXX.gro and 
		######then remove them so don't have any other files with this name in the current directory.
		######'''.format('do-order.py')))
		######return(False)

	# snapshots
	topofile = trajfile.replace('.xtc','.tpr')
	
	# (normalized) orientation of bilayer normal
	orientation_of_bilayer_normal = [float(XYZ[0]), float(XYZ[1]), float(XYZ[2])]
	norm = sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
	for i in range(3):
		orientation_of_bilayer_normal[i] /= norm
	stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
		orientation_of_bilayer_normal[0], \
		orientation_of_bilayer_normal[1], \
		orientation_of_bilayer_normal[2]  \
	))

	# output legend
	phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
	phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
	# PCs
	if   lipid_type == "DAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
	elif lipid_type == "DHPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
	elif lipid_type == "DLPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
	elif lipid_type == "DOPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
	elif lipid_type == "DEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
	elif lipid_type == "DPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
	elif lipid_type == "DSPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
	elif lipid_type == "POPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
	# PEs
	elif lipid_type == "DHPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
	elif lipid_type == "DLPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
	elif lipid_type == "DOPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
	elif lipid_type == "DSPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
	elif lipid_type == "DPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
	elif lipid_type == "POPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
	# PPCS
	elif lipid_type == "PPCS": bond_names = " NC3-PO4 PO4-AM1 AM1-AM2 AM1-C1A GL2-D1B C1A-C2A C2A-C3A C3A-C4A D1B-C2B C2B-C3B C3B-C4B\n"
	# output legend
	output_legend = "  Frame" + bond_names 

	# write the stuff
	stdout.write("\n " + output_legend)
	stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n")
	output = open('order.dat', 'w')
	output.write(output_legend)
	output.write(("-"*(len(output_legend) - 1)) + "\n")

	# Output all frame using trjconv 
	stdout.write("Output all coordinate files \n")
	command = "echo {0} | {5}trjconv -f {1} -s {2} {3} -sep -skip {4} -pbc whole -o frame_dump_.gro > /dev/null".format(lipid_type, trajfile, topofile, TimeRange, trajskip, GROMACS_LOC_prefixPath.replace('g_',''))
	print(command)
	sub.call(command, shell=True)

	# For each dumped frame
	stdout.write("Starting P2 calculation")
	order_parameters = []
	file_count = 0
	bonds = []
	while True:
		filename = "frame_dump_" + str(file_count) + ".gro"
		if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
				break
		
		stdout.write("Taking care of snapshot %s \n" % filename)

		# compute order parameter for each bond, for each snapshot
		current_order_parameters = []
		# bonds respectively involved in the head,
		#                             in the junction head-tail,
		#                             in each tail
		bonds = []

		for bond_name in bond_names.split():
			bonds.append(bond_name.split("-"))

		for bond in bonds:

			# parse .gro file, grep bead coordinates
			first, second = read_gro(filename, bond)

			# compute order parameter for each lipid
			order_parameter = 0.0
			for i in range(number_of_lipids):
				# vector between the two previous beads (orientation doesn't matter)
				vector = [0.0, 0.0, 0.0]
				for j in range(3):
					vector[j] = first[i][j] - second[i][j]
				norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
				# compute projection on the bilayer normal
				projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
				# order parameter
				order_parameter += projection**2/norm2

			# compute final averaged order parameter
			# store everything in lists
			current_order_parameters.append(0.5*(3.0*(order_parameter/number_of_lipids) - 1.0))
		order_parameters.append(current_order_parameters)

		# write results
		results = "%7i" % file_count
		for order_parameter in current_order_parameters:
			results += "%8.3f" % order_parameter
		stdout.write(" " + results + "\n")
		output.write(results + "\n")

		os.remove(filename)
		file_count += 1
	# End while loop

	stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
	stdout.write("Snapshots analysis done.%s\n" % (" "*56))
	stdout.write("Computing averages...\n")

	# average order parameter
	averaged_order_parameters = []
	for i in range(len(bonds)):
		sum = 0.0
		for j in range(len(order_parameters)):
			sum += order_parameters[j][i]
		averaged_order_parameters.append(sum/len(order_parameters))

# average order parameter2
	averaged_order_parameters2 = []
	for i in range(len(bonds)):
		sum = 0.0
		for j in range(len(order_parameters)):
			sum += (order_parameters[j][i] - averaged_order_parameters[i] )*(order_parameters[j][i]-averaged_order_parameters[i])
		averaged_order_parameters2.append(sqrt(sum)/(len(order_parameters)-1))
		
	# write results
	stdout.write("\n           " + bond_names)
	stdout.write(("-"*(len(output_legend) - 1)) + "\n")
	output.write(("-"*(len(output_legend) - 1)) + "\n")
	results = "average  "
	for order_parameter in averaged_order_parameters:
		results += "%8.3f" % order_parameter
	stdout.write(" " + results + "\n")
	output.write(results + "\n")
	stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")

	results = "er(uncor)"
	for order_parameter2 in averaged_order_parameters2:
		results += "%8.3f" % order_parameter2
	stdout.write(" " + results + "\n")
	output.write(results + "\n")
	stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")

	# Write abs average order parameters <Sn> (for carbon chains only)
	# WARNING this works with currenct lipids (all have defined x5 none carbon bonds) but for manually added lipids this might not be true
	ave_chain_s = 0
	for i in averaged_order_parameters[3:]: 
		ave_chain_s += abs(i)
	average_txt = "Abs average order parameters for carbon chains <Sn> = %8.3f \n\n" % (ave_chain_s / (len(averaged_order_parameters)-3))
	stdout.write(average_txt)
	output.write(average_txt)
	stdout.write("Results written in \"order.dat\".\n")

	output.close()
	return(True)
