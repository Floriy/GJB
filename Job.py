import subprocess as sub
import time
from shutil import copyfile
import numpy as np
import math
import Utility as ut
import os, sys, re
import glob
#from scipy.interpolate import griddata
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
from io import StringIO
#***********************************************************#
#***********************************************************#
#********************* PDB Files list **********************#
#***********************************************************#
#***********************************************************#

W_PDB = """
		CRYST1   04.000   04.000   04.000  90.00  90.00  90.00 P1           1
		ATOM      1  W   W   X   1      00.000  00.000  00.000 0.00  0.00
		END
		"""

OCO_PDB = """
		CRYST1   10.000   10.000   10.000  90.00  90.00 90.00 P 1           1
		ATOM      1  PC  OCO X   1       5.000   5.000   7.865 0.00  0.00          
		ATOM      2  C   OCO X   1       5.000   5.000   3.135 0.00  0.00          
		END
		"""

PW_PDB = """
		CRYST1   52.237   82.455   48.563  90.00  90.00  90.00 P1           1
		ATOM      1  W   PW  X   1       8.010  77.839  11.290 0.00  0.00       
		ATOM      2  WP  PW  X   1       7.810  77.479  12.630 0.00  0.00       
		ATOM      3  WM  PW  X   1       9.260  78.329  11.680 0.00  0.00       
		END
		"""

DSPC_PDB = """
		CRYST1  125.000  125.000  100.000  90.00  90.00  90.00 P 1           1
		ATOM      1  NC3 DSPCX   2     105.980  90.400  72.500  0.00  0.00
		ATOM      2  PO4 DSPCX   2     105.630  90.710  69.500  0.00  0.00
		ATOM      3  GL1 DSPCX   2     105.780  90.260  66.500  0.00  0.00
		ATOM      4  GL2 DSPCX   2     107.160  90.030  66.500  0.00  0.00
		ATOM      5  C1A DSPCX   2     106.230  90.050  63.500  0.00  0.00
		ATOM      6  C2A DSPCX   2     105.920  90.250  60.500  0.00  0.00
		ATOM      7  C3A DSPCX   2     106.390  90.500  57.500  0.00  0.00
		ATOM      8  C4A DSPCX   2     106.320  90.710  54.500  0.00  0.00
		ATOM      9  C5A DSPCX   2     106.050  90.770  51.500  0.00  0.00
		ATOM     10  C1B DSPCX   2     108.660  89.710  63.500  0.00  0.00
		ATOM     11  C2B DSPCX   2     108.450  89.500  60.500  0.00  0.00
		ATOM     12  C3B DSPCX   2     108.130  88.940  57.500  0.00  0.00
		ATOM     13  C4B DSPCX   2     108.740  88.990  54.500  0.00  0.00
		ATOM     14  C5B DSPCX   2     108.540  89.120  51.500  0.00  0.00
		END
		"""

DPPC_PDB = """
		TITLE     DPPC sim
		REMARK    THIS IS A SIMULATION BOX
		CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
		MODEL        1
		ATOM      1  NC3 DPPCX   1      13.368  14.653  25.970  1.00  0.00            
		ATOM      2  PO4 DPPCX   1      14.408  15.182  23.260  1.00  0.00            
		ATOM      3  GL1 DPPCX   1      13.608  14.862  19.630  1.00  0.00            
		ATOM      4  GL2 DPPCX   1      16.237  14.812  18.850  1.00  0.00            
		ATOM      5  C1A DPPCX   1      12.478  14.262  15.780  1.00  0.00            
		ATOM      6  C2A DPPCX   1      12.977  14.722  12.990  1.00  0.00            
		ATOM      7  C3A DPPCX   1      12.727  15.462  10.220  1.00  0.00            
		ATOM      8  C4A DPPCX   1      12.858  15.573   6.780  1.00  0.00            
		ATOM      9  C1B DPPCX   1      17.738  15.942  16.150  1.00  0.00            
		ATOM     10  C2B DPPCX   1      17.948  14.312  13.260  1.00  0.00            
		ATOM     11  C3B DPPCX   1      17.778  15.573  10.290  1.00  0.00            
		ATOM     12  C4B DPPCX   1      17.878  14.642   6.820  1.00  0.00            
		TER
		ENDMDL
		"""

DLPC_PDB = """
		CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
		ATOM      1  NC3 DLPCX   1       5.560   4.860  23.400  1.00  0.00            
		ATOM      2  PO4 DLPCX   1       5.220   6.260  20.220  1.00  0.00            
		ATOM      3  GL1 DLPCX   1       4.840   5.200  16.780  1.00  0.00            
		ATOM      4  GL2 DLPCX   1       7.930   5.720  16.400  1.00  0.00            
		ATOM      5  C1A DLPCX   1       4.150   4.880  13.220  1.00  0.00            
		ATOM      6  C2A DLPCX   1       3.880   5.100  10.250  1.00  0.00            
		ATOM      7  C3A DLPCX   1       3.950   5.290   6.840  1.00  0.00            
		ATOM      8  C1B DLPCX   1       9.590   5.700  13.540  1.00  0.00            
		ATOM      9  C2B DLPCX   1       8.850   6.380  10.430  1.00  0.00            
		ATOM     10  C3B DLPCX   1       9.110   5.310   6.910  1.00  0.00            
		END
		"""
		
		
SU_PDB = """
		CRYST1   04.000   04.000   04.000  90.00  90.00  90.00 P1           1
		ATOM      1  TEM TEMPS   1      00.000  00.000  00.000 0.00  0.00
		END
		"""


pdb_file_list = {'W':{'name':'water_single.pdb','content': W_PDB},
			   'OCO':{'name':'octanol_single.pdb','content': OCO_PDB},
			   'PW':{'name':'polwater_single.pdb','content': PW_PDB},
			   'DSPC':{ 'name':'dspc_single.pdb', 'content':DSPC_PDB},
			   'DPPC':{ 'name':'dppc_single.pdb', 'content':DPPC_PDB},
			   'DLPC':{ 'name':'dlpc_single.pdb', 'content':DPPC_PDB},
			   'SU':{ 'name':'su_single.pdb', 'content':SU_PDB}}

#***********************************************************#
#***********************************************************#
#***********************************************************#

class BaseProject(object):
	""" Main project class
	This class provides the basis for a job. Derived classes will have different fucntions based on base_project type.
	"""
	
	def __init__(self, sample, softwares, path_to_default):
		try:
			#Softwares
			self.softwares			= softwares
			self.path_to_default	= path_to_default
			self.protocol			= sample['PROTOCOL']
			
			self.job_id				= sample['JOBID']
			self.job_num			= sample['JOBNUM']
			self.job_seed			= sample['SEED']

			self.node				= sample['NODE']
			self.nbnodes			= sample['NBNODES']
			self.ppn				= sample['PPN']

			self.mdrun_opt			= sample['MDRUN_OPT']

			self.sample_molecules	= sample['SAMPLE']
			self.type				= sample['TYPE']
			
			self.geometry			= sample['GEOMETRY']
			self.dimensions			= {'LX': float(sample['GEOMETRY']['LX']),
										'LY': float(sample['GEOMETRY']['LY']), 
										'LZ': float(sample['GEOMETRY']['LZ'])}
			self.bottom_z			= 0.0
			
			self.sol_thickness		= None
			if self.type.startswith('MONO'):
				self.sol_thickness = float(sample['GEOMETRY']['SOL_THICKNESS'])
			
			self.lz_vaccum			= None
			if 'LzVac' in sample['GEOMETRY']:
				self.lz_vaccum	= float(sample['GEOMETRY']['LzVac'])
			
			self.fillmode			= None
			if 'FILLMODE' in sample['GEOMETRY']:
				
				self.fillmode = sample['GEOMETRY']['FILLMODE']
				
				if self.fillmode == 'PARTIAL-X':
					assert(False),"PARTIAL-X not implemented yet"
					
					if 'MIN' in sample['GEOMETRY']:
						
						self.minsolvent = float(sample['GEOMETRY']['MIN'])
					else:
						self.minsolvent = 0.0
						
					if 'MAX' in sample['GEOMETRY']:
						self.maxsolvent = float(sample['GEOMETRY']['MAX'])
					else:
						self.maxsolvent = float(self.dimensions['LX'])
						
				if self.fillmode == 'PARTIAL-Y':
					assert(False),"PARTIAL-Y not implemented yet"
					if 'MIN' in sample['GEOMETRY']:
						self.minsolvent = float(sample['GEOMETRY']['MIN'])
					else:
						self.minsolvent = 0.0
						
					if 'MAX' in sample['GEOMETRY']:
						self.maxsolvent = float(sample['GEOMETRY']['MAX'])
					else:
						self.maxsolvent = float(self.dimensions['LY'])
						
				if self.fillmode == 'PARTIAL-Z':
					print("PARTIAL-Z not tested with LIPIDS/DEFO/SU systems ! check bottom_z contraints before use ")
					
					if 'MIN' in sample['GEOMETRY']:
						self.minsolvent = float(sample['GEOMETRY']['MIN'])
					else:
						self.minsolvent = 0.0
						
					if 'MAX' in sample['GEOMETRY']:
						self.maxsolvent = float(sample['GEOMETRY']['MAX'])
					else:
						self.maxsolvent = float(self.dimensions['LZ'])
						
					if self.type.startswith('SOLVENT'):
						self.bottom_z = self.minsolvent          
						
			self.radius				= None
			if 'Radius' in sample['GEOMETRY']:
				self.radius		= sample['GEOMETRY']['Radius']
			
			#defo related
			self.defo				= None
			
			#su related
			self.su					= None
			self.nb_su				= None
			self.rugosity			= None
			# monolayer related
			self.mono				= None
			
			#wall related
			self.wall				= None
			
			#thermostat related
			self.thermostat			= None
			
			if 'DEFO' in sample:
				self.defo		= sample['DEFO']
				self.defo_dict	= {'DEFB' : {} }
			
			if 'SU' in sample:
				self.su			= sample['SU']
				
				if 'RUGOSITY' in sample['SU']:
					entries			= sample['SU']['RUGOSITY'].split('|')
					self.rugosity	= {}
					
					for e in entries:
						parameter_name	= e.split(':')[0]
						parameter_value	= e.split(':')[1]
						
						if ut.is_number(parameter_value):
							self.rugosity.update({parameter_name : float(parameter_value)})
						else:
							self.rugosity.update({parameter_name : parameter_value})
				
			if 'MONO' in sample:
				self.mono	= sample['MONO']
			
			if 'WALL' in sample:
				self.wall		= sample['WALL']
				
				
				if self.su is None:
					self.bottom_z	= 3.0
				
			if 'THERMOSTAT' in sample:
				self.thermostat	= sample['THERMOSTAT']

			self.protocol			= sample['PROTOCOL']
			
			self.system				= None
			self.shell				= 3.0
			
		except RuntimeError as e:
			print(e)
			
		self.lipid_list				= ['DSPC','DPPC','DLPC']
		self.solvent_list			= ['W','OCO','PW']
		self.substrate_list			= ['SSUP','SSUN','SSna','SSn0']
		
		self.packmol_seed			= None
		self.packmol_input			= None
		
		self.lipid_types			= None
		self.solvent_types			= None
		
		if sample['SEED'] == 'time':
			self.packmol_seed	= int(time.time())
			
		elif sample['SEED'] == 'jobnum':
			self.packmol_seed	= sample['JOBNUM'] * 123456789
			
		
		if self.packmol_seed is not None:
			self.packmol_input	= """
								# Packmol seed was set using {0}
								seed {1:d}
								tolerance 5.0
								nloop0   1000
								filetype pdb\n""".format(sample['SEED'], self.packmol_seed)
		else:
			self.packmol_input	= """
								# Packmol seed was set using default
								nloop0   1000
								tolerance 3.0
								filetype pdb\n"""
								
		self.packmol_input			+= """
								sidemax 1000000"""
	
	
	
	def from_input(self, inputs):
		""" Method to create sample from already existing configuration
		"""
		self.lipid_type = None
		self.solvent_type = None
		
		for lipid in self.lipid_list:
			if lipid in self.sample_molecules:
				self.lipid_type = lipid
		
		if self.lipid_type is not None:
			self.lipid_types = [self.lipid_type]
				
		self.solvent_types = []
		
		for sol in self.solvent_list:
			if sol in self.sample_molecules:
				self.solvent_types.append(sol)

		
		print("Looking for .gro file: {0}".format(inputs['GRO']))
		file_found = glob.glob(inputs['GRO'])
		
		assert(file_found), "{0} not found".format(inputs['GRO'])
		gro_file = file_found[0]
		
		path = ""
		if '/' in gro_file:
			path = '/'.join(gro_file.split('/')[:-1])
		
		gro_file = gro_file.split('/')[-1]
		print("Using: {0}".format(gro_file) )
		
		self.system = '_'.join(gro_file.split('_')[:-1])
		
		if 'TOP' not in inputs:
			top_file = self.system + '.top'
		else:
			top_file = inputs['TOP']
		
		if 'NDX' not in inputs:
			ndx_file = self.system + '.ndx'
		else:
			ndx_file = inputs['NDX']
		
		if 'TPR' not in inputs:
			tpr_file = gro_file.replace('.gro','.tpr')
		else:
			tpr_file = inputs['TPR']
			
		try:
			sub.call( "cp {0} ./{1}.gro".format( '{0}/{1}'.format(path, gro_file), self.system ), shell=True)
			sub.call( "cp {0} ./{1}.tpr".format( '{0}/{1}'.format(path, tpr_file), self.system ), shell=True)
			sub.call( "cp {0} .".format( '{0}/{1}'.format(path, ndx_file) ), shell=True)
		
		except IOError as e:
			print(e)
		
		read_index_cmd = """cat {0} | grep "\[" """.format(ndx_file)
		read_index_proc = sub.Popen(read_index_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
		read_index_out = read_index_proc.stdout.read()
		
		print(read_index_out.splitlines())
		nb_index = len(read_index_out.splitlines())
		
		#Retrieving the index
		self.nb_index = nb_index
		
		create_su		= False
		generate_su		= False
		add_monolayer	= False
		
		try:
			if self.defo:
				sub.call( "cp {0}/defos_topo.itp .".format(path), shell=True)
				sub.call( "cp {0}/defo_posres_DEFB_gen.itp .".format(path), shell=True)
				
			if self.su:
				su_posres = "{0}/su_posres_gen.itp".format(path)
				
				if 'GEN_SU' in inputs:
					generate_su		= True
				elif os.path.isfile(su_posres):
					sub.call( "cp {0}/su_posres_gen.itp .".format(path), shell=True)
				else:
					create_su	= True
			
			if self.mono:
				add_monolayer = True
			
		except IOError as e:
			print(e)
		
		translate = None
		
		# Get the dimensions of the box from the input gro file
		
		
		
		dimensions = str(ut.tail(file_found[0], 1)[0], 'utf-8').replace(repr('\n'),'')
		dimensions = dimensions.split()
		dimensions = [float(D) for D in dimensions]
		
		
		if 'TRANS' in inputs:
			translate = [float(T)/10. for T in inputs['TRANS'].strip().split()]
			
			translate_cmd = str("printf 'System\\n'  | {0} trjconv -f {1}.gro -s {1}.tpr -o {1}.gro -n {1}.ndx -pbc mol "
								"-trans {2} {3} {4}").format(self.softwares['GROMACS_LOC'], 
																self.system,
																translate[0], translate[1], translate[2],
																)
			sub.call(translate_cmd, shell=True)
		
		
		#Increase box size in each dimension
		if 'NEWBOX' in inputs:
			increment	= [float(N)/10. for N in inputs['NEWBOX'].strip().split()]
			newbox		= np.add(increment, dimensions)
			
			print("newbox :", newbox)
			new_box_cmd = str("printf 'System\\n'  | {0} trjconv -f {1}.gro -s {1}.tpr -o {1}.gro -n {1}.ndx -pbc mol "
								"-box {2} {3} {4}").format(self.softwares['GROMACS_LOC'], 
																self.system,
																newbox[0], newbox[1], newbox[2])
			
			dimensions	= newbox
			sub.call(new_box_cmd, shell=True)
		
		
		if 'THRESHOLD-Z' in inputs and 'THRESHOLD_ATOMS' in inputs: #TO MODIFY
			thresholds			= float(inputs['THRESHOLD-Z'])/10.
			thresholds_atoms	= inputs['THRESHOLD_ATOMS'].strip().split()
			#[float(T)/10. for T in inputs['THRESHOLD'].strip().split()]
			
			atoms		= ""
			nb_atoms	= 0
			modified_topology	= {}
			serials				= []
			with open("{0}.gro".format(self.system),'r') as gro_input:
				for line in gro_input:
					split_line = line.split()

					if len(split_line) < 3 or line.startswith('Generated'): continue
					
					# Add the pbc
					if len(split_line) == 3:
						atoms	+= line
						continue
					
					#Look at the z below the thickness
					#x		= float(line.split()[-6])
					#y		= float(line.split()[-5])
					z		= float(line.split()[-4])
					name	= line[5:9].strip()
					serial	= line[0:5].strip()
					
					# Dictionnary to update the topology
					if name not in modified_topology:
						modified_topology[name] = {'atoms': 0, 'molecules': 0}
					
					if name in thresholds_atoms:
						#print("{0} in {1}".format(name, thresholds_atoms))
						if z < thresholds:
							atoms += line
							modified_topology[name]['atoms'] += 1
							nb_atoms				+= 1
							
							if serial not in serials:
								serials.append(serial)
								modified_topology[name]['molecules']	+= 1
					else:
						atoms += line
						modified_topology[name]['atoms'] += 1
						nb_atoms				+= 1
						
						if serial not in serials:
							serials.append(serial)
							modified_topology[name]['molecules']	+= 1
					
			
			after_threshold_content = "AFTER_THRESHOLD-Z of {2}\n{0}\n{1}".format(nb_atoms, atoms, thresholds)
			
			system	= self.system.split('_')[0]
			for name in modified_topology:
				system	+= "_{0}{1}".format(modified_topology[name]['molecules'], name)
			
			with open('{0}.gro'.format(system),'w') as output_gro:
				output_gro.write(after_threshold_content)
			
			for lipid in self.lipid_list:
				if lipid in modified_topology:
					self.lipid_type = lipid
		
			if self.lipid_type is not None:
				self.lipid_types = [self.lipid_type]
					
			self.solvent_types = []
			
			for sol in self.solvent_list:
				if sol in modified_topology:
					self.solvent_types.append(sol)
					
			self.sample_molecules = {}
			for name in modified_topology:
				self.sample_molecules[name] = modified_topology[name]['molecules']
			
				
			grotopdb_cmd = "{0} editconf -f {1}.gro -o {1}.pdb -c no".format(self.softwares['GROMACS_LOC'], system)
			sub.call(grotopdb_cmd, shell=True)
			
			make_ndx_cmd	= """printf "q\\n" | {0} make_ndx -f {1}.gro -o {1}.ndx""".format(self.softwares['GROMACS_LOC'], system)
			sub.call(make_ndx_cmd, shell=True)
			
			list_of_ndx	= None
			with open("{0}.ndx".format(system)) as ndx:
				content = ndx.read()
				regex	= re.compile("\[ .*", re.MULTILINE)
				list_of_ndx = regex.findall(content)
				
			
			list_of_ndx	= [ndx[1:-1].strip() for ndx in list_of_ndx]
			nb_index	= len(list_of_ndx)
			
			index_for_bilayer	= []
			for lipid in self.lipid_list:
				if lipid in list_of_ndx:
					index_for_bilayer.append(list_of_ndx.index(lipid))
			
			index_for_defo	= []
			if "DEFB" in list_of_ndx:
				index_for_defo.append(list_of_ndx.index("DEFB"))
				
			index_for_substrate	= []
			for su in self.substrate_list:
				if su in list_of_ndx:
					index_for_substrate.append(list_of_ndx.index(su))
			
			add_to_ndx		= ""
			name_new_ndx	= ""
			
			
			if index_for_bilayer:
				index_for_bilayer	= [str(ndx) for ndx in index_for_bilayer]
				add_to_ndx		+= repr("{0}\n".format(" | ".join(index_for_bilayer)))
				name_new_ndx	+= repr("name {0} bilayer\n".format(nb_index))
				nb_index		+= 1
			
			if index_for_defo:
				index_for_defo	= [str(ndx) for ndx in index_for_defo]
				add_to_ndx		+= repr("{0}\n".format(" | ".join(index_for_defo)))
				name_new_ndx	+= repr("name {0} defo\n".format(nb_index))
				nb_index		+= 1
			
			if index_for_substrate:
				index_for_substrate	= [str(ndx) for ndx in index_for_substrate]
				add_to_ndx		+= repr("{0}\n".format(" | ".join(index_for_substrate)))
				name_new_ndx	+= repr("name {0} su\n".format(nb_index))
				nb_index		+= 1
				
			make_index_input	= add_to_ndx + name_new_ndx + repr("q\n")
			
			print(make_index_input)
			upating_ndx_cmd = """printf {2} | {0} make_ndx -f {1}.gro -n {1} -o {1}.ndx""".format(self.softwares['GROMACS_LOC'], system, make_index_input)
			
			sub.call(upating_ndx_cmd, shell=True)
			
			#name_to_remove	= thresholds_atoms + ["System"]
			#index_to_remove	= []
			
			#for name in name_to_remove:
				#for ndx in list_of_ndx:
					#if name in ndx:
						#index_to_remove.append(list_of_ndx.index(ndx))
			
			#new_index_input = ""
			##removing old indexes for atoms we removed
			#updating_index	= 0
			#for index in index_to_remove:
				#new_index_input	+= repr("del {0}\n".format(index-updating_index))
				#updating_index	+= 1
			

			##updating indexes for the modified atoms
			#new_index		= 0
			#for name in thresholds_atoms:
				#new_index_input	+= repr("r {0}\n".format(name))
				#new_index	+= 1
			
			
			##Retrieving the index
			#self.nb_index = len(list_of_ndx) - updating_index + new_index
			
			#group_index = [str(i) for i in range(0, self.nb_index, 1)]
			
			#if not create_su:
				#new_index_input	+= repr("{0}\nname {1} System\n".format(" | ".join(group_index), int(group_index[-1])+1))
			
			#new_index_input	+= repr("q\n")
			
			#print(new_index_input)
			#upating_ndx_cmd = """printf {3} | {0} make_ndx -f {1}.gro -n {2} -o {1}.ndx""".format(self.softwares['GROMACS_LOC'], system,
																				#self.system, new_index_input)
			
			
			#sub.call(upating_ndx_cmd, shell=True)
			
			self.system = system
			
		elif 'GEN_SU' in inputs:
			pass
		else:
			#if 'THRESHOLD-Z' not in inputs or 'THRESHOLD_ATOMS' not in inputs:
			if 'THRESHOLD-Z' in inputs:
				assert('THRESHOLD_ATOMS' in inputs), "You forgot to put THRESHOLD_ATOMS in your Parameters.csv"
			if 'THRESHOLD_ATOMS' in inputs:
				assert('THRESHOLD-Z' in inputs), "You forgot to put THRESHOLD-Z in your Parameters.csv"
			
			modified_topology	= {}
			serials				= []
			with open("{0}.gro".format(self.system),'r') as gro_input:
				for line in gro_input:
					split_line = line.split()
					
					if len(split_line) < 3 or line.startswith('Generated'): continue
					
					# Add the pbc
					if len(split_line) == 3: continue
					
					z		= float(line.split()[-4])
					name	= line[5:9].strip()
					serial	= line[0:5].strip()
					
					# Dictionnary to update the topology
					if name not in modified_topology:
						modified_topology[name] = {'atoms': 0, 'molecules': 0}
					

					modified_topology[name]['atoms'] += 1
					
					if serial not in serials:
						serials.append(serial)
						modified_topology[name]['molecules']	+= 1
				
			
			for lipid in self.lipid_list:
				if lipid in modified_topology:
					self.lipid_type = lipid
		
			if self.lipid_type is not None:
				self.lipid_types = [self.lipid_type]
					
			self.solvent_types = []
			
			for sol in self.solvent_list:
				if sol in modified_topology:
					self.solvent_types.append(sol)
					
			self.sample_molecules = {}
			for name in modified_topology:
				self.sample_molecules[name] = modified_topology[name]['molecules']
				
		
		if create_su:
			# pdb file for su
			assert(len(self.su['SuType']) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
			
			su_type = self.su['SuType'] + (4 - len(self.su['SuType']))*' '
			
			atom_type = None
			if self.su['SuType'].startswith('SU'):
				atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
			elif self.su['SuType'].startswith('S'):
				atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
				
			
			#if 'SCALE' in inputs:
				#scaleby = inputs['SCALE'].strip().split()
				#scalegro_cmd = "{0} editconf -f {1}.gro -o {1}.gro -scale {2} {3} {4}".format(self.softwares['GROMACS_LOC'], self.system,
																								#scale_by[0], scale_by[1], scale_by[2])
				#sub.call(scalegro_cmd, shell=True)
			
			grotopdb_cmd = "{0} editconf -f {1}.gro -o {1}.pdb -c no".format(self.softwares['GROMACS_LOC'], self.system)
			sub.call(grotopdb_cmd, shell=True)
			
			#self.nb_su = int( float(self.su['Density']) * dimensions[0] * dimensions[1] * float(self.su['Thickness']) / 10 )
			
			su_atoms = ""
			self.nb_su = 0
			total_thickness = 0.0
			if self.rugosity is not None:
				
				su_gro				= open(inputs['SU'],'r').readlines()
				read_from			= StringIO(''.join(su_gro[2:-1]))
				x_gro, y_gro, z_gro	= np.loadtxt(read_from, usecols=(3,4,5), unpack=True)
				
				#x_gro, y_gro = zip(*sorted(zip(x_gro, y_gro)))
				box_dim = str(ut.tail(inputs['SU'], 1)[0], 'utf-8')
				box_dim = box_dim.strip().split()
				box_dim = [ float(dim) for dim in box_dim]
				
				su_input = open(inputs['SU'],'r').readlines()
				
				if self.rugosity['type'].startswith('1dsine'):
					
					def sine1D(x, box_x, lx):
						to_radian	= 2. * math.pi
						freq		= to_radian * self.rugosity['lx']/10.
						height		= self.rugosity['z']/10.
						thickness	= float(self.su['Thickness']) / 10.
						result		= thickness + height * np.sin(freq*x/box_x)
						result[np.where(result < thickness)] = thickness
						
						return result
					
					x_height = sine1D(x_gro, box_dim[0], self.rugosity['lx'])
					for x,y,z,t,line in zip(x_gro, y_gro, z_gro, x_height, su_input[2:-1]):

						if z <= t:
							su_atoms += line
							self.nb_su += 1
				
					total_thickness = float(self.su['Thickness']) + self.rugosity['z']
				
				if self.rugosity['type'].startswith('1dbool'):
					
					def sine1D(x, box_x, lx):
						to_radian	= 2. * math.pi
						freq		= to_radian * self.rugosity['lx']/10.
						height		= self.rugosity['z']/10.
						thickness	= float(self.su['Thickness']) / 10.
						result		= thickness + height * np.sin(freq*x/box_x)
						result[np.where(result > thickness)] = thickness + height
						result[np.where(result < thickness)] = thickness
						
						return result
					
					x_height = sine1D(x_gro, box_dim[0], self.rugosity['lx'])
					for x,y,z,t,line in zip(x_gro, y_gro, z_gro, x_height, su_input[2:-1]):

						if z <= t:
							su_atoms += line
							self.nb_su += 1
							
					total_thickness = float(self.su['Thickness']) + self.rugosity['z']
				
				#if self.rugosity['type'].startswith('2dbool'):
					
					#def sine1D(x, box_x, lx):
						#to_radian	= 2. * math.pi
						#freq		= to_radian * self.rugosity['lx']/10.
						#height		= self.rugosity['z']/10.
						#thickness	= float(self.su['Thickness']) / 10.
						#result		= thickness + height * np.sin(freq*x/box_x)
						#result[np.where(result > thickness)] = thickness + height
						#result[np.where(result < thickness)] = thickness
						
						#return result
					
					#x_height = sine1D(x_gro, box_dim[0], self.rugosity['lx'])
					#for x,y,z,t,line in zip(x_gro, y_gro, z_gro, x_height, su_input[2:-1]):

						#if z <= t:
							#su_atoms += line
							#self.nb_su += 1
							
					#total_thickness = float(self.su['Thickness']) + self.rugosity['z']
				
			
			else:
				with open(inputs['SU'],'r') as su_input:
					for line in su_input:
						split_line = line.split()
						
						#print(split_line)
						if len(split_line) < 3: continue
						
						# Add the pbc
						if len(split_line) == 3:
							split_line[2] = float(self.su['Thickness'])/10.
							su_atoms += "  {0}  {1}  {2}".format(split_line[0], split_line[1], split_line[2])
							continue
						
						#Look at the z below the thickness
						z = float(line.split()[-4])
						if z < (float(self.su['Thickness']) / 10.):
							su_atoms += line
							self.nb_su += 1
				
				total_thickness = float(self.su['Thickness'])
			
			su_atoms = su_atoms.replace("TEMP", self.su['SuType'])
			su_atoms = su_atoms.replace("TEM", atom_type)
			
			substrate_content = "SUBSTRATE\n{0}\n{1}".format(self.nb_su, su_atoms)
			system = self.system + "_{0}{1}{2}".format(self.nb_su, self.su['SuType'], self.su['Version'])
			
			with open('substrate.gro','w') as substrate_file:
				substrate_file.write(substrate_content)
				
			grotopdb_cmd = "{0} editconf -f substrate.gro -o substrate.pdb -c no".format(self.softwares['GROMACS_LOC'])
			sub.call(grotopdb_cmd, shell=True)
			
			shell = 7.0
			self.packmol_input += """
								avoid_overlap yes
								
								output {0}.pdb
								
								# Input structure
								structure {1}.pdb
									fixed 0. 0. {2} 0.0 0.0 0.0
								end structure
								
								# Substrate
								structure substrate.pdb
									fixed 0. 0. 0. 0. 0. 0.
								end structure
								
								""".format(system, self.system, total_thickness+shell)
								
			
			#su_pdb_content = pdb_file_list['SU']['content'].replace("TEMP", su_type)
			#su_pdb_content = su_pdb_content.replace("TEM", atom_type)
			
			#with open(pdb_file_list['SU']['name'],'w') as su_pdb:
				#su_pdb.write(ut.RemoveUnwantedIndent(su_pdb_content))
								
			with open("packmol_{0}.input".format(self.system), 'w') as packmol_file:
				packmol_file.write(ut.RemoveUnwantedIndent(self.packmol_input))
								
			packmol_cmd = str("{0} < packmol_{1}.input > packmol_{1}.output "
					).format(self.softwares['PACKMOL'], self.system)
			
			sub.call(packmol_cmd, shell=True)
			
			
			pdbtogro_cmd = "{0} editconf -f {1}.pdb -o {1}.gro -box {2} {3} {4} -c no".format(self.softwares['GROMACS_LOC'], system,
																				dimensions[0], dimensions[1], 
																				dimensions[2] + (total_thickness + 2*shell)/10)
			sub.call(pdbtogro_cmd, shell=True)
			
			group_index = [str(i) for i in range(0, self.nb_index-1, 1)]
			
			make_ndx_input = """ "r {0}\\nname {1} su\\ndel 0\\ndel 0\\n{2}\\nname {3} System\\nq\\n" """.format(self.su['SuType'],
																									self.nb_index,
																									" | ".join(group_index),
																									int(group_index[-1])+1)
			
			add_su_ndx_cmd = """printf {3} | {0} make_ndx -f {1}.gro -n {2} -o {1}.ndx""".format(self.softwares['GROMACS_LOC'], system,
																				self.system, make_ndx_input)
			
			sub.call(add_su_ndx_cmd, shell=True)
		
			self.system = system
		
		if generate_su:
			""" Create a substrate for the given input"""
			su_density	= float(self.su['Density'])
			thickness	= float(self.su['Thickness'])
			
			assert(len(self.su['SuType']) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
			
			su_type = self.su['SuType'] + (4 - len(self.su['SuType']))*' '
			
			atom_type = None
			if self.su['SuType'].startswith('SU'):
				atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
			elif self.su['SuType'].startswith('S'):
				atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
			
			self.nb_su	= int(su_density * dimensions[0] * dimensions[1] * thickness/10.)
			
			system = "substrate_{0}{1}{2}T{3}A".format(self.nb_su, self.su['SuType'], self.su['Version'], int(thickness))
			
			shell = 1.0
			self.packmol_input += """
								
								output {0}.pdb
								
								structure {1}
									chain S
									number {2:g}
									inside box 0. 0. {3}  {4} {5} {6}
								end structure
								""".format(system, pdb_file_list['SU']['name'],
											self.nb_su, self.bottom_z + shell,
											dimensions[0]*10 - shell, dimensions[1]*10 - shell, thickness - shell)
			
			su_pdb_content = pdb_file_list['SU']['content'].replace("TEMP", su_type)
			su_pdb_content = su_pdb_content.replace("TEM", atom_type)
			
			with open(pdb_file_list['SU']['name'],'w') as su_pdb:
				su_pdb.write(ut.RemoveUnwantedIndent(su_pdb_content))
				
			with open("packmol_substrate.input".format(self.system), 'w') as packmol_file:
				packmol_file.write(ut.RemoveUnwantedIndent(self.packmol_input))
				
			packmol_cmd = str("{0} < packmol_substrate.input > packmol_substrate.output "
					).format(self.softwares['PACKMOL'], self.system)
			
			sub.call(packmol_cmd, shell=True)
			
			
			pdbtogro_cmd = "{0} editconf -f {1}.pdb -o {1}.gro -box {2} {3} {4} -c no".format(self.softwares['GROMACS_LOC'], system,
																				dimensions[0], dimensions[1], thickness/10. )
			
			sub.call(pdbtogro_cmd, shell=True)
			
			make_ndx_input = """ "r {0}\\nname 3 su\\nq\\n" """.format(self.su['SuType'])
			
			add_su_ndx_cmd = """printf {2} | {0} make_ndx -f {1}.gro -o {1}.ndx""".format(self.softwares['GROMACS_LOC'], system, 
																				 			make_ndx_input)
			
			sub.call(add_su_ndx_cmd, shell=True)
			
			self.system = system
			
		if add_monolayer:
			monolayer_z	= None
			
			if 'THRESHOLD-Z' in inputs:
				monolayer_z = thresholds
			else:
				monolayer_z = dimensions[2]
			
			grotopdb_cmd = "{0} editconf -f {1}.gro -o {1}.pdb -c no".format(self.softwares['GROMACS_LOC'], self.system)
			sub.call(grotopdb_cmd, shell=True)
			
			if 'APL' in self.mono:
				self.nb_lipid_monolayer = int( dimensions[0] * dimensions[1] / 100 / float(self.mono['APL']) )
			
			elif 'NBLIPIDS' in self.mono:
				self.nb_lipid_monolayer = int(self.mono['NBLIPIDS'])
			
			else:
				assert(False), "Your forgot to set the area per lipid (APL) or number of lipids in the monolayer to add"
			
			self.tmt = None
			self.dz = None
			
			if self.lipid_type == 'DSPC':
				self.tmt = 30.0
				self.dz = 7.0
			elif self.lipid_type in ('DPPC','DLPC'):
				self.tmt = 30.0
				self.dz = 10.0
			
			shell = 3.0
			defo_mono_packmol_input = ""
			if self.defo is not None:
				if self.defo['Height'] in ["followM", "box", "mono"]:
					defo_mono_packmol_input = """
								outside cylinder  {0} {1}  0. 0.  0.  1.  {2}  {3}
								""".format(self.dimensions[0]/2., self.dimensions[1]/2.,
										self.defo['Radius'], self.dimensions[2])
								
			if dimensions[2] < self.tmt + monolayer_z:
				dimensions[2] += self.tmt
			
			system += "_MONO_{1}{0}".format(self.lipid_type, self.nb_lipid_monolayer)
			self.packmol_input += """
								avoid_overlap yes
								
								output {0}.pdb
								
								# Input structure
								structure {1}.pdb
									fixed 0. 0. 0. 0.0 0.0 0.0
								end structure
								
								# Monolayer added
								structure {2}
									chain C
									resnumbers 3
									number {3:g}
									inside box 0. 0.  {4}  {5} {6} {7}
									constrain_rotation x 180 10
									constrain_rotation y 180 10
									{8}
								end structure
								
								""".format(system, self.system, pdb_file_list[self.lipid_type]['name'],
				   							self.nb_lipid_monolayer, monolayer_z,
				   							dimensions[0],dimensions[1],dimensions[2],defo_mono_packmol_input)
								
			
			#if 
								
			with open("packmol_{0}.input".format(system), 'w') as packmol_file:
				packmol_file.write(ut.RemoveUnwantedIndent(self.packmol_input))
				
			with open(pdb_file_list[self.lipid_type]['name'],'w') as lipid_pdb:
				lipid_pdb.write( ut.RemoveUnwantedIndent(pdb_file_list[self.lipid_type]['content']) )
								
			packmol_cmd = str("{0} < packmol_{1}.input > packmol_{1}.output "
					).format(self.softwares['PACKMOL'], system)
			
			sub.call(packmol_cmd, shell=True)
			
			
			pdbtogro_cmd = "{0} editconf -f {1}.pdb -o {1}.gro -box {2} {3} {4} -c no".format(self.softwares['GROMACS_LOC'], system,
																				dimensions[0], dimensions[1], dimensions[2])
			sub.call(pdbtogro_cmd, shell=True)
			
			group_index = [str(i) for i in range(0, self.nb_index-1, 1)]
			
			make_ndx_input = """ "r {0}\\nname mono{1} su\\ndel 0\\ndel 0\\n{2}\\nname {3} System\\nq\\n" """.format(self.lipid_type,
																									self.nb_index,
																									" | ".join(group_index),
																									int(group_index[-1])+1)
			
			add_su_ndx_cmd = """printf {3} | {0} make_ndx -f {1}.gro -n {2} -o {1}.ndx""".format(self.softwares['GROMACS_LOC'], system,
																				self.system, make_ndx_input)
			
			sub.call(add_su_ndx_cmd, shell=True)
		
			self.system = system
		
		
		self.output_file = "{0}.gro".format(self.system)
		self.index_file = "{0}.ndx".format(self.system)
		
		if create_su:
			self.create_topology()
		elif generate_su:
			self.create_topology()
		else:
			self.create_topology()
			sub.call( "cp {0} .".format( '{0}/{1}'.format(path, top_file) ), shell=True)
			sub.call( "cp {0} .".format( '{0}/martini*itp'.format(path) ), shell=True)
			
		
		
		return self.system, self.output_file, self.index_file
	
	
	def make(self):
		pass
	
	def write_packmol_input(self):
		pass
	
	def create_sample(self):
		pass
	
	def add_substrate(self):
		pass
	
	def check_pbc(self):
		""" Method to check the pbc vectors at the end of .gro files
		"""
		with open('check_pbc.sh','w') as check_pbc:
			check_pbc.write( str("#!/bin/bash +x\n"
							"read -r LX LY LZ <<< $(tail -n1 {0}.gro); "
							"""if [ "$LZ" == "" ]; then """
							"""echo "LZ was too big, correcting the vectors format"; """
							"LZ=$(cut -c 9- <<< $LY); LY=$(cut -c -8 <<< $LY); "
							"sed -i '$ d' {0}.gro; "
							"""echo "$LX $LY $LZ" >> {0}.gro; """
							"fi\n").format(self.system) )
		
		sub.call('chmod u+x check_pbc.sh; ./check_pbc.sh', shell=True)
	
	def write_check_pbc(self, input_file):
		""" Method to check the pbc vectors at the end of .gro files
		"""
		check_pbc_cmd = str("\n\n### Check pbc box vectors at the end of {0} ###\n"
						"if [ -a {0} ]; then "
						"read -r LX LY LZ <<< $(tail -n1 {0}); "
						"""if [ "$LZ" == "" ]; then """
						"""echo "LZ was too big, correcting the vectors format"; """
						"LZ=$(cut -c 9- <<< $LY); LY=$(cut -c -8 <<< $LY); "
						"sed -i '$ d' {0}; "
						"""echo "$LX $LY $LZ" >> {0}; """
						"fi;fi\n\n\n").format(input_file)
		
		return check_pbc_cmd
		 
	
	def add_substrate_in_run(self, sample, md_step, previous_cmd_files, prefix_gromacs_grompp, prefix_gromacs_mdrun):
		""" Method writing bash commands to include the substrate in the middle of the runs
		"""
		self.su = self.protocol[md_step]
		
		# pdb file for su
		assert(len(self.su['SuType']) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
		
		su_type = self.su['SuType'] + (4 - len(self.su['SuType']))*' '
		
		atom_type = None
		if self.su['SuType'].startswith('SU'):
			atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
			
		elif self.su['SuType'].startswith('S'):
			atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
		
		su_pdb_content = pdb_file_list['SU']['content'].replace("TEMP", su_type)
		su_pdb_content = su_pdb_content.replace("TEM", atom_type)
		
		with open(pdb_file_list['SU']['name'],'w') as su_pdb:
			su_pdb.write(ut.RemoveUnwantedIndent(su_pdb_content))
				
		density = self.su['Density']
		thickness = self.su['Thickness']
		
		file_input = previous_cmd_files['OUTPUT']
		file_output = file_input.replace('.gro','_{0}{1}d{2}.gro'.format(self.su['SuType'], self.su['Version'], self.su['Density']))
		
		cmd = "# BEGIN ####### INSERTING THE SUBSTRATE BELOW THE SAMPLE #######"
		cmd += "### Reading the box vectors\n"
		cmd += "read -r LX LY LZ<<<$(tail -n1 {0})\n\n".format(file_input)
		
		cmd += "density={0}\n".format(density)
		cmd += "thickness={0}\n".format(float(thickness)/10.)
		
		shell = 0.3
		cmd += "shell={0}\n".format(shell)
		
		cmd += "### Computing the number of substrate particles to add\n"
		cmd += str("""nb_su=$(awk -vX=$LX -vY=$LY -vZ=$LZ -vT=$thickness -vD=$density """
					"""'BEGIN{ printf "%.0f" ,X*Y*T*D }')\n\n""")
		
		cmd += "### Computing the new Z-vector for the box (new_LZ = LZ + thickness)\n"
		cmd += """nLZ=$(awk -vZ=$LZ -vT=$thickness -vS=$shell 'BEGIN{ print Z+T+2*S }')\n\n"""
		
		# Manipulating the sample 
		cmd += "### Removing pbc to manipulate it correcly\n"
		
		cmd += "printf {0} | {1} trjconv -f {2} -s {3} -n {4} -o {5} -pbc res\n".format(repr("0\n"), prefix_gromacs_grompp, file_input, 
																						file_input.replace('.gro','.tpr'),
																						self.index_file, 'fixed_pbc.gro')
		
		cmd += self.write_check_pbc('fixed_pbc.gro')
		
		cmd += "{0} editconf -f {1} -o {2} -pbc no\n".format(prefix_gromacs_grompp,
																'fixed_pbc.gro', 'temp_for_adding_su_1.gro')
		
		cmd += self.write_check_pbc('temp_for_adding_su_1.gro')
		
		cmd += "### Creating the new box and shifting all the previous atoms (0, 0, thickness)\n"
		cmd += """translate=$(awk -vT=$thickness -vS=$shell 'BEGIN{ print T+S }')\n\n"""
		cmd += "{0} editconf -f {1} -o {2} -box $LX $LY $nLZ -c no\n".format(prefix_gromacs_grompp,
																'temp_for_adding_su_1.gro', 'temp_for_adding_su_2.gro')
		
		cmd += self.write_check_pbc('temp_for_adding_su_2.gro')
		
		cmd += "{0} editconf -f {1} -o {2} -translate 0. 0. $translate\n\n".format(prefix_gromacs_grompp,
																			'temp_for_adding_su_2.gro',
																			'temp_for_adding_su_3.gro')
		
		cmd += self.write_check_pbc('temp_for_adding_su_3.gro')
		
		# Creating the substrate
		cmd += "### Creating the substrate with box size (LX, LY, thickness)\n"
		cmd += "{0} insert-molecules -ci su_single.pdb -nmol $nb_su -o sub.gro -box $LX $LY $thickness -radius {1}\n\n".format(prefix_gromacs_grompp,
																																self.su['Radius'])
		
		self.system += "_{0}{1}d{2}".format(self.su['SuType'], self.su['Version'], self.su['Density'])
		self.create_topology()
		
		# Energy minimisation of the substrate
		cmd += "### Writing topology for EM su\n"
		cmd += "printf " + repr( """#include "{0}"\n\n""".format(self.itp_file) ) + ">> su_insertion.top\n"
		cmd += "printf " + repr("[ system ]\nSUBSTRATE\n") + ">> su_insertion.top\n"
		cmd += "printf " + repr("[ molecules ]\n%s %d\n") +""" {0} $nb_su """.format(self.su['SuType']) + ">> su_insertion.top\n\n"
		
		
		#Modifying the rcut otherwise problem with gromacs
		try:
			path_to_EMmdp = self.path_to_default+'/SU/EM_su_insertion.mdp'
			new_rcut = float(self.su['Thickness'])/20.
			copyfile(path_to_EMmdp,'EM_su_insertion.mdp')
			
			"""
			#if float(self.su['Thickness'])/20. < 1.2:
				#with open(path_to_EMmdp,'r') as input_mdp:
					#with open('EM_su_insertion.mdp','w') as output_mdp:
						#copy_original_line = True
						
						#for i, line in enumerate(input_mdp):
							#if not line.startswith(';'):
								#for rcut in ('rcoulomb','rvdw',):
									#if rcut in line and rcut != '' and (line.partition(' ')[0] == rcut or line.partition('=')[0] == rcut):
										#output_mdp.write(rcut+'			= {0}\n'.format(new_rcut-0.1))
										#copy_original_line = False
										#continue
										
								#if copy_original_line:
									#output_mdp.write(line)
								
								#copy_original_line=True
			#else:
			"""
			
			
		except IOError as e:
			print(e)
		
		
		cmd += "### Energy minimisation for the substrate\n"
		cmd += str("{0}grompp -f EM_su_insertion.mdp -po EM_su_insertion.mdp -c {1} -p su_insertion.top -maxwarn 10 -o EM_su_insertion.tpr "
								"|& tee grompp_out/grompp_EM_su_insertion.output\n\n").format(prefix_gromacs_grompp, 'sub.gro',
																					previous_cmd_files['SYSTEM'])
								
		
		
		cmd += str("{0}mdrun -deffnm EM_su_insertion -c EM_su_insertion.gro "
					" |& tee mdrun_out/mdrun_EM_su_insertion.output \n\n").format(prefix_gromacs_mdrun)
		
		
		cmd += "### Removing pbc for energy minimlized su to manipulate it correcly\n"
		cmd += "{0} editconf -f {1} -o {2} -pbc no\n\n".format(prefix_gromacs_grompp,
																'EM_su_insertion.gro','EM_su_insertion_nopbc.gro' )
		
		#Inserting the substrate into the sample
		cmd += "### Inserting the substrate in the sample\n"
		cmd += str("{0} insert-molecules -f {1} -ci EM_su_insertion_nopbc.gro -nmol 1 -o {2} "
					"-rot none -ip position.dat \n\n").format(prefix_gromacs_grompp,
																'temp_for_adding_su_3.gro',
																file_output,
																self.index_file)
		
		cmd += self.write_check_pbc(file_output)
		
		cmd += "### Creating an .ndx with su group\n"
		file_index = self.index_file.replace( '.ndx', '_{0}{1}d{2}.ndx'.format(self.su['SuType'], self.su['Version'], self.su['Density']) )
		
		#self.nb_index += 1
		group_index = [str(i) for i in range(0, self.nb_index-1, 1)]
		creating_system = "del 0\ndel 0\n"
		creating_system += "{0}\nname {1} System".format(" | ".join(group_index), int(group_index[-1])+1)
		
		cmd += "printf "
		cmd += repr( "r {0}\nname {1} su\n{2}\nq\n".format(self.su['SuType'], self.nb_index, creating_system) )
		cmd += """ | {0} make_ndx -f {1} -n {2} -o {3} \n\n""".format(prefix_gromacs_grompp,
																	file_output, self.index_file,
																	file_index)
					
		with open('position.dat','w') as pos_dat:
			pos_dat.write("0. 0. 0.\n")
		
		self.index_file = file_index
		
		cmd += "### Appending the substrate particles at the end of the topology\n"
		cmd += "printf " + repr("%s %d\n") + """ {0} $nb_su """.format(self.su['SuType'])+ """>> {0}\n\n""".format(self.system+'.top')
		cmd += "# END ####### INSERTING THE SUBSTRATE BELOW THE SAMPLE #######\n"
					
		
		return cmd, file_output, file_index, self.system
		
	def Padd_substrate_in_run(self, sample, md_step, previous_cmd_files, prefix_gromacs_grompp, prefix_gromacs_mdrun):
		pass
	
	def add_wall_in_run(self, sample, md_step, previous_cmd_files, prefix_gromacs_grompp, prefix_gromacs_mdrun):
		""" Method writing bash commands to include walls in the middle of the runs
		"""
		file_input = previous_cmd_files['OUTPUT']
		file_output = file_input.replace('.gro','_WALL.gro')
		
		self.wall = self.protocol[md_step]
		cmd = "# BEGIN ####### INCLUDING WALL #######\n"
		cmd += "### Reading the box vectors\n"
		cmd += "read -r LX LY LZ<<<$(tail -n1 {0})\n\n".format(file_input)
		
		cmd += "### Get the index of atoms having a negative z\n"
		cmd += "printf {0} | {1} select -f {2} -on neg_z.ndx\n\n".format(repr('z < 0 \n'),prefix_gromacs_grompp, file_input)
		
		cmd += "    ### Removing pbc to manipulate it correcly\n"
		cmd += "    printf {0} | {1} trjconv -f {2} -s {2} -n {3} -o {4} -pbc nojump\n".format(repr("System\n"), prefix_gromacs_grompp, file_input,
																						 self.index_file, 'fixed_pbc_wall.gro')
		
		cmd += self.write_check_pbc('fixed_pbc_wall.gro')
		
		cmd += "### Loop on the index to find the minimum z\n"
		cmd += "declare -a indexes=( $(tail -n+2 neg_z.ndx) )\n"
		cmd += "declare -x min=0.0\n"
		cmd += "for i in ${indexes[@]}\ndo\n"
		cmd += "    line=( $(cat {0} | grep $i) )\n".format('fixed_pbc_wall.gro')
		cmd += str("""    min=$(awk -vZ=${line[5]} -vM=$min """
					"""'BEGIN{ printf "%f", Z<M?Z:M  }')\ndone\n\n""")
		
		# Manipulating the sample
		cmd += "echo $min\n"
		cmd += "### Creating the new box and shifting all the previous atoms (0, 0, min)\n"
		#cmd += """if (( $(echo "$min < 0.0" | bc -l) )); then\n"""
		cmd += "    offset=0.3\n"
		cmd += """    trans=$(echo "$offset - $min" | bc -l)\n"""
		cmd += """    nnLZ=$(awk -vZ=$LZ -vM=$offset 'BEGIN{ print Z+M }')\n"""
		cmd += """ echo $nnLZ \n"""
		
		
		cmd += "    {0} editconf -f {1} -o {2} -pbc no\n".format(prefix_gromacs_grompp,
																'fixed_pbc_wall.gro', 'temp_for_adding_wall_1.gro')
		cmd += self.write_check_pbc('temp_for_adding_wall_1.gro')
		
		cmd += "    {0} editconf -f {1} -o {2} -box $LX $LY $nnLZ -c no \n".format(prefix_gromacs_grompp,
																'temp_for_adding_wall_1.gro', file_output)
		
		
		#cmd += self.write_check_pbc('temp_for_adding_wall_2.gro')
		
		#cmd += "    {0} editconf -f {1} -o {2} -translate 0. 0. $trans\n\n".format(prefix_gromacs_grompp,
																			#'temp_for_adding_wall_2.gro',
																			#file_output)
		
		
		cmd += self.write_check_pbc(file_output)
		
		
		self.system += "_WALL"
		self.create_topology()
		
		if self.su is not None:
			cmd += "### Appending the substrate particles at the end of the topology used for wall\n"
			cmd += "printf " + repr("%s %d\n") + """ {0} $nb_su """.format(self.su['SuType'])+ """>> {0}\n\n""".format(self.system+'.top')
		
		cmd += "# END ####### INCLUDING WALL #######\n"
		
		return cmd, self.system, file_output
		
		
		
		
		
	def preparing_mdp(self, md_step):
		""" Method to prepare mdp parameters
		"""
		md_step_after_init = int(md_step)-1
		
		if self.defo is not None:
			preset_name = self.defo['defoProtocol'][md_step_after_init].strip(' ')
			defo_preset_for_md_step = self.defo['presets'][preset_name]
			
		if self.su is not None:
			md_step_after_init = int(md_step)-1
			preset_name = self.su['suProtocol'][md_step_after_init].strip(' ')
			su_preset_for_md_step = self.su['presets'][preset_name]
		
		#Showing the final parameters
		print("Writing step {0} ==========================".format(self.protocol[md_step]['stepType']))
		
		#Bool to check that the run is of NVT or NPT type
		NPT_or_NVT = False
		
		#Look if the run is NPT or NVT
		if self.protocol[md_step]['stepType'].startswith('NPT') or self.protocol[md_step]['stepType'].startswith('NVT'):
			NPT_or_NVT = True
			# if THERMOSTAT is in Parameters.csv the values set here will override the values set
			# for NPT and NVT runs
			if self.thermostat is not None:
				if 'ref-t' in self.thermostat:
					self.protocol[md_step].update( {'ref-t': self.thermostat['ref-t'] } )
		
		#Energy and Temperature groups
		energy_grps = ''
		T_coupling_grps = ''
		
		# Multiply tau-t and ref-t by the number of grps
		Tau_T_coupling_grps = ''
		ref_t_T_coupling_grps = ''
		
		auto_energy_grps = True
		auto_T_coupling_grps = True
		auto_tau_T_coupling_grps = True
		auto_ref_t_T_coupling_grps = True
		
		#Creates grps for all lipid found except if tau-t and ref-t contain multiple values
		#Exception
		if 'energygrps' in self.protocol[md_step]:
			if ' ' in str(self.protocol[md_step]['energygrps']).strip() or 'System' in str(self.protocol[md_step]['energygrps']).strip():
				energy_grps += str(self.protocol[md_step]['energygrps']) + ' '
				auto_energy_grps = False
				
		if 'tc-grps' in self.protocol[md_step]:
			if ' ' in str(self.protocol[md_step]['tc-grps']).strip() or 'System' in str(self.protocol[md_step]['tc-grps']).strip():
				T_coupling_grps += str(self.protocol[md_step]['tc-grps']) + ' '
				auto_T_coupling_grps = False
				
		if 'tau-t' in self.protocol[md_step]:
			if ' ' in str(self.protocol[md_step]['tau-t']).strip():
				Tau_T_coupling_grps += str(self.protocol[md_step]['tau-t']) + ' '
				auto_tau_T_coupling_grps = False
				
			elif not auto_T_coupling_grps:
				Tau_T_coupling_grps += str(self.protocol[md_step]['tau-t']) + ' '
				auto_tau_T_coupling_grps = False
				
		if 'ref-t' in self.protocol[md_step]:
			if ' ' in str(self.protocol[md_step]['ref-t']).strip():
				ref_t_T_coupling_grps += str(self.protocol[md_step]['ref-t']) + ' '
				auto_ref_t_T_coupling_grps = False
				
			elif not auto_T_coupling_grps:
				ref_t_T_coupling_grps += str(self.protocol[md_step]['ref-t']) + ' '
				auto_ref_t_T_coupling_grps = False
		
		#Automatic generation
		if auto_energy_grps or auto_T_coupling_grps or auto_tau_T_coupling_grps or auto_ref_t_T_coupling_grps:
			if self.lipid_types is not None:
				for lipid in self.lipid_types:
					if auto_energy_grps:
						energy_grps += lipid+' '
						
					if NPT_or_NVT:
						if auto_T_coupling_grps:
							T_coupling_grps += lipid+' '
							
						if auto_tau_T_coupling_grps:
							if 'tau-t' in self.protocol[md_step] and auto_tau_T_coupling_grps:
								Tau_T_coupling_grps += str(self.protocol[md_step]['tau-t'])+' '
								
							else:
								pass
							
						if auto_ref_t_T_coupling_grps:
							if 'ref-t' in self.protocol[md_step] and auto_ref_t_T_coupling_grps:
								ref_t_T_coupling_grps += str(self.protocol[md_step]['ref-t'])+' '
								
							else:
								pass
						
					
			#Creates grps for all solvents found
			if self.solvent_types is not None:
				for sol in self.solvent_types:
					if auto_energy_grps:
						energy_grps += sol+' '
						
					if NPT_or_NVT:
						if auto_T_coupling_grps:
							T_coupling_grps += sol+' '
							
						if 'tau-t' in self.protocol[md_step] and auto_tau_T_coupling_grps:
							Tau_T_coupling_grps += str(self.protocol[md_step]['tau-t'])+' '
							
						if 'ref-t' in self.protocol[md_step] and auto_ref_t_T_coupling_grps:
							ref_t_T_coupling_grps += str(self.protocol[md_step]['ref-t'])+' '
		
		#Updates the sample dictionnary with the values created above
		if(NPT_or_NVT):
			self.protocol[md_step].update({ 'tc-grps': T_coupling_grps })
			self.protocol[md_step].update({ 'tau-t': Tau_T_coupling_grps })
			self.protocol[md_step].update({ 'ref-t': ref_t_T_coupling_grps })
			
		#Creates grps for DEFO and modifies the parameters in .mdp if in the preset
		if self.defo is not None:
			self.defo_mdp_params = '\n\n ;;; Parameters for Defo ;;; \n\n'
			energy_grps += 'defo '
			
			if 'ref-t' in defo_preset_for_md_step or 'tau-t' in defo_preset_for_md_step: 
				T_coupling_grps += 'defo '
			
			for defo_param, defo_param_value in defo_preset_for_md_step.items():
				if defo_param not in ('posres','FCX','FCY','FCZ'):
					if defo_param in self.protocol[md_step]:
						current_value = str(self.protocol[md_step][defo_param])
						self.protocol[md_step][defo_param] = current_value + ' ' + str(defo_param_value) + ' '
					else:
						self.defo_mdp_params += defo_param + '		= '+ str(defo_param_value) + ' \n'
			
			if NPT_or_NVT:
				self.protocol[md_step].update({ 'tc-grps': T_coupling_grps })
		
		#Creates grps for SU and modifies the parameters in .mdp if in the preset
		if self.su is not None:
			self.su_mdp_params = '\n\n ;;; Parameters for Su ;;; \n\n'
			#if 'freezegrps' not in su_preset_for_md_step:
			energy_grps += 'su '
			
			if 'ref-t' in su_preset_for_md_step or 'tau-t' in su_preset_for_md_step: 
				T_coupling_grps += 'su '
			
			for su_param, su_param_value in su_preset_for_md_step.items():
				if su_param not in ('posres','FCX','FCY','FCZ'):
					if su_param in self.protocol[md_step]:
						current_value = str(self.protocol[md_step][su_param])
						self.protocol[md_step][su_param] = current_value + ' ' + str(su_param_value) + ' '
					
					else:
						self.su_mdp_params += su_param + '		= '+ str(su_param_value) + ' \n'
			
			if NPT_or_NVT:
				self.protocol[md_step].update({ 'tc-grps': T_coupling_grps })
				
		#Updates the Energy groups
		self.protocol[md_step].update({ 'energygrps': energy_grps })
		
		# End of preparation =======================================================
	
	
	
	
	def writing_mdp(self, md_step, default_mdp):
		""" Method to write mdp parameters to file
		"""
		# Beginning of writing .mdp ================================================
		
		#Looks for the parameters in Default MDP file and set the values
		#as chosen in Parameters.csv
		
		md_step_after_init = int(md_step)-1
		
		posres_define = ""
		with open(self.protocol[md_step]['stepType']+'.mdp','w') as output:
		
			#Position restraining for DEFO and Su if in preset
			if self.defo is not None:
				preset_name = self.defo['defoProtocol'][md_step_after_init].strip(' ')
				if preset_name:
					print("Parameters for defos during this step will be:")
					defo_preset_for_md_step = self.defo['presets'][preset_name]
					print("	{0} : {1}\n".format(preset_name, defo_preset_for_md_step))
					
					if 'posres' in defo_preset_for_md_step:
						if defo_preset_for_md_step['posres'] == 'on' :
							posres_maccro = "DEFO_POSRES_"+ self.protocol[md_step]['stepType']
							posres_define += " -D{0}".format(posres_maccro)
							
							nb_posres = len(glob.glob('defo*posres*gen*itp'))
							
							for defo_layer in self.defo_dict:
								#Modifying Posres for defo
								posres_gen_file = open('defo_posres_{0}_gen.itp'.format(defo_layer),'r')
								posres_file_name = 'defo_posres_'+defo_layer+'_'+self.protocol[md_step]['stepType']+'.itp'
								
								posres_file_defo = open(posres_file_name,'w')
								posres_defo = posres_gen_file.read()
								
								for FC in ('FCX','FCY','FCZ'):
									if FC in defo_preset_for_md_step:
										posres_defo = posres_defo.replace(FC, str(defo_preset_for_md_step[FC]))
									else:
										posres_defo = posres_defo.replace(FC, '0')
										
								posres_file_defo.write(posres_defo)
								
								#Including the posres file in the molecule
								defo_topology_file_name = "defos_topo.itp"
								defo_topology_file = open(defo_topology_file_name,'r')
								defo_topology = defo_topology_file.read()
								defo_topology_file.close()
								
								marker = ";PLACE_FOR_{0}_POSRES".format(defo_layer)
								preprocessor = """#ifdef {0}\n	#include "{1}"\n#endif\n;PLACE_FOR_{2}_POSRES """.format(posres_maccro, 
																														posres_file_name,
																														defo_layer)
								
								defo_topology = defo_topology.replace(marker, preprocessor)
								
								with open(defo_topology_file_name,'w') as defo_topology_file:
									defo_topology_file.write(defo_topology)
								
							###else:
								###posres_gen_file = open('defo_posres_gen.itp','r')
								###posres_file_name = 'defo_posres_'+self.protocol[md_step]['stepType']+'.itp'
								
								###posres_file_defo = open(posres_file_name,'w')
								###posres_defo = posres_gen_file.read()
								
								###for FC in ('FCX','FCY','FCZ'):
									###if FC in defo_preset_for_md_step:
										###posres_defo = posres_defo.replace(FC, str(defo_preset_for_md_step[FC]))
									###else:
										###posres_defo = posres_defo.replace(FC, '0')
										
								###posres_file_defo.write(posres_defo)
								
								####Including the posres file in the molecule
								###defo_topology_file_name = "defos_topo.itp"
								
								###with open(defo_topology_file_name,'a') as defo_topology_file:
									###defo_topology_file.write(ut.RemoveUnwantedIndent("""
																			
																			####ifdef {0}
																				####include "{1}"
																			####endif
																			
																			###""").format(posres_maccro, posres_file_name))
								
			if self.su is not None:
				preset_name = self.su['suProtocol'][md_step_after_init].strip(' ')
				
				if preset_name:
					print("Parameters for su during this step will be:")
					su_preset_for_md_step = self.su['presets'][preset_name]
					print("	{0} : {1}\n".format(preset_name, su_preset_for_md_step))
					
					if 'posres' in su_preset_for_md_step:
						if su_preset_for_md_step['posres'] == 'on' :
								
							posres_maccro = "SU_POSRES_" + self.protocol[md_step]['stepType']
							posres_define += " -D{0}".format(posres_maccro)
							
							posres_gen_file = open('su_posres_gen.itp','r')
							posres_file_name = 'su_posres_' + self.protocol[md_step]['stepType']+'.itp'
							
							
							posres_file_su = open(posres_file_name,'w')
							posres_su = posres_gen_file.read()
							
							for FC in ('FCX','FCY','FCZ'):
								if FC in su_preset_for_md_step:
									posres_su = posres_su.replace(FC, str(su_preset_for_md_step[FC]))
								else:
									posres_su = posres_su.replace(FC, '0')
									
							posres_file_su.write(posres_su)
							
							self.itp_file = "martini_v2.2"
							self.itp_file +="_{0}{1}".format(self.su['SuType'], self.su['Version'])
							if self.defo is not None:
								self.itp_file +="_DEF{0}".format(self.defo['Version'])
							if self.wall is not None:
								self.itp_file +="_WALL{0}".format(self.wall['Version'])
							self.itp_file += ".itp"
							
							su_topology_file_name = self.itp_file
							su_topology_file = open(su_topology_file_name,'r')
							su_topology = su_topology_file.read()
							su_topology_file.close()
							
							
							marker = ";POSITION_FOR_SU_POSRES"
							preprocessor = """#ifdef {0}\n	#include "{1}"\n#endif\n;POSITION_FOR_SU_POSRES """.format(posres_maccro,
																													posres_file_name)
							
							su_topology = su_topology.replace(marker, preprocessor)
							
							with open(su_topology_file_name,'w') as su_topology_file:
								su_topology_file.write(su_topology)
			
			
			if self.protocol[md_step]['stepType'].startswith('EM'):
				posres_define += ' -DPW_EM'
			else:
				posres_define += ' -DPW_NOTEM'
				
			output.write("define = {0}".format(posres_define))
			
			copy_original_line = True
			for i, line in enumerate(default_mdp):
				if not line.startswith(';'):
					for key in self.protocol[md_step]:
						if key in line and key != '' and (line.partition(' ')[0] == key or line.partition('=')[0] == key):
							output.write(key+'			= '+ str(self.protocol[md_step][key])+'\n')
							copy_original_line = False
							continue
						
					if copy_original_line:
						output.write(line)
						
					copy_original_line=True
					
		with open(self.protocol[md_step]['stepType']+'.mdp','a+') as output:
			output.write('\n ;Parameters not in default file : \n\n')
			
			# Open the file again to add the parameters not in default MDP
			COMPARE = open(self.protocol[md_step]['stepType']+'.mdp','r').read()
			
			for key in self.protocol[md_step]:
				if key not in COMPARE and key != 'stepType':
					output.write(key+'			= '+ str(self.protocol[md_step][key])+'\n')
			
			defo_parameters = ""
			# Writes parameters related to Defo at the end
			if self.defo is not None:
				defo_parameters += self.defo_mdp_params
				
			defo_su_parameters = ""
			# Writes parameters related to Su at the end
			if self.su is not None:
				parameters_not_to_write = ('posres','FCX','FCY','FCZ','ref-t','tau-t')
				
				if self.defo is not None:
					
					md_step_after_init = int(md_step)-1
					preset_name_su = self.su['suProtocol'][md_step_after_init].strip(' ')
					su_preset_for_md_step = self.su['presets'][preset_name_su]
					
					preset_name_def = self.defo['defoProtocol'][md_step_after_init].strip(' ')
					defo_preset_for_md_step = self.defo['presets'][preset_name_def]
					
					
					for i, line in enumerate(defo_parameters.split('\n')):
						if not line.startswith(';'):
							for key in su_preset_for_md_step:
								if key in line:
									previous_value = line.partition('=')[2]
									defo_su_parameters += "{0}	={1} {2}\n".format(key, previous_value, str(su_preset_for_md_step[key]))
									continue

					defo_su_parameters += "\n ;;; Parameters for Defo ;;;  \n"
					
					for key in defo_preset_for_md_step:
						if key not in parameters_not_to_write and key+'	' not in defo_su_parameters:
							defo_su_parameters += "{0}			= {1}\n".format(key, defo_preset_for_md_step[key])
							
					defo_su_parameters += "\n ;;; Parameters for Su ;;; \n"
					
					for key in su_preset_for_md_step:
						if key not in parameters_not_to_write and key+'	' not in defo_su_parameters:
							defo_su_parameters += "{0}			= {1}\n".format(key,su_preset_for_md_step[key])
							
				else:
					defo_su_parameters += self.su_mdp_params
			else:
				defo_su_parameters = defo_parameters
			
			output.write("\n ;;; PARAMETERS FOR SU AND DEFO ;;;  \n")
			output.write(defo_su_parameters)
			
			if self.wall is not None:
				output.write("\n ;;; Parameters for W ;;; \n")
				preset_name = self.wall['wallProtocol'][md_step_after_init].strip(' ')
				
				if preset_name:
					print("Parameters for wall during this step will be:")
					wall_preset_for_md_step = self.wall['presets'][preset_name]

					print("	{0} : {1}\n".format(preset_name, wall_preset_for_md_step))
					
					for key in wall_preset_for_md_step:
						output.write("{0}			= {1}\n".format(key, 
													wall_preset_for_md_step[key]))
		
		default_mdp.seek(0)
		
		
		
		
		
	def create_topology(self):
		""" Method for creating the topology of defo and su
		"""
		topo_headings = ["atomtypes","nonbond_params","moleculetype","atoms"]
		# out put file name
		self.itp_file = "martini_v2.2"
		
		self.include_su_topology = {}
		
		if self.su is not None:
			include_su_topology_file = """{0}/SU/SU_{1}.itp""".format(self.path_to_default, self.su['Version'])
		
			with open(include_su_topology_file,'r') as itpsu:
				for head in topo_headings:
					for heading_and_lines in ut.group_by_heading(itpsu, head):
						lines = []
						if head is not 'moleculetype':
							if head is not 'atoms':
								lines.extend([';;;;;;; SU_{0}\n'.format(self.su['Version'])])
						lines.extend(heading_and_lines[2:])
						self.include_su_topology.update({head:''.join(lines)})
			
			posres_su = ut.RemoveUnwantedIndent("""
											
											;POSITION_FOR_SU_POSRES\n\n
											
											""")
		
			if '[ moleculetype ]' in self.include_su_topology['moleculetype'] or '[moleculetype]' in self.include_su_topology['moleculetype']:
				self.include_su_topology['moleculetype'] = self.include_su_topology['moleculetype'].replace('[ moleculetype ]', posres_su + '[ moleculetype ]' )
				
			self.include_su_topology['moleculetype'] = '[ moleculetype ]\n' + self.include_su_topology['moleculetype'] + '\n' + posres_su
			
			self.itp_file +="_{0}{1}".format(self.su['SuType'], self.su['Version'])
			
		
		
		self.include_defo_topology = {}
		if self.defo is not None:
			include_defo_topology_file = """{0}/DEFO/DEFO_{1}.itp""".format(self.path_to_default, self.defo['Version'])
			
			with open(include_defo_topology_file,'r') as itpdefo:
				for head in topo_headings:
					for heading_and_lines in ut.group_by_heading(itpdefo, head):
						lines = []
						if head is not 'moleculetype':
							if head is not 'atoms':
								lines.extend([';;;;;;; DEFO_{0}\n'.format(self.defo['Version'])])
						lines.extend(heading_and_lines[2:])
						self.include_defo_topology.update({head:''.join(lines)})
			
			self.itp_file +="_DEF{0}".format(self.defo['Version'])
			
		self.include_wall_topology = {}
		if self.wall is not None:
			include_wall_topology_file = """{0}/WALL/WALL_{1}.itp""".format(self.path_to_default, self.wall['Version'])
			
			with open(include_wall_topology_file,'r') as itpdefo:
				for head in topo_headings:
					for heading_and_lines in ut.group_by_heading(itpdefo, head):
						lines = []
						if head is not 'moleculetype':
							if head is not 'atoms':
								lines.extend([';;;;;;; WALL_{0}\n'.format(self.wall['Version'])])
						lines.extend(heading_and_lines[2:])
						self.include_wall_topology.update({head:''.join(lines)})
			
			self.itp_file +="_WALL{0}".format(self.wall['Version'])
						
		
		self.itp_file += ".itp"
		if os.path.isfile(self.itp_file):
			os.remove(self.itp_file)
		
		
		
		with open(self.itp_file, 'a+') as output:
			temp_out = ""
			temp_out2 = ""
			with open(self.path_to_default + "/martini_v2.2.itp",'r') as def_martini:
				
				if self.su is not None:
					for line in def_martini:
						if 'nonbond_params' in line:
							temp_out += self.include_su_topology['atomtypes']
							temp_out += '\n[ nonbond_params ]\n'
							
						elif 'PLACE_FOR_SU' in line:
							temp_out += self.include_su_topology['nonbond_params']
							temp_out += '\n'
							
						else:
							temp_out += line
				else:
					temp_out = def_martini.read()
				
				
				if self.defo is not None:
					for line in temp_out.splitlines():
						if 'nonbond_params' in line:
							temp_out2 += self.include_defo_topology['atomtypes']
							temp_out2 += '\n[ nonbond_params ]\n'
							
						elif 'PLACE_FOR_DEFO' in line:
							temp_out2 += self.include_defo_topology['nonbond_params']
							if self.su is not None:
								atom_type = None
								if self.su['SuType'].startswith('SU'):
									atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
									
								elif self.su['SuType'].startswith('S'):
									atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
									
								temp_out2 += 'DEF  {0}  1     0.0 0.0 ;\n'.format(atom_type)
							temp_out2 += '\n'
							
						else:
							temp_out2 += line
							temp_out2 += '\n'
				else:
					temp_out2 = temp_out
							
				if self.wall is not None:
					for line in temp_out2.splitlines():
						if 'nonbond_params' in line:
							output.write(self.include_wall_topology['atomtypes'])
							output.write('\n[ nonbond_params ]\n')
							
						elif 'PLACE_FOR_WALL' in line:
							output.write(self.include_wall_topology['nonbond_params'])
							output.write('\n')
							
						else:
							output.write(line)
							output.write('\n')
				
				else:
					output.write(temp_out2)
				
				if self.su is not None:
					output.write("\n;;;;;;; SU_{0}\n".format(self.su['Version']))
					output.write(self.include_su_topology['moleculetype']+'\n')
					
					
		#=======================================================================
		# the topology file
		#=======================================================================
		topology = """#include "{0}"\n""".format(self.itp_file)
		
		if self.__class__.__name__ == "Membrane":
			topology += """#include "martini_v2.0_lipids.itp"\n"""
			sub.call("""cp {0}/martini_v2.0_lipids.itp  .""".format(self.path_to_default), shell= True)
			
		elif self.__class__.__name__  == "Solvent":
			topology += """#include "martini_v2.0_solvents.itp"\n"""
			sub.call("cp {0}/martini_v2.0_solvents.itp ./".format(self.path_to_default), shell= True)
			
		#Add the topology for the defos (bilayer and monolayer)
		if self.defo is not None:
			topology += ut.RemoveUnwantedIndent("""
				
				#include "defos_topo.itp"
				
				""")
			
		
		#Copy the topology files for martini forcefield
		
		if self.su is not None:
			sub.call("""cp {0}/SU/su_posres_gen.itp  .""".format(self.path_to_default), shell= True)
			
		with open(self.system+'.top','a') as topo_file:
			topo_file.write(topology)
			topo_file.write("\n[ system ]\n")

			topology = "{0}\n".format(self.system)
			topo_file.write(topology)

			topo_file.write("\n[ molecules ]\n")
			
			if self.lipid_types:
				for lipid in self.lipid_types:
					print(self.sample_molecules)
					lipid_number = self.sample_molecules[lipid]
					
					if self.mono is not None:
						lipid_number += float(self.mono['NbLipidsM'])
						
					topo_file.write( "{0} {1}\n".format(lipid, lipid_number) )
				
			
			for sol in self.solvent_types:
				topo_file.write("{0} {1}\n".format(sol, self.sample_molecules[sol]))
				
			if self.defo is not None:
				for defo in self.defo_dict:
					topo_file.write("{0} 1\n".format(defo))
				
			
			if self.su is not None:
				if self.nb_su is not None:
					topo_file.write("{0} {1}\n".format(self.su['SuType'], self.nb_su))
				else:
					for su in self.substrate_list:
						if su in self.sample_molecules:
							topo_file.write("{0} {1}\n".format(su, self.sample_molecules[su]))
	
	
	def pass_outputs(self):
		return self.system, self.output_file, self.index_file












class Membrane(BaseProject):
	""" Class for building membranes.
	This includes bilayer and trilayer with/without defos, substrate
	"""
	def __init__(self, sample, softwares, path_to_default):
		super(Membrane, self).__init__(sample, softwares, path_to_default)
		
		
		
		
		
	
	def make(self):
		self.nb_lipid_monolayer = {}
		self.nb_solvent_per_monolayer = {}
		self.lipid_type = None
		self.solvent_type = None
		
		self.finding_molecules()
		
		self.lipid_types = [self.lipid_type]
		self.solvent_types = [self.solvent_type]
		
		self.nb_sol_bottom = self.nb_solvent_per_monolayer[self.solvent_type]
		self.nb_sol_top = self.nb_solvent_per_monolayer[self.solvent_type]
		
		
		
		# Setting the total monolayer thickness and
		# the volume constraining the tails and heads
		self.tmt = None
		self.dz = None
		
		if self.lipid_type == 'DSPC':
			self.tmt = 30.0
			self.dz = 10.0
		elif self.lipid_type in ('DPPC','DLPC'):
			self.tmt = 30.0
			self.dz = 10.0
		else:
			assert(False), "Your Parameters.csv contains a lipid not yet defined in Prepare.py"
			
		# Setting the bilayer position
		self.height = None
		if 'HEIGHT' in self.geometry:
			self.height = self.geometry['HEIGHT']
		else:
			self.height = self.dimensions['LZ']/2.0
		
		self.m1_head_min, self.m1_head_max, self.m1_tail_min, self.m1_tail_max = (None,)*4
		self.m2_head_min, self.m2_head_max, self.m2_tail_min, self.m2_tail_max = (None,)*4
		self.m3_head_min, self.m3_head_max, self.m3_tail_min, self.m3_tail_max = (None,)*4
		
		self.update_tails_and_heads_positions()
		
		if self.type.startswith('BI'):
			self.adjust_bilayer_position()
			self.update_tails_and_heads_positions()
		
		self.system = "{0}_{1}{2}_{3}{4}".format(self.type, self.sample_molecules[self.lipid_type],
																	self.lipid_type, self.sample_molecules[self.solvent_type], 
																	self.solvent_type)
		#Adding the monolayer
		self.nb_sol_vapor = None
		if self.mono is not None: #self.type == 'TRILAYER':
			#if self.lz_vaccum is not None:
				#self.nb_sol_vapor = int( 0.007 * self.dimensions['LX'] * self.dimensions['LY'] * (self.lz_vaccum - self.dimensions['LZ'] - self.tmt) /1000)
			
			self.add_monolayer()
			
			self.system += "_{0}M{1}".format(self.sample_molecules[self.lipid_type], self.mono['NbLipidsM'])
			
		
		#Adding the defo
		self.nb_defo = 0
		self.defo_dict = {}
		
		if self.defo is not None:
			if self.type.startswith('BI'):
				self.nb_defo, self.defo_dict = self.add_defo_bilayer()
			elif self.type.startswith('MONO'):
				self.nb_defo, self.defo_dict = self.add_defo_monolayer()
			
			self.system += "_{0}DEF{1}".format(self.nb_defo, self.defo['Version'])
		
		#Adding the substrate
		self.nb_su = None
		
		if self.su is not None:
			self.add_substrate()
			self.system += "_{0}{1}{2}".format(self.nb_su, self.su['SuType'], self.su['Version'])
			
		if self.wall is not None:
			self.system += "_WALL"
		
		# Writing packmol_input while keeping track of the molecules in the sample
		self.nb_index = 2 # start at 2 because the index has already [system] and [other]
		
		if self.type.startswith('BI'):
			self.write_packmol_bilayer()
		elif self.type.startswith('MONO'):
			self.write_packmol_monolayer()
		
		if self.lz_vaccum is not None:
			self.dimensions['LZ'] = self.lz_vaccum
		
		# Calling packmol, vmd and gmx make_index
		if self.type.startswith('BI'):
			self.create_bilayer(self.nb_index)
		elif self.type.startswith('MONO'):
			self.create_monolayer(self.nb_index)
		
		# creating the topology for su and defo if in sample
		self.create_topology()
	
	
	def adjust_bilayer_position(self):
		""" Method checking that the bilayer does not overlap with box borders in Z direction
		"""
		lbz = self.dimensions['LZ']/2.0
		
		#Checking that the Bilayer lipids do not overlap with the Monolayer lipids and the substrate
		if (self.height - self.tmt) < 0.0:
			self.height -= self.height - self.tmt
			
		if (self.height + self.tmt) > self.dimensions['LZ']:
			self.height += self.dimensions['LZ'] - self.height - self.tmt
			
		#Checking that the Solvent is not too dense
		if self.height < lbz: #Too much Solvent below
			volume_lbz = self.dimensions['LX']*self.dimensions['LY']*(lbz - self.tmt)
			density_sol = self.nb_solvent_per_monolayer[self.solvent_type]/volume_lbz
			
			#volume of solvent to remove
			vol_to_remove = lbz - self.height
			vol_to_remove *= self.dimensions['LX']*self.dimensions['LY']
			
			#Nb of particles to remove
			nb_particles_to_relocate = int(density_sol * vol_to_remove)
			
			self.nb_sol_bottom -= nb_particles_to_relocate
			self.nb_sol_top += nb_particles_to_relocate
			
		if self.height > lbz: #Too much Solvent above
			volume_lbz = self.dimensions['LX']*self.dimensions['LY']
			volume_lbz *= self.dimensions['LZ'] - lbz - self.tmt
			
			density_sol = self.nb_solvent_per_monolayer[self.solvent_type]/volume_lbz
			
			vol_to_remove = self.height - lbz
			vol_to_remove *= self.dimensions['LX']*self.dimensions['LY']
			nb_particles_to_relocate = int(density_sol * vol_to_remove)
			
			self.nb_sol_bottom += nb_particles_to_relocate
			self.nb_sol_top -= nb_particles_to_relocate
	
	
	
	def update_tails_and_heads_positions(self):
		""" Method updating the positions of heads and tails
		"""
		
		if self.type.startswith('BI'):
			self.m1_head_min = self.height - self.tmt
			self.m1_head_max = self.height - self.tmt + self.dz
			self.m1_tail_min = self.height - self.dz
			self.m1_tail_max = self.height

			self.m2_head_min = self.height + self.tmt - self.dz
			self.m2_head_max = self.height + self.tmt
			self.m2_tail_min = self.height
			self.m2_tail_max = self.height + self.dz
			
		elif self.type.startswith('MONO'):
			self.m1_head_min = self.height - self.sol_thickness/2. - self.dz
			self.m1_head_max = self.height - self.sol_thickness/2.
			self.m1_tail_min = self.height - self.sol_thickness/2. - self.tmt
			self.m1_tail_max = self.height - self.sol_thickness/2. - self.tmt + self.dz

			self.m2_head_min = self.height + self.sol_thickness/2.
			self.m2_head_max = self.height + self.sol_thickness/2. + self.dz
			self.m2_tail_min = self.height + self.sol_thickness/2. + self.tmt - self.dz
			self.m2_tail_max = self.height + self.sol_thickness/2. + self.tmt
		
		if self.su is None and self.wall is not None:
			self.m1_head_min += self.bottom
			self.m1_head_max += self.bottom
			self.m1_tail_min += self.bottom
			self.m1_tail_max += self.bottom
			
			self.m2_head_min += self.bottom
			self.m2_head_max += self.bottom
			self.m2_tail_min += self.bottom
			self.m2_tail_max += self.bottom








	
	
	
	
	def finding_molecules(self):
		""" Method to find the molecules and associate their number from csv file.
		"""
		for lipid in self.lipid_list:
			if lipid in self.sample_molecules:
				self.nb_lipid_monolayer.update( {lipid : int(self.sample_molecules[lipid]/2)} )
				self.lipid_type = lipid
				
		for sol in self.solvent_list:
			if sol in self.sample_molecules:
				self.nb_solvent_per_monolayer.update( {sol:int(self.sample_molecules[sol]/2)} )
				self.solvent_type = sol
	
	
	
	
	def add_monolayer(self):
		""" Method adding the positions of heads and tails for monolayer
		"""
		if self.mono['LzM'] in ('top','default'):
			self.mono['LzM'] = self.dimensions['LZ']
		
		self.m3_head_min = self.mono['LzM'] 
		self.m3_head_max = self.mono['LzM'] + self.dz
		self.m3_tail_min = self.mono['LzM'] + self.tmt - self.dz
		self.m3_tail_max = self.mono['LzM'] + self.tmt
		
		if self.su is None and self.wall is not None:
			self.m3_head_min += self.bottom
			self.m3_head_max += self.bottom
			self.m3_tail_min += self.bottom
			self.m3_tail_max += self.bottom
		#add the monolayer at the top of box or chose the position and move the solvent ?
		
		
		
		
	
	
	
	def add_substrate(self):
		""" Method to add the substrate at the bottom of box
		"""
		thickness = float(self.su['Thickness'])
		# Shifting all the coordinates in Z
		self.height += thickness
		
		self.m1_head_min += thickness
		self.m1_head_max += thickness
		self.m1_tail_min += thickness
		self.m1_tail_max += thickness
		
		self.m2_head_min += thickness
		self.m2_head_max += thickness
		self.m2_tail_min += thickness
		self.m2_tail_max += thickness
		
		self.mono['LzM'] += thickness
		
		self.bottom_z += thickness
		
		if self.mono is not None:
			self.m3_head_min += thickness 
			self.m3_head_max += thickness
			self.m3_tail_min += thickness
			self.m3_tail_max += thickness
			
		
		self.nb_su = int( float(self.su['Density']) * self.dimensions['LX'] * self.dimensions['LY'] * thickness / 1000 )
		
		#Pushing the box limits
		self.dimensions['LZ'] += thickness
	
	
	
	
	def add_defo_bilayer(self):
		"""Method to add defo in the Membrane
		"""
		if self.defo['Height'] == 'box':
			#Set the Defo height from the substrate to the mono layer
			length_defo_bi = self.dimensions['LZ']
			
		elif self.defo['Height'] == 'bilayer':
			#Set the Defo height from the substrate to the top of the bilayer
			length_defo_bi = self.height + self.tmt + self.dz/4.0
			
		elif self.defo['Height'] == 'mono':
			#Set the Defo height from the substrate to the top of the monolayer
			length_defo_bi = self.dimensions['LZ'] + self.tmt + self.dz/4.0
			
		elif self.defo['Height'] == 'follow':
			#Set the Defo height from the bottom of the bilayer to its top
			length_defo_bi = 2*self.tmt + self.dz/2.0
			
		
		defo_per_layer = int(self.defo['DpL']) + 1
		
		#Set a dictionnary for automation in defos creation
		defo_dict = {"DEFB": None}
		
		nb_layers_bi = int(length_defo_bi/float(self.defo['DzDefo']))
		defo_dict['DEFB'] = {'length': length_defo_bi, 'nb layers': nb_layers_bi,'nb defo':0,'defo outside': [], 
										'defo inside': [], 'defo total': [], 'xyz file':'defo_bilayer.xyz',
										'format defo':"format_DEFB.vmd",'chain':'X', 'topo':'defo_bilayer_topo.itp'}
		
		if self.mono is not None:
			defo_dict["DEFM"] = {'length': 0.0, 'nb layers': 0,'nb defo':0,'defo outside': [], 
										'defo inside': [], 'defo total': [], 'xyz file':'defo_mono.xyz',
										'format defo':"format_DEFM.vmd",'chain':'Y', 'topo':'defo_monolayer_topo.itp'}
			
			L_defoMono = self.tmt + self.dz/4.0
			nb_layers_mono = int(L_defoMono/float(self.defo['DzDefo']))
			defo_dict['DEFM']['length'] = L_defoMono
			defo_dict['DEFM']['nb layers'] = int(nb_layers_bi/2)
			
		
		#Total number of defo
		defo_dict['DEFB']['nb defo'] = nb_layers_bi*defo_per_layer
		
		if self.mono is not None:
			defo_dict['DEFM']['nb defo'] = nb_layers_mono*defo_per_layer
		
		#Radius for the hole
		radius_defo = float(self.defo['Radius'])
		
		#Topology for defo
		topo_defo = """;;;;;; DEFOS FOR BILAYER\n"""
		
		
		#Total number of defo
		nb_defo = 0
		
		#Creating and Writing the defos configuration
		for defo_layer in defo_dict:
			
			posres_defo = """;;;;;;POSRES FOR DEFO_{0}\n""".format(defo_layer)
			posres_defo += ut.RemoveUnwantedIndent("""
						[position_restraints]
						;ai   funct   fcx     fcy     fcz
						
						""")
			
			xyz = defo_dict[defo_layer]['xyz file']
			
			if os.path.exists(xyz):
				os.remove(xyz)
			
			xyz_out = open(xyz,'a')
			
			#Writing the number of defo grains for each layer
			xyz_out.write( "{0}\n\n".format(defo_dict[defo_layer]['nb defo']) )
		
			defo_numb = 0

			#Creating the xyz file for the bilayer
			if defo_per_layer != 1:
				for i in range(0, defo_dict[defo_layer]['nb layers']):
					z_current_defo = float(self.defo['DzDefo'])*i
					
					for j in np.arange(0., 360., 360./(defo_per_layer-1)):
						angle = math.radians(j)
						xyz_out.write(ut.RemoveUnwantedIndent(
							"""
							DEF  {0}  {1}  {2} \n
							""".format(radius_defo*math.cos(angle), radius_defo*math.sin(angle), z_current_defo)
							))
						defo_numb += 1
						defo_dict['DEFB']['defo outside'].append(defo_numb)
					xyz_out.write(ut.RemoveUnwantedIndent(
							"""
							DEF  0.0 0.0 {0} \n
							""".format(z_current_defo)
							))
					defo_numb += 1
					defo_dict['DEFB']['defo inside'].append(defo_numb)
					
				
				defo_dict[defo_layer]['defo total'].extend(defo_dict[defo_layer]['defo inside'])
				defo_dict[defo_layer]['defo total'].extend(defo_dict[defo_layer]['defo outside'])
				defo_dict[defo_layer]['defo total'].sort()
				
			else:
				for i in range(0, defo_dict[defo_layer]['nb layers']):
					z_current_defo = float(self.defo['DzDefo'])*i
					xyz_out.write(ut.RemoveUnwantedIndent(
							"""
							DEF  0.0 0.0 {0} \n
							""".format(z_current_defo)
							))
					defo_numb += 1
					defo_dict[defo_layer]['defo total'].append(defo_numb)
			
			xyz_out.close()
			
			
		
		
			#Formating the defos
			format_defo = open(defo_dict[defo_layer]['format defo'],"w")
			format_defo.write(ut.RemoveUnwantedIndent(
				"""
				mol load xyz {0}
				set all [atomselect top "all"]
					
				$all set resname {1}
				$all set name DEF
				$all set type DEF
				$all set chain {2}
					
				package require pbctools
				pbc set {{0.5 0.5 {3}}}
					
				$all writepdb {4}
				unset all
					
				exit
				""".format(xyz, defo_layer, defo_dict[defo_layer]['chain'] , defo_dict[defo_layer]['length'], xyz.replace('xyz','pdb'))
				))
			format_defo.close()
			
			cmd = str("""{0} -dispdev text -e {1} > {2}""").format(self.softwares['VMD'], defo_dict[defo_layer]['format defo'],
																defo_dict[defo_layer]['format defo'].replace('vmd','log'))
			sub.call(cmd, shell=True)
			
								
							
			topo_defo += ut.RemoveUnwantedIndent("""
							[moleculetype]
							;molname    nrexcl
							{0} 1
							
							[atoms]
							;id     type     resnr    residu  atom    cgnr    charge
							
							""".format(defo_layer))
			
			for defnb in defo_dict[defo_layer]['defo total']:
				topo_defo += " {0}     DEF     1    {1}    DEF     {0}     0\n".format(defnb, defo_layer)
				posres_defo += " {0}    1     FCX   FCY   FCZ\n".format(defnb)
				
			#Checking for contrains
			if 'Interactions' in self.defo:
				"""
				The length a will depend on the circumradius : a = 2*R*math.sin(pi/(defo_per_layer-1))
				The angle theta will depend on the number of defo_per_layer theta = 360./(defo_per_layer-1)
				example with 4 defo per layer and one layer:
				      2          The bonds will be |  And the angles:
				    / | \          1-2   5-1       |     1-5-2
				   3--5--1         2-3   5-2       |     2-5-3
					\ | /          3-4   5-3       |     3-5-4
				      4            4-1   5-4       |     4-5-1
				"""
				if self.defo['Interactions'] == 'bond':
					topo_defo += self.bonds_topology(defo_dict[defo_layer])

				elif self.defo['Interactions'] == 'bond&angles':
					topo_defo += self.bonds_angles_topology( defo_dict[defo_layer])
					
				else:
					pass
				
			topo_defo += """\n;PLACE_FOR_{0}_POSRES\n\n\n""".format(defo_layer)
			
			posres_file_defo = open('defo_posres_{0}_gen.itp'.format(defo_layer),'w')
			posres_file_defo.write(posres_defo)
			posres_file_defo.close()
			
			nb_defo += len(defo_dict[defo_layer]['defo total'])
			
		
		
		topo_file_defo = open('defos_topo.itp','w')
		topo_file_defo.write(topo_defo)
		topo_file_defo.close()
		
		
		
		return nb_defo, defo_dict
		
	def add_defo_monolayer(self):
		"""Method to add defo in the Membrane
		"""
		if self.defo['Height'] == 'box':
			#Set the Defo height from the substrate to the mono layer
			length_defo_mono	= self.dimensions['LZ']
		
		elif self.defo['Height'] == 'mono':
			#Set the Defo height from the substrate to the top of the monolayer
			length_defo_mono	= 2*self.tmt + self.sol_thickness + 2*self.dz
			
		elif self.defo['Height'] == 'follow':
			#Set the Defo height from the bottom of the bilayer to its top
			length_defo_mono	= self.tmt + self.dz
			
		
		defo_per_layer	= int(self.defo['DpL']) + 1
		
		#Set a dictionnary for automation in defos creation
		defo_dict		= {"DEFT": None, "DEFB": None}
		
		nb_layers_mono	= int(length_defo_mono/float(self.defo['DzDefo']))
		
		#Total number of defo
		nb_defo_per_mono = nb_layers_mono * defo_per_layer
		
		defo_dict['DEFB'] = {'length': length_defo_mono, 'nb layers': nb_layers_mono,'nb defo':nb_defo_per_mono,'defo outside': [], 
										'defo inside': [], 'defo total': [], 'xyz file':'defo_mono_bot.xyz',
										'format defo':"format_DEFB.vmd",'chain':'X', 'topo':'defo_mono_topo.itp'}
		
		defo_dict['DEFT'] = {'length': length_defo_mono, 'nb layers': nb_layers_mono,'nb defo':nb_defo_per_mono,'defo outside': [], 
										'defo inside': [], 'defo total': [], 'xyz file':'defo_mono_top.xyz',
										'format defo':"format_DEFT.vmd",'chain':'Y', 'topo':'defo_mono_topo.itp'}
			
		
		
		if self.mono is not None:
			defo_dict['DEFM']['nb defo'] = nb_layers_mono*defo_per_layer
		
		#Radius for the hole
		radius_defo = float(self.defo['Radius'])
		
		#Topology for defo
		topo_defo = """;;;;;; DEFOS FOR BILAYER\n"""
		
		
		#Total number of defo
		nb_defo = 0
		
		#Creating and Writing the defos configuration
		for defo_layer in defo_dict:
			
			posres_defo = """;;;;;;POSRES FOR DEFO_{0}\n""".format(defo_layer)
			posres_defo += ut.RemoveUnwantedIndent("""
						[position_restraints]
						;ai   funct   fcx     fcy     fcz
						
						""")
			
			xyz = defo_dict[defo_layer]['xyz file']
			
			if os.path.exists(xyz):
				os.remove(xyz)
			
			xyz_out = open(xyz,'a')
			
			#Writing the number of defo grains for each layer
			xyz_out.write( "{0}\n\n".format(defo_dict[defo_layer]['nb defo']) )
		
			defo_numb = 0

			#Creating the xyz file for the bilayer
			if defo_per_layer != 1:
				for i in range(0, defo_dict[defo_layer]['nb layers']):
					z_current_defo = float(self.defo['DzDefo'])*i
					
					for j in np.arange(0., 360., 360./(defo_per_layer-1)):
						angle = math.radians(j)
						xyz_out.write(ut.RemoveUnwantedIndent(
							"""
							DEF  {0}  {1}  {2} \n
							""".format(radius_defo*math.cos(angle), radius_defo*math.sin(angle), z_current_defo)
							))
						defo_numb += 1
						defo_dict['DEFB']['defo outside'].append(defo_numb)
					xyz_out.write(ut.RemoveUnwantedIndent(
							"""
							DEF  0.0 0.0 {0} \n
							""".format(z_current_defo)
							))
					defo_numb += 1
					defo_dict['DEFB']['defo inside'].append(defo_numb)
					
				
				defo_dict[defo_layer]['defo total'].extend(defo_dict[defo_layer]['defo inside'])
				defo_dict[defo_layer]['defo total'].extend(defo_dict[defo_layer]['defo outside'])
				defo_dict[defo_layer]['defo total'].sort()
				
			else:
				for i in range(0, defo_dict[defo_layer]['nb layers']):
					z_current_defo = float(self.defo['DzDefo'])*i
					xyz_out.write(ut.RemoveUnwantedIndent(
							"""
							DEF  0.0 0.0 {0} \n
							""".format(z_current_defo)
							))
					defo_numb += 1
					defo_dict[defo_layer]['defo total'].append(defo_numb)
			
			xyz_out.close()
			
			
		
		
			#Formating the defos
			format_defo = open(defo_dict[defo_layer]['format defo'],"w")
			format_defo.write(ut.RemoveUnwantedIndent(
				"""
				mol load xyz {0}
				set all [atomselect top "all"]
					
				$all set resname {1}
				$all set name DEF
				$all set type DEF
				$all set chain {2}
					
				package require pbctools
				pbc set {{0.5 0.5 {3}}}
					
				$all writepdb {4}
				unset all
					
				exit
				""".format(xyz, defo_layer, defo_dict[defo_layer]['chain'] , defo_dict[defo_layer]['length'], xyz.replace('xyz','pdb'))
				))
			format_defo.close()
			
			cmd = str("""{0} -dispdev text -e {1} > {2}""").format(self.softwares['VMD'], defo_dict[defo_layer]['format defo'],
																defo_dict[defo_layer]['format defo'].replace('vmd','log'))
			sub.call(cmd, shell=True)
			
								
							
			topo_defo += ut.RemoveUnwantedIndent("""
							[moleculetype]
							;molname    nrexcl
							{0} 1
							
							[atoms]
							;id     type     resnr    residu  atom    cgnr    charge
							
							""".format(defo_layer))
			
			for defnb in defo_dict[defo_layer]['defo total']:
				topo_defo += " {0}     DEF     1    {1}    DEF     {0}     0\n".format(defnb, defo_layer)
				posres_defo += " {0}    1     FCX   FCY   FCZ\n".format(defnb)
				
			#Checking for contrains
			if 'Interactions' in self.defo:
				"""
				The length a will depend on the circumradius : a = 2*R*math.sin(pi/(defo_per_layer-1))
				The angle theta will depend on the number of defo_per_layer theta = 360./(defo_per_layer-1)
				example with 4 defo per layer and one layer:
				      2          The bonds will be |  And the angles:
				    / | \          1-2   5-1       |     1-5-2
				   3--5--1         2-3   5-2       |     2-5-3
					\ | /          3-4   5-3       |     3-5-4
				      4            4-1   5-4       |     4-5-1
				"""
				if self.defo['Interactions'] == 'bond':
					topo_defo += self.bonds_topology(defo_dict[defo_layer])

				elif self.defo['Interactions'] == 'bond&angles':
					topo_defo += self.bonds_angles_topology( defo_dict[defo_layer])
					
				else:
					pass
				
			topo_defo += """\n;PLACE_FOR_{0}_POSRES\n\n\n""".format(defo_layer)
			
			posres_file_defo = open('defo_posres_{0}_gen.itp'.format(defo_layer),'w')
			posres_file_defo.write(posres_defo)
			posres_file_defo.close()
			
			nb_defo += len(defo_dict[defo_layer]['defo total'])
			
		
		
		topo_file_defo = open('defos_topo.itp','w')
		topo_file_defo.write(topo_defo)
		topo_file_defo.close()
		
		
		
		return nb_defo, defo_dict
		
	
	
	
	
	
	def bonds_topology(self, defo_dict_layer):
		"""Method for writing bonds topology for the defo"""
		
		defo_per_layer = int(self.defo['DpL']) + 1
		
		topo_defo = ut.RemoveUnwantedIndent("""
					
					[bonds]
					; i j   funct   length  force.c.
					
					""")
		#Set the bonds in a layer
		# 0.2 prefactor as martini needs [nm]
		
		if int(self.defo['DpL']) != 0:
			topo_defo +="\n;Normal bonds  DEF1-DEF2\n"
			lenght_bond_plane = round(0.2* float(self.defo['Radius']) * math.sin( math.pi/float(defo_per_layer-1) ),2)
			
			for i in range(0, defo_dict_layer['nb layers']):
				j = i*(defo_per_layer - 1)
				k = j + defo_per_layer -2
				for DEF in defo_dict_layer['defo outside'][j:k+1]:
					nextDEF = DEF + 1
					if DEF == defo_dict_layer['defo outside'][k]:
						nextDEF = defo_dict_layer['defo outside'][j]
					topo_defo += """ {0} {1}   1       {2}    {3}\n""".format(DEF, nextDEF,
														lenght_bond_plane, self.defo['FbondP'])
			
			#Set the bonds inside the layer
			topo_defo +="\n;Normal bonds  DEFC-DEF\n"
			lenght_bond_plane_in = round(0.05*float(self.defo['Radius']), 2)
			
			for i,DEFC in zip(range(0, defo_dict_layer['nb layers']), defo_dict_layer['defo inside']):
				j = i*(defo_per_layer - 1)
				k = j + defo_per_layer - 2
				for DEF in defo_dict_layer['defo outside'][j:k+1]:
					topo_defo += """ {0} {1}   1       {2}    {3}\n""".format(DEFC, DEF,
																			lenght_bond_plane_in,
																	self.defo['FbondP'])
				
		#Set the bonds along normal to the layer
		#Distance in z [nm] (thus the 0.1 prefactor)
		topo_defo +="\n;Normal bonds  DEF-DEFAbove\n"
		lenght_bond_normal = round(0.1*float(self.defo['DzDefo']), 2)
		
		for DEF in defo_dict_layer['defo total']:
			DEFAbove = DEF + defo_per_layer
			if DEFAbove not in defo_dict_layer['defo total']:
				break
			topo_defo += """ {0} {1}   1       {2}    {3}\n""".format(DEF,DEFAbove,
															lenght_bond_normal, self.defo['FbondN'])
		return topo_defo
	
	
	
	def bonds_angles_topology(self, topo_defo, defo_dict_layer):
		"""Method for writing bonds and angles topology for the defo"""
		topo_defo = ut.RemoveUnwantedIndent("""
							
						[angles]
						; i j k  angle   force
						
						""")
		
		defo_per_layer = int(self.defo['DpL']) + 1
		
		if int(self.defo['DpL']) != 0:
			#Set the angles in plane
			topo_defo +="\n;In-plane angles  DEF-Center-nextDEF \n"
			angle = 360./(defo_per_layer-1)
			
			for i, DEFC in zip(range(0, defo_dict_layer['nb layers']), defo_dict_layer['defo inside']):
				j = i*(defo_per_layer-1)
				k = j + defo_per_layer-2
				for DEF in defo_dict_layer['defo outside'][j:k+1]:
					nextDEF = DEF + 1
					if DEF == defo_dict_layer['defo oustide'][k]:
						nextDEF = defo_dict_layer['defo oustide'][j]
					topo_defo += """ {0} {1} {2}   {3}   {4}  {5}\n""".format(DEF, DEFC, nextDEF,
																			self.defo['FtypeAngle'],
																			angle, 
																			self.defo['FangleP'])
			
			topo_defo +="\n;In-plane angles Center-DEF-nextDEF\n"
			outAngle = 0.5*(180. - angle)
			
			for i,DEFC in zip(range(0, defo_dict_layer['nb layers']), defo_dict_layer['defo inside']):
				j = i*(defo_per_layer-1)
				k = j + defo_per_layer-2
				for DEF in defo_dict_layer['defo oustide'][j:k+1]:
					nextDEF = DEF + 1
					if DEF == defo_dict_layer['defo oustide'][k]:
						nextDEF = defo_dict_layer['defo oustide'][j]
					topo_defo += """ {1} {0} {2}   {3}   {4}  {5}\n""".format(DEF, DEFC, nextDEF,
																				self.defo['FtypeAngle'],
																				outAngle,
																				self.defo['FangleP'])
			
			topo_defo +="\n;In-plane angles Center-nextDEF-DEF\n"
			for i,DEFC in zip(range(0, defo_dict_layer['nb layers']), defo_dict_layer['defo inside']):
				j = i*(defo_per_layer-1)
				k = j + defo_per_layer-2
				for DEF in defo_dict_layer['defo oustide'][j:k+1]:
					nextDEF = DEF + 1
					if DEF == defo_dict_layer['defo oustide'][k]:
						nextDEF = defo_dict_layer['defo oustide'][j]
					topo_defo += """ {1} {0} {2}   {3}   {4}  {5}\n""".format(DEF, DEFC, nextDEF,
																			self.defo['FtypeAngle'],
																			outAngle, 
																			self.defo['FangleP'])
			
			#Set the angles out-of-plane
			topo_defo +="\n;Out-of-plane angles Center-DEF-DEFAbove\n"
			for i,DEFC in zip(range(0, defo_dict_layer['nb layers']), defo_dict_layer['defo oustide']):
				j = i*(defo_per_layer-1)
				k = j + defo_per_layer-2
				
				for DEF in defo_dict_layer['defo oustide'][j:k+1]:
					DEFAbove = DEF + defo_per_layer
					
					if DEFAbove not in defo_dict_layer['defo total']:
						break
					
					topo_defo += """ {1} {0} {2}   {3}   90.  {4}\n""".format(DEF, DEFC, DEFAbove,
																				self.defo['FtypeAngle'],
																				self.defo['FangleN'])
		
			topo_defo +="\n;Out-of-plane angles DEF-Center-DEFCAbove\n"
			for i, DEFC in zip(range(0, defo_dict_layer['nb layers']), defo_dict_layer['defo total']):
				j = i*(defo_per_layer-1)
				k = j + defo_per_layer-2
				DEFCAbove = DEFC + defo_per_layer
				
				if DEFCAbove not in defo_dict_layer['defo inside']:
					break
				
				for DEF in defo_dict_layer['defo outside'][j:k+1]:
					topo_defo += """ {0} {1} {2}   {3}   90.  {4}\n""".format(DEF, DEFC, DEFCAbove,
																			self.defo['FtypeAngle'],
																				self.defo['FangleN'])
		
		return topo_defo
	
	
	
	
	def write_packmol_monolayer(self):
		""" Method to write packmol input for membrane creation"""
		shell = 3.0
		self.packmol_input += """
								output {0}.pdb""".format(self.system)
		
		
		defo_packmol_input = ""
		if self.defo is not None:
			defo_packmol_input = """outside cylinder  {0} {1}  0. 0.  0.  1.  {2}  {3}""".format(self.dimensions['LY']/2., self.dimensions['LY']/2.,
																										self.defo['Radius'], self.dimensions['LZ'])
		
		self.packmol_input += """
								
								#Bottom layer
								
								structure {1}
								chain A
								resnumbers 3
								number {2:g}
								inside box 0. 0. {3}  {4} {5} {6}
								constrain_rotation x 0 10
								constrain_rotation y 0 10
								{9}
								end structure
							
								#Top layer
								
								
								structure {1}
								chain B
								resnumbers 3
								number {2:g}
								inside box 0. 0.  {7}  {4} {5} {8}
								constrain_rotation x 180 10
								constrain_rotation y 180 10
								{9}
								end structure
									
							""".format(self.system, pdb_file_list[self.lipid_type]['name'],
										self.nb_lipid_monolayer[self.lipid_type],
										self.m1_tail_min, self.dimensions['LX'] - shell,
										self.dimensions['LY'] -shell, self.m1_head_max,
										self.m2_head_min, self.m2_tail_max, 
										defo_packmol_input)
		# 2 chains for the lipids
		self.nb_index += 1
			
		packmol_instruction_W_mol = "end structure"
		if self.solvent_type =='PW':
			packmol_instruction_W_mol = """
										atoms 2 3
										radius 0.2
										end atoms
									end structure"""
			
		self.packmol_input += """
								# Solvent
								structure {0}
									chain D
									number {1:g}
									inside box 0. 0. {2}  {3} {4} {5}
								{6}
						""".format(pdb_file_list[self.solvent_type]['name'], self.sample_molecules[self.solvent_type] ,
									self.m1_head_max,
									self.dimensions['LX'] - shell, self.dimensions['LY'] - shell,
									self.m2_head_min, packmol_instruction_W_mol)
			
		# A chain for solvent
		self.nb_index += 1
			
		if self.defo is not None:
			mid_top_layer		= self.height + self.tmt/2. + self.sol_thickness/2.
			mid_bottom_layer	= self.height - self.tmt/2. - self.sol_thickness/2.
			
			if self.defo['Height'] == 'follow':
				self.packmol_input += """
								# top layer defo
								structure {0}
									chain X
									number 1
									resnumbers 3
									fixed {1} {2} {3} 0.0 0.0 0.0
									radius 0.
									center
								end structure
								
								# bottom layer defo
								structure {4}
									chain Y
									number 1
									resnumbers 3
									fixed {1} {2} {5} 0.0 0.0 0.0
									radius 0.
									center
								end structure
							""".format(self.defo_dict['DEFT']['xyz file'].replace('xyz','pdb'),
										self.dimensions['LX']/2.0, self.dimensions['LY']/2.0,
										mid_top_layer, self.defo_dict['DEFB']['xyz file'].replace('xyz','pdb'),
										mid_bottom_layer)
				self.nb_index += 1
			
			else:
				self.packmol_input += """
								# defo
								structure {0}
									chain X
									number 1
									resnumbers 3
									fixed {1} {2} {3} 0.0 0.0 0.0
									radius 0.
								end structure
							""".format(self.defo_dict['DEFB']['xyz file'].replace('xyz','pdb').
										self.dimensions['LX']/2.0, self.dimensions['LY']/2.0,
										self.bottom_z)
				
			# for DEF
			self.nb_index += 1
						

		if self.su is not None:
			shell_su = 0.25
			self.packmol_input += """
								# Substrate
								structure {0}
									chain S
									number {1:g}
									inside box 0. 0. {5}  {2} {3} {4}
								end structure
							""".format(pdb_file_list['SU']['name'], self.nb_su, self.dimensions['LX'], self.dimensions['LY'], self.su['Thickness'],shell_su)
			
			# for SU
			self.nb_index += 1
					
			# pdb file for su
			assert(len(self.su['SuType']) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
			
			su_type = self.su['SuType'] + (4 - len(self.su['SuType']))*' '
			
			atom_type = None
			if self.su['SuType'].startswith('SU'):
				atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
			elif self.su['SuType'].startswith('S'):
				atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
			
			su_pdb_content = pdb_file_list['SU']['content'].replace("TEMP", su_type)
			su_pdb_content = su_pdb_content.replace("TEM", atom_type)
			
			with open(pdb_file_list['SU']['name'],'w') as su_pdb:
				su_pdb.write(ut.RemoveUnwantedIndent(su_pdb_content))
			
					
		with open('packmol_'+self.system+'.input','w') as pack_input:
			pack_input.write( ut.RemoveUnwantedIndent(self.packmol_input))
			
		with open(pdb_file_list[self.lipid_type]['name'],'w') as lipid_pdb:
			lipid_pdb.write( ut.RemoveUnwantedIndent(pdb_file_list[self.lipid_type]['content']) )
			
		with open(pdb_file_list[self.solvent_type]['name'],'w') as solvent_pdb:
			solvent_pdb.write( ut.RemoveUnwantedIndent(pdb_file_list[self.solvent_type]['content']) )
	
	
	
	
	def create_bilayer(self, index):
		"""Method calling packmol to create the sample"""
		packmol_cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """
					).format(self.softwares['PACKMOL'], self.system)
		
		sub.call(packmol_cmd, shell=True)
		
		#Convert the pdb to gro and put the pbc then check the box vectors at the end of the .gro
		editconf_cmd = str("{0} editconf -f {1}.pdb -o {1}.gro "
							"-box {2} {3} {4} -c no").format(self.softwares['GROMACS_LOC'], self.system,
														0.1*self.dimensions['LX'], 0.1*self.dimensions['LY'], 0.1*self.dimensions['LZ'])
		
		sub.call(editconf_cmd, shell=True, stdout=open('pdb2gro.txt','w'), stderr=open('err_pdb2gro.txt','w'))
		"""
		write_box = str("
			mol load pdb {0}.pdb
			set all [atomselect top "all"]

			package require pbctools
			pbc set {{{1} {2} {3}}}

			$all writepdb "{0}.withbox.pdb"
			unset all

			exit
			").format(self.system, self.dimensions['LX'], self.dimensions['LY'], self.dimensions['LZ'])

		with open('write_box.vmd','w+') as vmd_script:
			vmd_script.write(ut.RemoveUnwantedIndent(write_box) )
		
		vmd_cmd = "{0} -dispdev text -e write_box.vmd > write_box.log".format(self.softwares['VMD'])
		sub.call(vmd_cmd, shell=True)
		"""
		
		self.check_pbc()
		
		make_index = "chain A\nchain B\n"
		naming_index = "name {0} bottom{2}\nname {1} top{2}\n".format(self.nb_index, self.nb_index+1, self.lipid_type)
		count_index = 2
		
		
		if self.nb_sol_bottom != 0:
			make_index += "chain D\n"
			naming_index += "name {0} bottom{1}\n".format(self.nb_index+count_index, self.solvent_type)
			count_index += 1
			
		if self.nb_sol_top != 0:
			make_index += "chain E\n"
			naming_index += "name {0} top{1}\n".format(self.nb_index+count_index, self.solvent_type)
			count_index += 1
		
		if self.mono is not None: 
			make_index += "chain C\n"
			naming_index += "name {0} mono{1}\n".format(self.nb_index+count_index, self.lipid_type)
			count_index += 1
			
			make_index += "chain C\n"
			naming_index += "name {0} mono\n".format(self.nb_index+count_index, self.lipid_type)
			count_index += 1
		
		if self.nb_sol_vapor is not None:
			make_index += "chain V\n"
			naming_index += "name {1} vap{3}\n".format(self.nb_index+count_index, self.solvent_type)
			count_index += 1
		
		if self.su is not None:
			make_index += "chain S\n"
			naming_index += "name {0} su\n".format(self.nb_index+count_index)
			count_index += 1
		
		if self.defo is not None:
			if self.defo['Height'] == 'follow':
				make_index += "chain X\n"
				naming_index += "name {0} defo_bilayer\n".format(self.nb_index+count_index)
				
				if self.mono is not None:
					make_index += "chain Y\n"
					naming_index += "name {0} defo_mono\n".format(self.nb_index+count_index+1)
					naming_index += "{0} | {1}\nname {2} defo\n".format(self.nb_index+count_index, self.nb_index+count_index+1, self.nb_index+count_index+2)
					count_index += 3
					
				else:
					naming_index += "name {0} defo\n".format(self.nb_index+count_index)
					count_index += 1
							
				
			else:
				make_index += "chain X\n"
				naming_index += "name {0} defo\n".format(self.nb_index+count_index)
				count_index += 1
			
		
		#naming the bilayer
		naming_index += "{0} | {1}\nname {2} bilayer\nq\n".format(self.nb_index, self.nb_index+1, self.nb_index+count_index)
		count_index += 1
		
		self.nb_index += count_index
		# Writing the script to file
		make_index_str = make_index + naming_index
		
		with open('make_index.input','w+') as make_index_script:
			make_index_script.write(make_index_str)
		
		# Calling the script
		with open('make_index.output','w') as outfile:
			errfile = open('make_index.err','w')
			make_index_cmd = str("""{0}make_ndx -f {1}.pdb -o {1}.ndx < make_index.input"""
												).format(self.softwares['GROMACS_LOC'], self.system)
			sub.call(make_index_cmd, shell=True, stdout=outfile, stderr=errfile)
			
			errfile.close()
		
		self.output_file = "{0}.gro".format(self.system)
		self.index_file = "{0}.ndx".format(self.system)
		
		ApL = self.dimensions['LX'] * self.dimensions['LY']
		
		if self.defo is not None:
			ApL -= 3.141516 * float(self.defo['Radius'])**2
		
		ApL /= self.nb_lipid_monolayer[self.lipid_type]
		
		message = """
				================================
				================================
				Packmol finished initial input file
				{0} {1} || {2} {3} ||""".format(self.lipid_type, self.sample_molecules[self.lipid_type],
											self.solvent_type, self.sample_molecules[self.solvent_type])
		if self.mono is not None:
			message += " mono{0} {1} ||".format(self.lipid_type, self.mono['NbLipidsM'])
			
		if self.defo is not None:
			message += " DEF{0} {1} ".format(self.defo['Version'], self.nb_defo)
		
		if self.su is not None:
			message += " SU{0} {1} ".format(self.su['Version'], self.nb_su)
				
		message += """
				box sizes : {0}, {1}, {2}
				Area per lipid = {3} A**2
				===============================
				===============================
				""".format(self.dimensions['LX'], self.dimensions['LY'], self.dimensions['LZ'], ApL)
		
		print(ut.RemoveUnwantedIndent(message))
	
	
	
	def write_packmol_bilayer(self):
		""" Method to write packmol input for membrane creation"""
		self.packmol_input += """
							output {0}.pdb""".format(self.system)
		
		defo_packmol_input = ""
		if self.defo is not None:
			defo_packmol_input = """
									outside cylinder  {0} {1}  0. 0.  0.  1.  {2}  {3}
							""".format(self.dimensions['LY']/2., self.dimensions['LY']/2.,
										self.defo['Radius'], self.dimensions['LZ'])
		
		self.packmol_input += """
								# Bottom layer
								structure {1}
									chain A
									resnumbers 3
									number {2:g}
									inside box 0. 0. {3}  {4} {5} {6}
									constrain_rotation x 180 10
									constrain_rotation y 180 10
									{9}
								end structure
								
								# Top layer
								structure {1}
									chain B
									resnumbers 3
									number {2:g}
									inside box 0. 0.  {7}  {4} {5} {8}
									constrain_rotation x 0 10
									constrain_rotation y 0 10
									{9}
								end structure
									
							""".format(self.system, pdb_file_list[self.lipid_type]['name'],
										self.nb_lipid_monolayer[self.lipid_type],
										self.m1_head_min, self.dimensions['LX'],
										self.dimensions['LY'], self.m1_tail_max,
										self.m2_tail_min, self.m2_head_max, 
										defo_packmol_input)
		# 2 chains for the lipids
		self.nb_index += 1
		
		if self.mono is not None:
			defo_mono_packmol_input = ""
			if self.defo is not None:
				if self.defo['Height'] in ["follow", "box", "mono"]:
					defo_mono_packmol_input = defo_packmol_input
			
			self.packmol_input +="""
								# Monolayer
								structure {0}
									chain C
									resnumbers 3
									number {1:g}
									inside box 0. 0.  {2}  {3} {4} {5}
									constrain_rotation x 180 10
									constrain_rotation y 180 10
									{6}
								end structure
								
							""".format(pdb_file_list[self.lipid_type]['name'], self.mono['NbLipidsM'], self.mono['LzM'],
										self.dimensions['LX'], self.dimensions['LY'],
										self.mono['LzM'] + self.tmt,
										defo_mono_packmol_input)
			
		packmol_instruction_W_mol = ""
		if self.solvent_type =='PW':
			packmol_instruction_W_mol = """
									atoms 2 3
									radius 0.2
									end atoms"""
			
		if self.nb_sol_bottom:
			self.packmol_input += """
								# Bottom Solvent
								structure {0}
									chain D
									number {1:g}
									inside box 0. 0. {2}  {3} {4} {5}
									{6}
								end structure
							""".format(pdb_file_list[self.solvent_type]['name'], self.nb_sol_bottom,
										self.bottom_z,
										self.dimensions['LX'], self.dimensions['LY'],
										self.m1_head_min, packmol_instruction_W_mol)
			
			
		if self.nb_sol_top:
			self.packmol_input += """
								# Top solvent
								structure {0}
									chain E
									number {1:g}
									inside box 0. 0. {2}  {3} {4} {5}
									{6}
								end structure
							""".format(pdb_file_list[self.solvent_type]['name'], self.nb_sol_top,
										self.m2_head_max,
										self.dimensions['LX'], self.dimensions['LY'],
										self.dimensions['LZ'], packmol_instruction_W_mol)
			

		
		
		#if self.nb_sol_vapor:
			#self.packmol_input += """
							## Vapor solvent
							#structure {0}
								#chain V
								#number {1:g}
								#inside box 0. 0. {2}  {3} {4} {5}
								#{6}
							#end structure
							#""".format(pdb_file_list[self.solvent_type]['name'],self.nb_sol_vapor,
										#self.mono['LzM']+self.tmt,
										#self.dimensions['LX'], self.dimensions['LY'],
										#self.lz_vaccum, packmol_instruction_W_mol)
			
		# A chain for solvent
		self.nb_index += 1
			
		if self.defo is not None:
			if self.defo['Height'] == 'follow':
				self.packmol_input += """
								# Bilayer defo
								structure {0}
									chain X
									number 1
									resnumbers 3
									fixed {1} {2} {3} 0.0 0.0 0.0
									radius 0.
									center
								end structure
							""".format(self.defo_dict['DEFB']['xyz file'].replace('xyz','pdb'),
										self.dimensions['LX']/2.0, self.dimensions['LY']/2.0,
										self.height)
				
				if self.mono is not None:
					self.packmol_input += """
								# Monolayer defo
								structure {0}
									chain Y
									number 1
									resnumbers 3
									fixed {1} {2} {3} 0.0 0.0 0.0
									radius 0.
								end structure
							""".format(self.defo_dict['DEFM']['xyz file'].replace('xyz','pdb'),
										self.dimensions['LX']/2.0, self.dimensions['LY']/2.0,
										self.mono['LzM']-self.dz/8.0)
					self.nb_index += 1
				
				
			else:
				self.packmol_input += """
								# defo
								structure {0}
									chain X
									number 1
									resnumbers 3
									fixed {1} {2} {3} 0.0 0.0 0.0
									radius 0.
								end structure
							""".format(self.defo_dict['DEFB']['xyz file'].replace('xyz','pdb').
										self.dimensions['LX']/2.0, self.dimensions['LY']/2.0,
										self.bottom_z)
				
			# for DEF
			self.nb_index += 1
						
		if self.su is not None:
			shell_su = 0.25
			self.packmol_input += """
								# Substrate
								structure {0}
									chain S
									number {1:g}
									inside box 0. 0. {5}  {2} {3} {4}
								end structure
							""".format(pdb_file_list['SU']['name'], self.nb_su, self.dimensions['LX'], self.dimensions['LY'], self.su['Thickness'],shell_su)
			
			# for SU
			self.nb_index += 1
					
			# pdb file for su
			assert(len(self.su['SuType']) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
			
			su_type = self.su['SuType'] + (4 - len(self.su['SuType']))*' '
			
			atom_type = None
			if self.su['SuType'].startswith('SU'):
				atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
			elif self.su['SuType'].startswith('S'):
				atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
			
			su_pdb_content = pdb_file_list['SU']['content'].replace("TEMP", su_type)
			su_pdb_content = su_pdb_content.replace("TEM", atom_type)
			
			with open(pdb_file_list['SU']['name'],'w') as su_pdb:
				su_pdb.write(ut.RemoveUnwantedIndent(su_pdb_content))
			
					
		with open('packmol_'+self.system+'.input','w') as pack_input:
			pack_input.write( ut.RemoveUnwantedIndent(self.packmol_input))
			
		with open(pdb_file_list[self.lipid_type]['name'],'w') as lipid_pdb:
			lipid_pdb.write( ut.RemoveUnwantedIndent(pdb_file_list[self.lipid_type]['content']) )
			
		with open(pdb_file_list[self.solvent_type]['name'],'w') as solvent_pdb:
			solvent_pdb.write( ut.RemoveUnwantedIndent(pdb_file_list[self.solvent_type]['content']) )
	
	
	
	def create_monolayer(self, index):
		"""Method calling packmol to create the sample"""
		packmol_cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """
					).format(self.softwares['PACKMOL'], self.system)
		
		sub.call(packmol_cmd, shell=True)
		
		#Convert the pdb to gro and put the pbc then check the box vectors at the end of the .gro
		editconf_cmd = str("{0} editconf -f {1}.pdb -o {1}.gro "
							"-box {2} {3} {4} -c no").format(self.softwares['GROMACS_LOC'], self.system,
														0.1*self.dimensions['LX'], 0.1*self.dimensions['LY'], 0.1*self.dimensions['LZ'])
		
		sub.call(editconf_cmd, shell=True, stdout=open('pdb2gro.txt','w'), stderr=open('err_pdb2gro.txt','w'))
		
		self.check_pbc()
		
		make_index = "chain A\nchain B\n"
		naming_index = "name {0} bottom_mono\nname {1} top_mono\n".format(self.nb_index, self.nb_index+1, self.lipid_type)
		count_index = 2
		
		
		if self.su is not None:
			make_index += "chain S\n"
			naming_index += "name {0} su\n".format(self.nb_index+count_index)
			count_index += 1
		
		if self.defo is not None:
			if self.defo['Height'] == 'follow':
				make_index += "chain X\n"
				naming_index += "name {0} defo_top\n".format(self.nb_index + count_index)
				
				make_index += "chain Y\n"
				naming_index += "name {0} defo_bottom\n".format(self.nb_index + count_index + 1)
				
				naming_index += "{0} | {1}\nname {2} defo\n".format(self.nb_index + count_index, self.nb_index + count_index + 1, self.nb_index + count_index + 2)
				count_index += 3
					
			else:
				make_index += "chain X\n"
				naming_index += "name {0} defo\n".format(self.nb_index+count_index)
				count_index += 1
			
		
		# end of the make_index file
		naming_index += "q\n"
		
		self.nb_index += count_index
		# Writing the script to file
		make_index_str = make_index + naming_index
		
		with open('make_index.input','w+') as make_index_script:
			make_index_script.write(make_index_str)
		
		# Calling the script
		with open('make_index.output','w') as outfile:
			errfile = open('make_index.err','w')
			make_index_cmd = str("""{0}make_ndx -f {1}.pdb -o {1}.ndx < make_index.input"""
												).format(self.softwares['GROMACS_LOC'], self.system)
			sub.call(make_index_cmd, shell=True, stdout=outfile, stderr=errfile)
			
			errfile.close()
		
		self.output_file = "{0}.gro".format(self.system)
		self.index_file = "{0}.ndx".format(self.system)
		
		ApL = self.dimensions['LX'] * self.dimensions['LY']
		
		if self.defo is not None:
			ApL -= 3.141516 * float(self.defo['Radius'])**2
		
		ApL /= self.nb_lipid_monolayer[self.lipid_type]
		
		message = """
				================================
				================================
				Packmol finished initial input file
				{0} {1} || {2} {3} ||""".format(self.lipid_type, self.sample_molecules[self.lipid_type],
											self.solvent_type, self.sample_molecules[self.solvent_type])
			
		if self.defo is not None:
			message += " DEF{0} {1} ".format(self.defo['Version'], self.nb_defo)
		
		if self.su is not None:
			message += " SU{0} {1} ".format(self.su['Version'], self.nb_su)
				
		message += """
				box sizes : {0}, {1}, {2}
				Area per lipid = {3} A**2
				===============================
				===============================
				""".format(self.dimensions['LX'], self.dimensions['LY'], self.dimensions['LZ'], ApL)
		
		print(ut.RemoveUnwantedIndent(message))
	
	
	
	







class Solvent(BaseProject):
	""" Class for creating solvent samples
	"""
	def __init__(self, sample, softwares, path_to_default):
		super(Solvent, self).__init__(sample, softwares, path_to_default)
		
	def make(self):
		self.solvent_types = []
		
		self.LX = 0.0
		self.LY = 0.0
		self.LZ = 0.0
		
		if self.fillmode is not None:
			if self.fillmode == 'FULL':
				self.LX = self.dimensions['LX']
				self.LY = self.dimensions['LY']
				self.LZ = self.dimensions['LZ']
				
			elif self.fillmode == 'BOX':
				self.LX = self.LY = self.LZ = self.dimensions['LX']
				
			elif self.fillmode == 'SPHERE':
				self.center = (self.dimensions['LX']/2., self.dimensions['LY']/2., self.dimensions['LZ']/2.) 
				
			elif self.fillmode == 'HALF-X':
				self.LX = self.dimensions['LX']/2.0 
				self.LY = self.dimensions['LY']
				self.LZ = self.dimensions['LZ']
				
			elif self.fillmode == 'HALF-Y':
				self.LX = self.dimensions['LX']
				self.LY = self.dimensions['LY']/2.0 
				self.LZ = self.dimensions['LZ']
				
			elif self.fillmode == 'HALF-Z':
				self.LX = self.dimensions['LX']
				self.LY = self.dimensions['LY']
				self.LZ = self.dimensions['LZ']/2.0
			
			elif self.fillmode == 'PARTIAL-Z':
				self.LX = self.dimensions['LX']
				self.LY = self.dimensions['LY']
				#self.bottom_z= self.minsolvent
				self.LZ = self.maxsolvent
				
			elif self.fillmode == 'PARTIAL-X':
				assert(False),"PARTIAL-X not implemented yet"
				
			elif self.fillmode == 'PARTIAL-Y':
				assert(False),"PARTIAL-Y not implemented yet"	
				
			else:
				self.LX = self.dimensions['LX']
				self.LY = self.dimensions['LY']
				self.LZ = self.dimensions['LZ']
		else:
			self.LX = self.dimensions['LX']
			self.LY = self.dimensions['LY']
			self.LZ = self.dimensions['LZ']
				
		self.finding_molecules()
		
		self.system = self.system = "{0}".format(self.type)
		for sol in self.solvent_types:
			self.system += "_{0}{1}".format(self.sample_molecules[sol], sol)
		
		#Adding the substrate
		self.nb_su = None
		
		if self.su is not None:
			self.add_substrate()
			self.system += "_{0}{1}{2}".format(self.nb_su, self.su['SuType'], self.su['Version'])
		
		if self.wall is not None:
			self.system += "_WALL"
		
		self.nb_index = 2
		self.write_packmol_input()
		
		if self.lz_vaccum is not None:
			self.dimensions['LZ'] = self.lz_vaccum
		
		# Calling packmol, vmd and gmx make_index
		self.create_sample()
		
		# creating the topology for su and defo if in sample
		self.create_topology()
		
		
	
	def finding_molecules(self):
		""" Method to find the molecules and associate their number from csv file.
		"""
		for sol in self.solvent_list:
			if sol in self.sample_molecules:
				self.solvent_types.append(sol)
	
	
	
	
	def add_substrate(self):
		""" Method to add the substrate at the bottom of box
		"""
		# Shifting all the coordinates in Z
		self.bottom_z += float(self.su['Thickness'])
		#Pushing the box limits
		self.LZ += float(self.su['Thickness'])
		self.dimensions['LZ'] += float(self.su['Thickness'])
		
		self.nb_su = int( float(self.su['Density']) * self.dimensions['LX'] * self.dimensions['LY'] * float(self.su['Thickness']) / 1000 )
	
	
	
	def write_packmol_input(self):
		""" Method to write packmol input for membrane creation"""
		self.packmol_input += """
								output {0}.pdb
								""".format(self.system)
		
		shift = 0.
		nb_solvent = float(len(self.solvent_types))
		
		if self.fillmode == 'SPHERE':
			prev_radius = 0.0
			for chain, sol in enumerate(self.solvent_types):
				packmol_instruction_mol = ""
				
				if sol =='PW':
					packmol_instruction_mol = """
									atoms 2 3
									radius 0.2
									end atoms"""
				outside_of = ""
				if nb_solvent > 1:
					outside_of = "outside sphere {1} {2} {3} {4}".format(self.center[0],
																		self.center[1],
																		self.center[2],
																		prev_radius)
				self.packmol_input += """
								
								structure {0}
									chain {1}
									number {2:g}
									inside sphere {3} {4} {5} {6}
									{7}
									{8}
								end structure
								""".format(pdb_file_list[sol]['name'], chain,
											self.sample_molecules[sol], 
											self.center[0], self.center[1], self.center[2],
											self.radius/nb_solvent,
											outside_of,
											packmol_instruction_mol)
				prev_radius = self.radius
				shift += self.radius / nb_solvent
				
				self.nb_index += 1
		
		else:
			for chain, sol in enumerate(self.solvent_types):
				packmol_instruction_mol = ""
				
				if sol =='PW':
					packmol_instruction_mol = """
									atoms 2 3
									radius 0.2
									end atoms"""
					
				self.packmol_input += """
								
								structure {0}
									chain {1}
									number {2:g}
									inside box 0. 0. {3}  {4} {5} {6}
									{7}
								end structure
								""".format(pdb_file_list[sol]['name'], chain,
											self.sample_molecules[sol], self.bottom_z + shift,
											self.LX, self.LY, self.LZ/nb_solvent + shift,
											packmol_instruction_mol)
								
				shift += self.dimensions['LZ'] / nb_solvent
				
				self.nb_index += 1
		
		if self.su is not None:
			shell_su = 0.25            
			self.packmol_input += """
							# Substrate
							structure {0}
								chain S
								number {1:g}
								inside box 0. 0. {5}  {2} {3} {4}
							end structure
							""".format(pdb_file_list['SU']['name'], self.nb_su, self.dimensions['LX'],
										self.dimensions['LY'], self.su['Thickness'],shell_su)
			
			# for SU
			self.nb_index += 1
					
			# pdb file for su
			assert(len(self.su['SuType']) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
			
			su_type = self.su['SuType'] + (4 - len(self.su['SuType']))*' '
			
			atom_type = None
			if self.su['SuType'].startswith('SU'):
				atom_type = self.su['SuType'][-2:] + (3 - len(self.su['SuType'][-2:]))*' '
			elif self.su['SuType'].startswith('S'):
				atom_type = self.su['SuType'][-3:] + (3 - len(self.su['SuType'][-3:]))*' '
			
			su_pdb_content = pdb_file_list['SU']['content'].replace("TEMP", su_type)
			su_pdb_content = su_pdb_content.replace("TEM", atom_type)
			
			with open(pdb_file_list['SU']['name'],'w') as su_pdb:
				su_pdb.write(ut.RemoveUnwantedIndent(su_pdb_content))
			
					
		with open('packmol_'+self.system+'.input','w') as pack_input:
			pack_input.write( ut.RemoveUnwantedIndent(self.packmol_input) )
		
		for sol in self.solvent_types:
			with open(pdb_file_list[sol]['name'],'w') as solvent_pdb:
				solvent_pdb.write( ut.RemoveUnwantedIndent(pdb_file_list[sol]['content']) )
	
	
	
	
	
	def create_sample(self):
		packmol_cmd = str("{0} < packmol_{1}.input > "
							"packmol_{1}.output ").format(self.softwares['PACKMOL'], self.system)
		
		sub.call(packmol_cmd, shell=True)
		
		#Convert the pdb to gro and put the pbc then check the box vectors at the end of the .gro
		shell_wall = 0.0
		if self.wall is not None:
			shell_wall += 5.0
		editconf_cmd = str("{0} editconf -f {1}.pdb -o {1}.gro "
							"-box {2} {3} {4} -c no ").format(self.softwares['GROMACS_LOC'], self.system,
														0.1*self.dimensions['LX'], 0.1*self.dimensions['LY'], 0.1*(self.dimensions['LZ']+shell_wall))
		sub.call(editconf_cmd, shell=True, stdout=open('pdb2gro.txt','w'), stderr=open('err_pdb2gro.txt','w'))
		#-translate 0. 0. {5}
		"""
		write_box = str("
			mol load pdb {0}.pdb
			set all [atomselect top "all"]

			package require pbctools
			pbc set {{{1} {2} {3}}}

			$all writepdb "{0}.withbox.pdb"
			unset all

			exit
			").format(self.system, self.dimensions['LX'], self.dimensions['LY'], self.dimensions['LZ'])

		with open('write_box.vmd','w+') as vmd_script:
			vmd_script.write(ut.RemoveUnwantedIndent(write_box) )
		
		vmd_cmd = "{0} -dispdev text -e write_box.vmd > write_box.log".format(self.softwares['VMD'])
		sub.call(vmd_cmd, shell=True)
		"""
		
		self.check_pbc()
		
		make_index = ""
		naming_index = ""
		
		count_index = 0
		
		if self.su is not None:
			make_index += "chain S\n"
			naming_index += "name {0} su\n".format(self.nb_index)
			count_index += 1
		
		self.nb_index += count_index
		
		naming_index += "q\n"
		# Writing the script to file
		make_index_str = make_index + naming_index
		with open('make_index.input','w+') as make_index_script:
			make_index_script.write(make_index_str)
		
		# Calling the script
		with open('make_index.output','w') as outfile:
			errfile = open('make_index.err','w')
			make_index_cmd = str("""{0}make_ndx -f {1}.pdb -o {1}.ndx < make_index.input"""
												).format(self.softwares['GROMACS_LOC'], self.system)
			sub.call(make_index_cmd, shell=True, stdout=outfile, stderr=errfile)
			
			errfile.close()
		
		self.output_file = "{0}.gro".format(self.system)
		self.index_file = "{0}.ndx".format(self.system)
		
		message = """
				================================
				================================
				Packmol finished initial input file
				"""
		for sol in self.solvent_types:
			message += "{0} {1}".format(sol, self.sample_molecules[sol])
		
		message += """
				box size : {0}, {1}, {2}
				===============================
				===============================
				""".format(self.dimensions['LX'], self.dimensions['LY'], self.dimensions['LZ'])
		
		print(ut.RemoveUnwantedIndent(message))
