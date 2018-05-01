#!/usr/bin/python3
# -*- coding: utf-8 -*-


import csv		#Module to read csv files
import os
import re
import sys
import subprocess as sub
import copy
import shutil
import argparse
import glob

import Utility as ut
import Prepare

LipidsList = ['DSPC','DPPC','DLPC']

SolventsList = ['W','OCO','PW']

parameters_csv = {'SU': 'Parameters_su.csv', 'GADD_SU': 'Parameters_su.csv', 'PADD_SU': 'Parameters_su.csv',
				  'WALL':'Parameters_wall.csv', 'ADD_WALL': 'Parameters_wall.csv',
				  'DEFO': 'Parameters_defo.csv'}

# String storing the name of the project
project_name = ''
# String storing the path to default gromamcs files
path_to_default = ''
# Dictionnary with the path to the Softwares to use
Softwares = {}
# Dictionnary holding all information for each job
Jobs = {}
#Dictionnary holding the data for PBS
PBS = {}
#Dictionnary holding the data for TGCC
TGCC = {}

def main(argv=sys.argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s','--send', dest='send', type=str, required=True,
						help='Set the type for sending jobs: local, pbs, tgcc')
	
	parser.add_argument('-c','--create', dest='create', type=str, default='local',
						help='Set the type for creating the samples : local, pbs, tgcc')
	
	parser.add_argument('-p','--param', dest='parameter', type=str, default="Parameters.csv",
						help='Set the parameter file to use for jobs')

	cmd_param = parser.parse_args(sys.argv[1:])
	#**********************************************#
	#**********************************************#
	#**********************************************#
	#**********************************************#
	# Read the set of parameters for the different #
	# samples to study							   #
	#**********************************************#
	#**********************************************#
	#**********************************************#

	#Default name for the parameter file.
	parameter_file = cmd_param.parameter
	
	if cmd_param.send == 'pbs':
		with open('PBSinfo.csv') as pbs_info:
			reader = csv.reader(pbs_info, delimiter=',', skipinitialspace=True)
			for row in reader:
				info = row[0].strip()
				value = row[1].strip()
				PBS[info] = value
				
	elif cmd_param.send == 'tgcc':
		with open('TGCCinfo.csv') as pbs_info:
			reader = csv.reader(pbs_info, delimiter=',', skipinitialspace=True)
			for row in reader:
				info = row[0].strip()
				value = row[1].strip()
				TGCC[info] = value
	
	with open(parameter_file,'r') as Input_Params:
		inputParamsReader = csv.reader(Input_Params, delimiter=',', skipinitialspace=True)
		
		#Variable to keep track of the jobs ordering
		jobNumber = 0
		
		for row in inputParamsReader:
			row_name = row[0].strip()
			# Get the project name
			if row_name == 'PROJECT' :
				project_name = row[1]
				continue
			
			#Skip empty rows 
			test_row = list(filter(None, row))
			if not test_row: continue
			
			#Skip lines starting with '#' or empty lines
			if '#' in row_name or not test_row: continue
			
			# Each time a 'JOBID' is found a new job is created
			# it to separate different samples
			if row_name == 'JOBID':
				jobNumber += 1
				#Step of the run : EM, NVE, NVT, NPT ...
				stepNumber = 0
				Jobs.update( {jobNumber: {row[0].strip(' '): row[1].strip(' '), row[2].strip(' '): int(row[3]), row[4].strip(' '):row[5].strip(' ')} } )
				continue
			
			# Before reading information on each job the program stores the path to the Softwares and gromacs default files from Parameters.csv
			if jobNumber < 1 : 
				if row_name == 'GROMACS_LOC':
					Softwares.update({row[0].strip(' '): row[1].strip(' '), row[2].strip(' '): row[3].strip(' ')})
					continue
				
				if row_name == 'VMD':
					Softwares.update({row[0].strip(' '): row[1].strip(' ')})
					continue
				
				if row_name == 'PACKMOL':
					Softwares.update({row[0].strip(' '): row[1].strip(' ')})
					continue
				
				if row_name == 'GROMACS_DEFAULT':
					path_to_default = row[1].strip(' ')
			
			# Then if the program has reach a job, it will store all information available from Parameters.csv
			else:
				# Stores the node, number of nodes and processors per node.
				if row_name == 'NODE' :
					Jobs[jobNumber].update({row[0].strip(' '): row[1].strip(' '), row[2].strip(' '): row[3].strip(' '), row[4].strip(' '): row[5].strip(' ')})
					continue
				
				# Stores the options to pass to mdrun
				if row_name == 'MDRUN_OPT':
					Jobs[jobNumber].update({row[0].strip(' '): row[1].strip(' ')})
					continue
				
				# Stores the elements in the sample for the job
				if row_name == 'SAMPLE':
					Jobs[jobNumber]['SAMPLE'] = {}
					for element in range(1, len(row), 2):
						if row[element]:
							qtty = int(row[element+1])
							Jobs[jobNumber]['SAMPLE'].update( { row[element].strip(' '): qtty })
					continue
				
				# For a general thermostat parameters. 
				if row_name == 'THERMOSTAT':
					Jobs[jobNumber].update({ 'THERMOSTAT': {}})
					
					for param in range(1, len(row), 2):
						if row[param]:
							paramValue = row[param+1].strip(' ')
							
							if(ut.is_number(paramValue)):
								if('.' in ParameterValue or not float(ParameterValue).is_integer() ):
									Jobs[jobNumber]['BAROSTAT'].update( {row[param].strip(' '): float(paramValue) } )
									
								else:
									Jobs[jobNumber]['BAROSTAT'].update( {row[param].strip(' '): int(paramValue) } )
									
							else:
								Jobs[jobNumber]['BAROSTAT'].update( {row[param].strip(' '): paramValue })
								
						else:
							break
					continue
				
				# For general barostat parameters.
				if row_name == 'BAROSTAT':
					Jobs[jobNumber].update({ 'BAROSTAT': {}})
					
					for param in range(1, len(row), 2):
						if row[param]:
							paramValue = row[param+1].strip(' ')
							
							if(ut.is_number(paramValue)):
								if('.' in ParameterValue or not float(ParameterValue).is_integer() ):
									Jobs[jobNumber]['BAROSTAT'].update( {row[param].strip(' '): float(paramValue) } )
									
								else:
									Jobs[jobNumber]['BAROSTAT'].update( {row[param].strip(' '): int(paramValue) } )
									
							else:
								Jobs[jobNumber]['BAROSTAT'].update( {row[param].strip(' '): paramValue })
								
						else:
							break
					continue
				
				
				# Stores the type of the sample for the job
				if row_name == 'TYPE':
					Jobs[jobNumber].update({row[0].strip(' '): row[1].strip(' ')})
					continue
				
				# Stores the dimensions of the sample for the job
				if row_name == 'GEOMETRY' :
					Jobs[jobNumber]['GEOMETRY'] = {}
					for dimension in range(1, len(row), 2):
						if row[dimension]:
							dimValue = row[dimension+1].strip(' ')
							
							if ut.is_number(dimValue):
								Jobs[jobNumber]['GEOMETRY'].update( {row[dimension].strip(' '): float(dimValue) })
								
							else:
								Jobs[jobNumber]['GEOMETRY'].update( {row[dimension].strip(' '): dimValue })
								
						else:
							break
					
					continue
				
				# Stores the parameters for the monolayer of the sample for the job
				if row_name == 'MONO' :
					Jobs[jobNumber]['MONO'] = {}
					
					for mono_param in range(1, len(row), 2):
						if row[mono_param]:
							mono_value = row[mono_param+1].strip(' ')
							
							if ut.is_number(mono_value):
								Jobs[jobNumber]['MONO'].update( {row[mono_param].strip(' '): float(mono_value) })
								
							else:
								Jobs[jobNumber]['MONO'].update( {row[mono_param].strip(' '): mono_value })
								
						else:
							break
					
					continue
				
				"""
				## Stores information on the defo: defo per layer, Height, DzDefo, Radius, Version, TauDefo (thermostat), Tdefo (thermostat)
				#if row_name.startswith('DEFO'):
					#Jobs[jobNumber].update({ 'DEFO': {}})
					
					#for paramDefo in range(1, len(row), 2):
						#if row[paramDefo]:
							#paramDefoVal = row[paramDefo+1].strip(' ')
							#Jobs[jobNumber]['DEFO'].update( {row[paramDefo].strip(' '): paramDefoVal })
							
						#else:
							#break
						
					#with open(path_to_default+'/DEFO/Parameters_defo.csv','r') as Defo_Params:
						#defoParam = csv.reader(Defo_Params, delimiter=',', skipinitialspace=True)
						#Jobs[jobNumber]['DEFO'].update( { 'presets' : {} } )
						
						#for row in defoParam:
							##Skip empty rows 
							#test_row = list(filter(None, row))
							##Skip lines starting with '#' or empty lines
							#if '#' in row[0] or not test_row: continue
						
							#firstCol = row[0].strip(' ')
							#if firstCol:
								#name = firstCol
								##Creates a dictionary for each preset
								#Jobs[jobNumber]['DEFO']['presets'].update( { name : {}  })
								
								##Associates the parameters set in Parameters.csv to the correct protocol step dictionary
								#for index, param in enumerate(row[1:]):
									#if param:
										#Jobs[jobNumber]['DEFO']['presets'][name].update( { param.strip(' ') : index } )
							
							#else:
								#for param, paramValue in Jobs[jobNumber]['DEFO']['presets'][name].items():
									
									#paramValue = row[paramValue+1]
									
									#if paramValue:
										#if ut.is_number(paramValue):
											#if('.' in paramValue or not float(paramValue).is_integer() ):
												#Jobs[jobNumber]['DEFO']['presets'][name][param] = float(paramValue)
												
											#else:
												#Jobs[jobNumber]['DEFO']['presets'][name][param] = int(paramValue)
												
										#else:
											#Jobs[jobNumber]['DEFO']['presets'][name][param] = paramValue.strip(' ')
										
									
								
							
						
					
					continue
				"""
				
				# Stores information on the support: version, density, nblipids
				if row_name.startswith('SU') or row_name.startswith('WALL') or row_name.startswith('DEFO'):
					name = row[0].strip()
					Jobs[jobNumber].update({ name: {}})
					
					for param in range(1, len(row), 2):
						
						if row[param]:
							paramVal = row[param+1].strip(' ')
							Jobs[jobNumber][name].update( {row[param].strip(' '): paramVal })
							
						else:
							break
						
						
					path = "/{0}/{1}".format(name, parameters_csv[name])
					with open(path_to_default+path,'r') as Params:
						
						Param = csv.reader(Params, delimiter=',', skipinitialspace=True)
						Jobs[jobNumber][name].update( { 'presets' : {} } )
						
						for row in Param: 
							#Skip empty rows 
							test_row = list(filter(None, row))
							
							#Skip lines starting with '#' or empty lines
							if '#' in row[0] or not test_row: continue
						
							firstCol = row[0].strip(' ')
							if firstCol:
								preset_name = firstCol
								#Creates a dictionary for each preset
								Jobs[jobNumber][name]['presets'].update( { preset_name : {}  })
								
								#Associates the parameters set in Parameters.csv to the correct protocol step dictionary
								for index, param in enumerate(row[1:]):
									if param:
										if param.startswith('#'):
											continue
										Jobs[jobNumber][name]['presets'][preset_name].update( { param.strip(' ') : index } )
								
							else:
								for param, paramValue in Jobs[jobNumber][name]['presets'][preset_name].items():
									
									paramValue = row[paramValue+1]
									
									if paramValue:
										if ut.is_number(paramValue):
											if('.' in paramValue or not float(paramValue).is_integer() ):
												Jobs[jobNumber][name]['presets'][preset_name][param] = float(paramValue)
											else:
												Jobs[jobNumber][name]['presets'][preset_name][param] = int(paramValue)
										else:
											Jobs[jobNumber][name]['presets'][preset_name][param] = paramValue.strip(' ')
											
										
									
								
					continue
				
				if row_name.startswith('INPUT'):
					stepNb = str(stepNumber)
					Jobs[jobNumber]['PROTOCOL'].update( { stepNb: {'stepType' : row[0].strip(' ')}  })
					
					for dimension in range(1, len(row), 2):
						if row[dimension]:
							dimValue = row[dimension+1].strip(' ')
							
							if ut.is_number(dimValue):
								Jobs[jobNumber]['PROTOCOL'][stepNb].update( { row[dimension].strip(' '): float(dimValue) })
								
							else:
								Jobs[jobNumber]['PROTOCOL'][stepNb].update( { row[dimension].strip(' '): dimValue })
								
						else:
							break
						
					stepNumber += 1
					continue
				
				if row_name.startswith('TRANSLATE'):
					stepNb = str(stepNumber)
					Jobs[jobNumber]['PROTOCOL'].update( { stepNb: {'stepType' : row[0].strip(' ')}  })
					
					for dimension in range(1, len(row), 2):
						if row[dimension]:
							dimValue = row[dimension+1].strip(' ')
							
							if ut.is_number(dimValue):
								Jobs[jobNumber]['PROTOCOL'][stepNb].update( { row[dimension].strip(' '): float(dimValue) })
								
							else:
								Jobs[jobNumber]['PROTOCOL'][stepNb].update( { row[dimension].strip(' '): dimValue })
								
						else:
							break
						
					stepNumber += 1
					continue
				"""
				## Stores information on the support: version, density, nblipids
				#if row_name.startswith('WALL'):
					#Jobs[jobNumber].update({ 'WALL': {}})
					
					#for paramWall in range(1, len(row), 2):
						
						#if row[paramWall]:
							#paramWallVal = row[paramWall+1].strip(' ')
							#Jobs[jobNumber]['WALL'].update( {row[paramWall].strip(' '): paramWallVal })
							
						#else:
							#break
						
					#with open(path_to_default+'/WALL/Parameters_wall.csv','r') as Wall_Params:
						
						#wallParam = csv.reader(Wall_Params, delimiter=',', skipinitialspace=True)
						#Jobs[jobNumber]['WALL'].update( { 'presets' : {} } )
						
						#for row in wallParam: 
							##Skip empty rows 
							#test_row = list(filter(None, row))
							
							##Skip lines starting with '#' or empty lines
							#if '#' in row[0] or not test_row: continue
						
							#firstCol = row[0].strip(' ')
							#if firstCol:
								#name = firstCol
								##Creates a dictionary for each preset
								#Jobs[jobNumber]['WALL']['presets'].update( { name : {}  })
								
								##Associates the parameters set in Parameters.csv to the correct protocol step dictionary
								#for index, param in enumerate(row[1:]):
									#if param:
										#Jobs[jobNumber]['WALL']['presets'][name].update( { param.strip(' ') : index } )
								
							#else:
								#for param, paramValue in Jobs[jobNumber]['WALL']['presets'][name].items():
									
									#paramValue = row[paramValue+1]
									
									#if paramValue:
										#if ut.is_number(paramValue):
											#if('.' in paramValue or not float(paramValue).is_integer() ):
												#Jobs[jobNumber]['WALL']['presets'][name][param] = float(paramValue)
											#else:
												#Jobs[jobNumber]['WALL']['presets'][name][param] = int(paramValue)
										#else:
											#Jobs[jobNumber]['WALL']['presets'][name][param] = paramValue.strip(' ')
					#continue
					"""
				
				
				if row_name.startswith('GADD_SU') or row_name.startswith('ADD_WALL'):
					name = row[0].strip()
					
					stepNb = str(stepNumber)
					#Creates a dictionary for the step with an associated protocol step number to keep track of the order of execution
					Jobs[jobNumber]['PROTOCOL'].update( { stepNb: {'stepType' : row[0].strip(' ')}  })
					
					for param in range(1, len(row), 2):
						
						if row[param]:
							if row[param].strip() == 'PROTOCOL':
								Protocol = ['']*(stepNumber)
								if '+' in row[param+1]:
									Protocol.extend(row[param+1].split('+'))
								else:
									Protocol.extend([row[param+1].strip(' ')])
									
								if name.startswith('GADD_SU') or name.startswith('PADD_SU'):
									Jobs[jobNumber]['PROTOCOL'][stepNb].update( {'suProtocol' : Protocol} )
									
								elif name.startswith('ADD_WALL'):
									Jobs[jobNumber]['PROTOCOL'][stepNb].update( {'wallProtocol' : Protocol} )
								
							else:
								paramVal = row[param+1].strip(' ')
								Jobs[jobNumber]['PROTOCOL'][stepNb].update( {row[param].strip(' '): paramVal })
							
						else:
							break
					
					folder = name
					
					if folder.startswith('GADD_SU') or folder.startswith('PADD_SU'):
						folder = 'SU'
					elif folder.startswith('ADD_WALL'):
						folder = 'WALL'
						
					path = "/{0}/{1}".format(folder, parameters_csv[name])
					with open(path_to_default+path,'r') as Params:
						
						Param = csv.reader(Params, delimiter=',', skipinitialspace=True)
						Jobs[jobNumber]['PROTOCOL'][stepNb].update( { 'presets' : {} } )
						
						for row in Param: 
							#Skip empty rows 
							test_row = list(filter(None, row))
							
							#Skip lines starting with '#' or empty lines
							if '#' in row[0] or not test_row: continue
						
							firstCol = row[0].strip(' ')
							if firstCol:
								preset_name = firstCol
								#Creates a dictionary for each preset
								Jobs[jobNumber]['PROTOCOL'][stepNb]['presets'].update( { preset_name : {}  })
								
								#Associates the parameters set in Parameters.csv to the correct protocol step dictionary
								for index, param in enumerate(row[1:]):
									if param:
										if param.startswith('#'):
											continue
										Jobs[jobNumber]['PROTOCOL'][stepNb]['presets'][preset_name].update( { param.strip(' ') : index } )
								
								
							else:
								for param, paramValue in Jobs[jobNumber]['PROTOCOL'][stepNb]['presets'][preset_name].items():
									
									paramValue = row[paramValue+1]
									
									if paramValue:
										if ut.is_number(paramValue):
											if('.' in paramValue or not float(paramValue).is_integer() ):
												Jobs[jobNumber]['PROTOCOL'][stepNb]['presets'][preset_name][param] = float(paramValue)
											else:
												Jobs[jobNumber]['PROTOCOL'][stepNb]['presets'][preset_name][param] = int(paramValue)
										else:
											Jobs[jobNumber]['PROTOCOL'][stepNb]['presets'][preset_name][param] = paramValue.strip(' ')
										
									
								
							
					#Increases the protocol step number
					stepNumber += 1
					
					continue
				
				
				if row_name.startswith('PROTOCOL'):
					Jobs[jobNumber].update( {'PROTOCOL': {}} )
					
					for index, protocol in enumerate(row[1:]):
						if 'DEFO' in Jobs[jobNumber] and protocol == 'DEFO':
							defoProtocol = row[index+2].split('+')
							Jobs[jobNumber]['DEFO'].update( {'defoProtocol' : defoProtocol} )
							
						if 'SU' in Jobs[jobNumber] and protocol == 'SU':
							suProtocol = row[index+2].split('+')
							Jobs[jobNumber]['SU'].update( {'suProtocol' : suProtocol} )
						
						if 'WALL' in Jobs[jobNumber] and protocol == 'WALL':
							wallProtocol = row[index+2].split('+')
							Jobs[jobNumber]['WALL'].update( {'wallProtocol' : wallProtocol} )
						
						if protocol == 'RUNTYPE':
							# Either CREATE or PROD for creating and production resp.
							runProtocol = ['CREATE']
							runProtocol.extend(row[index+2].split('+'))
							Jobs[jobNumber]['RUNTYPE'] = runProtocol
						
						
					if 'RUNTYPE' not in Jobs[jobNumber]:
						runProtocol = ['CREATE']
						runProtocol.extend(['PROD']*100)
						Jobs[jobNumber]['RUNTYPE'] =  runProtocol
					
					continue
				
				
				# Job protocol. associates each parameter with an index
				if(row_name.startswith('COPY') or row_name.startswith('INIT') or row_name.startswith('EM') or row_name.startswith('NVE') or row_name.startswith('NVT') or row_name.startswith('NPT')):
					
					stepNb = str(stepNumber)
					#Creates a dictionary for the step with an associated protocol step number to keep track of the order of execution
					Jobs[jobNumber]['PROTOCOL'].update( { stepNb: {'stepType' : row[0].strip(' ')}  })

					#Associates the parameters set in Parameters.csv to the correct protocol step dictionary
					for index, param in enumerate(row[1:]):
						if(param):
							Jobs[jobNumber]['PROTOCOL'][stepNb].update( { param.strip(' ') : index + 1 })

					#Increases the protocol step number
					stepNumber += 1
					
					continue
				
				#Job protocol. Associates each parameter with its value based on the index
				if(row_name.strip() == ''):
					stepNb = str(stepNumber-1)
					for param, paramValue in Jobs[jobNumber]['PROTOCOL'][stepNb].items():
						if param != 'stepType':
							paramValue = row[paramValue]
							
							if(paramValue):
								if ut.is_number(paramValue):
									if '.' in paramValue or not float(paramValue).is_integer() :
										Jobs[jobNumber]['PROTOCOL'][stepNb][param] = float(paramValue)
									else:
										Jobs[jobNumber]['PROTOCOL'][stepNb][param] = int(paramValue)
								else:
									Jobs[jobNumber]['PROTOCOL'][stepNb][param] = paramValue.strip(' ')
				
				if row_name.strip() == 'END_OF_PROJECT':
					break

	#**********************************************#
	#**********************************************#
	#**********************************************#
	#**********************************************#
	#**********************************************#
	#**********************************************#
	
	#Finding the software version. If not found ask the user.

	print('Gromacs local :', Softwares['GROMACS_LOC'])

	software_version = re.sub("[^0-9]","",Softwares['GROMACS_LOC'].split('/')[-2])

	if not software_version:
		software_version = input('No version detected in path. Please select the version you want between 4, 5 or 2016: ')

	print('Gromacs version :', software_version)

	# Checking that the path to default files is provided. If not ask the user for a path.
	if not path_to_default:
		path_to_default = input('Path to default gromacs files was not found. Please enter the path to the GROMACS default files directory: ')
	print('Gromacs Default :', path_to_default)

	#Variables to set the correct path for gromacs
	if software_version.startswith('4') :
		Softwares['GROMACS_REM'] = Softwares['GROMACS_REM']+'/'
		Softwares['GROMACS_LOC'] = Softwares['GROMACS_LOC']+'/'
	else:
		Softwares['GROMACS_REM'] = Softwares['GROMACS_REM']+'/gmx '
		Softwares['GROMACS_LOC'] = Softwares['GROMACS_LOC']+'/gmx '
	#Creates the base directory
	if(not os.path.isdir(project_name)):
		os.makedirs(project_name, exist_ok=True)
	else:
		shutil.rmtree(project_name)
		os.makedirs(project_name, exist_ok=True)

	#Dictionnary to perform sample copy
	Init_for_copy={}


	#**********************************************#
	#*****      Generating MDrun files        *****#
	#**********************************************#
	for job_number in range(1, len(Jobs)+1):
		current_job = Jobs[job_number]
		name = current_job['JOBID']
		print("===============\n===============\nSample "+name+"\n===============\n===============\n")
		
		#Creates the sample directory
		directory = "{0}/{1}".format(project_name, name)
		if not os.path.isdir(directory):
			os.mkdir(directory)

		script_file_create = open("{0}/create.sh".format(directory),'w')
		script_file_create.truncate()
		
		script_file_prod = open("{0}/run.sh".format(directory),'w')
		script_file_prod.truncate()
		
		create_to = ""
		prefix_gromacs_mdrun_create = ""
		prefix_gromacs_grompp_create = ""
		
		copy_to = ""
		prefix_gromacs_mdrun_prod = ""
		prefix_gromacs_grompp_prod = ""
		
		
		if cmd_param.send == 'local':
			copy_to = """
							#!/bin/bash +x
							mkdir -p ../{0}_OUTPUT
							cp * ../{0}_OUTPUT/
							cd ../{0}_OUTPUT/

							""".format(name)
							
			prefix_gromacs_grompp_prod = Softwares['GROMACS_LOC']
			prefix_gromacs_mdrun_prod = Softwares['GROMACS_LOC']
			
		elif cmd_param.send == 'pbs':
			copy_to = """
						#!/bin/bash +x
						mkdir -p /scratch/{2}/gromacs/{0}
						cd ..
						LOCALDIR="$(pwd)"
						cp -r {1} /scratch/{2}/gromacs/{0}
						cd /scratch/{2}/gromacs/{0}/{1}

						""".format(project_name, name, PBS['username'] )
			
			prefix_gromacs_grompp_prod = Softwares['GROMACS_REM']
			prefix_gromacs_mdrun_prod = Softwares['GROMACS_REM']
			
		elif cmd_param.send == 'tgcc':
			copy_to = """
						#!/bin/bash +x
						
						# ADDAPT AUTOMATICALLY TO GROMACS VERSION ?
						# OR AT LEAST VERIFY THAT THE TWO VERSION INFORMATION
						# IN {4} and TGCCinfo.csv ARE CONSISTENT  ?
						
						module load {3}
						
						""".format(project_name, name, TGCC['username'], TGCC['GMXversion'], parameter_file)
						
			prefix_gromacs_grompp_prod = "gmx_mpi "
			prefix_gromacs_mdrun_prod = "ccc_mprun gmx_mpi "
			
		# Samples creation #
		if cmd_param.create == 'local':
			create_to = """
							#!/bin/bash +x

							""".format(name)
							
			prefix_gromacs_grompp_create = Softwares['GROMACS_LOC']
			prefix_gromacs_mdrun_create = Softwares['GROMACS_LOC']
			
		elif cmd_param.create == 'pbs':
			create_to = """
						#!/bin/bash +x
						mkdir -p /scratch/{2}/gromacs/{0}
						cd ..
						LOCALDIR="$(pwd)"
						cp -r {1} /scratch/{2}/gromacs/{0}
						cd /scratch/{2}/gromacs/{0}/{1}

						""".format(project_name, name, PBS['username'] )
			
			prefix_gromacs_grompp_create = Softwares['GROMACS_REM']
			prefix_gromacs_mdrun_create = Softwares['GROMACS_REM']
			
		elif cmd_param.create == 'tgcc':
			create_to = """
						#!/bin/bash +x
						
						# ADDAPT AUTOMATICALLY TO GROMACS VERSION ?
						# OR AT LEAST VERIFY THAT THE TWO VERSION INFORMATION
						# IN {4} and TGCCinfo.csv ARE CONSISTENT  ?
						
						module load {3}
						
						""".format(project_name, name, TGCC['username'], TGCC['GMXversion'], parameter_file)
						
			prefix_gromacs_grompp_create = "gmx_mpi "
			prefix_gromacs_mdrun_create = "ccc_mprun gmx_mpi "
		
		script_file_create.write(ut.RemoveUnwantedIndent(create_to))
		script_file_prod.write(ut.RemoveUnwantedIndent(copy_to))
		
		previous_cmd_files={}
		
		#**********************************************#
		# Read the default files for EM, NVT, NVP, NVE #
		#**********************************************#
		#**********************************************#
		# Lipids files #
		if( current_job['TYPE'] == 'BILAYER' or current_job['TYPE'] == 'TRILAYER'):
			if software_version.startswith('4') :
				EMdefault = open(path_to_default+'/mdp_lipids/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(path_to_default+'/mdp_lipids/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(path_to_default+'/mdp_lipids/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(path_to_default+'/mdp_lipids/GROMACS4/NVEdefault.mdp','r')
			
			else:
				EMdefault = open(path_to_default+'/mdp_lipids/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(path_to_default+'/mdp_lipids/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(path_to_default+'/mdp_lipids/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(path_to_default+'/mdp_lipids/GROMACS2016/NVEdefault.mdp','r')
		
		# Solvent files #
		if( current_job['TYPE'] == 'SOLVENT'):
			if software_version.startswith('4'):
				EMdefault = open(path_to_default+'/mdp_sol/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(path_to_default+'/mdp_sol/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(path_to_default+'/mdp_sol/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(path_to_default+'/mdp_sol/GROMACS4/NVEdefault.mdp','r')

				
			else:
				EMdefault = open(path_to_default+'/mdp_sol/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(path_to_default+'/mdp_sol/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(path_to_default+'/mdp_sol/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(path_to_default+'/mdp_sol/GROMACS2016/NVEdefault.mdp','r')
		#**********************************************#
		#**********************************************#
		
		
		
		Sample = None
		#Go into the corresponding directory
		with ut.cd(project_name+'/'+name):
			output_filename=''
			script_file_prod = open('run.sh','a+')
			script_file_prod.write('mkdir -p mdrun_out\n')
			script_file_prod.write('mkdir -p grompp_out\n')
			
			script_file_create = open('create.sh','a+')
			script_file_create.write('mkdir -p create_mdrun_out\n')
			script_file_create.write('mkdir -p create_grompp_out\n')
			
			so_far_file = """{0}_SoFar.txt""".format(name)
			so_far_file = open(so_far_file, 'w')
			
			#Executes the steps
			for md_step in range(0, len(current_job['PROTOCOL'])):
				
				md_step = str(md_step)
				
				
				#============================================================
				#					STEP IS INITIALIZATION
				#============================================================
				if current_job['PROTOCOL'][md_step]['stepType'].startswith('INIT'):
					if current_job['TYPE'] in ['BILAYER','TRILAYER']:
						Sample = Prepare.Membrane(current_job, Softwares, path_to_default)
						
						
					if current_job['TYPE'] in ['SOLVENT']:
						Sample = Prepare.Solvent(current_job, Softwares, path_to_default)
						
					
					Sample.make()
					
					so_far_file.write(str("""{0} done \n""").format(current_job['PROTOCOL'][md_step]['stepType']))
					
					previous_cmd_files = {'SYSTEM': None, 'OUTPUT': None, 'INDEX': None}
					previous_cmd_files['SYSTEM'], previous_cmd_files['OUTPUT'], previous_cmd_files['INDEX'] = Sample.pass_outputs()
				
					continue
				
				#============================================================
				#					STEP IS INPUT (need .gro, .ndx, .top, .itp, .tpr)
				#============================================================
				elif current_job['PROTOCOL'][md_step]['stepType'].startswith('INPUT'):
					""" This step will copy the file you put in Parameters.csv.
						If a defo is in the file you are copying you need to specify the version using:
						DEFO, Version, version_number, ...
						You can also change the bonds parameter and Mass if wanted.
						the protocol will be define in the same way as INPUT
					"""
					
					if current_job['TYPE'] in ['BILAYER','TRILAYER']:
						Sample = Prepare.Membrane(current_job, Softwares, path_to_default)
						
					if current_job['TYPE'] in ['SOLVENT']:
						Sample = Prepare.Solvent(current_job, Softwares, path_to_default)
					
					previous_cmd_files = {'SYSTEM': None, 'OUTPUT': None, 'INDEX': None}
					previous_cmd_files['SYSTEM'], previous_cmd_files['OUTPUT'], previous_cmd_files['INDEX'] = Sample.from_input(current_job['PROTOCOL'][md_step])
				
					continue
				
				#============================================================
				#					STEP IS SAMPLE COPY
				#============================================================
				
				elif current_job['PROTOCOL'][md_step]['stepType'].startswith('COPY') :
					assert(0), "Copy not yet implemented"
					"""
					CopyFrom = Jobs[ Jobs[job_number]['PROTOCOL'][step]['samplenumber'] ]['JOBID'] 
					try:
						CopiedJob = Prepare.CopySample(Jobs, current_job, md_step, path_to_default)
						if CopiedJob == -1:
							raise ValueError("Please Check your "+parameter_file)
					except ValueError as error:
						print(error)
						sys.exit(1)
						
					so_far_file.write( "{0} {1} done\n".format( current_job['PROTOCOL'][md_step]['stepType'], CopyFrom) )
					previous_cmd_files = copy.deepcopy(CopiedJob)
					
					continue
					"""
				#============================================================
				#					STEP IS TRANSLATE
				#============================================================
				elif current_job['PROTOCOL'][md_step]['stepType'].startswith('TRANSLATE'):
					TX = 0.1*current_job['PROTOCOL'][md_step]['X']
					TY = 0.1*current_job['PROTOCOL'][md_step]['Y']
					TZ = 0.1*current_job['PROTOCOL'][md_step]['Z']
					
					output_file = previous_cmd_files['OUTPUT'].replace('.gro','_TRANSX{0}Y{1}Z{2}.gro'.format(TX, TY, TZ) )
					
					prefix_gromacs_grompp = None
					prefix_gromacs_mdrun = None
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							prefix_gromacs_grompp = prefix_gromacs_grompp_create
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_create
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							prefix_gromacs_grompp = prefix_gromacs_grompp_prod
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
					else:
						prefix_gromacs_grompp = prefix_gromacs_grompp_prod
						prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
						
					assert(prefix_gromacs_grompp is not None or prefix_gromacs_mdrun is not None), str("There was a problem in translating the "
																								"sample (prefix for gromacs)")
						
					cmd = "    printf {0} | {1} trjconv -f {2} -s {2} -o {3} -trans {4} {5} {6} -pbc atom\n\n".format(repr("System\n"), 
																													prefix_gromacs_grompp,
																													previous_cmd_files['OUTPUT'],
																													output_file,
																													TX, TY, TZ)
					
					sofar_cmd = "echo 'Translating the sample done' >> {0}_SoFar.txt\n\n\n".format(name)
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							script_file_create.write(cmd)
							script_file_create.write(sofar_cmd)
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							script_file_prod.write(cmd)
							script_file_prod.write(sofar_cmd)
					else:
						script_file_prod.write(cmd)
						script_file_prod.write(sofar_cmd)
						
					previous_cmd_files['OUTPUT'] = output_file
					
					continue
				
				#============================================================
				#					STEP IS ADDING SUBSTRATE WITH GROMACS
				#============================================================
				elif current_job['PROTOCOL'][md_step]['stepType'].startswith('GADD_SU'):
					file_output = None
					file_index = None
					system = None
					
					prefix_gromacs_grompp = None
					prefix_gromacs_mdrun = None
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							prefix_gromacs_grompp = prefix_gromacs_grompp_create
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_create
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							prefix_gromacs_grompp = prefix_gromacs_grompp_prod
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
					else:
						prefix_gromacs_grompp = prefix_gromacs_grompp_prod
						prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
						
					assert(prefix_gromacs_grompp is not None or prefix_gromacs_mdrun is not None), str("There was a problem in adding the "
																								"substrate (prefix for gromacs)")
						
					cmd,  file_output, file_index, system = Sample.add_substrate_in_run(current_job, md_step,
																							previous_cmd_files, 
																							prefix_gromacs_grompp,
																							prefix_gromacs_mdrun)
					
					sofar_cmd = "echo 'Inserting substrate done' >> {0}_SoFar.txt\n\n\n".format(name)
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							script_file_create.write(cmd)
							script_file_create.write(sofar_cmd)
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							script_file_prod.write(cmd)
							script_file_prod.write(sofar_cmd)
					else:
						script_file_prod.write(cmd)
						script_file_prod.write(sofar_cmd)
					
					assert(file_output is not None or file_index is not None or system is not None), "There was a problem in adding the substrate"
					
					previous_cmd_files['OUTPUT'] = file_output
					previous_cmd_files['INDEX'] = file_index
					previous_cmd_files['SYSTEM'] = system
					continue
				
				#============================================================
				#					STEP IS ADDING WALL
				#============================================================
				elif current_job['PROTOCOL'][md_step]['stepType'].startswith('ADD_WALL') :
					file_output = None
					system = None
					
					prefix_gromacs_grompp = None
					prefix_gromacs_mdrun = None
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							prefix_gromacs_grompp = prefix_gromacs_grompp_create
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_create
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							prefix_gromacs_grompp = prefix_gromacs_grompp_prod
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
					else:
						prefix_gromacs_grompp = prefix_gromacs_grompp_prod
						prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
						
					assert(prefix_gromacs_grompp is not None or prefix_gromacs_mdrun is not None), "There was a problem in adding the wall (prefix for gromacs)"
						
					cmd, system , file_output = Sample.add_wall_in_run(current_job, md_step, previous_cmd_files, 
																		prefix_gromacs_grompp, prefix_gromacs_mdrun)
					
					sofar_cmd = "echo 'Adding wall done' >> {0}_SoFar.txt\n\n\n".format(name)
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							script_file_create.write(cmd)
							script_file_create.write(sofar_cmd)
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							script_file_prod.write(cmd)
							script_file_prod.write(sofar_cmd)
					else:
						script_file_prod.write(cmd)
						script_file_prod.write(sofar_cmd)
					
					assert(file_output is not None or system is not None), "There was a problem in adding the wall"
					
					previous_cmd_files['OUTPUT'] = file_output
					previous_cmd_files['SYSTEM'] = system
					continue
				#============================================================
				#						STEP IS EM
				#============================================================
				elif current_job['PROTOCOL'][md_step]['stepType'].startswith('EM'):
					prefix_gromacs_grompp = None
					prefix_gromacs_mdrun = None
					grompp_out = None
					mdrun_out = None
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							prefix_gromacs_grompp = prefix_gromacs_grompp_create
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_create
							grompp_out = "create_grompp_out"
							mdrun_out = "create_mdrun_out"
							
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							prefix_gromacs_grompp = prefix_gromacs_grompp_prod
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
							grompp_out = "grompp_out"
							mdrun_out = "mdrun_out"
					else:
						prefix_gromacs_grompp = prefix_gromacs_grompp_prod
						prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
						grompp_out = "grompp_out"
						mdrun_out = "mdrun_out"
						
					assert(prefix_gromacs_grompp is not None or prefix_gromacs_mdrun is not None), str("There was a problem in EM prep "
																								"(prefix for gromacs)")
					
					#Output for Pulling data
					pulling_output = ''
					if 'DEFO' in current_job:
						pulling_output = "-pf {0}_pf.xvg -px {0}_px.xvg".format( current_job['PROTOCOL'][md_step]['stepType'])
					
					Sample.preparing_mdp(md_step)
					Sample.writing_mdp(md_step, EMdefault)

					cmd = str("#Energy Minimization using parameters in {0}.mdp:\n#Input : "
								"{1}\n#Output: {2}_{3}-{0}.gro \n\n").format(current_job['PROTOCOL'][md_step]['stepType'],
																				previous_cmd_files['OUTPUT'],previous_cmd_files['SYSTEM'],
																				current_job['JOBNUM'])

					cmd += str("{0}grompp -f {1}.mdp -po {2}-{1}.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} "
								"|& tee {6}/grompp_{5}_{1}.output\n\n").format(prefix_gromacs_grompp,
																					current_job['PROTOCOL'][md_step]['stepType'],
																					previous_cmd_files['SYSTEM'],previous_cmd_files['OUTPUT'], 
																					previous_cmd_files['INDEX'],current_job['JOBNUM'],
																					grompp_out)

					cmd += str("{0}mdrun -deffnm {2}_{3}-{1} -c {2}_{3}-{1}.gro {4}"
								" |& tee {5}/mdrun_{3}_{1}.output \n\n").format(prefix_gromacs_mdrun,
																				current_job['PROTOCOL'][md_step]['stepType'],
																				previous_cmd_files['SYSTEM'], current_job['JOBNUM'],
																				pulling_output, mdrun_out)
								
					
					
					cmd += "echo '{0} done' >> {1}_SoFar.txt\n\n\n".format(current_job['PROTOCOL'][md_step]['stepType'], name)

					output_filename = "{1}_{2}-{0}.gro".format(current_job['PROTOCOL'][md_step]['stepType'], previous_cmd_files['SYSTEM'], current_job['JOBNUM'])
					
					
					cmd += str("### Correcting the box vectors if LZ is too large ###\n"
								"if [ -a {0} ]; then read -r LX LY LZ <<< $(tail -n1 {0}); "
								"""if [ "$LZ" == "" ]; then """
								"""echo "LZ was too big, correcting the box vectors format"; """
								"LZ=$(cut -c 9- <<< $LY); LY=$(cut -c -8 <<< $LY); "
								"sed -i '$ d' {0}; "
								"""echo "$LX $LY $LZ" >> {0}; """
								"fi;fi\n\n\n").format(output_filename)
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							script_file_create.write(cmd)
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							script_file_prod.write(cmd)
					else:
						script_file_prod.write(cmd)

					previous_cmd_files['OUTPUT'] = output_filename
					
					continue

				#============================================================
				#						STEP IS NVE, NVT or NPT
				#============================================================
				else:
					prefix_gromacs_grompp = None
					prefix_gromacs_mdrun = None
					grompp_out = None
					mdrun_out = None
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							prefix_gromacs_grompp = prefix_gromacs_grompp_create
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_create
							grompp_out = "create_grompp_out"
							mdrun_out = "create_mdrun_out"
							
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							prefix_gromacs_grompp = prefix_gromacs_grompp_prod
							prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
							grompp_out = "grompp_out"
							mdrun_out = "mdrun_out"
					else:
						prefix_gromacs_grompp = prefix_gromacs_grompp_prod
						prefix_gromacs_mdrun = prefix_gromacs_mdrun_prod
						grompp_out = "grompp_out"
						mdrun_out = "mdrun_out"
						
					assert(prefix_gromacs_grompp is not None or prefix_gromacs_mdrun is not None), str("There was a problem in NVE, NVT, NPT prep "
																								"(prefix for gromacs)")
					
					cmd = ""
					
					if current_job['PROTOCOL'][md_step]['stepType'].startswith('NVE'):
						cmd = str("""#NVE using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}.gro \n\n""").format(current_job['PROTOCOL'][md_step]['stepType'],
																														previous_cmd_files['OUTPUT'], previous_cmd_files['SYSTEM'], current_job['JOBNUM'])
						Sample.preparing_mdp(md_step)
						Sample.writing_mdp(md_step, NVEdefault)
						
					elif current_job['PROTOCOL'][md_step]['stepType'].startswith('NVT'):
						cmd = str("""#NVT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}.gro \n\n""").format(current_job['PROTOCOL'][md_step]['stepType'],
																														previous_cmd_files['OUTPUT'], previous_cmd_files['SYSTEM'], current_job['JOBNUM'])
						Sample.preparing_mdp(md_step)
						Sample.writing_mdp(md_step, NVTdefault)
						
					
					elif current_job['PROTOCOL'][md_step]['stepType'].startswith('NPT'):
						cmd = str("""#NPT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}.gro \n\n""").format(current_job['PROTOCOL'][md_step]['stepType'],
																														previous_cmd_files['OUTPUT'], previous_cmd_files['SYSTEM'], current_job['JOBNUM'])
						Sample.preparing_mdp(md_step)
						Sample.writing_mdp(md_step, NPTdefault)
					
					
					
					#Output for Pulling data
					pulling_output = ''
					if 'DEFO' in current_job:
						pulling_output = "-pf {0}_pf.xvg -px {0}_px.xvg".format( current_job['PROTOCOL'][md_step]['stepType'])

					cmd += str("{0}grompp -f {1}.mdp -po {2}-{1}.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} "
								"|& tee {6}/grompp_{5}_{1}.output\n\n""").format(prefix_gromacs_grompp,
																						current_job['PROTOCOL'][md_step]['stepType'], 
																						previous_cmd_files['SYSTEM'], previous_cmd_files['OUTPUT'], previous_cmd_files['INDEX'],
																						current_job['JOBNUM'],
																						grompp_out)

					cmd += str("{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}.gro {5} "
								"|& tee {6}/mdrun_{3}_{1}.output \n\n").format(prefix_gromacs_mdrun, 
																					current_job['PROTOCOL'][md_step]['stepType'], 
																					previous_cmd_files['SYSTEM'], current_job['JOBNUM'], current_job['MDRUN_OPT'],
																					pulling_output,
																					mdrun_out)

					cmd += "echo '{0} done' >> {1}_SoFar.txt\n\n\n".format(current_job['PROTOCOL'][md_step]['stepType'],name)
					

					output_filename = "{1}_{2}-{0}.gro".format(current_job['PROTOCOL'][md_step]['stepType'], previous_cmd_files['SYSTEM'], current_job['JOBNUM'])
					
					cmd += str("### Correcting the box vectors if LZ is too large ###\n"
								"if [ -a {0} ]; then read -r LX LY LZ <<< $(tail -n1 {0}); "
								"""if [ "$LZ" == "" ]; then """
								"""echo "LZ was too big, correcting the box vectors format"; """
								"LZ=$(cut -c 9- <<< $LY); LY=$(cut -c -8 <<< $LY); "
								"sed -i '$ d' {0}; "
								"""echo "$LX $LY $LZ" >> {0}; """
								"fi;fi\n\n\n").format(output_filename)
					
					if 'RUNTYPE' in current_job:
						if current_job['RUNTYPE'][int(md_step)] == 'CREATE':
							script_file_create.write(cmd)
							
						elif current_job['RUNTYPE'][int(md_step)] == 'PROD':
							script_file_prod.write(cmd)
					else:
						script_file_prod.write(cmd)
					
					previous_cmd_files['OUTPUT'] = output_filename
					
					continue
			
			end_of_run = ""
			end_of_create = ""
			if cmd_param.send== 'local':
				end_of_run = """
							echo "End of run for {1}"
							""".format(project_name, name)
			
			elif cmd_param.send == 'pbs':
				end_of_run = """
							echo "End of run for {1}"
							cp -r /scratch/{2}/gromacs/{0}/{1} ${{LOCALDIR}}/{1}_OUTPUT
							rm -r /scratch/{2}/gromacs/{0}
							""".format(project_name, name, PBS['username'])
				pbs_name = "{0}_{1}_{2}".format(software_version, project_name, name)
				pbs_file_content = """
						#! /bin/sh -x
						#PBS -N {0}
						#PBS -j oe
						#PBS -m abe
						#PBS -M {1}
						#PBS -q {2}
						#PBS -l nodes={3}:ppn={4}
						#PBS -l walltime=10000:0:0

						NPROCS=`cat $PBS_NODEFILE | wc -l`
						NODES=`uniq $PBS_NODEFILE | wc -l`
						JOBINFO="$PBS_O_WORKDIR/$PBS_JOBNAME-$PBS_JOBID.jobinfo"


						printf "Time =  `date`\\n" > $JOBINFO
						printf "PBS work directory = $PBS_O_WORKDIR\\n" >> $JOBINFO
						printf "PBS queue = $PBS_O_QUEUE\\n" >> $JOBINFO
						printf "PBS job ID = $PBS_JOBID\\n" >> $JOBINFO
						printf "PBS job name = $PBS_JOBNAME\\n" >> $JOBINFO
						printf "This job will run on $NPROCS processors\\n" >> $JOBINFO
						printf "List of nodes in $PBS_NODEFILE\\n" >> $JOBINFO
						uniq $PBS_NODEFILE >> $JOBINFO
						NODES=`uniq $PBS_NODEFILE`

						cd $PBS_O_WORKDIR
						echo  "files in $PBS_O_WORKDIR"
						ls -ltr
						echo "============================="

						echo "Run GMX"
						date

						./run.sh  $PBS_JOBID >>  $PBS_O_WORKDIR/$PBS_JOBNAME-$PBS_JOBID.out

						echo "End of GMX"
						date

						echo "============================="
						""".format(pbs_name, PBS['mail'], current_job['NODE'], current_job['NBNODES'], current_job['PPN'])
								
				with open(name+'.pbs','w') as pbs_file:
					pbs_file.write( ut.RemoveUnwantedIndent(pbs_file_content) )
							
			elif cmd_param.send == 'tgcc':
				
				tgcc_file_content = """
						#! /bin/sh -x
						#MSUB -q {2}      # queue = standard, test, long (in {6}, key "node")
						#MSUB -n {3}      # total number of cores, 16 cores per nodes (32 logical cores ) (in Parameters.cvs, key "ppn")
						#MSUB -T {5}      # times in seconds (in TGCCinfo.csv, key "time_s")
						#MSUB -r {0}	  # job name (automatically generated)
						#MSUB -A {4} 	  # group for allocation (gen7662) (in TgccInfo.csv, key "group")
						#MSUB -V          # transfer the variable of ENVIRONNEMENT
						#MSUB -@ {7}:begin,end
						#MSUB -oe %I.eo   # input and output of JOB
						
						# the number of nodes is deduced automatically from the number of cores (fixed nb core/node = 16 on curie noeud fin)
						# ccc_mprun gives this information to gmx_mpi (seems to work !?)
						
						echo "===================== BEGIN JOB $SLURM_JOBID =============================== "
						
						JOBINFO="${{SLURM_SUBMIT_DIR}}/${{SLURM_JOB_NAME}}-${{SLURM_JOBID}}.jobinfo"
						printf "Time =  `date`\\n" > ${{JOBINFO}}
						printf "SLURM submit directory = ${{SLURM_SUBMIT_DIR}}\\n" >> ${{JOBINFO}}
						printf "TGCC queue = {2}, user {4}, max time = {5} seconds \\n" >> ${{JOBINFO}}
						printf "SLURM job ID = ${{SLURM_JOBID}} \\n" >> ${{JOBINFO}}
						printf "SLURM job name = ${{SLURM_JOB_NAME}} \\n" >> ${{JOBINFO}}
						printf "This job will run on {3} processors\\n" >> ${{JOBINFO}}
						printf "List of nodes : ${{SLURM_NODEID}} \\n\\n" >> ${{JOBINFO}}
						
						
						export MYTMPDIR="${{SCRATCHDIR}}/JOB_${{SLURM_JOBID}}"
						export OUTPUTDIR="${{SLURM_SUBMIT_DIR}}/JOB_${{SLURM_JOBID}}_OUTPUT"
						
						mkdir -p ${{MYTMPDIR}}
						mkdir -p ${{OUTPUTDIR}}
						
						rsync ${{SLURM_SUBMIT_DIR}}/* ${{MYTMPDIR}} 
						cd ${{MYTMPDIR}} 
						./run.sh  >> ${{OUTPUTDIR}}/JOB_${{SLURM_JOBID}}.out
						rsync -r ./* ${{OUTPUTDIR}}/.
						cd ${{SLURM_SUBMIT_DIR}}
						rm -rf ${{MYTMPDIR}}
						
						echo "===================== END  JOB $SLURM_JOBID =============================== "
						""".format(name, TGCC['mail'], TGCC['queue'], current_job['PPN'],TGCC['group'], TGCC['time_s'], parameter_file, TGCC['mail'])
				
				with open(name+'.ccc_msub','w') as tgcc_file:
					tgcc_file.write( ut.RemoveUnwantedIndent(tgcc_file_content) )
					
			
			if cmd_param.create == 'local':
				end_of_create = """
							echo "End of create for {1}"
							""".format(project_name, name)
			
			elif cmd_param.create == 'pbs':
				end_of_create = """
							echo "End of create for {1}"
							cp -r /scratch/{2}/gromacs/{0}/{1} ${{LOCALDIR}}/{1}
							rm -r /scratch/{2}/gromacs/{0}
							""".format(project_name, name, PBS['username'])
				
				pbs_name = "{0}_{1}_{2}".format(software_version, project_name, name)
				pbs_file_content = """
						#! /bin/sh -x
						#PBS -N {0}
						#PBS -j oe
						#PBS -m abe
						#PBS -M {1}
						#PBS -q {2}
						#PBS -l nodes={3}:ppn={4}
						#PBS -l walltime=10000:0:0

						NPROCS=`cat $PBS_NODEFILE | wc -l`
						NODES=`uniq $PBS_NODEFILE | wc -l`
						JOBINFO="$PBS_O_WORKDIR/$PBS_JOBNAME-$PBS_JOBID.jobinfo"


						printf "Time =  `date`\\n" > $JOBINFO
						printf "PBS work directory = $PBS_O_WORKDIR\\n" >> $JOBINFO
						printf "PBS queue = $PBS_O_QUEUE\\n" >> $JOBINFO
						printf "PBS job ID = $PBS_JOBID\\n" >> $JOBINFO
						printf "PBS job name = $PBS_JOBNAME\\n" >> $JOBINFO
						printf "This job will run on $NPROCS processors\\n" >> $JOBINFO
						printf "List of nodes in $PBS_NODEFILE\\n" >> $JOBINFO
						uniq $PBS_NODEFILE >> $JOBINFO
						NODES=`uniq $PBS_NODEFILE`

						cd $PBS_O_WORKDIR
						echo  "files in $PBS_O_WORKDIR"
						ls -ltr
						echo "============================="

						echo "Run GMX"
						date

						./run.sh  $PBS_JOBID >>  $PBS_O_WORKDIR/$PBS_JOBNAME-$PBS_JOBID.out

						echo "End of GMX"
						date

						echo "============================="
						""".format(pbs_name, PBS['mail'], current_job['NODE'], current_job['NBNODES'], current_job['PPN'])
								
				with open(name+'.create_pbs','w') as pbs_file:
					pbs_file.write( ut.RemoveUnwantedIndent(pbs_file_content) )
							
			elif cmd_param.create == 'tgcc':
				
				tgcc_file_content = """
						#! /bin/sh -x
						#MSUB -q {2}      # queue = standard, test, long (in {6}, key "node")
						#MSUB -n {3}      # total number of cores, 16 cores per nodes (32 logical cores ) (in Parameters.cvs, key "ppn")
						#MSUB -T {5}      # times in seconds (in TGCCinfo.csv, key "time_s")
						#MSUB -r {0}_created	  # job name (automatically generated)
						#MSUB -A {4} 	  # group for allocation (gen7662) (in TgccInfo.csv, key "group")
						#MSUB -V          # transfer the variable of ENVIRONNEMENT
						#MSUB -@ {7}:begin,end
						#MSUB -oe %I.eo   # input and output of JOB
						
						# the number of nodes is deduced automatically from the number of cores (fixed nb core/node = 16 on curie noeud fin)
						# ccc_mprun gives this information to gmx_mpi (seems to work !?)
						
						echo "===================== BEGIN JOB $SLURM_JOBID =============================== "
						
						JOBINFO="${{SLURM_SUBMIT_DIR}}/${{SLURM_JOB_NAME}}-${{SLURM_JOBID}}.jobinfo"
						printf "Time =  `date`\\n" > ${{JOBINFO}}
						printf "SLURM submit directory = ${{SLURM_SUBMIT_DIR}}\\n" >> ${{JOBINFO}}
						printf "TGCC queue = {2}, user {4}, max time = {5} seconds \\n" >> ${{JOBINFO}}
						printf "SLURM job ID = ${{SLURM_JOBID}} \\n" >> ${{JOBINFO}}
						printf "SLURM job name = ${{SLURM_JOB_NAME}} \\n" >> ${{JOBINFO}}
						printf "This job will run on {3} processors\\n" >> ${{JOBINFO}}
						printf "List of nodes : ${{SLURM_NODEID}} \\n\\n" >> ${{JOBINFO}}
						
						
						export MYTMPDIR="${{SCRATCHDIR}}/JOB_${{SLURM_JOBID}}"
						export OUTPUTDIR="${{SLURM_SUBMIT_DIR}}"
						
						mkdir -p ${{MYTMPDIR}}
						
						rsync ${{SLURM_SUBMIT_DIR}}/* ${{MYTMPDIR}} 
						cd ${{MYTMPDIR}} 
						./create.sh  >> ${{OUTPUTDIR}}/JOB_${{SLURM_JOBID}}_create.out
						rsync -r ./* ${{OUTPUTDIR}}/.
						cd ${{SLURM_SUBMIT_DIR}}
						rm -rf ${{MYTMPDIR}}
						ccc_msub {0}.ccc_msub
						
						echo "===================== END  JOB $SLURM_JOBID =============================== "
						""".format(name, TGCC['mail'], TGCC['queue'], current_job['PPN'],TGCC['group'], TGCC['time_s'], parameter_file, TGCC['mail'])
				
				with open(name+'.create_ccc_msub','w') as tgcc_file:
					tgcc_file.write( ut.RemoveUnwantedIndent(tgcc_file_content) )
			
			
			
			
			
					
			script_file_create.write(ut.RemoveUnwantedIndent(end_of_create))
			script_file_prod.write(ut.RemoveUnwantedIndent(end_of_run))

			sub.call("""chmod a+x run.sh ; if [ -e create.sh ] ; then chmod a+x create.sh;fi """,shell=True)
			so_far_file.close()
			
		script_file_prod.close()
		script_file_create.close()
		EMdefault.close()
		NVTdefault.close()
		NPTdefault.close()
		NVEdefault.close()

	#**********************************************#
	# Creating scripts for starting jobs  *****#
	#**********************************************#
	with ut.cd(project_name):
		#Creates a script for sending all samples to queue
		script_to_send = open('send_jobs.sh','w')
		script_to_send.truncate()
		script_to_send.write("""#!/bin/bash\n\n""")
		script_to_send.close()
		
		script_to_create = open('create_jobs.sh','w')
		script_to_create.truncate()
		script_to_create.write("""#!/bin/bash\n\n""")
		script_to_create.close()

		script_to_send = open('send_jobs.sh','a+')
		script_to_create = open('create_jobs.sh','a+')

		for job_number in Jobs:
			current_job = Jobs[job_number]
			name = current_job['JOBID']
			
			send_cmd = ""
			create_cmd = ""
			if cmd_param.send == 'local':
				send_cmd = "cd {0} \n ./run.sh \ncd .. \n".format(name)
				
			elif cmd_param.send == 'pbs':
				send_cmd = "cd {0} \n qsub {0}.pbs \n cd .. \n".format(name)
				
			elif cmd_param.send == 'tgcc':
				send_cmd = "cd {0} \n ccc_msub {0}.ccc_msub \n cd .. \n".format(name)
			
			if cmd_param.create == 'local':
				create_cmd = "cd {0} \n ./create.sh \ncd .. \n".format(name)
				
			elif cmd_param.create == 'pbs':
				create_cmd = "cd {0} \n qsub {0}.create_pbs \n cd .. \n".format(name)
				
			elif cmd_param.create == 'tgcc':
				create_cmd = "cd {0} \n ccc_msub {0}.create_ccc_msub \n cd .. \n".format(name)
				
			script_to_send.write(send_cmd)
			script_to_create.write(create_cmd)
		
		see_jobs_cmd = ""
		see_create_cmd = ""
		if cmd_param.send == 'pbs':
			see_jobs_cmd = "qstat -n -u {0}\necho ' To follow your JOBS, type qstat -n -u {0} ' \n".format(PBS['username'])
		
		elif cmd_param.send == 'tgcc':
			see_jobs_cmd = "ccc_mpp -u {0}\n\necho ' To follow your JOBS, type ccc_mpp -u {0} ' \n".format(TGCC['username'])
		
		if cmd_param.create == 'pbs':
			see_create_cmd = "qstat -n -u {0}\necho ' To follow your JOBS, type qstat -n -u {0} ' \n".format(PBS['username'])
		
		elif cmd_param.create == 'tgcc':
			see_create_cmd = "ccc_mpp -u {0}\n\necho ' To follow your JOBS, type ccc_mpp -u {0} ' \n".format(TGCC['username'])
		
		script_to_send.write(see_jobs_cmd)
		script_to_create.write(see_create_cmd)
		script_to_send.close()
		script_to_create.close()
		sub.call("""chmod a+x send_jobs.sh; chmod a+x create_jobs.sh""", shell=True)

	#**********************************************#
	# Copy Parameters.csv to project dir  *****#
	#**********************************************#
	sub.call("cp {0} {1}".format(parameter_file, project_name), shell=True)
	
if __name__ == "__main__":
	main(sys.argv)
