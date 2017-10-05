#!/usr/bin/python3
# -*- coding: utf-8 -*-


import csv		#Module to read csv files
import os
import re
import sys
import subprocess as sub
import copy
import shutil

import Prepare
import Utility


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

LipidsList = ['DSPC','DPPC','DLPC']
SolventsList = ['W','OCO','PW']

# String storing the name of the project
ProjectName = ''
# String storing the path to default gromamcs files
PathToDefault = ''
# Dictionnary with the path to the softwares to use
Softwares = {}
# Dictionnary holding all information for each job
Jobs = {}
#Dictionnary holding the data for PBS
PBS = {}
#Dictionnary holding the data for TGCC
TGCC = {}

#**********************************************#
#**********************************************#
#**********************************************#
#**********************************************#
# Read the set of parameters for the different #
# samples to study							   #
#**********************************************#
#**********************************************#
#**********************************************#


with open('Parameters.csv','r') as Input_Params:
	inputParamsReader = csv.reader(Input_Params, delimiter='|', skipinitialspace=True)
	
	#Variable to keep track of the jobs ordering
	jobNumber = 0
	
	for row in inputParamsReader:
		
		# Get the project name
		if row[0] == 'PROJECT' :
			ProjectName = row[1]
			continue
		
		#Skip empty rows 
		test_row = list(filter(None, row))
		if not test_row: continue
		
		#Skip lines starting with '#' or empty lines
		if '#' in row[0] or not test_row: continue
		
		# Each time a 'JOBID' is found a new job is created
		# it to separate different samples
		if row[0] == 'JOBID':
			jobNumber += 1
			#Step of the run : EM, NVE, NVT, NPT ...
			stepNumber = 0
			Jobs.update( {jobNumber: {row[0].strip(' '): row[1].strip(' '), row[2].strip(' '): int(row[3])} } )
			continue
		
		# Before reading information on each job the program stores the path to the softwares and gromacs default files from Parameters.csv
		if jobNumber < 1 : 
			if row[0] == 'GROMACS_LOC':
				Softwares.update({row[0].strip(' '): row[1].strip(' '), row[2].strip(' '): row[3].strip(' ')})
				continue
			if row[0] == 'VMD':
				Softwares.update({row[0].strip(' '): row[1].strip(' ')})
				continue
			if row[0] == 'PACKMOL':
				Softwares.update({row[0].strip(' '): row[1].strip(' ')})
				continue
			if row[0] == 'GROMACS_DEFAULT':
				PathToDefault = row[1].strip(' ')
		
		# Then if the program has reach a job, it will store all information available from Parameters.csv
		else:
			# Stores the node, number of nodes and processors per node.
			if row[0] == 'node' :
				Jobs[jobNumber].update({row[0].strip(' '): row[1].strip(' '), row[2].strip(' '): row[3].strip(' '), row[4].strip(' '): row[5].strip(' ')})
				continue
			
			# Stores the options to pass to mdrun
			if row[0] == 'MDRUN_OPT':
				Jobs[jobNumber].update({row[0].strip(' '): row[1].strip(' ')})
				continue
			
			# Stores the elements in the sample for the job
			if row[0] == 'SAMPLE':
				for element in range(1, len(row), 2):
					if row[element]:
						qtty = int(row[element+1])
						Jobs[jobNumber].update( { row[element].strip(' '): qtty })
				continue
			
			# For a general thermostat parameters. 
			if row[0] == 'THERMOSTAT':
				Jobs[jobNumber].update({ 'THERMOSTAT': {}})
				for param in range(1, len(row), 2):
					if row[param]:
						paramValue = row[param+1].strip(' ')
						if(is_number(paramValue)):
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
			if row[0] == 'BAROSTAT':
				Jobs[jobNumber].update({ 'BAROSTAT': {}})
				for param in range(1, len(row), 2):
					if row[param]:
						paramValue = row[param+1].strip(' ')
						if(is_number(paramValue)):
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
			if row[0] == 'TYPE':
				Jobs[jobNumber].update({row[0].strip(' '): row[1].strip(' ')})
				continue
			
			# Stores the dimensions of the sample for the job
			if row[0] == 'LX' :
				for dimension in range(0, len(row), 2):
					if row[dimension]:
						dimValue = row[dimension+1].strip(' ')
						if is_number(dimValue):
							Jobs[jobNumber].update( {row[dimension].strip(' '): float(dimValue) })
						else:
							Jobs[jobNumber].update( {row[dimension].strip(' '): dimValue })
					else:
						break
				
				continue
			
			# Stores information on the defo: defo per layer, Height, DzDefo, Radius, Version, TauDefo (thermostat), Tdefo (thermostat)
			if row[0].startswith('DEFO'):
				Jobs[jobNumber].update({ 'DEFO': {}})
				
				for paramDefo in range(1, len(row), 2):
					if row[paramDefo]:
						paramDefoVal = row[paramDefo+1].strip(' ')
						Jobs[jobNumber]['DEFO'].update( {row[paramDefo].strip(' '): paramDefoVal })
					else:
						break
					
				with open(PathToDefault+'/DEFO/Parameters_defo.csv','r') as Defo_Params:
					defoParam = csv.reader(Defo_Params, delimiter='|', skipinitialspace=True)
					Jobs[jobNumber]['DEFO'].update( { 'presets' : {} } )
					
					for row in defoParam: 
						firstCol = row[0].strip(' ')
						if firstCol:
							name = firstCol
							#Creates a dictionary for each preset
							Jobs[jobNumber]['DEFO']['presets'].update( { name : {}  })
							
							#Associates the parameters set in Parameters.csv to the correct protocol step dictionary
							for index, param in enumerate(row[1:]):
								if param:
									Jobs[jobNumber]['DEFO']['presets'][name].update( { param.strip(' ') : index } )
						
						else:
							for param, paramValue in Jobs[jobNumber]['DEFO']['presets'][name].items():
								
								paramValue = row[paramValue+1]
								
								if paramValue:
									if is_number(paramValue):
										if('.' in paramValue or not float(paramValue).is_integer() ):
											Jobs[jobNumber]['DEFO']['presets'][name][param] = float(paramValue)
										else:
											Jobs[jobNumber]['DEFO']['presets'][name][param] = int(paramValue)
									else:
										Jobs[jobNumber]['DEFO']['presets'][name][param] = paramValue.strip(' ')
				continue
			
			# Stores information on the support: version, density, nblipids
			if row[0].startswith('SU'):
				Jobs[jobNumber].update({ 'SU': {}})
				
				for paramSu in range(1, len(row), 2):
					if row[paramSu]:
						paramSuVal = row[paramSu+1].strip(' ')
						Jobs[jobNumber]['SU'].update( {row[paramSu].strip(' '): paramSuVal })
					else:
						break
					
				with open(PathToDefault+'/SU/Parameters_su.csv','r') as Su_Params:
					suParam = csv.reader(Su_Params, delimiter='|', skipinitialspace=True)
					Jobs[jobNumber]['SU'].update( { 'presets' : {} } )
					for row in suParam: 
						firstCol = row[0].strip(' ')
						if firstCol:
							name = firstCol
							#Creates a dictionary for each preset
							Jobs[jobNumber]['SU']['presets'].update( { name : {}  })
							
							#Associates the parameters set in Parameters.csv to the correct protocol step dictionary
							for index, param in enumerate(row[1:]):
								if param:
									Jobs[jobNumber]['SU']['presets'][name].update( { param.strip(' ') : index } )
							
						else:
							for param, paramValue in Jobs[jobNumber]['SU']['presets'][name].items():
								
								paramValue = row[paramValue+1]
								
								if paramValue:
									if is_number(paramValue):
										if('.' in paramValue or not float(paramValue).is_integer() ):
											Jobs[jobNumber]['SU']['presets'][name][param] = float(paramValue)
										else:
											Jobs[jobNumber]['SU']['presets'][name][param] = int(paramValue)
									else:
										Jobs[jobNumber]['SU']['presets'][name][param] = paramValue.strip(' ')
				continue
			
			
			if row[0].startswith('PROTOCOL'):
				Jobs[jobNumber].update( {'PROTOCOL': {}} )
				
				for index, protocol in enumerate(row[1:]):
					if 'DEFO' in Jobs[jobNumber] and protocol == 'DEFO':
						defoProtocol = row[index+2].strip(',')
						Jobs[jobNumber]['DEFO'].update( {'defoProtocol' : defoProtocol.split(',')} )
						
					if 'SU' in Jobs[jobNumber] and protocol == 'SU':
						suProtocol = row[index+2].strip(',')
						Jobs[jobNumber]['SU'].update( {'suProtocol' : suProtocol.split(',')} )
					
				continue
			
			
			# Job protocol. associates each parameter with an index
			if(row[0].startswith('COPY') or row[0].startswith('INIT') or row[0].startswith('EM') or row[0].startswith('NVE') or row[0].startswith('NVT') or row[0].startswith('NPT')):
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
			if(row[0] == ''):
				stepNb = str(stepNumber-1)
				for param, paramValue in Jobs[jobNumber]['PROTOCOL'][stepNb].items():
					if param != 'stepType':
						paramValue = row[paramValue]
						
						if(paramValue):
							if is_number(paramValue):
								if '.' in paramValue or not float(paramValue).is_integer() :
									Jobs[jobNumber]['PROTOCOL'][stepNb][param] = float(paramValue)
								else:
									Jobs[jobNumber]['PROTOCOL'][stepNb][param] = int(paramValue)
							else:
								Jobs[jobNumber]['PROTOCOL'][stepNb][param] = paramValue.strip(' ')

#**********************************************#
#**********************************************#
#**********************************************#
#**********************************************#
#**********************************************#
#**********************************************#

#Finding the software version. If not found ask the user.

print('Gromacs local :', Softwares['GROMACS_LOC'])

SoftwareVersion = re.sub("[^0-9]","",Softwares['GROMACS_LOC'].split('/')[-2])

if not SoftwareVersion:
	SoftwareVersion = input('No version detected in path. Please select the version you want between 4, 5 or 2016: ')

print('Gromacs version :', SoftwareVersion)

# Checking that the path to default files is provided. If not ask the user for a path.
if not PathToDefault:
	PathToDefault = input('Path to default gromacs files was not found. Please enter the path to the GROMACS default files directory: ')
print('Gromacs Default :', PathToDefault)

#Variables to set the correct path for gromacs
GROMACS_REM_prefixPath = ''
GROMACS_LOC_prefixPath = ''

#Creates the base directory
if(not os.path.isdir(ProjectName)):
	os.makedirs(ProjectName, exist_ok=True)
else:
	shutil.rmtree(ProjectName)
	os.makedirs(ProjectName, exist_ok=True)

#Dictionnary to perform sample copy
InitForCopy={}



#*******************************************************************************#
#*******************************************************************************#
#***********************								    ************************#
#********************** Generates files for use with TGCC ***********************#
#***********************								    ************************#
#*******************************************************************************#
#*******************************************************************************#

### MEMORY SHEET FOR TGCC SLURM, TO INCLUDE INTO the ccc_msub file
### ccc_msub to submit a job
### ccc_mpinfo displays info on utilisation of queues
### ccc_mqinfo displays info on queue definition
### ccc_mpp -u user to display all jobs from a user
### ccc_mpeek gives information about a job during its run.
### ccc_mdel kills jobs
### ccc_myproject give global info on time consumption of the group

#Slurm Variable Name 	Description 	                     Example values            	                    PBS/Torque analog
#$SLURM_JOBID       	Job ID             	                 5741192                	                    $PBS_JOBID
#$SLURM_JOB_NAME    	Job Name           	                 myjob                  	                    $PBS_JOBNAME
#$SLURM_SUBMIT_DIR  	Submit Directory           	         /lustre/payerle/work 	                        $PBS_O_WORKDIR
#$SLURM_JOB_NODELIST  	Nodes assigned to job    	         compute-b24-[1-3,5-9],compute-b25-[1,4,8] 	    cat $PBS_NODEFILE
#$SLURM_SUBMIT_HOST 	Host submitted from 	             login-1.deepthought2.umd.edu 	                $PBS_O_HOST
#$SLURM_JOB_NUM_NODES 	Number of nodes allocated to job 	  2 	                                        $PBS_NUM_NODES
#$SLURM_CPUS_ON_NODE 	Number of cores/node 	              8,3 	                                        $PBS_NUM_PPN
#$SLURM_NTASKS  	    Total number of cores for job??? 	  11 	                                        $PBS_NP
#$SLURM_NODEID 	        Index to node running on
#relative to nodes assigned to job 	                           0 	                                        $PBS_O_NODENUM
#$PBS_O_VNODENUM 	    Index to core running on
#within node 	  4 	$SLURM_LOCALID
#$SLURM_PROCID 	        Index to task relative to job 	     0 	                                             $PBS_O_TASKNUM - 1

#####################EXAMPLE FROM TGCC WEB SITE

##!/bin/bash
##MSUB -r MyJob_Para          # Request name
##MSUB -n 32                  # Number of tasks to use
##MSUB -T 1800                # Elapsed time limit in seconds
##MSUB -o example_%I.o        # Standard output. %I is the job id
##MSUB -e example_%I.e        # Error output. %I is the job id
##MSUB -A paxxxx              # Project ID

#set -x
#cd ${BRIDGE_MSUB_PWD}
#ccc_mprun ./a.out
##or
## ccc_mprun -n 32 ./a.out
##or
## ccc_mprun -n ${BRIDGE_MSUB_NPROC} ./a.out
## BRIDGE_MSUB_NPROC represents the number of tasks
#####################EXAMPLE FROM TGCC WEB SITE



if(sys.argv[1] == '--tgcc'):
	#Read informations provided in TGCCinfo.csv if the option --tgcc is used
	with open('TGCCinfo.csv','r') as TGCCinfo:
		Reader = csv.reader(TGCCinfo, delimiter='|',skipinitialspace=True)
		for row in Reader:
			if Reader.line_num == 1:
				continue
			else:
				TGCC.update( {row[0]: row[1] } )
				
	
	#**********************************************#
	#*****      Generating MDrun files        *****#
	#**********************************************#
	for sampleNumber in range(1, len(Jobs)+1):
		name = Jobs[sampleNumber]['JOBID']
		print("===============\n===============\nSample : "+name+"\n===============\n===============\n")
		#Creates the sample directory
		if(not os.path.isdir(ProjectName+'/'+name)):
			os.mkdir(ProjectName+'/'+name)

		ScriptFile = open(ProjectName+'/'+name+'/run.sh','w')
		ScriptFile.truncate()
		CopyToScratch = str("""
						#!/bin/bash +x
						
						# ADDAPT AUTOMATICALLY TO GROMACS VERSION ?
						# OR AT LEAST VERIFY THAT THE TWO VERSION INFORMATION
						# IN Parameters.csv and TGCCinfo.csv ARE CONSISTENT  ?
						
						module load {3}
						alias gmx="gmx_mpi"
						
						""").format(ProjectName,name,TGCC['username'],TGCC['GMXversion'])
		ScriptFile.write(Utility.RemoveUnwantedIndent(CopyToScratch))

		PrevCmdFiles={}
		
		#**********************************************#
		# Read the default files for EM, NVT, NVP, NVE #
		#**********************************************#
		#**********************************************#
		# Lipids files #
		if( Jobs[sampleNumber]['TYPE'] == 'BILAYER' or Jobs[sampleNumber]['TYPE'] == 'TRILAYER'):
			if SoftwareVersion.startswith('4'):
				EMdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/'
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/'
				
			else:
				EMdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/gmx '
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/gmx '
		
		# Solvent files #
		if( Jobs[sampleNumber]['TYPE'] == 'SOLVENT'):
			if SoftwareVersion.startswith('4'):
				EMdefault = open(PathToDefault+'/mdp_sol/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/'
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/'
				
			else:
				EMdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/gmx '
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/gmx '
		#**********************************************#
		#**********************************************#
		
		
		
		#**********************************************#
		#Go into the corresponding directory
		#**********************************************#
		with Utility.cd(ProjectName+'/'+name):
			OutputFilename=''
			ScriptFile = open('run.sh','a+')
			ScriptFile.write('mkdir -p mdrun_out\n')
			ScriptFile.write('mkdir -p grompp_out\n')
			SoFarFile = str("""{0}_SoFar.txt""").format(name)
			SoFar = open(SoFarFile, 'w')
			
			#**********************************************#
			#Prepare the steps **************************#
			#**********************************************#
			for step in range(0, len(Jobs[sampleNumber]['PROTOCOL'])):

				step = str(step)
				
				
				#============================================================
				#					STEP IS INITIALIZATION
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('INIT')):
					if(Jobs[sampleNumber]['TYPE'] == 'BILAYER'):
						InitForCopy = Prepare.InitBilayer(Jobs[sampleNumber], Softwares, GROMACS_LOC_prefixPath, PathToDefault)
					if(Jobs[sampleNumber]['TYPE'] == 'SOLVENT'):
						InitForCopy = Prepare.InitSolvent(Jobs[sampleNumber], Softwares, GROMACS_LOC_prefixPath, PathToDefault)
					SoFar.write(str("""{0} done""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType']))
					PrevCmdFiles = copy.deepcopy(InitForCopy)
					continue

				#============================================================
				#					STEP IS SAMPLE COPY
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('COPY')):
#? ERROR ?				CopyFrom = Jobs[ 'SAMPLE'+str(Jobs[sampleNumber]['PROTOCOL'][step]['samplenumber']) ]['JOBID']
					CopyFrom = Jobs[ str(Jobs[sampleNumber]['PROTOCOL'][step]['samplenumber']) ]['JOBID']
					Prepare.CopySample(CopyFrom)
					SoFar.write(  str("""{0} {1} done""").format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], CopyFrom) )
					PrevCmdFiles = copy.deepcopy(InitForCopy)
					continue

				#============================================================
				#						STEP IS EM
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('EM')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, EMdefault, SoftwareVersion)

					Comment = str("""#Energy Minimization using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""gmx_mpi grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""ccc_mprun gmx_mpi mdrun {4} -deffnm {2}_{3}-{1} -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

				#============================================================
				#						STEP IS NVE
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NVE')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, NVEdefault, SoftwareVersion)

					Comment = str("""#NVE using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""gmx_mpi grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""ccc_mdrun gmx_mpi mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

				#============================================================
				#						STEP IS NVT
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NVT')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, NVTdefault, SoftwareVersion)

					Comment = str("""#NVT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""gmx_mpi grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""ccc_mdrun gmx_mpi mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])

					PrevCmdFiles['OUTPUT'] = OutputFilename

					continue

				#============================================================
				#						STEP IS NPT
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NPT')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, NPTdefault, SoftwareVersion)

					Comment = str("""#NPT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""gmx_mpi grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""ccc_mdrun gmx_mpi mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue
			
			#**********************************************#
			# End for prepare the steps **************#
			#**********************************************#
			
			#**********************************************#
			# Creating TGCC files  **************#
			#**********************************************#
			sub.call("""chmod a+x run.sh""",shell=True)
			SoFar.close()
			TGCCfile = str("""
						#! /bin/sh -x
						#MSUB -q {2}      # queue = standard, test, long (in Parameters.csv, key "node")
						#MSUB -n {3}      # total number of cores, 16 cores per nodes (32 logical cores ) (in Parameters.cvs, key "ppn")
						#MSUB -T {5}      # times in seconds (in TGCCinfo.csv, key "time_s")
						#MSUB -r {0}	  # job name (automatically generated)
						#MSUB -A {4} 	  # group for allocation (gen7662) (in TgccInfo.csv, key "group")
						#MSUB -V          # transfer the variable of ENVIRONNEMENT
						#MSUB -m be       # email at begin and at end
						#MSUB -M {1}      # email address
						#MSUB -oe %I.eo   # input and output of JOB
						
						# the number of nodes is deduced automatically from the number of cores (fixed nb core/node = 16 on curie noeud fin)
						# ccc_mprun gives this information to gmx_mpi (seems to work !?)
						
						echo "===================== BEGIN JOB $SLURM_JOBID =============================== "
						
						JOBINFO=$SLURM_SUBMIT_DIR/$SLURM_JOB_NAME-$SLURM_JOBID.jobinfo"
						printf "Time =  `date`\\n" > $JOBINFO
						printf "SLURM submit directory = $SLURM_SUBMIT_DIR\\n" >> $JOBINFO
						printf "TGCC queue = {2}, user {4}, max time = {5} seconds \\n" >> $JOBINFO
						printf "SLURM job ID = $SLURM_JOBID \\n" >> $JOBINFO
						printf "SLURM job name = $SLURM_JOB_NAME \\n" >> $JOBINFO
						printf "This job will run on {3} processors\\n" >> $JOBINFO
						printf "List of nodes : $SLURM_NODEID \\n\\n" >> $JOBINFO
						
						
						export MYTMPDIR="$SCRATCHDIR/JOB_$SLURM_JOBID"
						export OUTPUTDIR="$SLURM_SUBMIT_DIR/JOB_$SLURM_JOBID_OUTPUT
						
						mkdir -p $MYTMPDIR
						mkdir -p $OUTPUTDIR
						
						rsync $SLURM_SUBMIT_DIR/* $MYTMPDIR 
						cd $MYTMPDIR ./run.sh  > $OUTPUTDIR/JOB_$SLURM_JOBID.out
						rsync -r ./* $OUTPUTDIR/.
						cd $SLURM_SUBMIT_DIR
						rm -rf $MYTMPDIR
						
						echo "===================== END  JOB $SLURM_JOBID =============================== "
						""").format(ProjectName+'_'+name+'__'+PrevCmdFiles['SYSTEM']+'__'+SoftwareVersion, TGCC['mail'], TGCC['queue'], Jobs[sampleNumber]['ppn'],TGCC['group'],TGCC['time_s'])
			f = open(name+'.ccc_msub','w')
			f.write( Utility.RemoveUnwantedIndent(TGCCfile) )
			f.close()
			
			#**********************************************#
			# End Creating TGCC files  **************#
			#**********************************************#
			#
			## ALREADY DONE IN THE MSUB SCRIPT
			#
			#CopyToScratch = str("""
							#echo "End of run for {1}"
							#cp -r /scratch/{2}/gromacs/{0}/{1} ${{LOCALDIR}}/{1}_OUTPUT
							#rm -r /scratch/{2}/gromacs/{0}
							#""").format(ProjectName, name, TGCC['username'])
			#ScriptFile.write(Utility.RemoveUnwantedIndent(CopyToScratch) )
		ScriptFile.close()
		
		EMdefault.close()
		NVTdefault.close()
		NPTdefault.close()
		NVEdefault.close()
	
	#**********************************************#
	# Creating scripts for starting jobs  *****#
	#**********************************************#
	with Utility.cd(ProjectName):
		#Creates a script for sending all samples to queue
		ScriptForJobs = open('SendToSLURM.sh','w')
		ScriptForJobs.truncate()
		ScriptForJobs.write("""#!/bin/bash\n\n""")
		ScriptForJobs.close()

		ScriptForJobs = open('SendToSLURM.sh','a+')

		for sampleNumber in Jobs:
			name = Jobs[sampleNumber]['JOBID']
			qsubCmd = str("""cd {0} \n ccc_msub {0}.ccc_msub \n cd .. \n""").format(name)
			ScriptForJobs.write(qsubCmd)

		ScriptForJobs.write( str("""ccc_mpp -u {0}""").format(TGCC['username']) )
		ScriptForJobs.close()
		sub.call("""chmod a+x SendToSLURM.sh""",shell=True)
	
	#**********************************************#
	# Copy Parameters.csv to project dir  *****#
	#**********************************************#
	sub.call('cp Parameters.csv '+ ProjectName, shell=True)



#*******************************************************************************#
#*******************************************************************************#
#***********************								************************#
#********************** Generates files for use with PBS ***********************#
#***********************								************************#
#*******************************************************************************#
#*******************************************************************************#

if(sys.argv[1] == '--pbs'):
	#Read informations provided in PBSinfo.csv if the option --pbs is used
	with open('PBSinfo.csv','r') as PBSinfo:
		Reader = csv.reader(PBSinfo, delimiter='|',skipinitialspace=True)
		for row in Reader:
			if Reader.line_num == 1:
				continue
			else:
				PBS.update( {row[0]: row[1] } )
				
	
	#**********************************************#
	#*****      Generating MDrun files        *****#
	#**********************************************#
	for sampleNumber in range(1, len(Jobs)+1):
		name = Jobs[sampleNumber]['JOBID']
		print("===============\n===============\nSample : "+name+"\n===============\n===============\n")
		#Creates the sample directory
		if(not os.path.isdir(ProjectName+'/'+name)):
			os.mkdir(ProjectName+'/'+name)

		ScriptFile = open(ProjectName+'/'+name+'/run.sh','w')
		ScriptFile.truncate()
		CopyToScratch = str("""
						#!/bin/bash
						mkdir -p /scratch/{2}/gromacs/{0}
						cd ..
						LOCALDIR="$(pwd)"
						cp -r {1} /scratch/{2}/gromacs/{0}
						cd /scratch/{2}/gromacs/{0}/{1}

						""").format(ProjectName,name,PBS['username'] )
		ScriptFile.write(Utility.RemoveUnwantedIndent(CopyToScratch))

		PrevCmdFiles={}
		
		#**********************************************#
		# Read the default files for EM, NVT, NVP, NVE #
		#**********************************************#
		#**********************************************#
		# Lipids files #
		if( Jobs[sampleNumber]['TYPE'] == 'BILAYER' or Jobs[sampleNumber]['TYPE'] == 'TRILAYER'):
			if SoftwareVersion.startswith('4'):
				EMdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/'
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/'
				
			else:
				EMdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/gmx '
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/gmx '
		
		# Solvent files #
		if( Jobs[sampleNumber]['TYPE'] == 'SOLVENT'):
			if SoftwareVersion.startswith('4'):
				EMdefault = open(PathToDefault+'/mdp_sol/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/'
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/'
				
			else:
				EMdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/gmx '
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/gmx '
		#**********************************************#
		#**********************************************#
		
		
		
		#**********************************************#
		#Go into the corresponding directory
		#**********************************************#
		with Utility.cd(ProjectName+'/'+name):
			OutputFilename=''
			ScriptFile = open('run.sh','a+')
			ScriptFile.write('mkdir -p mdrun_out\n')
			ScriptFile.write('mkdir -p grompp_out\n')
			SoFarFile = str("""{0}_SoFar.txt""").format(name)
			SoFar = open(SoFarFile, 'w')
			
			#**********************************************#
			#Prepare the steps **************************#
			#**********************************************#
			for step in range(0, len(Jobs[sampleNumber]['PROTOCOL'])):

				step = str(step)
				
				
				#============================================================
				#					STEP IS INITIALIZATION
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('INIT')):
					if(Jobs[sampleNumber]['TYPE'] == 'BILAYER'):
						InitForCopy = Prepare.InitBilayer(Jobs[sampleNumber], Softwares, GROMACS_LOC_prefixPath, PathToDefault)
					if(Jobs[sampleNumber]['TYPE'] == 'SOLVENT'):
						InitForCopy = Prepare.InitSolvent(Jobs[sampleNumber], Softwares, GROMACS_LOC_prefixPath, PathToDefault)
					SoFar.write(str("""{0} done""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType']))
					PrevCmdFiles = copy.deepcopy(InitForCopy)
					continue

				#============================================================
				#					STEP IS SAMPLE COPY
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('COPY')):
# TRY FOR CLAIRE ? ERROR					CopyFrom = Jobs[ str(Jobs[sampleNumber]['PROTOCOL'][step]['samplenumber']) ]['JOBID']
					CopyFrom = Jobs[ 'SAMPLE'+str(Jobs[sampleNumber]['PROTOCOL'][step]['samplenumber']) ]['JOBID']
					Prepare.CopySample(CopyFrom)
					SoFar.write(  str("""{0} {1} done""").format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], CopyFrom) )
					PrevCmdFiles = copy.deepcopy(InitForCopy)
					continue

				#============================================================
				#						STEP IS EM
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('EM')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, EMdefault, SoftwareVersion)

					Comment = str("""#Energy Minimization using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1} -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

				#============================================================
				#						STEP IS NVE
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NVE')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, NVEdefault, SoftwareVersion)

					Comment = str("""#NVE using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

				#============================================================
				#						STEP IS NVT
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NVT')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, NVTdefault, SoftwareVersion)

					Comment = str("""#NVT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])

					PrevCmdFiles['OUTPUT'] = OutputFilename

					continue

				#============================================================
				#						STEP IS NPT
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NPT')):

					Prepare.WriteMDP(Jobs[sampleNumber], step, NPTdefault, SoftwareVersion)

					Comment = str("""#NPT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} &> grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro &> mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_REM_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'])
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue
			
			#**********************************************#
			# End for prepare the steps **************#
			#**********************************************#
			
			#**********************************************#
			# Creating PBS files  **************#
			#**********************************************#
			sub.call("""chmod a+x run.sh""",shell=True)
			SoFar.close()
			PBSfile = str("""
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

						./run.sh  $PBS_JOBID >  $PBS_O_WORKDIR/$PBS_JOBNAME-$PBS_JOBID.out

						echo "End of GMX"
						date

						echo "============================="
						""").format(ProjectName+'_'+name+'__'+PrevCmdFiles['SYSTEM']+'__'+SoftwareVersion, PBS['mail'], Jobs[sampleNumber]['node'], Jobs[sampleNumber]['nbnodes'], Jobs[sampleNumber]['ppn'])
			f = open(name+'.pbs','w')
			f.write( Utility.RemoveUnwantedIndent(PBSfile) )
			f.close()
			#**********************************************#
			# End Creating PBS files  **************#
			#**********************************************#
			
			CopyToScratch = str("""
							echo "End of run for {1}"
							cp -r /scratch/{2}/gromacs/{0}/{1} ${{LOCALDIR}}/{1}_OUTPUT
							rm -r /scratch/{2}/gromacs/{0}
							""").format(ProjectName, name, PBS['username'])
			ScriptFile.write(Utility.RemoveUnwantedIndent(CopyToScratch) )
		ScriptFile.close()
		
		EMdefault.close()
		NVTdefault.close()
		NPTdefault.close()
		NVEdefault.close()
	
	#**********************************************#
	# Creating scripts for starting jobs  *****#
	#**********************************************#
	with Utility.cd(ProjectName):
		#Creates a script for sending all samples to queue
		ScriptForJobs = open('SendToPBS.sh','w')
		ScriptForJobs.truncate()
		ScriptForJobs.write("""#!/bin/bash\n\n""")
		ScriptForJobs.close()

		ScriptForJobs = open('SendToPBS.sh','a+')

		for sampleNumber in Jobs:
			name = Jobs[sampleNumber]['JOBID']
			qsubCmd = str("""cd {0} \n qsub {0}.pbs \n cd .. \n""").format(name)
			ScriptForJobs.write(qsubCmd)

		ScriptForJobs.write( str("""qstat -n -u {0}""").format(PBS['username']) )
		ScriptForJobs.close()
		sub.call("""chmod a+x SendToPBS.sh""",shell=True)
	
	#**********************************************#
	# Copy Parameters.csv to project dir  *****#
	#**********************************************#
	sub.call('cp Parameters.csv '+ ProjectName, shell=True)







#*********************************************************************************#
#*********************************************************************************#
#***********************                         *********************************#
#********************** Generates files for local ********************************#
#***********************                         *********************************#
#*********************************************************************************#
#*********************************************************************************#

if(sys.argv[1] == '--local'):
	#**********************************************#
	#*****      Generating MDrun files        *****#
	#**********************************************#
	for i in range(1,len(Jobs)+1):
		sampleNumber = str("""SAMPLE{0}""").format(i)
		sampleNumber = i
		name = Jobs[sampleNumber]['JOBID']
		print("===============\n===============\nSample "+name+"\n===============\n===============\n")
		#Creates the sample directory
		if(not os.path.isdir(ProjectName+'/'+name)):
			os.mkdir(ProjectName+'/'+name)

		ScriptFile = open(ProjectName+'/'+name+'/run.sh','w')
		ScriptFile.truncate()
		CopyToOutput = str("""
						#!/bin/bash
						mkdir -p ../{0}_OUTPUT
						cp * ../{0}_OUTPUT/
						cd ../{0}_OUTPUT/

						""").format(name)
		ScriptFile.write(Utility.RemoveUnwantedIndent(CopyToOutput))
		
		PrevCmdFiles={}
		
		#**********************************************#
		# Read the default files for EM, NVT, NVP, NVE #
		#**********************************************#
		#**********************************************#
		# Lipids files #
		if( Jobs[sampleNumber]['TYPE'] == 'BILAYER' or Jobs[sampleNumber]['TYPE'] == 'TRILAYER'):
			if SoftwareVersion.startswith('4') :
				EMdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_lipids/GROMACS4/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/'
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/'
				
			else:
				EMdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_lipids/GROMACS2016/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/gmx '
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/gmx '
		
		# Solvent files #
		if( Jobs[sampleNumber]['TYPE'] == 'SOLVENT'):
			if SoftwareVersion.startswith('4'):
				EMdefault = open(PathToDefault+'/mdp_sol/GROMACS4/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_sol/GROMACS4/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/'
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/'
				
			else:
				EMdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/EMdefault.mdp','r')
				NVTdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NVTdefault.mdp','r')
				NPTdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NPTdefault.mdp','r')
				NVEdefault = open(PathToDefault+'/mdp_sol/GROMACS2016/NVEdefault.mdp','r')

				GROMACS_REM_prefixPath = Softwares['GROMACS_REM']+'/gmx '
				GROMACS_LOC_prefixPath = Softwares['GROMACS_LOC']+'/gmx '
		#**********************************************#
		#**********************************************#
			
		#Go into the corresponding directory
		with Utility.cd(ProjectName+'/'+name):
			OutputFilename=''
			ScriptFile = open('run.sh','a+')
			ScriptFile.write('mkdir -p mdrun_out\n')
			ScriptFile.write('mkdir -p grompp_out\n')
			SoFarFile = str("""{0}_SoFar.txt""").format(name)
			SoFar = open(SoFarFile, 'w')
			
			#Executes the steps
			for step in range(0, len(Jobs[sampleNumber]['PROTOCOL'])):

				step = str(step)

				#============================================================
				#					STEP IS INITIALIZATION
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('INIT')):
					if(Jobs[sampleNumber]['TYPE'] == 'BILAYER'):
						InitForCopy = Prepare.InitBilayer(Jobs[sampleNumber], Softwares, GROMACS_LOC_prefixPath, PathToDefault)
					if(Jobs[sampleNumber]['TYPE'] == 'SOLVENT'):
						InitForCopy = Prepare.InitSolvent(Jobs[sampleNumber], Softwares, GROMACS_LOC_prefixPath, PathToDefault)
					SoFar.write(str("""{0} done \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType']))
					PrevCmdFiles = copy.deepcopy(InitForCopy)
					continue

				#============================================================
				#					STEP IS SAMPLE COPY
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('COPY')):

					CopyFrom = Jobs[ 'SAMPLE'+str(Jobs[sampleNumber]['PROTOCOL'][step]['samplenumber']) ]['JOBID']
					Prepare.CopySample(CopyFrom)
					SoFar.write(  str("""{0} {1} done \n\n""").format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], CopyFrom) )
					PrevCmdFiles = copy.deepcopy(InitForCopy)
					continue

				#============================================================
				#						STEP IS EM
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('EM')):
					
					#Output for Pulling data
					PullingOutput = ''
					if 'DEFO' in Jobs[sampleNumber]:
						PullingOutput = "-pf {0}_pf.xvg -px {0}_px.xvg".format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'])
					
					Prepare.WriteMDP(Jobs[sampleNumber], step, EMdefault, SoftwareVersion)

					Comment = str("""#Energy Minimization using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} |& tee grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun -deffnm {2}_{3}-{1} -c {2}_{3}-{1}_out.gro {4} |& tee mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], PullingOutput)
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

				#============================================================
				#						STEP IS NVE
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NVE')):
					
					#Output for Pulling data
					PullingOutput = ''
					if 'DEFO' in Jobs[sampleNumber]:
						PullingOutput = "-pf {0}_pf.xvg -px {0}_px.xvg".format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'])

					Prepare.WriteMDP(Jobs[sampleNumber], step, NVEdefault, SoftwareVersion)

					Comment = str("""#NVE using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} |& tee grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro {5} |& tee mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'], PullingOutput)
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

				#============================================================
				#						STEP IS NVT
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NVT')):
					
					#Output for Pulling data
					PullingOutput = ''
					if 'DEFO' in Jobs[sampleNumber]:
						PullingOutput = "-pf {0}_pf.xvg -px {0}_px.xvg".format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'])

					Prepare.WriteMDP(Jobs[sampleNumber], step, NVTdefault, SoftwareVersion)

					Comment = str("""#NVT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} |& tee grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro {5} |& tee mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'], PullingOutput)
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])

					PrevCmdFiles['OUTPUT'] = OutputFilename

					continue

				#============================================================
				#						STEP IS NPT
				#============================================================
				if(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'].startswith('NPT')):
					
					#Output for Pulling data
					PullingOutput = ''
					if 'DEFO' in Jobs[sampleNumber]:
						PullingOutput = "-pf {0}_pf.xvg -px {0}_px.xvg".format( Jobs[sampleNumber]['PROTOCOL'][step]['stepType'])

					Prepare.WriteMDP(Jobs[sampleNumber], step, NPTdefault, SoftwareVersion)

					Comment = str("""#NPT using parameters in {0}.mdp:\n#Input : {1}\n#Output: {2}_{3}-{0}_out.gro \n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(Comment)

					cmd = str("""{0}grompp -f {1}.mdp -po {2}-{1}_out.mdp -c {3} -p {2}.top -maxwarn 10 -o {2}_{5}-{1}.tpr -n {4} |& tee grompp_out/grompp_{5}_{1}.output\n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], PrevCmdFiles['OUTPUT'], PrevCmdFiles['INDEX'],Jobs[sampleNumber]['JOBNUM'])
					ScriptFile.write(cmd)

					cmd = str("""{0}mdrun {4} -deffnm {2}_{3}-{1}  -c {2}_{3}-{1}_out.gro {5} |& tee mdrun_out/mdrun_{3}_{1}.output \n\n""").format(GROMACS_LOC_prefixPath, Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'], Jobs[sampleNumber]['MDRUN_OPT'], PullingOutput)
					ScriptFile.write(cmd)

					cmd = str("""echo '{0} done' >> {1}_SoFar.txt\n\n\n""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'],name)
					ScriptFile.write(cmd)

					OutputFilename = str("""{1}_{2}-{0}_out.gro""").format(Jobs[sampleNumber]['PROTOCOL'][step]['stepType'], PrevCmdFiles['SYSTEM'], Jobs[sampleNumber]['JOBNUM'])
					PrevCmdFiles['OUTPUT'] = OutputFilename
					continue

			CopyToOutput = str("""
							echo "End of run for {1}"
							""").format(ProjectName,name)
			ScriptFile.write(Utility.RemoveUnwantedIndent(CopyToOutput) )

			sub.call("""chmod a+x run.sh""",shell=True)
			SoFar.close()
			
		ScriptFile.close()
		EMdefault.close()
		NVTdefault.close()
		NPTdefault.close()
		NVEdefault.close()

	#**********************************************#
	# Creating scripts for starting jobs  *****#
	#**********************************************#
	with Utility.cd(ProjectName):
		#Creates a script for sending all samples to queue
		ScriptToStart = open('Start.sh','w')
		ScriptToStart.truncate()
		ScriptToStart.write("""#!/bin/bash\n\n""")
		ScriptToStart.close()

		ScriptToStart = open('Start.sh','a+')

		for sampleNumber in Jobs:
			name = Jobs[sampleNumber]['JOBID']
			startCmd = str("""cd {0} \n ./run.sh \ncd .. \n""").format(name)
			ScriptToStart.write(startCmd)
			
		ScriptToStart.close()
		sub.call("""chmod a+x Start.sh""",shell=True)

	#**********************************************#
	# Copy Parameters.csv to project dir  *****#
	#**********************************************#
	sub.call('cp Parameters.csv '+ProjectName,shell=True)

