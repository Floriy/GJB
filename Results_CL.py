#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 1st param : directory to study
#
import os
import glob
import sys
import csv
import shutil
import subprocess as sub
import numpy as np
import copy
import matplotlib as mpl
import Utility

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import readline

import matplotlib.pyplot as plt
#*******************************************#
# Variables to change if files are named    #
# Differently                               #

suffix = '_OUTPUT'
mdrun_out = 'mdrun_out'

SolventsList = ['W','OCO']
LipidsList = ['DSPC','DPPC','DLPC']
#											#
#*******************************************#
#*******************************************#
# HTML related variables					#
#											#

HTMLheader = str("""
				<!DOCTYPE html>
				<html>
				<meta name="viewport" content="width=device-width, initial-scale=1.0">

				<head>
				</head>

				<style>
				.sum {
					font-family: arial, sans-serif;
					border-collapse: collapse;
					width: 100%;
				}

				.sum td, th {
					border: 1px solid #000000;
					text-align: center;
					padding: 8px;
				}

				.sum .th-run {
					background-color: #7592B7;
				}
				.sum .th-smpl {
					background-color: #C9E0FF;
				}
				</style>

				<body>

	""")
HTMLfoot = str("""
		</table>
		</body>
		</meta>
		</html>
		""")


#*******************************************#
#*******************************************#
#*******************************************#
#*******************************************#
#				Script begins here			#
#											#
#*******************************************#
#*******************************************#
#*******************************************#
#*******************************************#
gmx_path = sys.argv[2]
OrderParam = sys.argv[3] #if number use this number for P2 computation if -- P2 not computed
ProjectName = sys.argv[4] #Name of the main directory
Parameters = sys.argv[5:] #Name of the runs to study

#List for word completion in the interactive part
WordsForCompletion = ['all','-run','-avg','exit','-sameplot']

#Dictionary containing the averages and Standard Error
Averages = {}

#Dictionary for the performances
Perf = {}

#Types of the systems
Types = {}

RunsToStudy = {Parameters[i]:Parameters[i+1] for i in range(0,len(Parameters),2)}

##Variable holding the correct path for different gromacs versions
#GROMACS_LOC_prefixPath = ''
print(RunsToStudy)


#Go into each sample directory and calls g_energy for each specified run
with Utility.cd(ProjectName):

	#Finds Gromacs path#
	GROMACS_PATH = gmx_path
        TEST_PATH = GROMACS_PATH.uppercase()
	if(TEST_PATH[TEST_PATH.rfind('S')+1] == '4'):
		GROMACS_LOC_prefixPath = GROMACS_PATH + '/g_'
	if( TEST_PATH.find('2016') != -1):
		GROMACS_LOC_prefixPath = GROMACS_PATH + '/gmx '
	
	##########print(GROMACS_LOC_prefixPath)
	#Create HTML file for summary table
	HTMLout = open('Analysis_summary.html','a+')
	Utility.deleteContent(HTMLout)
	HTMLout.write(Utility.RemoveUnwantedIndent(HTMLheader))
	HTMLout.write(Utility.RemoveUnwantedIndent("""
												<h1>{0}</h1>
												<table class="sum">
												""".format(ProjectName)) )

	#Creates the base directory
	if(not os.path.isdir("Analysis")):
		os.mkdir("Analysis")
                #Creates the base directory
	        if(not os.path.isdir("Analysis/Density")):
		    os.mkdir("Analysis/Density")

	else:
		shutil.rmtree('Analysis')
		os.mkdir("Analysis")

	#List of the samples studied
	SampleDir = [ name for name in os.listdir(os.getcwd()) if os.path.isdir(name) if name.endswith(suffix)]

	#Data to plot and analyse
	#Contains a list of the quantities in the .xvg file with their associated column number
	DATA = {}

	#Dictionary for stocking data summary on each specified run
	Runsmry = {}
	for Run in RunsToStudy:

		if Run.startswith('NPT'):
			HTMLRunhead = str("""
			<!-- Information on {0} -->
				<tr>
					<th class="th-run" colspan="8">{0}</th>
				</tr>
			<!-- Samples -->
			""".format(Run))
			Runsmry.update( {Run : Utility.RemoveUnwantedIndent(HTMLRunhead)} )
			continue

		if Run.startswith('NVE'):
			HTMLRunhead = str("""
			<!-- Information on {0} -->
				<tr>
					<th class="th-run" colspan="8">{0}</th>
				</tr>
			<!-- Samples -->
			""".format(Run))
			Runsmry.update( {Run : Utility.RemoveUnwantedIndent(HTMLRunhead)} )
			continue
        
	for Smpl in SampleDir:
		samplename = Smpl.split(suffix)[0]
		Perf.update( {samplename : {}} )
		DATA.update( {samplename: {}} )
		
		with Utility.cd(Smpl):
			
			molecules = {}
			Topo = glob.glob('*.top')
			
			if Topo: 
				f_avg = open(Topo[0],'r') #*************************
				for heading_and_lines in Utility.group_by_heading( f_avg, '[ molecules ]'):
					lines = []
					lines.extend(heading_and_lines[2:])
					for mol in lines:
						print(mol)
						Molname = mol.split(' ')[0]
						Molnumber = int(mol.split(' ')[1].split('\n')[0])
						molecules.update({Molname : Molnumber})
			DATA[samplename].update( { 'Mol': molecules } )
			
			
			
			for RTS in RunsToStudy:
				print("""Looking for {0} in {1}""".format(RTS,Smpl))
				RunsFound = [name for name in os.listdir(os.getcwd()) if name.endswith('edr') if (RTS+'.edr') in name]
				print(RunsFound)
				for RF in RunsFound:
					Perf[samplename].update( {RF.rsplit('-')[1].split('.edr')[0] : {}} )
				#####print(RunsFound)
				Types.update( {Smpl.split(suffix)[0] : (''.join(RunsFound)).split('_')[0]} )
				######print(Types)
				VariablesForOutput = ''
				
				for name in os.listdir(os.getcwd()):
					if (name.endswith('gro')):
						if ((RTS+'.gro') in name):
							GroFound = name
				
				print("""Found final gro file {0} in {1}""".format(GroFound,Smpl))
				f_gro = open(GroFound,'r')
				for i, x in enumerate(f_gro):
					if (i == 1):
						ParticleNum = x
				f_gro.close()
				DATA[samplename]['ParticleNum'] = int(ParticleNum)
				#print("""Found {0} Particles """.format(ParticleNum))
				
				for i in range(1,1000):
					VariablesForOutput += 'echo {0}\\n;'.format(i)
				VariablesForOutput += '\\n;\\n;'
				
							
				for RF in RunsFound:
					Input = RF
					Output = Smpl + '_' + RF.replace('edr','xvg')
					RunName = RF.rsplit('-')[1].split('.edr')[0]
					print(RunsToStudy[RunName])
					cmd = """({0}) | {1}energy {4} -f {2} -o ../Analysis/{3}""".format(VariablesForOutput, GROMACS_LOC_prefixPath, Input, Output, RunsToStudy[RunName])
					g_energyProc = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE,shell=True)
					g_energyOut = g_energyProc.stdout.read()

					for line in str(g_energyOut):
						if "Potential" in line:
							print(line)

					OutputAvg = Output.replace('xvg','avg')

					f_out = open('../Analysis/'+OutputAvg,'wb+')
					f_out.write(g_energyOut)
					f_out.close()
					if OrderParam != '--':
						if 'NPT' in RF and 'BILAYER' in RF:
							trajfiles = [trjf for trjf in os.listdir(os.getcwd()) if trjf == RF.replace('edr','xtc')]
							if trajfiles:
								RunName = trajfiles[0].rsplit('-')[1].split('.xtc')[0]
								for lipid in LipidsList:
									if lipid in DATA[samplename]['Mol']:
										Nblipids = DATA[samplename]['Mol'][lipid]
										Utility.do_order(trajfiles[0],RunsToStudy[RunName], OrderParam ,[0,0,1] , Nblipids, lipid, GROMACS_LOC_prefixPath)
										shutil.copyfile('order.dat', """../Analysis/{0}__{1}_order.dat""".format(samplename+'_'+RunName, lipid))
								os.remove('order.dat')
							
			
			#for RF in RunsFound:
					
						
				# Get the performance
				with Utility.cd(mdrun_out):
					MDrunOut = [name for name in os.listdir(os.getcwd())]

					for MDout in MDrunOut:
							if RTS in MDout:
								if RTS in Perf[samplename]:
									temp = (Utility.tail(MDout, 10)[6]).decode("utf-8")
									temp = temp
									temp = temp.split()
									if "Performance:" not in temp:
										print("Something went wrong with this one: see {0}".format(MDout))
										continue
									Perf[samplename][RTS].update( { 'ns/day' : temp[1]} )
									Perf[samplename][RTS].update( { 'hour/ns' : temp[2]} )

	############print(Perf)
	with Utility.cd("Analysis"): 
		XVGfiles = [ name for name in os.listdir(os.getcwd()) if name.endswith('.xvg') if not name.startswith('#')]
		OrderDATfiles = [ name for name in os.listdir(os.getcwd()) if name.endswith('order.dat') if not name.startswith('#')]
		
		#Initialize the DATA dictionary with the name of the samples studied
		for XVG in XVGfiles:
			SmplName = XVG[:XVG.find(suffix)]
			

		#Fill the DATA dictionary with the quantities and their associated columns in the .xvg files
		for XVG in XVGfiles:

			SmplName = XVG[:XVG.find(suffix)]
			#add the sample name to the list of WordsForCompletion
			WordsForCompletion.append(SmplName)
			ColumnIndex = 1

			with open(XVG) as data:
				Runtype = XVG[XVG.rfind('-')+1:XVG.rfind('.xvg')] #NVE, NPT, NVT
				Type = XVG[XVG.find(suffix)+len(suffix)+1 : XVG.find('_', XVG.find(suffix)+len(suffix)+1)]


				#add the run type to the list of WordsForCompletion
				WordsForCompletion.append(Runtype)
				#####print(SmplName+' '+Runtype)
				if DATA[SmplName] == None:
					DATA[SmplName] = {Runtype : {}}
					#####print('NONE')
				else:
					DATA[SmplName].update( { Runtype: {} } )
					#####print('Already created')

				# Create a temporary file with the comments removed
				TempOut = open('temp','a+')
				# Fill in DATA for Time
				DATA[SmplName][Runtype].update( { 'Time':0  } )
				#add the Qtty name to the list of WordsForCompletion
				WordsForCompletion.append('Time')

				for line in data:

					if line.startswith('#') or line.startswith('@'):
						#Extract the legend to get the Column/quantity match
						if 's{0}'.format(ColumnIndex-1) in line:
							Qtty = line[  line.find('"')+1:line.rfind('"')  ]
							Qtty = Qtty.replace(' ','-')
							DATA[SmplName][Runtype].update( { Qtty : ColumnIndex} )
							#add the Qtty name to the list of WordsForCompletion
							WordsForCompletion.append(Qtty)
							ColumnIndex += 1
						continue
					TempOut.write(line)
				TempOut.close()

			for Quantity in DATA[SmplName][Runtype]:
				######print(SmplName+' '+Runtype+' __ '+Quantity+' __ '+ str(DATA[SmplName][Runtype][Quantity])+' in '+XVG)
				DATA[SmplName][Runtype][Quantity] = np.loadtxt( 'temp', dtype='float', usecols= (int(DATA[SmplName][Runtype][Quantity]),) )
				#######print(DATA[SmplName])
			os.remove('temp')
		
		#Plot P2 the order parameter
		if OrderParam != '--':
			for ODF in OrderDATfiles:
				P2fig = plt.figure(num= 1, figsize= [12.80,7.20],dpi=100, facecolor='w')
				G_P2 = P2fig.add_subplot(111)
				
				Xaxislabels = []
				Yaxis = []
				
				with open(ODF) as data:
					Xaxislabels = data.readline().split()[1:]
					for line in data:
						if 'average' in line and 'order' not in line:
							Yaxis = line.split()[1:]
							print(line)
							
				Xaxis = [x for x in range(0,2*len(Xaxislabels),2)]
				print(Yaxis)
				G_P2.plot(Xaxis, Yaxis)
				G_P2.set_yscale('linear')
				G_P2.set_xticks(Xaxis)
				G_P2.set_xticklabels(Xaxislabels, rotation='vertical', fontsize=16)
				G_P2.set_xlabel('Bond')
				G_P2.set_ylabel('$P_{2}$',rotation=0, fontsize=16)
				G_P2.yaxis.set_label_coords(0, 1.01, transform=None)
				plt.tight_layout()
				plt.savefig(ODF.replace('dat','eps'),format='eps', dpi=100)
				
				P2fig.clf()




	#Compute the area per lipid and surfacetension
	for data in DATA:
		for Run in DATA[data]:
			if Types[data] == 'BILAYER':
				if Run.find('NPT') != -1:
					for Mol in DATA[data]['Mol']:
						if Mol in LipidsList:
							if 'Box-X' in DATA[data][Run] and 'Box-Y' in DATA[data][Run]:
								# 2*Box-X*Box-Y/nb(DSPC)
								ApL = 2.0*DATA[data][Run]['Box-X']*DATA[data][Run]['Box-Y']/DATA[data]['Mol'][Mol]

								DATA[data][Run].update( { 'Area-per-Lipid' : ApL } )
								#Appends the Area-per-Lipid to the list of keywords
								WordsForCompletion.append('Area-per-Lipid')

							if 'Box-Z' in DATA[data][Run] and 'Pres-ZZ' in DATA[data][Run] and 'Pres-YY' in DATA[data][Run] and 'Pres-XX' in DATA[data][Run]:
								# Lz*(Pzz - (Pxx+Pyy)/2)
								SurfTens = DATA[data][Run]['Box-Z']/4.0 * (  DATA[data][Run]['Pres-ZZ'] - (DATA[data][Run]['Pres-XX']+DATA[data][Run]['Pres-YY'])/2.0  )

								DATA[data][Run].update( {  'SurfTens' : SurfTens  } )
								#Appends the SurfTens to the list of keywords
								WordsForCompletion.append('SurfTens')

	#################print(Perf)
	#################print('\n\n\n')
	#################print(DATA)

	# Creates the html file summary with the fluctuations ratio, avg temperature, ...
	with Utility.cd('Analysis'):
		Avgfiles = [name for name in os.listdir(os.getcwd()) if name.endswith('avg')]
		Avgfiles.sort()

		for AVGF in Avgfiles:
			f_in = open(AVGF,'r')
			Runtype = AVGF.rsplit('-',1)[1].split('.avg')[0]
			Smpl = AVGF.split(suffix)[0]
			Type = AVGF.split(suffix+'_')[1].split('_')[0]
			######print(Type)
			#######print(Runtype)
			#######if Runtype in DATA[Smpl]: print("""{0} in {1}""".format(Runtype,Smpl))

			if Runtype.startswith('NVE') and Perf[Smpl][Runtype]:
				RMSDPot = 0.0
				DriftPot = 0.0
				RMSDTot = 0.0
				DriftTot = 0.0
				Temp = 0.0
				Temperr = 0.0
				DriftTemp = 0.0
				for line in f_in:
					if "Potential" in line :
						RMSDPot = float(line.split()[3])
						DriftPot = float(line.split()[4])
						continue
					if "Total Energy" in line :
						RMSDTot = float(line.split()[4])
						DriftTot = float(line.split()[5])
						continue
					if "Temperature" in line :
						Temp = float(line.split()[1])
						Temperr = float(line.split()[2])
						DriftTemp = float(line.split()[4])
						continue
					if "through" in line :
						BeginAvg = float(line.split('through')[0].rsplit('[')[1])
						EndAvg = float(line.split('through')[1].split('ps')[0])
						continue
					if "Kinetic En." in line :
						EkinAvg = float(line.split()[2])
						continue
				NVESAMPLE = str("""
								<tr>
									<th class="th-smpl" rowspan="3">{0}</th>
									<th class="th-smpl" colspan="2">T(K)</th>
									<th class="th-smpl" colspan="2">
										<math>
											<mfrac bevelled="True">
												<mrow>
													<msub>
														<mi>&Delta;E</mi>
														<mtext>tot</mtext>
													</msub>
												</mrow>
													<msub>
														<mi>&Delta;E</mi>
														<mtext>pot</mtext>
													</msub>
											</mfrac>
										</math>
									</th>
								<th class="th-smpl" rowspan="2" >Time for Avg (ps)</th>
								<th class="th-smpl" colspan="2" >Performance</th>
								</tr>
								<tr>
									<th class="th-smpl" >Average</th>
									<th class="th-smpl" >Drift</th>

									<th class="th-smpl" >RMSD ratio</th>
									<th class="th-smpl" >Drift ( in Ekin/microsec.) </th>

									<th class="th-smpl" >ns/day</th>
									<th class="th-smpl" >hour/ns</th>
								</tr>
								<tr>
									<td>{1} <mo>&pm;</mo> {2}</td>
									<td>{3}</td>
									<td>{4}</td>
									<td>{5}</td>
									<td>{6} - {7}</td>
									<td>{8}</td>
									<td>{9}</td>
								</tr>
								#""".format(Smpl,round(Temp,1),round(Temperr,1), round(DriftTemp,3), round(RMSDTot/RMSDPot,2), round(DriftTot/(EndAvg-BeginAvg)*1000000/EkinAvg,3), BeginAvg, EndAvg, Perf[Smpl][Runtype]['ns/day'], Perf[Smpl][Runtype]['hour/ns']))
				Runsmry[Runtype] += Utility.RemoveUnwantedIndent(NVESAMPLE)
				continue


			if Runtype.startswith('NPT') and Perf[Smpl][Runtype] and (Type == 'BILAYER' or Type == 'FREE_BILAYER'):
				Temp = 0.0
				Temperr = 0.0
				PZZ = 0.0
				PZZerr = 0.0
				for line in f_in:
					if "Temp" in line :
						Temp = float(line.split()[1])
						Temperr = float(line.split()[2])
						continue
					if "Pres-ZZ" in line :
						PZZ = float(line.split()[1])
						PZZerr = float(line.split()[2])
						continue
					if "through" in line :
						BeginAvg = float(line.split('through')[0].rsplit('[')[1])
						EndAvg = float(line.split('through')[1].split('ps')[0])
						continue
				NPTSAMPLE = str("""
									<tr>
										<th class="th-smpl" rowspan="3">{0}</th>
										<th class="th-smpl" rowspan="2" >T(K)</th>
										<th class="th-smpl" rowspan="2" >
											<math>
												<msub>
													<mn>P</mn>
													<mn>zz</mn>
												</msub>
											</math>
										</th>
									<th class="th-smpl" rowspan="2" >Area per lipid (<math><msup><mtext>nm</mtext> <mn>2</mn></msup></math>)</th>
									<th class="th-smpl" rowspan="2" ><math><mi>&gamma;</mi></math> (bar.nm)</th>
									<th class="th-smpl" rowspan="2" >Time for Avg (ps)</th>
									<th class="th-smpl" colspan="2" >Performance</th>
								</tr>
								<tr>
										<th class="th-smpl" >ns/day</th>
										<th class="th-smpl" >hour/ns</th>
								</tr>
								<tr>
									<td>{1} <mo>&pm;</mo> {2}</td>
									<td>{3} <mo>&pm;</mo> {4}</td>
									<td>{5}</td>
									<td>{6}</td>
									<td>{7} - {8}</td>
									<td>{9}</td>
									<td>{10}</td>
								</tr>
								""".format(Smpl, round(Temp,1), round(Temperr,2), round(PZZ,2), round(PZZerr,2), round(np.mean(DATA[Smpl][Runtype]['Area-per-Lipid']), 2), round(np.mean(DATA[Smpl][Runtype]['SurfTens']), 2), BeginAvg, EndAvg, Perf[Smpl][Runtype]['ns/day'], Perf[Smpl][Runtype]['hour/ns'] ) )
				Runsmry[Runtype] += Utility.RemoveUnwantedIndent(NPTSAMPLE)
				continue


			if Runtype.startswith('NPT') and Perf[Smpl][Runtype] and Type == 'SOLVENT' :
				Temp = 0.0
				Temperr = 0.0
				Potential = 0.0
				for line in f_in:
					if "Temp" in line :
						Temp = float(line.split()[1])
						Temperr = float(line.split()[2])
						continue
					if "Potential" in line :
						Potential = float(line.split()[1])
						continue
					if "through" in line :
						BeginAvg = float(line.split('through')[0].rsplit('[')[1])
						EndAvg = float(line.split('through')[1].split('ps')[0])
						continue
				NPTSAMPLE = str("""
									<tr>
										<th class="th-smpl" rowspan="3">{0}</th>
										<th class="th-smpl" rowspan="2" >T(K)</th>
										<th class="th-smpl" rowspan="2" >
											<math>
												<msub>
													<mn>E</mn>
													<mn>pot</mn>
												</msub>
											</math>
										</th>
									<th class="th-smpl" rowspan="2" >Time for Avg(ps)</th>
									<th class="th-smpl" colspan="2" >Performance</th>
								</tr>
								<tr>
										<th class="th-smpl" >ns/day</th>
										<th class="th-smpl" >hour/ns</th>
								</tr>
								<tr>
									<td>{1} <mo>&pm;</mo> {2}</td>
									<td>{3}</td>
									<td>{4} - {5}</td>
									<td>{6}</td>
									<td>{7}</td>
								</tr>
								""".format(Smpl, round(Temp,1), round(Temperr,2), round(Potential,2), BeginAvg, EndAvg, Perf[Smpl][Runtype]['ns/day'], Perf[Smpl][Runtype]['hour/ns'] ) )
				Runsmry[Runtype] += Utility.RemoveUnwantedIndent(NPTSAMPLE)
				continue

	for key in Runsmry: HTMLout.write(Runsmry[key])
	HTMLout.write(Utility.RemoveUnwantedIndent(HTMLfoot))
	HTMLout.close()

	#Removes the duplicates in WordsForCompletion
	list(set(WordsForCompletion))


	#***************************************************#
	#***************************************************#
	#			BEGIN INTERACTIVITY						#
	#***************************************************#
	#***************************************************#
	# Variables for plotting :
	# maximum y value for the plot
	ymin = None
	ymax = None
	#Condition to stay in menu
	Stay = True
	#
	if(sys.argv[1] == "--interactive"):

		# Word completer
		completer = Utility.Completer(WordsForCompletion)
		readline.parse_and_bind("tab: complete")
		readline.set_completer_delims(' ')
		readline.set_completer(completer.complete)

		while Stay:
			OnScreenMessage = str("""
					#******************************************************************************************************#
					#******************************************************************************************************#
					To analyze all samples use :
						all
					To analyze specified samples write their names as follows:
						<sample1> <sample2> ... <instructions>
					To select which MD run to analyze:
						-run RunName
					To plot quantities use :
						2D : -x <Qtty1> <Qtty2> ... -y <Qtty3> <Qtty4> ...
						3D : -x <Qtty1> <Qtty2> ... -y <Qtty3> <Qtty4> ... -z <Qtty4> <Qtty6> ... (not implemented)
					To compute Averages use:
						-avg <Qtty1> <Qtty2> ... option[-b <TimeToBegin> -e <TimeToEnd>] ...
					#******************************************************************************************************#
					#******************************************************************************************************#
					""")
			print(Utility.RemoveUnwantedIndent(OnScreenMessage))

			#*** Print the quantities to analyse ***#
			for Smpl in DATA:
				print(Smpl+" :")
				for Run in DATA[Smpl]:
					if Run != 'mol':
						print("  "+Run+" :", end='')
						for Qtty in DATA[Smpl][Run]:
							print(Qtty, end=' | ')
						print('')
				print()

			#*** User keyboard input ***#
			KbrdCmd = []
			while True:
				line = input()
				if line:
					KbrdCmd.append(line)
				else:
					break

			#*** Reads the keyboard input and execute the commands ***#
			for cmd in KbrdCmd:

				#If user input is 'exit', exit the program
				if "exit" in cmd:
					Stay = False
					print("Leaving")
					break

				# **********************************************************#
				#If user uses -x -y, plot graph in 3D of the chosen qtties  #
				# **********************************************************#
				if "-z" in cmd:
					print("-z in command -> to implement")

				# **********************************************************#
				#If user uses -x -y, plot graph in 2D of the chosen qtties  #
				# **********************************************************#
				if "-x" in cmd or "-y" in cmd:
					X = cmd.split("-x")[1].split("-y")[0].strip().split(' ')
					Y = cmd.rsplit("-y")[1].strip().split(' ')

					#If user uses -sameplot option plot all the qtties on the same plot  #
					if '-sameplot' in cmd:

						for x in X:
							ymin = None
							ymax = None
							plotname = ''
							fig = plt.figure()
							ax = fig.add_subplot(111)

							for y in Y:

								plotname += y + ', '
								#If user input starts with all, repeat this command for each sample
								if "all" in cmd:
									#Reads the name of the run from the cmd
									Run = cmd.split("-run")[1].rsplit(' ')[1].strip()

									for data in DATA:
										######print(data+' with '+Run)
										if Run in DATA[data]:
											ax.plot(DATA[data][Run][x], DATA[data][Run][y],linewidth='2',label=data+' '+y)
											ax.set_xlabel(x)

											# Set ymin and ymax to the global extremum amongst all sets
											if( (ymin is None) or (ymax is None)):
												ymin = np.amin( DATA[data][Run][y] )
												ymax = np.amax( DATA[data][Run][y] )
											if(ymin > np.amin( DATA[data][Run][y] )): ymin = np.amin( DATA[data][Run][y] )
											if(ymax < np.amax( DATA[data][Run][y] )): ymax = np.amax( DATA[data][Run][y] )
								else:
									#Reads the name of the run from the cmd
									Run = cmd.split("-run")[1].rsplit(' ')[1].strip()
									Smpls = cmd.split("-run")[0].strip().split(' ')
									for smpl in Smpls:
										#####print(smpl+' with '+Run)
										if Run in DATA[smpl]:
											ax.plot(DATA[smpl][Run][x], DATA[smpl][Run][y],linewidth='2',label=smpl+' '+y)
											ax.set_xlabel(x)

											# Set ymin and ymax to the global extremum amongst all sets
											if( (ymin is None) or (ymax is None)):
												ymin = np.amin( DATA[smpl][Run][y] )
												ymax = np.amax( DATA[smpl][Run][y] )
											if(ymin > np.amin( DATA[smpl][Run][y] )): ymin = np.amin( DATA[smpl][Run][y] )
											if(ymax < np.amax( DATA[smpl][Run][y] )): ymax = np.amax( DATA[smpl][Run][y] )


							plotname += 'as a function of '+ x
							ax.set_title(plotname)
							ax.set_ylim( [ymin, ymax] )
							ax.legend()
							plt.show()
							#plt.savefig("all__{0}_asfunc_{1}.jpg".format(y,x))

					else:
						for x in X:
							ymin = None
							ymax = None
							for y in Y:

								fig = plt.figure()
								ax = fig.add_subplot(111)

								#If user input starts with all, repeat this command for each sample
								if "all" in cmd:
									#Reads the name of the run from the cmd
									Run = cmd.split("-run")[1].rsplit(' ')[1].strip()

									for data in DATA:
										##########print(data+' with '+Run)
										if Run in DATA[data]:
											ax.plot(DATA[data][Run][x], DATA[data][Run][y],linewidth='2',label=data)
											ax.set_title(y + ' as a function of ' + x + ' for ' + Run)
											ax.set_xlabel(x)
											ax.set_ylabel(y)

											# Set ymin and ymax to the global extremum amongst all sets
											if( (ymin is None) or (ymax is None)):
												ymin = np.amin( DATA[data][Run][y] )
												ymax = np.amax( DATA[data][Run][y] )
											if(ymin > np.amin( DATA[data][Run][y] )): ymin = np.amin( DATA[data][Run][y] )
											if(ymax < np.amax( DATA[data][Run][y] )): ymax = np.amax( DATA[data][Run][y] )

								else:
									#Reads the name of the run from the cmd
									Run = cmd.split("-run")[1].rsplit(' ')[1].strip()
									Smpls = cmd.split("-run")[0].strip().split(' ')
									for smpl in Smpls:
										#####print(smpl+' with '+Run)
										if Run in DATA[smpl]:
											ax.plot(DATA[smpl][Run][x], DATA[smpl][Run][y],linewidth='2',label=smpl)
											ax.set_title(y + ' as a function of ' + x + ' for ' + Run)
											ax.set_xlabel(x)
											ax.set_ylabel(y)

											# Set ymin and ymax to the global extremum amongst all sets
											if( (ymin is None) or (ymax is None)):
												ymin = np.amin( DATA[smpl][Run][y] )
												ymax = np.amax( DATA[smpl][Run][y] )
											if(ymin > np.amin( DATA[smpl][Run][y] )): ymin = np.amin( DATA[smpl][Run][y] )
											if(ymax < np.amax( DATA[smpl][Run][y] )): ymax = np.amax( DATA[smpl][Run][y] )

								ax.set_ylim( [ymin, ymax] )
								ax.legend()
								plt.show()
								#plt.savefig("all__{0}_asfunc_{1}.jpg".format(y,x))


				# **********************************************************#
				#If user uses -avrg, compute the average and Standard Error #
				# **********************************************************#
				if "-avg" in cmd:
					#Get the name of the run to study
					Run = cmd.split("-run")[1].split("-avg")[0].strip()
					Qtties = cmd.split("-avg")[1].strip().split(' ')

					# Case where no beginning nor end is provided ***************

					#Analyze all samples
					if 'all' in cmd:
						for Qtty in Qtties:
							for data in DATA:
								Averages.update( {data : {}} )
								if Run in DATA[data]:
									Averages[data].update( {Run : {}})
									temp = Utility.ComputeAverage(DATA[data][Run][Qtty])

									Avg = np.mean(temp[:,0])
									Var = temp[:,1]/temp[:,2]
									Varerr = np.sqrt(2*(temp[:,1])**2/(temp[:,2])**3)

									x = np.linspace(0,len(Var)-1,len(Varerr))

									EcT = np.sqrt(Var)
									EcTerror = np.sqrt(Varerr)

									plt.figure(0)
									plt.xlabel("Number of blockings")
									plt.xlim(0,len(x) +1)
									plt.ylabel('$\sigma$', fontsize=16)
									plt.errorbar(x, EcT, yerr=EcTerror,fmt='-o', label = data)
						plt.title("$\sigma_{{{0}}}$ as a function of the number of blocking operation".format(Qtty))
						plt.legend()
						plt.show()

					#Analyze the provided samples
					else:
						Smpls = cmd.split("-run")[0].strip().split(' ')
						for smpl in Smpls:
							Averages.update( {smpl : {}} )
							if Run in DATA[smpl]:
								Averages[smpl].update( {Run : {}})
								for Qtty in Qtties:
									temp = Utility.ComputeAverage(DATA[Smpl][Run][Qtty])

									Avg = np.mean(temp[:,0])
									Var = temp[:,1]/temp[:,2]
									Varerr = np.sqrt(2*(temp[:,1])**2/(temp[:,2])**3)

									x = np.linspace(0,len(Var)-1,len(Varerr))

									EcT = np.sqrt(Var)
									EcTerror = np.sqrt(Varerr)

									plt.figure(0)
									plt.xlabel("Number of blockings")
									plt.xlim(0,len(x) +1)
									plt.ylabel('$\sigma$', fontsize=16)
									plt.errorbar(x, EcT, yerr=EcTerror,fmt='-o', label = data)
									plt.legend()
									plt.show()

