import subprocess as sub
import time
import numpy as np
import math
import Utility
import os
import glob

#List of Lipids and Solvents:
LipidsList = ['DSPC','DPPC','DLPC']
SolventsList = ['W','OCO','PW']

#***********************************************************#
#***********************************************************#
#********************* PDB Files list **********************#
#***********************************************************#
#***********************************************************#
PDBfileList = {'W':'water_single.pdb','OCO':'octanol_single.pdb',
			   'PW':'polwater_single.pdb','DSPC':'dspc_single.pdb',
			   'DPPC':'dppc_single.pdb','DLPC':'dlpc_single.pdb','SU':'su_single.pdb'}

W_PDB = """
		CRYST1   04.000   04.000   04.000  90.00  90.00  90.00 P1           1
		ATOM      1  W   W   X   1      00.000  00.000  00.000 0.00  0.00
		END
		"""

OCO_PDB = """
		CRYST1   10.000   10.000   10.000  90.00  90.00 90.00 P 1           1
		ATOM      1  PC  OCO     1      5.000  5.000  7.865 0.00  0.00          
		ATOM      2  C   OCO     1      5.000  5.000  3.135 0.00  0.00          
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
		ATOM      1  TEMPSU  S   1      00.000  00.000  00.000 0.00  0.00
		END
		"""
		
#***********************************************************#
#***********************************************************#
#***********************************************************#




















def InitBilayer(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault):
	packmolSeed = 0
	if Sample['SEED'] == 'time':
		packmolSeed = int(time.time())
	if Sample['SEED'] == 'jobnum':
		packmolSeed = Sample['JOBNUM'] * 123456789
		
	if 'DEFO' in Sample and not 'SU' in Sample:
		return InitBilayerWithHole(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed)

	if 'SU' in Sample and not 'DEFO' in Sample: #Bilayer + Wall:
		return InitBilayerWithWall(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed)

	if 'DEFO' in Sample and 'SU' in Sample: #Bilayer + Hole + Wall:
		return InitBilayerWithHoleAndWall(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed)

	else: #Default for Free Bilayer
		return InitFreeBilayer(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed)
		
		
		
		
def InitFreeBilayer(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed):
	# packmol constraints are not strict !
	# put the molecules in a smaller box to
	# avoid initial crash because of PBC
	# the SHELL parameter defines the shell in which
	# packmol is not supposed to put particles
	SHELL = 3.0
	
	LXS = Sample['LX']-SHELL
	LYS = Sample['LY']-SHELL
	LZS = Sample['LZ']-SHELL

	LZ2 = LZS/2.0

	#Number of lipid and water per layer
	NLM = {}
	NSOLM = {}
	LipidType = ''
	
	for L in LipidsList:
		if(L in Sample):
			NLM.update( {L:int(Sample[L]/2)} )
			LipidType = L
	for Sol in SolventsList:
		print(Sol,' is in sample:', Sol in Sample)
		if Sol in Sample:
			print(Sol)
			NSOLM.update( {Sol:int(Sample[Sol]/2)} )
			SolventType = Sol
	
	#If the number of lipids per Monolayer is high
	# we cut the sample till we get a reasonnable amount
	# of lipis per monolayer for packmol
	NbMult = 0
	NeedToExpand = False
	while NLM[LipidType] > 512:
		NLM[LipidType] /= 2
		NSOLM[SolventType] /= 2
		NbMult += 1
		#diviser boîte
		NeedToExpand = True
		print(NLM[LipidType])
	print(LipidType)
	#==================================================================================
	# creating the initial bilayer (lipids + water) using packmol
	#==================================================================================

	System = str('{0}_{1}{2}_{3}{4}').format(Sample['TYPE'],Sample[LipidType], LipidType, Sample[SolventType], SolventType)
	
	#=======================================================================
	# Setting the bilayers  ================================================
	#=======================================================================
	
	TMT = 0.0
	DZ = 0.0
	# DSPC bilayer =========================================================
	if LipidType == "DSPC":
		# Geometry to preprare the bilayer with packmol
		# total monolayer thickness (Angstrom)
		TMT = 30.0
		# deltaz : thickness in Z direction for
		# the volume constraining the heads and tails beads
		DZ = 7.0
	
	# DPPC bilayer =========================================================
	if LipidType == "DPPC":
		TMT = 30
		DZ = 10
	
	# DLPC bilayer =========================================================
	if LipidType == "DLPC":
		TMT = 30
		DZ = 10

	M1headMIN = LZ2 - TMT
	M1headMAX = LZ2 - TMT + DZ
	M1tailMIN = LZ2 - DZ
	M1tailMAX = LZ2

	M2headMIN = LZ2 + TMT - DZ
	M2headMAX = LZ2 + TMT
	M2tailMIN = LZ2
	M2tailMAX = LZ2 + DZ
	
	#=======================================================================
	# Setting the bilayers for packmol =====================================
	#=======================================================================
	if packmolSeed:
		PackmolInput = """
					#Packmol seed was set using {0}
					seed {1}
				
				""".format(Sample['SEED'], packmolSeed)
	else:
		PackmolInput = """
					#Packmol Seed was set using default
				
							"""
	#DSPC BILAYER
	if LipidType == "DSPC":
		PackmolInput += """
						#
						# Lipid bilayer with water, perpendicular to z axis
						#

						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.

						tolerance 3.0

						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)

						filetype pdb

						# The output pdb file

						output {0}.pdb

						structure {14}
							chain A
							resnumbers 3
							number {1:g}
							inside box 0. 0. {2}  {3} {4} {5}
							atoms 1
								below plane 0. 0. 1. {6}
							end atoms
							atoms 9 14
								over plane 0. 0. 1. {7}
							end atoms
						end structure

						structure {14}
							chain B
							resnumbers 3
							number {1:g}
							inside box 0. 0.  {8}  {3} {4} {9}
							atoms 1
								over plane 0. 0. 1. {10}
							end atoms
							atoms 9 14
								below plane 0. 0. 1. {11}
							end atoms
						end structure

						""".format(System, NLM[LipidType], M1headMIN, LXS, LYS, M1tailMAX, M1headMAX, M1tailMIN, M2tailMIN, M2headMAX, M2headMIN, M2tailMAX, NSOLM[SolventType], LZS, PDBfileList[LipidType], PDBfileList[SolventType])
	
	# DPPC bilayer, DLPC bilayer =========================================================
	if LipidType in {'DPPC','DLPC'}:
		PackmolInput += """
						#
						# Lipid bilayer with water, perpendicular to z axis
						#

						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.

						tolerance 3.0

						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)

						filetype pdb

						# The output pdb file

						output {0}.pdb

						structure {14}
							chain A
							resnumbers 3
							number {1:g}
							inside box 0. 0. {2}  {3} {4} {5}
							constrain_rotation x 180 10
							constrain_rotation y 180 10
						end structure

						structure {14}
							chain B
							resnumbers 3
							number {1:g}
							inside box 0. 0.  {8}  {3} {4} {9}
							constrain_rotation x 0 10
							constrain_rotation y 0 10
						end structure

						""".format(System, NLM[LipidType], M1headMIN, LXS, LYS, M1tailMAX, M1headMAX, M1tailMIN, M2tailMIN, M2headMAX, M2headMIN, M2tailMAX, NSOLM[SolventType], LZS, PDBfileList[LipidType], PDBfileList[SolventType])
							
	#=======================================================================
	#=======================================================================
	#=======================================================================
	
	#=======================================================================
	# Setting the solvent for packmol ============================
	#=======================================================================
	if SolventType == 'PW':
		PackmolInput += """
		
						structure {0}
							chain C
							number {1:g}
							inside box 0. 0. 0.  {2} {3} {4}
							atoms 2 3
								radius 0.2
							end atoms
						end structure


						structure {0}
							chain D
							number {1:g}
							inside box 0. 0. {5}  {2} {3} {6}
							atoms 2 3
								radius 0.2
							end atoms
						end structure
						""".format(PDBfileList[SolventType], NSOLM[SolventType], LXS, LYS,  M1headMIN,
							M2headMAX, LZS)
	else:
		PackmolInput += """
						
						structure {0}
							chain C
							number {1:g}
							inside box 0. 0. 0.  {2} {3} {4}
						end structure
						
						structure {0}
							chain D
							number {1:g}
							inside box 0. 0. {5}  {2} {3} {6}
						end structure
						""".format(PDBfileList[SolventType], NSOLM[SolventType], LXS, LYS,  M1headMIN,
							M2headMAX, LZS)

	f = open('packmol_'+System+'.input','w')
	f.write(Utility.RemoveUnwantedIndent(PackmolInput))
	f.close()
	
	# Lipids input pdb =====================================================
	if( LipidType == "DSPC"):
		f = open('dspc_single.pdb','w')
		f.write(Utility.RemoveUnwantedIndent(DSPC_PDB))
		f.close()
	if( LipidType == "DPPC"):
		f = open('dppc_single.pdb','w')
		f.write(Utility.RemoveUnwantedIndent(DPPC_PDB))
		f.close()
	if( LipidType == "DLPC"):
		f = open('dlpc_single.pdb','w')
		f.write(Utility.RemoveUnwantedIndent(DLPC_PDB))
		f.close()
		
	# Solvents input pdb ===================================================
	if( SolventType == 'W'):
		f = open('water_single.pdb','w')
		f.write(Utility.RemoveUnwantedIndent(W_PDB))
		f.close()
	if( SolventType == 'PW'):
		f = open('polwater_single.pdb','w')
		f.write(Utility.RemoveUnwantedIndent(PW_PDB))
		f.close()
	if( SolventType == 'OCO'):
		f = open('octanol_single.pdb','w')
		f.write(Utility.RemoveUnwantedIndent(OCO_PDB))
		f.close()
	

	#===============================================================================
	#lauching packmol
	#===============================================================================

	cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """).format(Softwares['PACKMOL'],System)
	sub.call(cmd, shell=True)

	
	#===============================================================================
	# ensure the right box in pdb file
	#===============================================================================

	WriteBox = str("""
			mol load pdb {0}.pdb
			set all [atomselect top "all"]

			package require pbctools
			pbc set {{{1} {2} {3}}}

			$all writepdb "{0}.withbox.pdb"
			unset all

			exit
			""").format(System, Sample['LX'], Sample['LY'], Sample['LZ'])

	f = open('write_box.vmd','w+')
	f.write(Utility.RemoveUnwantedIndent(WriteBox) )
	f.close()

	## ======================================================================
	cmd = str("""{0} -dispdev text -e write_box.vmd > write_box.log""" 
				).format(Softwares['VMD'])
	sub.call(cmd, shell=True)

	
	#===============================================================================
	# For big samples
	#===============================================================================
	if NeedToExpand:
		cmd = str("""{0}genconf -f {1}.withbox.pdb -nbox {2} {2} 1 -o {1}.withbox.pdb""").format(GROMACS_LOC_prefixPath, System, NbMult)
		sub.call(cmd, shell=True)
	## ======================================================================
	ApL = Sample['LX']*Sample['LY']/NLM[LipidType]
	print(Utility.RemoveUnwantedIndent(str("""
			================================
			================================
			Packmol finished initial input file
			{0} {1}, {2} {3}
			box sizes : {4}, {5}, {6}
			Area per lipid = {7} A**2
			===============================
			===============================
			""").format(LipidType, Sample[LipidType],SolventType, Sample[SolventType], Sample['LX'],Sample['LY'],Sample['LZ'],ApL)))

	#==================================================================================
	# the topology file
	#==================================================================================

	Topology = str("""
				#include "martini_v2.2.itp" ; modified with polarisable water
				#include "martini_v2.0_lipids.itp"
				""")
	#Copy the topology files for martini forcefield
	sub.call("""cp {0}/martini_v2.0_lipids.itp {0}/martini_v2.2.itp ./""".format(PathToDefault), shell= True)
	f = open(System+'.top','w')
	f.write(Utility.RemoveUnwantedIndent(Topology))

	f = open(System+'.top','a')
	f.write("""\n[ system ]\n""")

	Topology = str("""{0} {1}\n""").format(LipidType, Sample['TYPE'])
	f.write(Topology)

	###add a for loop for multiple types (Later)
	f.write("""\n[ molecules ]\n""")
	Topology = str("""{0} {1}\n{2} {3}""").format(LipidType, Sample[LipidType],SolventType, Sample[SolventType])
	f.write(Topology)
	f.close()

	#==================================================================================
	# the index file
	#==================================================================================

	cmd = str("""echo q | {0}make_ndx -f {1}.withbox.pdb -o {1}.ndx""").format(GROMACS_LOC_prefixPath,System)
	sub.call(cmd, shell=True)

	#==================================================================================
	# Output the files for other steps
	#==================================================================================
	Output = str("""{0}.withbox.pdb""").format(System)
	Index = str("""{0}.ndx""").format(System)
	return { 'SYSTEM': System, 'OUTPUT': Output, 'INDEX':Index}



def InitBilayerWithWall(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed):
	THICKNESS = float(Sample['SU']['Thickness'])
	DENSITY = float(Sample['SU']['Density'])
	NBLIPIDS_MONO = int(Sample['SU']['NbLipidsM'])
	SU_VERSION = Sample['SU']['Version']
	SU_TYPE = Sample['SU']['SuType']
	LZM = float(Sample['SU']['LzM'])
	
	
	SHELL = 3.0
	LXS = Sample['LX'] - SHELL
	LYS = Sample['LY'] - SHELL

	

	#Number of lipid and water per layer
	NLM = {}
	NSOLM = {}
	LipidType = ''
	
	for L in LipidsList:
		if(L in Sample):
			NLM.update( {L:int(Sample[L]/2)} )
			LipidType = L
			
	for Sol in SolventsList:
		if Sol in Sample:
			NSOLM.update( {Sol:int(Sample[Sol]/2)} )
			SolventType = Sol
			
	#Computing the number of SU using Density in CG/nm³ thus the divided by 1000 for Volume
	nbSu = int( DENSITY * Sample['LX'] * Sample['LY'] * THICKNESS /1000 )
	
	#=======================================================================
	# Setting the bilayers  ================================================
	#=======================================================================
	
	TMT = 0.0
	DZ = 0.0
	# DSPC bilayer =========================================================
	if LipidType == "DSPC":
		# Geometry to preprare the bilayer with packmol
		# total monolayer thickness (Angstrom)
		TMT = 30.0
		# deltaz : thickness in Z direction for
		# the volume constraining the heads and tails beads
		DZ = 7.0
	
	# DPPC bilayer =========================================================
	if LipidType == "DPPC":
		TMT = 30.0
		DZ = 10
	
	# DLPC bilayer =========================================================
	if LipidType == "DLPC":
		TMT = 30.0
		DZ = 10
	
	assert(TMT > 0.0 and DZ > 0.0),"Your Parameters.csv contains a lipid not yet defined in Prepare.py"
	#Set the box height with the monolayer + shell for PBC
	LZ = LZM + TMT + SHELL
	LBZ = (LZM-THICKNESS)/2.0 + THICKNESS
	LZ2 = 0.0
	NbParticlesToRelocate = 0
	NbSolTop = NSOLM[SolventType]
	NbSolBottom = NSOLM[SolventType]
	if 'BilayerHeight' in Sample['SU']:
		LZ2 = float(Sample['SU']['BilayerHeight'])
		#Checking that the Bilayer lipids do not overlap with the Monolayer lipids and the substrate
		if (LZ2 - TMT) < THICKNESS:
			LZ2 -= LZ2 -TMT - THICKNESS
		if (LZ2 +TMT) > LZM:
			LZ2 += LZM - LZ2 - TMT
		#Checking that the Solvent is not too dense
		if LZ2 < LBZ: #Too much Solvent below
			DensitySol = NSOLM[SolventType]/Sample['LX']/Sample['LY']/(LBZ - TMT - THICKNESS)
			VolToRemove = LBZ - LZ2
			VolToRemove *= LXS*LYS
			NbParticlesToRelocate = int(DensitySol * VolToRemove)
			
			NbSolBottom -= NbParticlesToRelocate
			NbSolTop += NbParticlesToRelocate
			
		if LZ2 > LBZ: #Too much Solvent above
			DensitySol = NSOLM[SolventType]/LXS/LYS/(LZM - LBZ - TMT)
			VolToRemove = LZ2 - LBZ
			VolToRemove *= LXS*LYS
			NbParticlesToRelocate = int(DensitySol * VolToRemove)
			
			NbSolBottom += NbParticlesToRelocate
			NbSolTop -= NbParticlesToRelocate
	else:
		LZ2 = LBZ
	
	#Set the Geometry
	M1headMIN = LZ2 - TMT
	M1headMAX = LZ2 - TMT + DZ
	M1tailMIN = LZ2 - DZ
	M1tailMAX = LZ2

	M2headMIN = LZ2 + TMT - DZ
	M2headMAX = LZ2 + TMT
	M2tailMIN = LZ2
	M2tailMAX = LZ2 + DZ
	
	#=======================================================================
	#=======================================================================
	#=======================================================================

	System = str('{0}_{1}{2}_{3}M{2}_{4}{5}_{7}{6}{8}').format(Sample['TYPE'], Sample[LipidType],
											LipidType, NBLIPIDS_MONO , 
											Sample[SolventType], SolventType, SU_TYPE, nbSu, SU_VERSION)
	
	#=======================================================================
	# Setting the bilayers for packmol =====================================
	#=======================================================================
	if packmolSeed:
		PackmolInput = """
					#Packmol seed was set using {0}
					seed {1}
				
				""".format(Sample['SEED'], packmolSeed)
	else:
		PackmolInput = """
					#Packmol Seed was set using default
				
							"""
							
	# DLPC bilayer =========================================================
	if LipidType == "DSPC":
		PackmolInput += """
					#
					# Lipid bilayer with water, perpendicular to z axis
					#
					
					# Every atom from diferent molecules will be far from each other at
					# least 3.0 Anstroms at the solution.
					
					tolerance 3.0
					
					# Coordinate file types will be in pdb format (keyword not required for
					# pdb file format, but required for tinker, xyz or moldy)
					
					filetype pdb
					
					# do not avoid the fixed molecules
					#avoid_overlap no
					
					# The output pdb file
					
					output {0}.pdb
					
					structure {12}
						chain A
						resnumbers 3
						number {1:g}
						inside box 0. 0. {2}  {3} {4} {5}
						atoms 1
							below plane 0. 0. 1. {6}
						end atoms
						atoms 9 14
							over plane 0. 0. 1. {7}
						end atoms
					end structure
					
					structure {12}
						chain B
						resnumbers 3
						number {1:g}
						inside box 0. 0.  {8}  {3} {4} {9}
						atoms 1
							over plane 0. 0. 1. {10}
						end atoms
						atoms 9 14
							below plane 0. 0. 1. {11}
						end atoms
					end structure
					
					""".format(System, NLM[LipidType], M1headMIN, LXS,
								LYS, M1tailMAX, M1headMAX, M1tailMIN, 
								M2tailMIN, M2headMAX, M2headMIN, M2tailMAX, 
								PDBfileList[LipidType])
					
		PackmolInput += """
					
					structure {0}
						chain C
						resnumbers 3
						number {1:g}
						inside box 0. 0.  {2}  {3} {4} {5}
						atoms 1
							below plane 0. 0. 1. {6}
						end atoms
						atoms 9 14
							over plane 0. 0. 1. {7}
						end atoms
					end structure
					
					""".format(PDBfileList[LipidType], NBLIPIDS_MONO, LZM, LXS, LYS, TMT+LZM, LZM+DZ, TMT+LZM-DZ)
	
	# DPPC bilayer, DLPC bilayer =========================================================
	if LipidType == 'DPPC' or LipidType == 'DLPC':
		PackmolInput += """
						#
						# Lipid bilayer with water, perpendicular to z axis
						#

						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.

						tolerance 3.0

						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)

						filetype pdb

						# The output pdb file

						output {0}.pdb

						structure {8}
							chain A
							resnumbers 3
							number {1:g}
							inside box 0. 0. {2}  {3} {4} {5}
							constrain_rotation x 180 10
							constrain_rotation y 180 10
						end structure

						structure {8}
							chain B
							resnumbers 3
							number {1:g}
							inside box 0. 0.  {6}  {3} {4} {7}
							constrain_rotation x 0 10
							constrain_rotation y 0 10
						end structure

						""".format(System, NLM[LipidType], M1headMIN, LXS, 
									LYS, M1tailMAX, M2tailMIN, M2headMAX,
									PDBfileList[LipidType])
		
		PackmolInput += """
						
						structure {0}
							chain C
							resnumbers 3
							number {1:g}
							inside box 0. 0.  {2}  {3} {4} {5}
							constrain_rotation x 180 10
							constrain_rotation y 180 10
						end structure
						
						""".format(PDBfileList[LipidType], NBLIPIDS_MONO, LZM, LXS, LYS, TMT+LZM)
	
	#=======================================================================
	#=======================================================================
	#=======================================================================

	#=======================================================================
	# Setting the defo with solvent for packmol ============================
	#=======================================================================
	#Computing the number of Solvent vapor
	NbSolVapor = int( 0.007 * Sample['LX'] * Sample['LY'] * (Sample['LZ']-LZM-TMT) /1000)
	
	if SolventType == 'PW':
		if NbSolBottom:
			PackmolInput += """
							
							structure {0}
								chain D
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType],NbSolBottom,THICKNESS,LXS,LYS,M1headMIN)
			
		if NbSolTop:
			PackmolInput += """
							
							structure {0}
								chain E
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType],NbSolTop,M2headMAX,LXS,LYS,LZM)
		PackmolInput += """
							
							structure {0}
								chain V
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType],NbSolVapor,LZM+TMT,LXS,LYS,Sample['LZ'])
	else:
		if NbSolBottom:
			PackmolInput += """
							
							structure {0}
								chain D
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
							end structure
							""".format(PDBfileList[SolventType],NbSolBottom,THICKNESS,LXS,LYS,M1headMIN)
		if NbSolTop:
			PackmolInput += """
							
							structure {0}
								chain E
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
							end structure
							""".format(PDBfileList[SolventType],NbSolTop,M2headMAX,LXS,LYS,LZM)
		PackmolInput += """
						
						structure {0}
							chain V
							number {1:g}
							inside box 0. 0. {2}  {3} {4} {5}
						end structure
						""".format(PDBfileList[SolventType],NbSolVapor,LZM+TMT,LXS,LYS,Sample['LZ'])
					
	
		
	#=======================================================================
	# Adding the substrat ==================================================
	#=======================================================================
	PackmolInput += """
					
					structure {0}
						chain S
						number {1:g}
						inside box 0. 0. 0.  {2} {3} {4}
					end structure
					""".format(PDBfileList['SU'], nbSu, LXS, LYS, THICKNESS)
	
	f = open('packmol_'+System+'.input','w')
	f.write(Utility.RemoveUnwantedIndent(PackmolInput))
	f.close()
	
	# Lipids input pdb =====================================================
	if( LipidType == "DSPC"):
		f = open(PDBfileList['DSPC'],'w')
		f.write(Utility.RemoveUnwantedIndent(DSPC_PDB))
		f.close()
	if( LipidType == "DPPC"):
		f = open(PDBfileList['DPPC'],'w')
		f.write(Utility.RemoveUnwantedIndent(DPPC_PDB))
		f.close()
	if( LipidType == "DLPC"):
		f = open(PDBfileList['DLPC'],'w')
		f.write(Utility.RemoveUnwantedIndent(DLPC_PDB))
		f.close()
		
	# Solvents input pdb ===================================================
	if( SolventType == 'W'):
		f = open(PDBfileList['W'],'w')
		f.write(Utility.RemoveUnwantedIndent(W_PDB))
		f.close()
	if( SolventType == 'PW'):
		f = open(PDBfileList['PW'],'w')
		f.write(Utility.RemoveUnwantedIndent(PW_PDB))
		f.close()
	if( SolventType == 'OCO'):
		f = open(PDBfileList['OCO'],'w')
		f.write(Utility.RemoveUnwantedIndent(OCO_PDB))
		f.close()
		
	# pdb file for su
	suType = SU_TYPE
	assert(len(SU_TYPE) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
	suType = SU_TYPE + (4 - len(SU_TYPE))*' '
	suPdb = SU_PDB.replace("TEMP", suType)
	f = open(PDBfileList['SU'],'w')
	f.write(Utility.RemoveUnwantedIndent(suPdb))
	f.close()
		
		
	#=======================================================================
	#lauching packmol
	#=======================================================================

	cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """).format(Softwares['PACKMOL'], System)
	sub.call(cmd, shell=True)

	#=======================================================================
	# ensure the right box in pdb file
	#=======================================================================

	WriteBox = str("""
			mol load pdb {0}.pdb
			set all [atomselect top "all"]

			package require pbctools
			pbc set {{{1} {2} {3}}}

			$all writepdb "{0}.withbox.pdb"
			unset all

			exit
			""").format(System, Sample['LX'], Sample['LY'], Sample['LZ'])

	f = open('write_box.vmd','w+')
	f.write(Utility.RemoveUnwantedIndent(WriteBox) )
	f.close()

	## ======================================================================
	cmd = str("""{0} -dispdev text -e write_box.vmd > write_box.log"""
											).format(Softwares['VMD'])
	sub.call(cmd, shell=True)
	
	MakeIndex = str("""
			chain A
			chain B
			chain C
			chain D
			chain E
			chain V
			chain S
			name 5 bottom{0}
			name 6 top{0}
			name 7 mono{0}
			name 8 bottom{1}
			name 9 top{1}
			name 10 vap{1}
			name 11 su
			5 | 6
			name 12 bi{0}
			q
			
			""").format(LipidType, SolventType)

	f = open('make_index.input','w+')
	f.write(Utility.RemoveUnwantedIndent(MakeIndex) )
	f.close()
	
	cmd = str("""{0}make_ndx -f {1}.withbox.pdb -o {1}.ndx < make_index.input""").format(GROMACS_LOC_prefixPath,
																					System)
	sub.call(cmd, shell=True)
	
	
	#=======================================================================
	# the topology file
	#=======================================================================
	TopologySu(Sample, PathToDefault)
	
	Topology = """
				#include "martini_v2.2_{1}_{0}.itp"
				#include "martini_v2.0_lipids.itp"
				
				""".format(Sample['SU']['Version'], SU_TYPE)
	
	#Copy the topology files for martini forcefield
	sub.call("""cp {0}/martini_v2.0_lipids.itp  .""".format(PathToDefault), shell= True)
	sub.call("""cp {0}/SU/su_posres.itp  .""".format(PathToDefault), shell= True)
	f = open(System+'.top','w')
	f.write(Utility.RemoveUnwantedIndent(Topology))

	f = open(System+'.top','a')
	f.write("""\n[ system ]\n""")

	Topology = str("""{0}_{1}_WITH_{3}_{2}\n""").format(LipidType, Sample['TYPE'], Sample['SU']['Version'],
													SU_TYPE)
	f.write(Topology)

	f.write("""\n[ molecules ]\n""")
	Topology = str("""{0} {1}\n{2} {3}\n{4} {5}\n""").format(LipidType, Sample[LipidType]+NBLIPIDS_MONO,
														SolventType, Sample[SolventType], SU_TYPE, nbSu)
	f.write(Topology)
	f.close()
		
	#==================================================================================
	# Output the files for other steps
	#==================================================================================
	Output = str("""{0}.withbox.pdb""").format(System)
	Index = str("""{0}.ndx""").format(System)
	
	ApL = (Sample['LX']*Sample['LY'])/NLM[LipidType]
	print(Utility.RemoveUnwantedIndent("""
			================================
			================================
			Packmol finished initial input file
			{0} {1}, {2} {3} with SU{4} {5} and mono{0} {6}
			box sizes : {7}, {8}, {9}
			Area per lipid = {10} A**2
			===============================
			===============================
			""".format(LipidType, Sample[LipidType],SolventType, Sample[SolventType],
				Sample['SU']['Version'],nbSu,NBLIPIDS_MONO,
				Sample['LX'],Sample['LY'],Sample['LZ'],ApL)))
	return { 'SYSTEM': System, 'OUTPUT': Output, 'INDEX':Index}



def InitBilayerWithHole(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed):
	#Bilayer + Hole
		SHELL = 3.0
		
		LXS = Sample['LX']-SHELL
		LYS = Sample['LY']-SHELL
		LZS = Sample['LZ']-SHELL

		LZ2 = LZS/2.0
		
		#Number of lipid and water per layer
		NLM = {}
		NSOLM = {}
		
		LipidType = ''
		for L in LipidsList:
			if(L in Sample):
				NLM.update( {L:int(Sample[L]/2)} )
				LipidType = L
		for Sol in SolventsList:
			if(Sol in Sample):
				NSOLM.update( {Sol:int(Sample[Sol]/2)} )
				SolventType = Sol
		
		
		#=======================================================================
		# Setting the bilayers  ================================================
		#=======================================================================
		
		TMT = 0.0
		DZ = 0.0
		# DSPC bilayer =========================================================
		if LipidType == "DSPC":
			# Geometry to preprare the bilayer with packmol
			# total monolayer thickness (Angstrom)
			TMT = 30.0
			# deltaz : thickness in Z direction for
			# the volume constraining the heads and tails beads
			DZ = 7.0
		
		# DPPC bilayer =========================================================
		if LipidType == "DPPC":
			TMT = 30
			DZ = 10
		
		# DLPC bilayer =========================================================
		if LipidType == "DLPC":
			TMT = 30
			DZ = 10

		M1headMIN = LZ2 - TMT
		M1headMAX = LZ2 - TMT + DZ
		M1tailMIN = LZ2 - DZ
		M1tailMAX = LZ2

		M2headMIN = LZ2 + TMT - DZ
		M2headMAX = LZ2 + TMT
		M2tailMIN = LZ2
		M2tailMAX = LZ2 + DZ
		
		
		#=======================================================================
		#=======================================================================
		#=======================================================================
		
		#=======================================================================
		# Creating the defo ====================================================
		#=======================================================================
		
		if Sample['DEFO']['Height'] == 'box':
			L_defo = LZS
		if Sample['DEFO']['Height'] == 'bilayer':
			L_defo = LZ2+TMT+DZ/4.0
		if Sample['DEFO']['Height'] == 'follow':
			L_defo = 2*TMT+DZ/2.0
		NbLayers = int(L_defo/float(Sample['DEFO']['DzDefo']))
		
		#Defo per Layer
		defoPerLayer = int(Sample['DEFO']['DpL']) + 1
		#Total number of defo
		nbDefo = defoPerLayer*NbLayers
		#Radius for the hole
		radiusDefo = float(Sample['DEFO']['Radius'])
		
		#Number of Solvent inside the hole
		if(0): #Set to 1 to insert solvents in the Defo
			NbSolvIn = 0
			if SolventType == 'W':
				# Density x Volume
				NbSolvIn = int(8.26 * float(2*TMT*radiusDefo*radiusDefo*math.pi/1000.))
			if SolventType == 'OCO':
				# Density x Volume
				NbSolvIn = int(8.26 * float(2*TMT*radiusDefo*radiusDefo*math.pi/1000.)/2.)
			if SolventType == 'PW':
				# Density x Volume
				NbSolvIn = int(8.26 * float(2*TMT*radiusDefo*radiusDefo*math.pi/1000.)/3.)
		
		#Creating and Writing the defos configuration
		DefoXYZ_filename = 'defo.xyz'
		if os.path.exists(DefoXYZ_filename):
			os.remove(DefoXYZ_filename)
		XYZout = open(DefoXYZ_filename,"a")
		XYZout.write(str(nbDefo)+'\n\n')
		
		for i in range(0, NbLayers):
			ZcurrentDefo = float(Sample['DEFO']['DzDefo'])*i
			for i in np.arange(0., 360., 360./(defoPerLayer-1)):
				angle = math.radians(i)
				XYZout.write(Utility.RemoveUnwantedIndent(
					"""
					DEF  {0}  {1}  {2} \n
					""".format(radiusDefo*math.cos(angle), radiusDefo*math.sin(angle), ZcurrentDefo)
					))
			XYZout.write(Utility.RemoveUnwantedIndent(
					"""
					DEF  0.0 0.0 {0} \n
					""".format(ZcurrentDefo)
					))
		XYZout.close()
		
		#Formating the defos
		formatDEFO = open("format_DEFO.vmd","w")
		formatDEFO.write(Utility.RemoveUnwantedIndent(
			"""
			mol load xyz {0}
			set all [atomselect top "all"]
				
			$all set resname DEFO
			$all set name DEF
			$all set type DEF
			$all set chain X
				
			package require pbctools
			pbc set {{0.5 0.5 {1}}}
				
			$all writepdb {2}
			unset all
				
			exit
			""".format(DefoXYZ_filename, L_defo, DefoXYZ_filename.replace('xyz','pdb'))
			))
		formatDEFO.close()
		
		cmd = str("""{0} -dispdev text -e format_DEFO.vmd > format_DEFO.log""").format(Softwares['VMD'])
		sub.call(cmd, shell=True)
		
		#=======================================================================
		#=======================================================================
		#=======================================================================
		
		#=======================================================================
		# Setting the bilayers for packmol =====================================
		#=======================================================================
		
		#Setting the name of the system
		System = """{0}_{1}{2}_{3}{4}_{5}DEFO{6}""".format(Sample['TYPE'], Sample[LipidType],LipidType,
													Sample[SolventType], SolventType, nbDefo, 
													Sample['DEFO']['Version'])
		
		if packmolSeed:
			PackmolInput = """
						#Packmol seed was set using {0}
						seed {1}
					
					""".format(Sample['SEED'], packmolSeed)
		else:
			PackmolInput = """
						#Packmol Seed was set using default
					
							"""
		# DLPC bilayer =========================================================
		if LipidType == "DSPC":
			PackmolInput += """
						#
						# Lipid bilayer with water, perpendicular to z axis
						#
						
						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.
						
						tolerance 3.0
						
						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)
						
						filetype pdb
						
						# do not avoid the fixed molecules
						#avoid_overlap no
						
						# The output pdb file
						
						output {0}.pdb
						
						structure {1}
							chain A
							resnumbers 3
							number {2:g}
							inside box 0. 0. {3}  {4} {5} {6}
							atoms 1
								below plane 0. 0. 1. {7}
							end atoms
							atoms 9 14
								over plane 0. 0. 1. {8}
							end atoms
							outside cylinder  {13} {14}  {3}  0.  0.  1.  {15}  {16}
						end structure
						
						structure {1}
							chain B
							resnumbers 3
							number {2:g}
							inside box 0. 0.  {9}  {4} {5} {10}
							atoms 1
								over plane 0. 0. 1. {11}
							end atoms
							atoms 9 14
								below plane 0. 0. 1. {12}
							end atoms
							outside cylinder  {13} {14}  {3}  0.  0.  1.  {15}  {16}
						end structure
						
						""".format(System, PDBfileList[LipidType], NLM[LipidType], M1headMIN, 
									LXS, LYS, M1tailMAX, M1headMAX,
									M1tailMIN, M2tailMIN,M2headMAX, M2headMIN,
									M2tailMAX, LXS/2.0, LYS/2.0, radiusDefo,
									L_defo)
		
		# DPPC bilayer, DLPC bilayer =========================================================
		if LipidType in {'DPPC','DLPC'}:
			PackmolInput += """
						#
						# Lipid bilayer with water, perpendicular to z axis
						#
						
						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.
						
						tolerance 3.0
						
						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)
						
						filetype pdb
						
						# do not avoid the fixed molecules
						#avoid_overlap no
						
						# The output pdb file
						
						output {0}.pdb
						
						structure {1}
							chain A
							resnumbers 3
							number {2:g}
							inside box 0. 0. {3}  {4} {5} {6}
							constrain_rotation x 180 10
							constrain_rotation y 180 10
							outside cylinder  {9} {10}  0.  0.  0.  1.  {11}  {12}
						end structure
						
						structure {1}
							chain B
							resnumbers 3
							number {2:g}
							inside box 0. 0.  {7}  {4} {5} {8}
							constrain_rotation x 0 10
							constrain_rotation y 0 10
							outside cylinder  {9} {10}  0.  0.  0.  1.  {11}  {12}
						end structure
						
						""".format(System, PDBfileList[LipidType], NLM[LipidType], M1headMIN, LXS, LYS, M1tailMAX,
							M2tailMIN, M2headMAX, LXS/2.0, LYS/2.0, radiusDefo, L_defo)
		
		#=======================================================================
		#=======================================================================
		#=======================================================================
		
		#=======================================================================
		# Setting the defo with solvent for packmol ============================
		#=======================================================================
		
		if SolventType == 'PW':
			PackmolInput += """
			
							structure {0}
								chain C
								number {1:g}
								inside box 0. 0. 0.  {2} {3} {4}
								atoms 2 3
									radius 0.2
								end atoms
							end structure


							structure {0}
								chain D
								number {1:g}
								inside box 0. 0. {5}  {2} {3} {6}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType], NSOLM[SolventType], LXS, LYS,  M1headMIN,
								M2headMAX, LZS)
		else:
			PackmolInput += """
							
							structure {0}
								chain C
								number {1:g}
								inside box 0. 0. 0.  {2} {3} {4}
							end structure
							
							structure {0}
								chain D
								number {1:g}
								inside box 0. 0. {5}  {2} {3} {6}
							end structure
							""".format(PDBfileList[SolventType], NSOLM[SolventType], LXS, LYS,  M1headMIN,
								M2headMAX, LZS)
		if Sample['DEFO']['Height'] == 'follow':
			PackmolInput += """
							
							structure {0}
								chain X
								number 1
								center
								resnumbers 3
								fixed {1} {2} {3} 0.0 0.0 0.0
								radius 0.
							end structure
							""".format(DefoXYZ_filename.replace('xyz','pdb'), 
									LXS/2.0, LYS/2.0, LZ2)
		else:
			PackmolInput += """
							
							structure {0}
								chain X
								number 1
								resnumbers 3
								fixed {1} {2} {3} 0.0 0.0 0.0
								radius 0.
							end structure
							""".format(DefoXYZ_filename.replace('xyz','pdb'), 
									LXS/2.0, LYS/2.0, 0.0)
						
		f = open('packmol_'+System+'.input','w')
		f.write(Utility.RemoveUnwantedIndent(PackmolInput))
		f.close()
		
		# Lipids input pdb =====================================================
		if( LipidType == "DSPC"):
			f = open('dspc_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(DSPC_PDB))
			f.close()
		if( LipidType == "DPPC"):
			f = open('dppc_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(DPPC_PDB))
			f.close()
		if( LipidType == "DLPC"):
			f = open('dlpc_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(DLPC_PDB))
			f.close()
			
		# Solvents input pdb ===================================================
		if( SolventType == 'W'):
			f = open('water_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(W_PDB))
			f.close()
		if( SolventType == 'PW'):
			f = open('polwater_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(PW_PDB))
			f.close()
		if( SolventType == 'OCO'):
			f = open('octanol_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(OCO_PDB))
			f.close()
		
		#=======================================================================
		#lauching packmol
		#=======================================================================

		cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """).format(Softwares['PACKMOL'],System)
		sub.call(cmd, shell=True)

		#=======================================================================
		# ensure the right box in pdb file
		#=======================================================================

		WriteBox = str("""
				mol load pdb {0}.pdb
				set all [atomselect top "all"]

				package require pbctools
				pbc set {{{1} {2} {3}}}

				$all writepdb "{0}.withbox.pdb"
				unset all

				exit
				""").format(System, Sample['LX'], Sample['LY'], Sample['LZ'])

		f = open('write_box.vmd','w+')
		f.write(Utility.RemoveUnwantedIndent(WriteBox) )
		f.close()

		## ======================================================================
		cmd = str("""{0} -dispdev text -e write_box.vmd > write_box.log"""
												).format(Softwares['VMD'])
		sub.call(cmd, shell=True)
		
		MakeIndex = str("""
				chain A
				chain B
				chain C
				chain D
				chain X
				name 5 bottom{0}
				name 6 top{0}
				name 7 bottom{1}
				name 8 top{1}
				name 9 defo
				q
				
				""").format(LipidType,SolventType)

		f = open('make_index.input','w+')
		f.write(Utility.RemoveUnwantedIndent(MakeIndex) )
		f.close()
		
		cmd = str("""{0}make_ndx -f {1}.withbox.pdb -o {1}.ndx < make_index.input""").format(GROMACS_LOC_prefixPath,
																					   System)
		sub.call(cmd, shell=True)
		
		#=======================================================================
		# the topology file
		#=======================================================================
		TopologyDefo(Sample, PathToDefault)
		
		Topology = str("""
					#include "martini_v2.2_DEFO_{0}.itp"
					#include "martini_v2.0_lipids.itp"
					""".format(Sample['DEFO']['Version']))
		
		#Copy the topology files for martini forcefield
		sub.call("""cp {0}/martini_v2.0_lipids.itp  .""".format(PathToDefault), shell= True)
		sub.call("""cp {0}/DEFO/defo_posres.itp  .""".format(PathToDefault), shell= True)
		f = open(System+'.top','w')
		f.write(Utility.RemoveUnwantedIndent(Topology))

		f = open(System+'.top','a')
		f.write("""\n[ system ]\n""")

		Topology = str("""{0}_{1}_WITH_DEFO_{2}\n""").format(LipidType, Sample['TYPE'], Sample['DEFO']['Version'])
		f.write(Topology)

		f.write("""\n[ molecules ]\n""")
		Topology = str("""{0} {1}\n{2} {3}\n{4} {5}\n""").format(LipidType, Sample[LipidType], SolventType, 
														   Sample[SolventType],'DEFO',nbDefo)
		f.write(Topology)
		f.close()
		
		Output = str("""{0}.withbox.pdb""").format(System)
		Index = str("""{0}.ndx""").format(System)
		
		
		ApL = (Sample['LX']*Sample['LY'] - 3.141516*radiusDefo*radiusDefo)/NLM[LipidType]
		print(Utility.RemoveUnwantedIndent("""
				================================
				================================
				Packmol finished initial input file
				{0} {1}, {2} {3} with DEFO{4} {5}
				box sizes : {6}, {7}, {8}
				Area per lipid = {9} A**2
				===============================
				===============================
				""".format(LipidType, Sample[LipidType],SolventType, Sample[SolventType],
					Sample['DEFO']['Version'],nbDefo,
					Sample['LX'],Sample['LY'],Sample['LZ'],ApL)))
	
		return { 'SYSTEM': System, 'OUTPUT': Output, 'INDEX':Index}
	
	

def InitBilayerWithHoleAndWall(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault, packmolSeed):
	THICKNESS = float(Sample['SU']['Thickness'])
	DENSITY = float(Sample['SU']['Density'])
	NBLIPIDS_MONO = int(Sample['SU']['NbLipidsM'])
	SU_VERSION = Sample['SU']['Version']
	SU_TYPE = Sample['SU']['SuType']
	LZM = float(Sample['SU']['LzM'])
	
	
	SHELL = 3.0
	LXS = Sample['LX'] - SHELL
	LYS = Sample['LY'] - SHELL
	
	#Number of lipid and water per layer
	NLM = {}
	NSOLM = {}
	LipidType = ''
	
	for L in LipidsList:
		if(L in Sample):
			NLM.update( {L:int(Sample[L]/2)} )
			LipidType = L
			
	for Sol in SolventsList:
		if Sol in Sample:
			NSOLM.update( {Sol:int(Sample[Sol]/2)} )
			SolventType = Sol
	
	#=======================================================================
	# Setting the bilayers  ================================================
	#=======================================================================
	
	TMT = 0.0
	DZ = 0.0
	# DSPC bilayer =========================================================
	if LipidType == "DSPC":
		# Geometry to preprare the bilayer with packmol
		# total monolayer thickness (Angstrom)
		TMT = 30.0
		# deltaz : thickness in Z direction for
		# the volume constraining the heads and tails beads
		DZ = 7.0
	
	# DPPC bilayer =========================================================
	if LipidType == "DPPC":
		TMT = 30.0
		DZ = 10
	
	# DLPC bilayer =========================================================
	if LipidType == "DLPC":
		TMT = 30.0
		DZ = 10
	
	assert(TMT > 0.0 and DZ > 0.0),"Your Parameters.csv contains a lipid not yet defined in Prepare.py"
	#Set the box height with the monolayer + shell for PBC
	LBZ = (LZM-THICKNESS)/2.0 + THICKNESS
	LZ2 = 0.0
	NbParticlesToRelocate = 0
	NbSolTop = NSOLM[SolventType]
	NbSolBottom = NSOLM[SolventType]
	if 'BilayerHeight' in Sample['SU']:
		LZ2 = float(Sample['SU']['BilayerHeight'])
		#Checking that the Bilayer lipids do not overlap with the Monolayer lipids and the substrate
		if (LZ2 - TMT) < THICKNESS:
			LZ2 -= LZ2 -TMT - THICKNESS
		if (LZ2 + TMT) > LZM:
			LZ2 += LZM - LZ2 - TMT
		#Checking that the Solvent is not too dense
		if LZ2 < LBZ: #Too much Solvent below
			DensitySol = NSOLM[SolventType]/LXS/LYS/(LBZ - TMT - THICKNESS)
			VolToRemove = LBZ - LZ2
			VolToRemove *= LXS*LYS
			NbParticlesToRelocate = int(DensitySol * VolToRemove)
			
			NbSolBottom -= NbParticlesToRelocate
			NbSolTop += NbParticlesToRelocate
			
		if LZ2 > LBZ: #Too much Solvent above
			DensitySol = NSOLM[SolventType]/LXS/LYS/(LZM - LBZ - TMT)
			VolToRemove = LZ2 - LBZ
			VolToRemove *= LXS*LYS
			NbParticlesToRelocate = int(DensitySol * VolToRemove)
			
			NbSolBottom += NbParticlesToRelocate
			NbSolTop -= NbParticlesToRelocate
	else:
		LZ2 = LBZ
			
	
	M1headMIN = LZ2 - TMT
	M1headMAX = LZ2 - TMT + DZ
	M1tailMIN = LZ2 - DZ
	M1tailMAX = LZ2

	M2headMIN = LZ2 + TMT - DZ
	M2headMAX = LZ2 + TMT
	M2tailMIN = LZ2
	M2tailMAX = LZ2 + DZ
	
	#=======================================================================
	# Creating the defo ====================================================
	#=======================================================================
	
	if Sample['DEFO']['Height'] == 'box':
		#Set the Defo height from the substrate to the mono layer
		L_defo = Sample['LZ'] - THICKNESS
	if Sample['DEFO']['Height'] == 'bilayer':
		#Set the Defo height from the substrate to the top of the bilayer
		L_defo = LZ2 + TMT - THICKNESS + DZ/4.0
	if Sample['DEFO']['Height'] == 'mono':
		#Set the Defo height from the substrate to the top of the monolayer
		L_defo = LZM + TMT - THICKNESS + DZ/4.0
	if Sample['DEFO']['Height'] == 'follow':
		#Set the Defo height from the bottom of the bilayer to its top
		L_defo = TMT+DZ/4.0
	
	NbLayers = int(L_defo/float(Sample['DEFO']['DzDefo']))
	
	#Defo per Layer
	defoPerLayer = int(Sample['DEFO']['DpL']) + 1
	#Total number of defo
	nbDefo = defoPerLayer*NbLayers
	#Radius for the hole
	radiusDefo = float(Sample['DEFO']['Radius'])
	#Number of Solvent inside the hole
	if(0): #Set to 1 to insert solvents in the Defo
		NbSolvIn = 0
		if SolventType == 'W':
			# Density x Volume
			NbSolvIn = int(8.26 * float(2*TMT*radiusDefo*radiusDefo*math.pi/1000.))
		if SolventType == 'OCO':
			# Density x Volume
			NbSolvIn = int(8.26 * float(2*TMT*radiusDefo*radiusDefo*math.pi/1000.)/2.)
		if SolventType == 'PW':
			# Density x Volume
			NbSolvIn = int(8.26 * float(2*TMT*radiusDefo*radiusDefo*math.pi/1000.)/3.)
	
	#Creating and Writing the defos configuration
	DefoXYZ_filename = 'defo.xyz'
	if os.path.exists(DefoXYZ_filename):
		os.remove(DefoXYZ_filename)
	XYZout = open(DefoXYZ_filename,"a")
	XYZout.write(str(nbDefo)+'\n\n')
	
	for i in range(0, NbLayers):
		ZcurrentDefo = float(Sample['DEFO']['DzDefo'])*i
		for i in np.arange(0., 360., 360./(defoPerLayer-1)):
			angle = math.radians(i)
			XYZout.write(Utility.RemoveUnwantedIndent(
				"""
				DEF  {0}  {1}  {2} \n
				""".format(radiusDefo*math.cos(angle), radiusDefo*math.sin(angle), ZcurrentDefo)
				))
		XYZout.write(Utility.RemoveUnwantedIndent(
				"""
				DEF  0.0 0.0 {0} \n
				""".format(ZcurrentDefo)
				))
	XYZout.close()
	
	#Formating the defos
	formatDEFO = open("format_DEFO.vmd","w")
	formatDEFO.write(Utility.RemoveUnwantedIndent(
		"""
		mol load xyz {0}
		set all [atomselect top "all"]
			
		$all set resname DEFO
		$all set name DEF
		$all set type DEF
		$all set chain X
			
		package require pbctools
		pbc set {{0.5 0.5 {1}}}
			
		$all writepdb {2}
		unset all
			
		exit
		""".format(DefoXYZ_filename, L_defo, DefoXYZ_filename.replace('xyz','pdb'))
		))
	formatDEFO.close()
	
	cmd = str("""{0} -dispdev text -e format_DEFO.vmd > format_DEFO.log""").format(Softwares['VMD'])
	sub.call(cmd, shell=True)
	
	#=======================================================================
	# Number of su particles ===============================================
	#=======================================================================
	#Computing the number of SU Density is in CG/nm³ thus the divided by 1000 for Volume
	nbSu = int( DENSITY * Sample['LX'] * Sample['LY'] * THICKNESS /1000 )
	
	#=======================================================================
	#=======================================================================
	#=======================================================================
	
	#=======================================================================
	#=======================================================================
	#=======================================================================

	System = str('{0}_{1}{2}_{3}M{2}_{4}{5}_{6}DEFO{7}_{8}{10}{9}').format(Sample['TYPE'], Sample[LipidType],
											LipidType, NBLIPIDS_MONO , 
											Sample[SolventType], SolventType,
											nbDefo, Sample['DEFO']['Version'], nbSu, Sample['SU']['Version'],SU_TYPE)
	
	#=======================================================================
	# Setting the bilayers for packmol =====================================
	#=======================================================================
	if packmolSeed:
		PackmolInput = """
					#Packmol seed was set using {0}
					seed {1}
					
					""".format(Sample['SEED'], packmolSeed)
	else:
		PackmolInput = """
					#Packmol seed was set using default
					
					"""
	# DLPC bilayer =========================================================
	if LipidType == "DSPC":
		PackmolInput += """
					#
					# Lipid bilayer with water, perpendicular to z axis
					#
					
					# Every atom from diferent molecules will be far from each other at
					# least 3.0 Anstroms at the solution.
					
					tolerance 3.0
					
					# Coordinate file types will be in pdb format (keyword not required for
					# pdb file format, but required for tinker, xyz or moldy)
					
					filetype pdb
					
					# do not avoid the fixed molecules
					#avoid_overlap no
					
					# The output pdb file
					
					output {0}.pdb
					
					structure {1}
						chain A
						resnumbers 3
						number {2:g}
						inside box 0. 0. {3}  {4} {5} {6}
						atoms 1
							below plane 0. 0. 1. {7}
						end atoms
						atoms 9 14
							over plane 0. 0. 1. {8}
						end atoms
						outside cylinder  {13} {14}  {3}  0.  0.  1.  {15}  {16}
					end structure
					
					structure {1}
						chain B
						resnumbers 3
						number {2:g}
						inside box 0. 0.  {9}  {4} {5} {10}
						atoms 1
							over plane 0. 0. 1. {11}
						end atoms
						atoms 9 14
							below plane 0. 0. 1. {12}
						end atoms
						outside cylinder  {13} {14}  {3}  0.  0.  1.  {15}  {16}
					end structure
					
					""".format(System, PDBfileList[LipidType], NLM[LipidType], M1headMIN, LXS,
								LYS, M1tailMAX, M1headMAX, M1tailMIN, 
								M2tailMIN, M2headMAX, M2headMIN, M2tailMAX, 
								LXS/2.0, LYS/2.0, radiusDefo, L_defo)
		if Sample['DEFO']['Height'] in {'follow','mono','box'}:
			PackmolInput += """
						
						structure {0}
							chain C
							resnumbers 3
							number {1:g}
							inside box 0. 0.  {2}  {3} {4} {5}
							atoms 1
								below plane 0. 0. 1. {6}
							end atoms
							atoms 9 14
								over plane 0. 0. 1. {7}
							end atoms
							outside cylinder  {8} {9}  {2}  0.  0.  1.  {10}  {11}
						end structure
						
						""".format(PDBfileList[LipidType], NBLIPIDS_MONO, LZM, LXS, LYS, TMT+LZM, LZM+DZ, TMT+LZM-DZ,
									LXS/2.0, LYS/2.0, radiusDefo, L_defo)
		else:
			PackmolInput += """
						
						structure {0}
							chain C
							resnumbers 3
							number {1:g}
							inside box 0. 0.  {2}  {3} {4} {5}
							atoms 1
								below plane 0. 0. 1. {6}
							end atoms
							atoms 9 14
								over plane 0. 0. 1. {7}
							end atoms
						end structure
						
						""".format(PDBfileList[LipidType], NBLIPIDS_MONO, LZM, LXS, LYS, TMT+LZM, LZM+DZ, TMT+LZM-DZ)
	
	# DPPC bilayer, DLPC bilayer =========================================================
	if LipidType in {'DPPC','DLPC'}:
		PackmolInput += """
						#
						# Lipid bilayer with water, perpendicular to z axis
						#

						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.

						tolerance 3.0

						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)

						filetype pdb

						# The output pdb file

						output {0}.pdb

						structure {1}
							chain A
							resnumbers 3
							number {2:g}
							inside box 0. 0. {3}  {4} {5} {6}
							constrain_rotation x 180 10
							constrain_rotation y 180 10
							outside cylinder  {9} {10}  {3}  0.  0.  1.  {11}  {12}
						end structure

						structure {1}
							chain B
							resnumbers 3
							number {2:g}
							inside box 0. 0.  {7}  {4} {5} {8}
							constrain_rotation x 0 10
							constrain_rotation y 0 10
							outside cylinder  {9} {10}  {3}  0.  0.  1.  {11}  {12}
						end structure

						""".format(System, PDBfileList[LipidType], NLM[LipidType], M1headMIN, LXS, 
									LYS, M1tailMAX, M2tailMIN, M2headMAX,
									LXS/2.0, LYS/2.0, radiusDefo, L_defo
									)
		if Sample['DEFO']['Height'] in {'follow','mono','box'}:
			PackmolInput += """
							
							structure {0}
								chain C
								resnumbers 3
								number {1:g}
								inside box 0. 0.  {2}  {3} {4} {5}
								constrain_rotation x 180 10
								constrain_rotation y 180 10
								outside cylinder  {6} {7}  {2}  0.  0.  1.  {8}  {9}
							end structure
							
							""".format(PDBfileList[LipidType], NBLIPIDS_MONO, LZM, LXS, LYS, TMT+LZM, 
										LXS/2.0, LYS/2.0, radiusDefo, L_defo/2.0)
		else:
			PackmolInput += """
							
							structure {0}
								chain C
								resnumbers 3
								number {1:g}
								inside box 0. 0.  {2}  {3} {4} {5}
								constrain_rotation x 180 10
								constrain_rotation y 180 10
							end structure
							
							""".format(PDBfileList[LipidType], NBLIPIDS_MONO, LZM, LXS, LYS, TMT+LZM)
							
	#=======================================================================
	#=======================================================================
	#=======================================================================
	
	
					
	#=======================================================================
	# Setting the defo with solvent for packmol ============================
	#=======================================================================
	#Computing the number of Solvent vapor
	NbSolVapor = int( 0.007 * Sample['LX'] * Sample['LY'] * (Sample['LZ']-LZM-TMT) /1000)
	
	if SolventType == 'PW':
		if NbSolBottom:
			PackmolInput += """
							
							structure {0}
								chain D
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType],NbSolBottom,THICKNESS,LXS,LYS,M1headMIN)
			
		if NbSolTop:
			PackmolInput += """
							
							structure {0}
								chain E
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType],NbSolTop,M2headMAX,LXS,LYS,LZM)
							
		PackmolInput += """
							
							structure {0}
								chain V
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""".format(PDBfileList[SolventType],NbSolVapor,LZM+TMT,LXS,LYS,Sample['LZ'])
	else:
		if NbSolBottom:
			PackmolInput += """
							
							structure {0}
								chain D
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
							end structure
							""".format(PDBfileList[SolventType],NbSolBottom,THICKNESS,LXS,LYS,M1headMIN)
		if NbSolTop:
			PackmolInput += """
							
							structure {0}
								chain E
								number {1:g}
								inside box 0. 0. {2}  {3} {4} {5}
							end structure
							""".format(PDBfileList[SolventType],NbSolTop,M2headMAX,LXS,LYS,LZM)
		PackmolInput += """
						
						structure {0}
							chain V
							number {1:g}
							inside box 0. 0. {2}  {3} {4} {5}
						end structure
						""".format(PDBfileList[SolventType],NbSolVapor,LZM+TMT,LXS,LYS,Sample['LZ'])
					
	if Sample['DEFO']['Height'] == 'follow':
		PackmolInput += """
					structure {0}
						chain X
						number 1
						resnumbers 3
						fixed {1} {2} {3} 0.0 0.0 0.0
						radius 0.
					end structure
					
					structure {0}
						chain X
						number 1
						resnumbers 3
						fixed {1} {2} {4} 0.0 0.0 0.0
						radius 0.
					end structure
					""".format(DefoXYZ_filename.replace('xyz','pdb'),LXS/2.0, LYS/2.0, LZ2-TMT-DZ/4.0,LZ2)
		PackmolInput += """
					structure {0}
						chain Y
						number 1
						resnumbers 3
						fixed {1} {2} {3} 0.0 0.0 0.0
						radius 0.
					end structure
					""".format(DefoXYZ_filename.replace('xyz','pdb'),LXS/2.0, LYS/2.0, LZM-DZ/4.0)
			
						
	else:
		PackmolInput += """
					structure {0}
						chain X
						number 1
						resnumbers 3
						fixed {1} {2} {3} 0.0 0.0 0.0
						radius 0.
					end structure
					""".format(DefoXYZ_filename.replace('xyz','pdb'),LXS/2.0, LYS/2.0, THICKNESS)
					
	
	
	#=======================================================================
	# Adding the substrat ==================================================
	#=======================================================================
	PackmolInput += """
					
					structure {0}
						chain S
						number {1:g}
						inside box 0. 0. 0.  {2} {3} {4}
					end structure
					""".format(PDBfileList['SU'], nbSu, LXS, LYS, THICKNESS)
	
	
	f = open('packmol_'+System+'.input','w')
	f.write(Utility.RemoveUnwantedIndent(PackmolInput))
	f.close()
	
	# Lipids input pdb =====================================================
	if( LipidType == "DSPC"):
		f = open(PDBfileList['DSPC'],'w')
		f.write(Utility.RemoveUnwantedIndent(DSPC_PDB))
		f.close()
	if( LipidType == "DPPC"):
		f = open(PDBfileList['DPPC'],'w')
		f.write(Utility.RemoveUnwantedIndent(DPPC_PDB))
		f.close()
	if( LipidType == "DLPC"):
		f = open(PDBfileList['DLPC'],'w')
		f.write(Utility.RemoveUnwantedIndent(DLPC_PDB))
		f.close()
		
	# Solvents input pdb ===================================================
	if( SolventType == 'W'):
		f = open(PDBfileList['W'],'w')
		f.write(Utility.RemoveUnwantedIndent(W_PDB))
		f.close()
	if( SolventType == 'PW'):
		f = open(PDBfileList['PW'],'w')
		f.write(Utility.RemoveUnwantedIndent(PW_PDB))
		f.close()
	if( SolventType == 'OCO'):
		f = open(PDBfileList['OCO'],'w')
		f.write(Utility.RemoveUnwantedIndent(OCO_PDB))
		f.close()
		
	# pdb file for su
	suType = SU_TYPE
	assert(len(SU_TYPE) <= 4),"The name of the SU should be 4 letters max (Otherwise problem with PDB file format)" 
	suType = SU_TYPE + (4 - len(SU_TYPE))*' '
	suPdb = SU_PDB.replace("TEMP", suType)
	f = open(PDBfileList['SU'],'w')
	f.write(Utility.RemoveUnwantedIndent(suPdb))
	f.close()
		
		
	#=======================================================================
	#lauching packmol
	#=======================================================================

	cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """
							).format(Softwares['PACKMOL'], System)
	sub.call(cmd, shell=True)

	#=======================================================================
	# ensure the right box in pdb file
	#=======================================================================

	WriteBox = str("""
			mol load pdb {0}.pdb
			set all [atomselect top "all"]

			package require pbctools
			pbc set {{{1} {2} {3}}}

			$all writepdb "{0}.withbox.pdb"
			unset all

			exit
			""").format(System, Sample['LX'], Sample['LY'], Sample['LZ'])

	f = open('write_box.vmd','w+')
	f.write(Utility.RemoveUnwantedIndent(WriteBox) )
	f.close()

	## ======================================================================
	cmd = str("""{0} -dispdev text -e write_box.vmd > write_box.log"""
											).format(Softwares['VMD'])
	sub.call(cmd, shell=True)
	if Sample['DEFO']['Height'] == 'follow':
		MakeIndex = """
			chain A
			chain B
			chain C
			chain D
			chain E
			chain V
			chain S
			chain X
			chain Y
			name 6 bottom{0}
			name 7 top{0}
			name 8 mono{0}
			name 9 bottom{1}
			name 10 top{1}
			name 11 vap{1}
			name 12 su
			name 13 defoBi
			name 14 defoMono
			6 | 7
			name 15 bi{0}
			q
			
			""".format(LipidType, SolventType)
	else:
		MakeIndex = """
			chain A
			chain B
			chain C
			chain D
			chain E
			chain V
			chain S
			chain X
			name 6 bottom{0}
			name 7 top{0}
			name 8 mono{0}
			name 9 bottom{1}
			name 10 top{1}
			name 11 vap{1}
			name 12 su
			name 13 defo
			6 | 7
			name 14 bi{0}
			q
			
			""".format(LipidType, SolventType)

	f = open('make_index.input','w+')
	f.write(Utility.RemoveUnwantedIndent(MakeIndex) )
	f.close()
	
	cmd = str("""{0}make_ndx -f {1}.withbox.pdb -o {1}.ndx < make_index.input"""
										).format(GROMACS_LOC_prefixPath, System)
	sub.call(cmd, shell=True)
	
	
	#=======================================================================
	# the topology file
	#=======================================================================
	TopologyDefoSu(Sample, PathToDefault)
	
	Topology = """
					#include "martini_v2.2_{2}_{0}_DEFO_{1}.itp"
					#include "martini_v2.0_lipids.itp"
					
				""".format( Sample['SU']['Version'], Sample['DEFO']['Version'], SU_TYPE )
	
	#Copy the topology files for martini forcefield
	sub.call("""cp {0}/martini_v2.0_lipids.itp  .""".format(PathToDefault), shell= True)
	sub.call("""cp {0}/DEFO/defo_posres.itp  .""".format(PathToDefault), shell= True)
	sub.call("""cp {0}/SU/su_posres.itp  .""".format(PathToDefault), shell= True)
	f = open(System+'.top','w')
	f.write(Utility.RemoveUnwantedIndent(Topology))

	f = open(System+'.top','a')
	f.write("""\n[ system ]\n""")

	Topology = str("""{0}_{1}_WITH_DEFO_{2}_AND_{4}_{3}\n""").format(LipidType, Sample['TYPE'], Sample['DEFO']['Version'], Sample['SU']['Version'], SU_TYPE)
	f.write(Topology)

	f.write("""\n[ molecules ]\n""")
	Topology = str("""{0} {1}\n{2} {3}\n{4} {5}\n{6} {7}\n""").format(LipidType, Sample[LipidType]+NBLIPIDS_MONO,SolventType, Sample[SolventType],'DEFO', nbDefo, SU_TYPE, nbSu)
	f.write(Topology)
	f.close()
		
	#==================================================================================
	# Output the files for other steps
	#==================================================================================
	Output = str("""{0}.withbox.pdb""").format(System)
	Index = str("""{0}.ndx""").format(System)
	
	ApL = (Sample['LX']*Sample['LY'] - 3.141516*radiusDefo*radiusDefo)/NLM[LipidType]
	print(Utility.RemoveUnwantedIndent("""
			================================
			================================
			Packmol finished initial input file
			{0} {1}, {2} {3} with DEFO{4} {5}, SU{6} {7} and mono{0} {8}
			box sizes : {9}, {10}, {11}
			Area per lipid = {12} A**2
			===============================
			===============================
			""".format(LipidType, Sample[LipidType],SolventType, Sample[SolventType],
				Sample['DEFO']['Version'],nbDefo,
				Sample['SU']['Version'],nbSu,
				NBLIPIDS_MONO,
				Sample['LX'],Sample['LY'],Sample['LZ'],ApL)))
	return { 'SYSTEM': System, 'OUTPUT': Output, 'INDEX':Index}

























def InitSolvent(Sample, Softwares, GROMACS_LOC_prefixPath, PathToDefault):
	# packmol constraints are not strict !
	# put the molecules in a smaller box to
	# avoid initial crash because of PBC
	# the SHELL parameter defines the shell in which
	# packmol is not supposed to put particles
	
	
	SHELL = 3.0
	
	LXS = 0.0
	LYS = 0.0
	LZS = 0.0
	
	if 'FILLMODE' in Sample:
		if Sample['FILLMODE'] == 'FULL':
			LXS = Sample['LX']-SHELL
			LYS = Sample['LY']-SHELL
			LZS = Sample['LZ']-SHELL
		if Sample['FILLMODE'] == 'BOX':
			LXS = LYS = LZS = Sample['LX']-SHELL
		if Sample['FILLMODE'] == 'HALF-X':
			LXS = Sample['LX']/2.0-SHELL
			LYS = Sample['LY']-SHELL
			LZS = Sample['LZ']-SHELL
		if Sample['FILLMODE'] == 'HALF-Y':
			LXS = Sample['LX']-SHELL
			LYS = Sample['LY']/2.0-SHELL
			LZS = Sample['LZ']-SHELL
		if Sample['FILLMODE'] == 'HALF-Z':
			LXS = Sample['LX']-SHELL
			LYS = Sample['LY']-SHELL
			LZS = Sample['LZ']/2.0-SHELL
	else:
		LXS = Sample['LX']-SHELL
		LYS = Sample['LY']-SHELL
		LZS = Sample['LZ']-SHELL
	
	Solvents = []
	NbSol = 0
	for Sol in SolventsList:
		if(Sol in Sample):
			Solvents.append(Sol)
			NbSol += 1
			
	#==================================================================================
	# Only one solvent case
	#==================================================================================
	if NbSol == 1:
		SolventType = Solvents[0]
		#==================================================================================
		# creating the initial bilayer (lipids + water) using packmol
		#==================================================================================
		
		System = str('{0}_{1}{2}').format(Sample['TYPE'], Sample[SolventType], SolventType)
		if SolventType == 'PW':
			# creating input for packmol
			PackmolInput = str("""
							#
							#  solvent
							#
							
							# Every atom from diferent molecules will be far from each other at
							# least 3.0 Anstroms at the solution.
							
							tolerance 3.0
							
							# Coordinate file types will be in pdb format (keyword not required for
							# pdb file format, but required for tinker, xyz or moldy)
							
							filetype pdb
							
							# The output pdb file
							
							output {0}.pdb
							
							structure {5}
								number {1:g}
								inside box 0. 0. 0.  {2} {3} {4}
								atoms 2 3
									radius 0.2
								end atoms
							end structure
							""").format(System, Sample[SolventType], LXS, LYS, LZS, PDBfileList[SolventType])
			
			f = open('packmol_'+System+'.input','w')
			f.write(Utility.RemoveUnwantedIndent(PackmolInput))
			f.close()
		else:
			# creating input for packmol
			PackmolInput = str("""
							#
							#  solvent
							#
							
							# Every atom from diferent molecules will be far from each other at
							# least 3.0 Anstroms at the solution.
							
							tolerance 3.0
							
							# Coordinate file types will be in pdb format (keyword not required for
							# pdb file format, but required for tinker, xyz or moldy)
							
							filetype pdb
							
							# The output pdb file
							
							output {0}.pdb
							
							structure {5}
								number {1:g}
								inside box 0. 0. 0.  {2} {3} {4}
							end structure
							""").format(System, Sample[SolventType], LXS, LYS, LZS, PDBfileList[SolventType])
			
			f = open('packmol_'+System+'.input','w')
			f.write(Utility.RemoveUnwantedIndent(PackmolInput))
			f.close()
		
		# water input pdb
		if( SolventType == 'W'):
			f = open('water_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(W_PDB))
			f.close()
		
		if( SolventType == 'OCO'):
			f = open('octanol_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(OCO_PDB))
			f.close()
		
		if( SolventType == 'PW'):
			f = open('polwater_single.pdb','w')
			f.write(Utility.RemoveUnwantedIndent(PW_PDB))
			f.close()

		
		#==================================================================================
		#lauching packmol
		#==================================================================================
		
		cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """).format(Softwares['PACKMOL'],System)
		sub.call(cmd, shell=True)
		
		#==================================================================================
		# ensure the right box in pdb file
		#==================================================================================
		
		WriteBox = str("""
				mol load pdb {0}.pdb
				set all [atomselect top "all"]
		
				package require pbctools
				pbc set {{{1} {2} {3}}}
				
				$all writepdb "{0}.withbox.pdb"
				unset all
				
				exit
				""").format(System, Sample['LX'], Sample['LY'], Sample['LZ'])
		
		f = open('write_box.vmd','w+')
		f.write(Utility.RemoveUnwantedIndent(WriteBox) )
		f.close()
		
		## ======================================================================
		cmd = str("""{0} -dispdev text -e write_box.vmd > write_box.log""").format(Softwares['VMD'])
		sub.call(cmd, shell=True)
		
		## ======================================================================
		print(Utility.RemoveUnwantedIndent(str("""
				================================
				================================
				Packmol finished initial input file
				{0} {1}
				box sizes : {2}, {3}, {4}
				===============================
				===============================
				""").format(SolventType, Sample[SolventType], Sample['LX'],Sample['LY'],Sample['LZ']) ))
		
		#==================================================================================
		# the topology file
		#==================================================================================
		
		Topology = str("""
					#include "martini_v2.2.itp" ; modified with polarisable water
					#include "martini_v2.0_solvents.itp"
					""")
		#Copy the topology files for martini forcefield
		sub.call("""cp {0}/martini_v2.0_solvents.itp {0}/martini_v2.2.itp ./""".format(PathToDefault), shell= True)
		f = open(System+'.top','w')
		f.write(Utility.RemoveUnwantedIndent(Topology))
		
		f = open(System+'.top','a')
		f.write(Utility.RemoveUnwantedIndent("""
						
						[ system ]
						{0}
						
						""".format(System)))
		
		###add a for loop for multiple types (Later)
		f.write("""\n[ molecules ]\n""")
		Topology = str("""{0} {1}""").format(SolventType, Sample[SolventType])
		f.write(Topology)
		f.close()
		
		#==================================================================================
		# the index file
		#==================================================================================
		cmd = str("""echo q | {0}make_ndx -f {1}.withbox.pdb -o {1}.ndx""").format(GROMACS_LOC_prefixPath,System)
		sub.call(cmd, shell=True)
		
		#==================================================================================
		# Output the files for other steps
		#==================================================================================
		Output = str("""{0}.withbox.pdb""").format(System)
		Index = str("""{0}.ndx""").format(System)
		return { 'SYSTEM': System, 'OUTPUT': Output, 'INDEX':Index}
	
	#==================================================================================
	# More than one solvent case
	#==================================================================================
	if NbSol > 1:
		#==================================================================================
		# creating the initial bilayer (lipids + water) using packmol
		#==================================================================================
		Qt = ''
		for Sol in Solvents:
			Qt += "{0}{1}_".format(Sample[Sol], Sol)
		#Remove last underscore
		li = Qt.rsplit('_',1)
		Qt = ('').join(li)
		#**********************
		System = str('{0}_{1}').format(Sample['TYPE'], Qt)
		
		# creating input for packmol
		f = open('packmol_'+System+'.input','a+')
		PackmolInput = str("""
						#
						#  solvent
						#
						
						# Every atom from diferent molecules will be far from each other at
						# least 3.0 Anstroms at the solution.
						
						tolerance 3.0
						
						# Coordinate file types will be in pdb format (keyword not required for
						# pdb file format, but required for tinker, xyz or moldy)
						
						filetype pdb
						
						# The output pdb file
						
						output {0}.pdb
						
						
						""").format(System)
		f.write(Utility.RemoveUnwantedIndent(PackmolInput))
		Shift = 0.
		for Sol in Solvents:
			print(Sol)
			if Sol == 'PW':
				PackmolInput = str("""
						
						structure {5}
							resnumbers 3
							number {0:g}
							inside box 0. 0. {4} {1} {2} {3}
							atoms 2 3
								radius 0.2
							end atoms
						end structure
						
						""").format(Sample[Sol], LXS, LYS, (LZS/NbSol)+Shift, Shift, PDBfileList[Sol])
				f.write(Utility.RemoveUnwantedIndent(PackmolInput))
			else:
				PackmolInput = str("""
						
						structure {5}
							resnumbers 3
							number {0:g}
							inside box 0. 0. {4} {1} {2} {3}
						end structure
						
						""").format(Sample[Sol], LXS, LYS, (LZS/NbSol)+Shift, Shift, PDBfileList[Sol])
				f.write(Utility.RemoveUnwantedIndent(PackmolInput))
				Shift += LZS/NbSol
				
		f.close()
		
		for Sol in Solvents:
			if( Sol == 'W'):
				f = open('water_single.pdb','w')
				f.write(Utility.RemoveUnwantedIndent(W_PDB))
				f.close()
		
			if( Sol == 'OCO'):
				f = open('octanol_single.pdb','w')
				f.write(Utility.RemoveUnwantedIndent(OCO_PDB))
				f.close()
			
			if( Sol == 'PW'):
				f = open('polwater_single.pdb','w')
				f.write(Utility.RemoveUnwantedIndent(PW_PDB))
				f.close()
		
		
		#==================================================================================
		#lauching packmol
		#==================================================================================
		
		cmd = str("""{0} < packmol_{1}.input > packmol_{1}.output """).format(Softwares['PACKMOL'],System)
		sub.call(cmd, shell=True)
		
		#==================================================================================
		# ensure the right box in pdb file
		#==================================================================================
		
		WriteBox = str("""
				mol load pdb {0}.pdb
				set all [atomselect top "all"]
		
				package require pbctools
				pbc set {{{1} {2} {3}}}
				
				$all writepdb "{0}.withbox.pdb"
				unset all
				
				exit
				""").format(System, Sample['LX'], Sample['LY'], Sample['LZ'])
		
		f = open('write_box.vmd','w+')
		f.write(Utility.RemoveUnwantedIndent(WriteBox) )
		f.close()
		
		## ======================================================================
		cmd = str("""{0} -dispdev text -e write_box.vmd > write_box.log""").format(Softwares['VMD'])
		sub.call(cmd, shell=True)
		
		## ======================================================================
		print(Utility.RemoveUnwantedIndent("""
				================================
				================================
				Packmol finished initial input file"""))
		for Sol in Solvents:
			print("""{0} {1}""".format(Sol, Sample[Sol]))
		print(Utility.RemoveUnwantedIndent("""
				box sizes : {0}, {1}, {2}
				===============================
				===============================
				""").format(Sample['LX'],Sample['LY'],Sample['LZ']) )
		
		#==================================================================================
		# the topology file
		#==================================================================================
		
		Topology = str("""
					#include "martini_v2.2.itp" ; modified with polarisable water
					#include "martini_v2.0_solvents.itp"
					""")
		#Copy the topology files for martini forcefield
		sub.call("""cp {0}/martini_v2.0_solvents.itp {0}/martini_v2.2.itp ./""".format(PathToDefault), shell= True)
		f = open(System+'.top','w')
		f.write(Utility.RemoveUnwantedIndent(Topology))
		
		f = open(System+'.top','a')
		f.write(Utility.RemoveUnwantedIndent("""
						
						[ system ]
						{0}
						
						""".format(System)))
		
		f.write("""\n[ molecules ]\n""")
		for Sol in Solvents:
			Topology = str("""{0} {1} \n""").format(Sol, Sample[Sol])
			f.write(Topology)
		f.close()
		
		#==================================================================================
		# the index file
		#==================================================================================
		cmd = str("""echo q | {0}make_ndx -f {1}.withbox.pdb -o {1}.ndx""").format(GROMACS_LOC_prefixPath,System)
		sub.call(cmd, shell=True)
		
		#==================================================================================
		# Output the files for other steps
		#==================================================================================
		Output = str("""{0}.withbox.pdb""").format(System)
		Index = str("""{0}.ndx""").format(System)
		return { 'SYSTEM': System, 'OUTPUT': Output, 'INDEX':Index}
































def TopologyDefoSu(Sample, PathToDefault):
	Headings = ["atomtypes","nonbond_params","moleculetype","atoms"]
	IsuTopology = {}
	IdefoTopology = {}
	DefoAtomtype = ""
	DefoNonBondParams = ""
	DefoMolType = ""
	DefoAtoms = ""
	
	IncludeSuTopologyFile = PathToDefault+"""/SU/SU_{0}.itp""".format(Sample['SU']['Version'])
	
	IncludeDefoTopologyFile = PathToDefault+"""/DEFO/DEFO_{0}.itp""".format(Sample['DEFO']['Version'])
	
	with open(IncludeSuTopologyFile,'r') as ITPSu:
		for Head in Headings:
			for heading_and_lines in Utility.group_by_heading(ITPSu, Head):
				lines = []
				if Head is not 'moleculetype':
					if Head is not 'atoms':
						lines.extend([';;;;;;; SU_{0}\n'.format(Sample['SU']['Version'])])
				lines.extend(heading_and_lines[2:])
				IsuTopology.update({Head:''.join(lines)})
	
	PosResSu = Utility.RemoveUnwantedIndent("""
										#ifdef SU_POSRES
										#include "su_posres.itp"
										#endif
										
										""")
	
	if '[ moleculetype ]' in IsuTopology['moleculetype'] or '[moleculetype]' in IsuTopology['moleculetype']:
		IsuTopology['moleculetype'] = IsuTopology['moleculetype'].replace('[ moleculetype ]', PosResSu+'[ moleculetype ]' )
	IsuTopology['moleculetype'] = '[ moleculetype ]\n'+IsuTopology['moleculetype']+'\n'+PosResSu
			
		
	with open(IncludeDefoTopologyFile,'r') as ITPDefo:
		for Head in Headings:
			for heading_and_lines in Utility.group_by_heading(ITPDefo, Head):
				lines = []
				if Head is not 'moleculetype':
					if Head is not 'atoms':
						lines.extend([';;;;;;; DEFO_{0}\n'.format(Sample['DEFO']['Version'])])
				lines.extend(heading_and_lines[2:])
				IdefoTopology.update({Head:''.join(lines)})
	
	if os.path.isfile("""martini_v2.2_{2}_{0}_DEFO_{1}.itp""".format(Sample['SU']['Version'], Sample['DEFO']['Version'],Sample['SU']['SuType'])):
		os.remove("""martini_v2.2_{2}_{0}_DEFO_{1}.itp""".format(Sample['SU']['Version'], Sample['DEFO']['Version'], Sample['SU']['SuType']))
	
	TempOut = ''
	
	Output = open("""martini_v2.2_{2}_{0}_DEFO_{1}.itp""".format(Sample['SU']['Version'], Sample['DEFO']['Version'], Sample['SU']['SuType']),'a+')
	
	with open(PathToDefault+"""/martini_v2.2.itp""",'r') as DefMartini:
		
		for line in DefMartini:
			if 'nonbond_params' in line:
				TempOut += IsuTopology['atomtypes']
				TempOut += '\n'
			if 'PLACE_FOR_SU' in line:
				TempOut += IsuTopology['nonbond_params']
				TempOut += '\n'
			else:
				TempOut += line
		
		
		for line in TempOut.splitlines():
			if 'nonbond_params' in line:
				Output.write(IdefoTopology['atomtypes'])
				Output.write('\n')
			if 'PLACE_FOR_DEFO' in line:
				Output.write(IdefoTopology['nonbond_params'])
				Output.write('\n')
			else:
				Output.write(line)
				Output.write('\n')
	
	
	Output.write("\n;;;;;;; DEFO_{0}\n".format(Sample['DEFO']['Version']))
	Output.write("[ moleculetype ]\n"+IdefoTopology['moleculetype']+'\n')
	Output.write("[ atoms ]\n"+IdefoTopology['atoms'])
	Output.write(Utility.RemoveUnwantedIndent("""
		
											#ifdef DEFO_POSRES
											#include "defo_posres.itp"
											#endif
											
											"""))
	
	Output.write("\n;;;;;;; SU_{0}\n".format(Sample['SU']['Version']))
	Output.write(IsuTopology['moleculetype']+'\n')
	Output.close()


def TopologyDefo(Sample, PathToDefault):
	Headings = ["atomtypes","nonbond_params","moleculetype","atoms"]
	Itopology = {}
	DefoAtomtype = ""
	DefoNonBondParams = ""
	DefoMolType = ""
	DefoAtoms = ""
	
	IncludeTopologyFile = PathToDefault+"""/DEFO/DEFO_{0}.itp""".format(Sample['DEFO']['Version'])
	
	with open(IncludeTopologyFile,'r') as ITPDefo:
		for Head in Headings:
			for heading_and_lines in Utility.group_by_heading( ITPDefo, Head):
				lines = []
				if Head is not 'moleculetype':
					if Head is not 'atoms':
						lines.extend([';;;;;;; DEFO_{0}\n'.format(Sample['DEFO']['Version'])])
				lines.extend(heading_and_lines[2:])
				Itopology.update({Head:''.join(lines)})
	
	if os.path.isfile("""martini_v2.2_DEFO_{0}.itp""".format(Sample['DEFO']['Version'])):
		os.remove("""martini_v2.2_DEFO_v{0}.itp""".format(Sample['DEFO']['Version']))
			
	Output = open("""martini_v2.2_DEFO_{0}.itp""".format(Sample['DEFO']['Version']),'a')
	with open(PathToDefault+"""/martini_v2.2.itp""",'r') as DefMartini:
		
		for line in DefMartini:
			if 'nonbond_params' in line:
				Output.write(Itopology['atomtypes'])
				#Output.write('[ nonbond_params ]')
			if 'PLACE_FOR_DEFO' in line:
				Output.write(Itopology['nonbond_params'])
				Output.write('')
			else:
				Output.write(line)
		
	Output.write(";;;;;;; DEFO_{0}\n".format(Sample['DEFO']['Version']))
	Output.write("[ moleculetype ]\n"+Itopology['moleculetype']+'\n')
	Output.write("[ atoms ]\n"+Itopology['atoms'])
	Output.write(Utility.RemoveUnwantedIndent("""
											#ifdef DEFO_POSRES
											#include "defo_posres.itp"
											#endif
											"""))
	Output.close()


def TopologySu(Sample, PathToDefault):
	Headings = ["atomtypes","nonbond_params","moleculetype","atoms"]
	IsuTopology = {}
	DefoAtomtype = ""
	DefoNonBondParams = ""
	DefoMolType = ""
	DefoAtoms = ""
	
	IncludeSuTopologyFile = PathToDefault+"""/SU/SU_{0}.itp""".format(Sample['SU']['Version'])
	
	with open(IncludeSuTopologyFile,'r') as ITPSu:
		for Head in Headings:
			for heading_and_lines in Utility.group_by_heading( ITPSu, Head):
				lines = []
				if Head is not 'moleculetype':
					if Head is not 'atoms':
						lines.extend([';;;;;;; SU_{0}\n'.format(Sample['SU']['Version'])])
				lines.extend(heading_and_lines[2:])
				IsuTopology.update({Head:''.join(lines)})
	
	PosResSu = Utility.RemoveUnwantedIndent("""
										#ifdef SU_POSRES
										#include "su_posres.itp"
										#endif
										
										""")
	
	if '[ moleculetype ]' in IsuTopology['moleculetype'] or '[moleculetype]' in IsuTopology['moleculetype']:
		IsuTopology['moleculetype'] = IsuTopology['moleculetype'].replace('[ moleculetype ]', PosResSu+'[ moleculetype ]' )
	IsuTopology['moleculetype'] = '[ moleculetype ]\n'+IsuTopology['moleculetype']+'\n'+PosResSu
	
	if os.path.isfile("""martini_v2.2_{1}_{0}.itp""".format(Sample['SU']['Version'], Sample['SU']['SuType'])):
		os.remove(""""martini_v2.2_{1}_{0}.itp""".format(Sample['SU']['Version'], Sample['SU']['SuType']))
			
	Output = open("""martini_v2.2_{1}_{0}.itp""".format(Sample['SU']['Version'], Sample['SU']['SuType']),'a')
	with open(PathToDefault+"""/martini_v2.2.itp""",'r') as DefMartini:
		
		for line in DefMartini:
			if 'nonbond_params' in line:
				Output.write(IsuTopology['atomtypes'])
				#Output.write('[ nonbond_params ]')
			if 'PLACE_FOR_SU' in line:
				Output.write(IsuTopology['nonbond_params'])
				Output.write('')
			else:
				Output.write(line)
		
	Output.write(";;;;;;; SU_{0}\n".format(Sample['SU']['Version']))
	Output.write(IsuTopology['moleculetype']+'\n')
	Output.close()














































def CopySample(Jobs, currentJob, step, PathToDefault):
	
	SampleToCopy = Jobs[ currentJob['PROTOCOL'][step]['samplenumber'] ]
	SampleToCopyName = Jobs[ currentJob['PROTOCOL'][step]['samplenumber'] ]['JOBID']
	currentJobName = currentJob['JOBID']
	
	copyMethod = currentJob['PROTOCOL'][step]['method']
	
	PDBfilepath = glob.glob('../'+SampleToCopyName+'/*.withbox.pdb')[0]
	TOPfilepath = glob.glob('../'+SampleToCopyName+'/*.top')[0]
	NDXfilepath = glob.glob('../'+SampleToCopyName+'/*.ndx')[0]
	
	if copyMethod:
		if copyMethod == "structure":
			if 'DEFO' in SampleToCopy or 'DEFO' in currentJob:
				if not ('DEFO' in SampleToCopy and 'DEFO' in currentJob):
					if 'SU' in SampleToCopy:
						print("DEFO is in {0} but not in {1} while using structure copy method.".format(SampleToCopyName, currentJobName))
					else:
						print("DEFO is in {0} but not in {1} while using structure copy method.".format(currentJobName, SampleToCopyName))
					return -1
				
			if 'SU' in SampleToCopy or 'SU' in currentJob:
				if not ('SU' in SampleToCopy and 'SU' in currentJob):
					if 'SU' in SampleToCopy:
						print("SU is in {0} but not in {1} while using structure copy method.".format(SampleToCopyName, currentJobName))
					else:
						print("SU is in {0} but not in {1} while using structure copy method.".format(currentJobName, SampleToCopyName))
					return -1
			
			if 'DEFO' in currentJob and 'SU' in currentJob:
				defoVersionBefore = SampleToCopy['DEFO']['Version']
				defoVersionCurrent = currentJob['DEFO']['Version']
				
				suVersionBefore = SampleToCopy['SU']['Version']
				suTypeBefore = SampleToCopy['SU']['SuType']
				
				suVersionCurrent = currentJob['SU']['Version']
				suTypeCurrent = currentJob['SU']['SuType']
				
				#Add columns for pdb format
				nbCols = len(suTypeBefore) - len(suTypeCurrent)
				if nbCols < 0:
					suTypeBefore += (-nbCols)*' '
				elif nbCols > 0:
					suTypeCurrent += nbCols*' '
					
				# If copy pdb file when modifying SuType, the pdb file should be changed with the correct SuType
				# The name should also be changed with the correct version and SuType, the number of SU being the same
				
				PDBfile = open(PDBfilepath,'r')
				PDBfileCurrent = PDBfile.read().replace(suTypeBefore, suTypeCurrent)
				PDBfile.close()
				
				
				PDBfilenameBefore = PDBfilepath.split('/')[-1]
				
				fixedPart = PDBfilenameBefore.split('_')[:-3]
				defoPart = PDBfilenameBefore.split('_')[-2]
				suPart = PDBfilenameBefore.split('_')[-1]
				
				suPart = suPart.replace(suTypeBefore.strip(' '), suTypeCurrent.strip(' '))
				suPart = suPart.replace(suVersionBefore, suVersionCurrent)
				
				defoPart = defoPart.replace(defoVersionBefore, defoVersionCurrent)
				
				PDBfilenameCurrent = fixedPart
				PDBfilenameCurrent.extend([defoPart, suPart])
				PDBfilenameCurrent = '_'.join(PDBfilenameCurrent)
				
				PDBfile = open(PDBfilenameCurrent, 'w')
				PDBfile.write(PDBfileCurrent)
				PDBfile.close()
				
				#Topology ##########################################################
				TOPfile = open(TOPfilepath, 'r')
				
				suTypeNVersionBefore = suTypeBefore+'_'+suVersionBefore
				suTypeNVersionCurrent = suTypeCurrent.strip(' ')+'_'+suVersionCurrent
				
				defoTypeNVersionBefore = 'DEFO_'+defoVersionBefore
				defoTypeNVersionCurrent = 'DEFO_'+defoVersionCurrent
				
				#Replacing for new .top file
				#For headers and [system]
				TOPfileCurrent = TOPfile.read().replace(suTypeNVersionBefore, suTypeNVersionCurrent)
				
				TOPfileCurrent = TOPfileCurrent.replace(defoTypeNVersionBefore, defoTypeNVersionCurrent)
				
				#For [molecules]
				TOPfileCurrent = TOPfileCurrent.replace(suTypeBefore, suTypeCurrent.strip(' '))
				
				TOPfile.close()
				TOPfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.top')
				TOPfile = open(TOPfilenameCurrent,'w')
				TOPfile.write(TOPfileCurrent)
				TOPfile.close()
				
				#Topology ITP ##########################################################
				TopologyDefoSu(currentJob, PathToDefault)
				
				#Index #########################################################
				NDXfile = open(NDXfilepath,'r')
				NDXfileCurrent = NDXfile.read()
				NDXfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.ndx')
				NDXfile.close()
				
				NDXfile = open(NDXfilenameCurrent,'w')
				NDXfile.write(NDXfileCurrent)
				NDXfile.close()
				
				#POSRES ########################################################
				sub.call('cp ../'+SampleToCopyName+'/*posres.itp .', shell= True)
				#LIPIDS ITP ####################################################
				sub.call('cp ../'+SampleToCopyName+'/*lipids.itp .', shell= True)
			
			
			elif 'DEFO' not in currentJob and 'SU' in currentJob:
				suVersionBefore = SampleToCopy['SU']['Version']
				suTypeBefore = SampleToCopy['SU']['SuType']
				
				suVersionCurrent = currentJob['SU']['Version']
				suTypeCurrent = currentJob['SU']['SuType']
				
				#Add columns for pdb format
				nbCols = len(suTypeBefore) - len(suTypeCurrent)
				if nbCols < 0:
					suTypeBefore += (-nbCols)*' '
				elif nbCols > 0:
					suTypeCurrent += nbCols*' '
				
				PDBfile = open(PDBfilepath,'r')
				PDBfileCurrent = PDBfile.read().replace(suTypeBefore, suTypeCurrent)
				PDBfile.close()
				
				PDBfilenameBefore = PDBfilepath.split('/')[-1]
				
				fixedPart = PDBfilenameBefore.split('_')[:-2]
				suPart = PDBfilenameBefore.split('_')[-1]
				
				suPart = suPart.replace(suTypeBefore.strip(' '), suTypeCurrent.strip(' '))
				suPart = suPart.replace(suVersionBefore, suVersionCurrent)
				
				PDBfilenameCurrent = fixedPart
				PDBfilenameCurrent.append(suPart)
				PDBfilenameCurrent = '_'.join(PDBfilenameCurrent)
				
				PDBfile = open(PDBfilenameCurrent, 'w')
				PDBfile.write(PDBfileCurrent)
				PDBfile.close()
				
				#Topology ##########################################################
				TOPfile = open(TOPfilepath, 'r')
				
				suTypeNVersionBefore = suTypeBefore.strip(' ')+'_'+suVersionBefore
				suTypeNVersionCurrent = suTypeCurrent.strip(' ')+'_'+suVersionCurrent
				
				#Replacing for new .top file
				#For headers and [system]
				TOPfileCurrent = TOPfile.read().replace(suTypeNVersionBefore, suTypeNVersionCurrent)
				
				#For [molecules]
				TOPfileCurrent = TOPfileCurrent.replace(suTypeBefore, suTypeCurrent.strip(' '))
				
				TOPfile.close()
				TOPfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.top')
				TOPfile = open(TOPfilenameCurrent,'w')
				TOPfile.write(TOPfileCurrent)
				TOPfile.close()
				
				#Topology ITP ##########################################################
				TopologySu(currentJob, PathToDefault)
				
				#Index #########################################################
				NDXfile = open(NDXfilepath,'r')
				NDXfileCurrent = NDXfile.read()
				NDXfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.ndx')
				NDXfile.close()
				
				NDXfile = open(NDXfilenameCurrent,'w')
				NDXfile.write(NDXfileCurrent)
				NDXfile.close()
				
				#POSRES ########################################################
				sub.call('cp ../'+SampleToCopyName+'/*posres.itp .', shell= True)
				#LIPIDS ITP ####################################################
				sub.call('cp ../'+SampleToCopyName+'/*lipids.itp .', shell= True)
			
			
			elif 'DEFO' in currentJob and 'SU' not in currentJob:
				defoVersionBefore = SampleToCopy['DEFO']['Version']
				defoVersionCurrent = currentJob['DEFO']['Version']
				
				# If copy pdb file when modifying SuType, the pdb file should be changed with the correct SuType
				# The name should also be changed with the correct version and SuType, the number of SU being the same
				
				PDBfile = open(PDBfilepath,'r')
				PDBfileCurrent = PDBfile.read()
				PDBfile.close()
				
				
				PDBfilenameBefore = PDBfilepath.split('/')[-1]
				
				fixedPart = PDBfilenameBefore.split('_')[:-2]
				defoPart = PDBfilenameBefore.split('_')[-1]
				
				defoPart = defoPart.replace(defoVersionBefore, defoVersionCurrent)
				
				PDBfilenameCurrent = fixedPart
				PDBfilenameCurrent.append(defoPart)
				PDBfilenameCurrent = '_'.join(PDBfilenameCurrent)
				
				PDBfile = open(PDBfilenameCurrent, 'w')
				PDBfile.write(PDBfileCurrent)
				PDBfile.close()
				
				#Topology ##########################################################
				TOPfile = open(TOPfilepath, 'r')
				
				defoTypeNVersionBefore = 'DEFO_'+defoVersionBefore
				defoTypeNVersionCurrent = 'DEFO_'+defoVersionCurrent
				
				#Replacing for new .top file
				#For headers and [system]
				TOPfileCurrent = TOPfile.read().replace(defoTypeNVersionBefore, defoTypeNVersionCurrent)
				
				TOPfile.close()
				TOPfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.top')
				
				TOPfile = open(TOPfilenameCurrent,'w')
				TOPfile.write(TOPfileCurrent)
				TOPfile.close()
				
				#Topology ITP ##########################################################
				TopologyDefo(currentJob, PathToDefault)
				
				#Index #########################################################
				NDXfile = open(NDXfilepath,'r')
				NDXfileCurrent = NDXfile.read()
				NDXfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.ndx')
				NDXfile.close()
				
				NDXfile = open(NDXfilenameCurrent,'w')
				NDXfile.write(NDXfileCurrent)
				NDXfile.close()
				
				#POSRES ########################################################
				sub.call('cp ../'+SampleToCopyName+'/*posres.itp .', shell= True)
				#LIPIDS ITP ####################################################
				sub.call('cp ../'+SampleToCopyName+'/*lipids.itp .', shell= True)
			
			
			else:
				sub.call('cp '+PDBfilepath+' .', shell= True)
				sub.call('cp '+TOPfilepath+' .', shell= True)
				sub.call('cp '+NDXfilepath+' .', shell= True)
				sub.call('cp ../'+SampleToCopyName+'/*.itp .', shell= True)
				
				PDBfilenameCurrent = PDBfilepath.split('/')[-1]
				NDXfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.ndx')
				
			
				
				
		elif copyMethod == "all" :
			sub.call('cp '+PDBfilepath+' .', shell= True)
			sub.call('cp '+TOPfilepath+' .', shell= True)
			sub.call('cp '+NDXfilepath+' .', shell= True)
			sub.call('cp ../'+SampleToCopyName+'/*.itp .', shell= True)
			
			PDBfilenameCurrent = PDBfilepath.split('/')[-1]
			NDXfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.ndx')
	
	
	
	else:
		sub.call('cp '+PDBfilepath+' .', shell= True)
		sub.call('cp '+TOPfilepath+' .', shell= True)
		sub.call('cp '+NDXfilepath+' .', shell= True)
		sub.call('cp ../'+SampleToCopyName+'/*.itp .', shell= True)
		
		PDBfilenameCurrent = PDBfilepath.split('/')[-1]
		NDXfilenameCurrent = PDBfilenameCurrent.replace('.withbox.pdb','.ndx')
	
	System = PDBfilenameCurrent.strip('.withbox.pdb')
	return { 'SYSTEM': System, 'OUTPUT': PDBfilenameCurrent, 'INDEX':NDXfilenameCurrent}



def WriteMDP(Sample, step, defaultMDP, Version):
	
	# Get the parameters for the step
	stepMD = Sample['PROTOCOL'][step]
	
	if 'DEFO' in Sample:
		stepAfterInit = int(step)-1
		presetName = Sample['DEFO']['defoProtocol'][stepAfterInit].strip(' ')
		defoPresetForStep = Sample['DEFO']['presets'][presetName]
		
	if 'SU' in Sample:
		stepAfterInit = int(step)-1
		presetName = Sample['SU']['suProtocol'][stepAfterInit].strip(' ')
		suPresetForStep = Sample['SU']['presets'][presetName]
		SU_TYPE = Sample['SU']['SuType']
	
	#Showing the final parameters
	print("Writing step {0} ==========================".format(stepMD['stepType']))
	######for param in stepMD:
		######if param != 'stepType':
			######print(param +' = '+str(stepMD[param]))
	
	# Preparation with Parameters.csv before writing .mdp files ================
	
	#Bool to check that the run is of NVT or NPT type
	NPT_or_NVT = False
	
	#Look if the run is NPT or NVT
	if stepMD['stepType'].startswith('NPT') or stepMD['stepType'].startswith('NVT'):
		NPT_or_NVT = True
		# if THERMOSTAT is in Parameters.csv the values set here will override the values set
		# for NPT and NVT runs
		if 'THERMOSTAT' in Sample:
			if 'ref-t' in Sample['THERMOSTAT']:
				Sample['PROTOCOL'][step].update( {'ref-t': Sample['THERMOSTAT']['ref-t'] } )
	
	#Energy and Temperature groups
	Egrps = ''
	Tcgrps = ''
	# Multiply tau-t and ref-t by the number of grps
	Tau_Tcgrps = ''
	T_Tcgrps = ''
	autoEGrps = True
	autoTcGrps = True
	autoTau_TcGrps = True
	autoT_TcGrps = True
	#Creates grps for all lipid found except if tau-t and ref-t contain multiple values
	#Exception
	if 'energygrps' in Sample['PROTOCOL'][step]:
		if ' ' in str(Sample['PROTOCOL'][step]['energygrps']):
			Egrps += str(Sample['PROTOCOL'][step]['energygrps'])
			autoEGrps = False
			
	if 'tc-grps' in Sample['PROTOCOL'][step]:
		if ' ' in str(Sample['PROTOCOL'][step]['tc-grps']):
			Tcgrps += str(Sample['PROTOCOL'][step]['tc-grps'])
			autoTcGrps = False
			
	if 'tau-t' in Sample['PROTOCOL'][step]:
		if ' ' in str(Sample['PROTOCOL'][step]['tau-t']):
			Tau_Tcgrps += str(Sample['PROTOCOL'][step]['tau-t'])
			autoTau_TcGrps = False
			
	if 'ref-t' in Sample['PROTOCOL'][step]:
		if ' ' in str(Sample['PROTOCOL'][step]['ref-t']):
			T_Tcgrps += str(Sample['PROTOCOL'][step]['ref-t'])
			autoT_TcGrps = False
			
	#Automatic generation
	if autoEGrps or autoTcGrps or autoTau_TcGrps or autoT_TcGrps:
		for lipid in LipidsList:
			if lipid in Sample:
				if autoEGrps:
					Egrps += lipid+' '
				if(NPT_or_NVT):
					if autoTcGrps:
						Tcgrps += lipid+' '
					if autoTau_TcGrps:
						if 'tau-t' in Sample['PROTOCOL'][step] and autoTau_TcGrps:
							Tau_Tcgrps += str(Sample['PROTOCOL'][step]['tau-t'])+' '
						else:
							pass
					if autoT_TcGrps:
						if 'ref-t' in Sample['PROTOCOL'][step] and autoT_TcGrps:
							T_Tcgrps += str(Sample['PROTOCOL'][step]['ref-t'])+' '
						else:
							pass
					
		#Creates grps for all solvents found
		for sol in SolventsList:
			if sol in Sample:
				if autoEGrps:
					Egrps += sol+' '
				if NPT_or_NVT:
					if autoTcGrps:
						Tcgrps += sol+' '
						if 'tau-t' in Sample['PROTOCOL'][step] and autoTau_TcGrps:
							Tau_Tcgrps += str(Sample['PROTOCOL'][step]['tau-t'])+' '
						if 'ref-t' in Sample['PROTOCOL'][step] and autoT_TcGrps:
							T_Tcgrps += str(Sample['PROTOCOL'][step]['ref-t'])+' '
	
	#Updates the sample dictionnary with the values created above
	if(NPT_or_NVT):
		Sample['PROTOCOL'][step].update({ 'tc-grps': Tcgrps })
		Sample['PROTOCOL'][step].update({ 'tau-t': Tau_Tcgrps })
		Sample['PROTOCOL'][step].update({ 'ref-t': T_Tcgrps })
		
	#Creates grps for DEFO and modifies the parameters in .mdp if in the preset
	if 'DEFO' in Sample:
		defoMdpParams = '\n\n ;;; Parameters for Defo ;;; \n\n'
		Egrps += 'DEFO'
		if 'ref-t' in defoPresetForStep or 'tau-t' in defoPresetForStep: 
			Tcgrps += ' DEFO'
		
		for defoParam, defoParamValue in defoPresetForStep.items():
			if defoParam != 'posres':
				if defoParam in Sample['PROTOCOL'][step]:
					currValue = str(Sample['PROTOCOL'][step][defoParam])
					Sample['PROTOCOL'][step][defoParam] = currValue + ' ' + str(defoParamValue) + ' '
				else:
					defoMdpParams += defoParam + '		= '+ str(defoParamValue) + ' \n'
		
		if(NPT_or_NVT):
			Sample['PROTOCOL'][step].update({ 'tc-grps': Tcgrps })
	
	#Creates grps for SU and modifies the parameters in .mdp if in the preset
	if 'SU' in Sample:
		suMdpParams = '\n\n ;;; Parameters for Su ;;; \n\n'
		Egrps += ' '+ SU_TYPE
		if 'ref-t' in suPresetForStep or 'tau-t' in suPresetForStep: 
			Tcgrps += ' '+ SU_TYPE
		
		for suParam, suParamValue in suPresetForStep.items():
			if suParam != 'posres':
				if suParam in Sample['PROTOCOL'][step]:
					currValue = str(Sample['PROTOCOL'][step][suParam])
					Sample['PROTOCOL'][step][suParam] = currValue + ' ' + str(suParamValue) + ' '
				
				else:
					suMdpParams += suParam + '		= '+ str(suParamValue) + ' \n'
		suMdpParams += 'wall-atomtype' + '		= '+ str(Sample['SU']['SuType']) + ' \n'
		if(NPT_or_NVT):
			Sample['PROTOCOL'][step].update({ 'tc-grps': Tcgrps })
			
	#Updates the Energy groups
	Sample['PROTOCOL'][step].update({ 'energygrps': Egrps })
	
	# End of preparation =======================================================
	
	
	
	# Beginning of writing .mdp ================================================
	
	#Looks for the parameters in Default MDP file and set the values
	#as chosen in Parameters.csv
	OUTPUT = open(stepMD['stepType']+'.mdp','w')
	
	#Position restraining for DEFO and Su if in preset
	
	if 'DEFO' in Sample and 'SU' not in Sample:
		if 'posres'  in defoPresetForStep:
			if defoPresetForStep['posres'] == 'on': OUTPUT.write("define = -DDEFO_POSRES")
			
	elif 'SU' in Sample and 'DEFO' not in Sample:
		if 'posres' in suPresetForStep:
			if suPresetForStep['posres'] == 'on': OUTPUT.write("define = -DSU_POSRES")
			
	elif 'SU' in Sample and 'DEFO' in Sample:
		if 'posres' in suPresetForStep and 'posres' in defoPresetForStep:
			if suPresetForStep['posres'] == 'on' and defoPresetForStep['posres'] == 'on' :
				OUTPUT.write("define = -DSU_POSRES -DDEFO_POSRES")
				
			elif suPresetForStep['posres'] == 'off' and defoPresetForStep['posres'] == 'on' :
				OUTPUT.write("define = -DDEFO_POSRES")
				
			elif suPresetForStep['posres'] == 'on' and defoPresetForStep['posres'] == 'off' :
				OUTPUT.write("define = -DSU_POSRES")
				
				
	
		
	CopyOriginalLine = True
	for i, line in enumerate(defaultMDP):
		if not line.startswith(';'):
			for key in stepMD:
				if key in line and key != '' and (line.partition(' ')[0] == key or line.partition('=')[0] == key):
					OUTPUT.write(key+'			= '+ str(stepMD[key])+'\n')
					CopyOriginalLine = False
					continue
			if(CopyOriginalLine):
				OUTPUT.write(line)
			CopyOriginalLine=True

	OUTPUT = open(stepMD['stepType']+'.mdp','a+')
	OUTPUT.write('\n ;Parameters not in default file : \n\n')
	
	# Open the file again to add the parameters not in default MDP
	COMPARE = open(stepMD['stepType']+'.mdp','r').read()
	for key in stepMD:
		if key not in COMPARE and key != 'stepType':
			OUTPUT.write(key+'			= '+ str(stepMD[key])+'\n')
	
	# Writes parameters related to Defo at the end
	if 'DEFO' in Sample:
		OUTPUT.write(defoMdpParams)
	
	# Automatic filling for defo (not done yet as it is not 
	# really useful ) 
	if(0): #Currently off as not yet implemented
		if 'DEFO' in Sample and Sample['TYPE'] == 'BILAYER':
			Lipid = ''
			for lipid in LipidsList:
				if lipid in Sample:
					Lipid = lipid
			if Sample['DEFO']['Height'] == 'follow':
				if 'SU' in Sample:
					if Version.startswith('4'):
						#finding current lipid used
						OUTPUT.write('pull-group1		= {0}\n'.format(Lipid))
					else:
						#finding current lipid used
						pullParams = """
									pull-ngroups			= 2
									pull-group1-name			= defoBi
									pull-group1-pbcatom			= 0
									
									pull-group2-name			= defoMono
									pull-group2-pbcatom = 0
									
									pull-ncoords			= 2
									pull-coord1-groups
									pull-coord2-groups			= 0 2
									""".format()
						OUTPUT.write('pull-group1-name		= {0}\n'.format(Lipid))
				
			else:
				pass
		
	# Writes parameters related to Su at the end
	if 'SU' in Sample:
		OUTPUT.write(suMdpParams)
	
	
	defaultMDP.seek(0)

	OUTPUT.close()
