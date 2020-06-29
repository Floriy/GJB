#!/home/cloison/Softwares/ANACONDA/anaconda3/bin/python
# -*- coding: utf-8 -*-

#
# example of usage :
# > ./multiply_substrate_nXnYnZ.py -i subd16_no_surftens_defo_DSPC.gro -nx 2 -ny 2 -o subd16_no_surftens_defo_DSPC_mult_2_2_1.gro -gmx /home/cloison/Softwares/GROMACS/2016/bin
#

import subprocess as sub
import argparse

def main():     
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input', dest='inputfilename', type=str, required=True,
						help='input gro filename')
	
	parser.add_argument('-o','--output', dest='outputfilename', type=str, default='out.gro',
						help='output gro filename')
	
	parser.add_argument('-nx','--numberX', dest='nx', type=int, default="1",
						help='number of time the substrate is reproduced in x direction')
	parser.add_argument('-ny','--numberY', dest='ny', type=int, default="1",
						help='number of time the substrate is reproduced in y direction')
	parser.add_argument('-nz','--numberZ', dest='nz', type=int, default="1",
						help='number of time the substrate is reproduced in z direction')
	parser.add_argument('-cfn','--commandFileName', dest='cfn', type=str, default="gmxcommand.sh",
						help='number of time the substrate is reproduced in z direction')
	parser.add_argument('-gmx','--gmxpath', dest='gmxpath', type=str, default="/home/cloison/Softwares/GROMACS/2016/bin",
						help='number of time the substrate is reproduced in z direction')
	parser.add_argument('-clean','--clean', dest='clean', type=bool, default= True,
						help='delete intermediate files after work')
   
   
	# OPEN GMX COMMAND FILE AND WRITE PREAMBLE
	args = parser.parse_args()      
	script_file_run = open(args.cfn,'w')
	sub.call("""chmod a+x {0}""".format(args.cfn), shell=True)
	script_file_run.write("#!/bin/bash -x \n")
	script_file_run.write("source {0}/GMXRC \n".format(args.gmxpath))
	
	
	# READ BOX SIZE AND PARTICLE NUMBER
	boxsize = sub.check_output("""tail -1 {0}""".format(args.inputfilename), shell=True)
	LX, LY, LZ = [float(L) for L in  boxsize.split()]
	newLX = LX * args.nx
	newLY = LY * args.ny
	newLZ = LZ * args.nz
	newboxsize = f"{newLX} {newLY} {newLZ}"
	
	# CREATE NEW SYSTEM WITH RIGH BOX SIZE
	outbox=args.inputfilename.replace('.gro','')+f"_newbox.gro"
	script_file_run.write(f"gmx editconf -f {args.inputfilename} -box {newboxsize} -o {outbox} -noc \n")
	
	# UT THE INITIAL SYSTEM
	inputfile=open(args.inputfilename, "r")
	inputfile_lines=inputfile.readlines()
	inputfile.close()
	
	particlenumber = int(inputfile_lines[1])
	print(f"particlenumber = {particlenumber}")
	
	
	# CREATE MULTIPLIED SYSTEMS IN THE GMX INPUT AND EXECUTE IT
	partindex = 1
	for ix in range(0,args.nx):
		for iy in range(0,args.ny):
			for iz in range(0,args.nz):
				if ix+iy+iz is not 0:
					print(f"translate box of vector {ix},{iy},{iz}")
					out=args.inputfilename.replace('.gro','')+f"_trans_{ix}_{iy}_{iz}.gro"
					vector=f"{ix*LX} {iy*LY} {iz*LZ}"
					partindex = partindex + particlenumber
					script_file_run.write(f"gmx editconf -f {outbox} -translate {vector} -resnr {partindex} -o {out}\n")
	    
	script_file_run.close()
	sub.call(f"./{args.cfn} \n", shell = True)
	
	### CREATE THE FINAL GRO FILE BY CONCATENATING THE OLD SYSTEMS
	
	# ADD THE INITIAL SYSTEM IN THE FINAL GRO
	outfile=open(args.outputfilename, "w")
	outfile.write(inputfile_lines[0])
	outfile.write(f"{partindex + particlenumber-1}\n") 
	for line in inputfile_lines[2:particlenumber+2]:
		outfile.write(line)
	
	# ADD THE MULTIPLIED SYSTEMS IN THE FINAL GRO
	for ix in range(0,args.nx):
		for iy in range(0,args.ny):
			for iz in range(0,args.nz):
				if ix+iy+iz is not 0:
					out=args.inputfilename.replace('.gro','')+f"_trans_{ix}_{iy}_{iz}.gro"
					nextfile=open(out, "r")
					for line in nextfile.readlines()[2:particlenumber+2]:
						outfile.write(line)
					nextfile.close()
					
	outfile.write(f"{newboxsize}\n")
	script_file_run.close() 
	
	# CLEANN AFTER WORK
	if args.clean :
		for ix in range(0,args.nx):
			for iy in range(0,args.ny):
				for iz in range(0,args.nz):
					if ix+iy+iz is not 0:
						out=args.inputfilename.replace('.gro','')+f"_trans_{ix}_{iy}_{iz}.gro"       
						sub.call(f"rm {out}\n", shell = True)
		sub.call(f"rm {outbox}\n", shell = True)
		sub.call(f"rm {args.cfn}\n", shell = True)
    
if __name__ == "__main__":
	main()
	
	
