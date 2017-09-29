# GJB
###################################################
# USAGE for creating an input serie of gromacs simulations 
# - MARTINI lipidic bilayers or trilayers (DSPC, DPPC, DLPC) + solvent
# - one or two solvants only (W, )
# - eventually with vaccum
# - eventually with default pores, or support, or wall
# (once the Parameters.csv  and PbsInfo.csv have been created) 
###################################################
TO DO INPUTS FOR LYNX
> Script.py --PBS 

TO DO INPUTS LOCALLY
> Script.py --missing

###################################################
# FORMAT for Parameters.cvs
################################
Basic format Rules : 
- comment with # 
- separator with | 



=> DEFO 

Some problem with the fact that the parameters_defo.cvs has the priority of Parameters.cvs, 
the variables are added to T323K_v1.0 instead of replaced by..


=> SUPPORT 


=> WALL


=> TRILAYER



# TOOLS for ANALYSIS
##########################

Results.py =>  Analysis each simulation of a serie, depending on the type of simulation, and create a Summary over all simulations of this serie.

USAGE 

