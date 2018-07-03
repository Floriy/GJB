#!/usr/bin/python3
# -*- coding: utf-8 -*-


import argparse #Argument parser module !
import sys
import glob
import os
import shutil
import subprocess as sub
import traceback
import math, cmath
from collections import OrderedDict
import gc #garbage collector
import time
import Utility as ut
from io import StringIO

#global variables
extensions = ['ndx','gro','xtc']

#def check_existing_dir
def depth(d, level=1):
	if not isinstance(d, dict) or not d:
		return level
	
	return max(depth(d[k], level + 1) for k in d)

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

#function for padding
def padding_grid(grid_x, grid_y, box_x, box_y, pad_x, pad_y):
	padx_begin = np.subtract(grid_x[-pad_x-1:-1], box_x) 
	padx_end = np.add(grid_x[1:pad_x+1], box_x)
	
	pady_begin = [np.subtract(grid_y[0][-pad_y-1:-1], box_y)]
	pady_end = [np.add(grid_y[0][1:pad_y+1], box_y)]

	x_right_min = np.amin(grid_x[-pad_x-1:-1])
	x_left_max = np.amax(grid_x[1:pad_x+1])
	
	y_up_min = np.amin(grid_y[0][-pad_x-1:-1])
	y_down_max = np.amax(grid_y[0][1:pad_x+1])
	
	# padd the x axis
	grid_x = np.insert(grid_x, [0], padx_begin, axis=0)
	grid_x_el_len = len(grid_x)
	grid_x = np.insert(grid_x, [grid_x_el_len], padx_end, axis=0)
	
	# padd the y axis
	grid_y = np.insert(grid_y,[0], pady_begin, axis=1)
	grid_y_el_len = grid_y[0].size
	grid_y = np.insert(grid_y, [grid_y_el_len], pady_end, axis=1)
	
	
	grid_x = np.array([ np.pad(X, pad_width=(pad_x, pad_y), mode='edge') for X in grid_x ])
	grid_y = np.insert(grid_y, [len(grid_y)], grid_y[:pad_x+pad_y], axis=0)
	
	
	#print(x_right_min)
	#print(x_left_max)
	
	#print(y_up_min)
	#print(y_down_max)
	return grid_x, grid_y, x_right_min, x_left_max, y_up_min, y_down_max

def padding3d_grid(grid_x, grid_y, grid_z, box_x, box_y, box_z, pad_x, pad_y, pad_z):
	padx_begin = np.subtract(grid_x[-pad_x-1:-1], box_x) 
	padx_end = np.add(grid_x[1:pad_x+1], box_x)
	
	pady_begin = [np.subtract(grid_y[0][-pad_y-1:-1], box_y)]
	pady_end = [np.add(grid_y[0][1:pad_y+1], box_y)]
	
	padz_begin = [np.subtract(grid_z[0][-pad_z-1:-1], box_z)]
	padz_end = [np.add(grid_z[0][1:pad_z+1], box_z)]

	x_right_min = np.amin(grid_x[-pad_x-1:-1])
	x_left_max = np.amax(grid_x[1:pad_x+1])
	
	y_up_min = np.amin(grid_y[0][-pad_y-1:-1])
	y_down_max = np.amax(grid_y[0][1:pad_y+1])
	
	z_up_min = np.amin(grid_z[0][-pad_z-1:-1])
	z_down_max = np.amax(grid_z[0][1:pad_z+1])
	
	# padd the x axis
	grid_x = np.insert(grid_x, [0], padx_begin, axis=0)
	grid_x_el_len = len(grid_x)
	grid_x = np.insert(grid_x, [grid_x_el_len], padx_end, axis=0)
	
	# padd the y axis
	grid_y = np.insert(grid_y,[0], pady_begin, axis=1)
	grid_y_el_len = grid_y[0].size
	grid_y = np.insert(grid_y, [grid_y_el_len], pady_end, axis=1)
	
	# padd the y axis
	grid_z = np.insert(grid_z,[0], padz_begin, axis=2)
	grid_z_el_len = grid_z[0].size
	grid_z = np.insert(grid_z, [grid_z_el_len], padz_end, axis=2)
	
	
	grid_x = np.array([ np.pad(X, pad_width=(pad_x, pad_y), mode='edge') for X in grid_x ])
	grid_y = np.insert(grid_y, [len(grid_y)], grid_y[:pad_x+pad_y], axis=0)
	grid_z = np.insert(grid_z, [len(grid_z)], grid_z[:pad_x+pad_y], axis=1)
	
	
	#print(x_right_min)
	#print(x_left_max)
	
	#print(y_up_min)
	#print(y_down_max)
	return grid_x, grid_y, grid_z, x_right_min, x_left_max, y_up_min, y_down_max, z_up_min, z_down_max, 

def padding_values(dataframe, pad_x, pad_y, box_x, box_y, x_right_min, x_left_max, y_up_min, y_down_max):
	
	x_sorted_df = dataframe.sort_values('X coords', ascending=True)
	y_sorted_df = dataframe.sort_values('Y coords', ascending=True)
	
	padded = dataframe
	
	concatenate = [dataframe]
	for leaflet in ['lower leaflet', 'upper leaflet']:
		xpbc_right = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['X coords'] < x_left_max),:]
		xpbc_left = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['X coords'] > x_right_min),:]
		
		ypbc_up = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['Y coords'] < y_down_max),:]
		ypbc_down = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['Y coords'] > y_up_min),:]
		
		xpbc_right.loc[:,'X coords'] += box_x
		xpbc_left.loc[:,'X coords'] -= box_x
		ypbc_up.loc[:,'Y coords'] += box_y
		ypbc_down.loc[:,'Y coords'] -= box_y
		
		concatenate.extend([xpbc_right, xpbc_left,ypbc_up,ypbc_down])
	
	padded = pd.concat(concatenate, ignore_index=True)
	return padded

def padding3d_values(dataframe, pad_x, pad_y, box_x, box_y, x_right_min, x_left_max, y_up_min, y_down_max):
	
	x_sorted_df = dataframe.sort_values('X coords', ascending=True)
	y_sorted_df = dataframe.sort_values('Y coords', ascending=True)
	
	padded = dataframe
	
	concatenate = [dataframe]
	for leaflet in ['lower leaflet', 'upper leaflet']:
		xpbc_right = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['X coords'] < x_left_max),:]
		xpbc_left = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['X coords'] > x_right_min),:]
		
		ypbc_up = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['Y coords'] < y_down_max),:]
		ypbc_down = dataframe.loc[(dataframe['leaflet'] == leaflet) & (dataframe['Y coords'] > y_up_min),:]
		
		xpbc_right.loc[:,'X coords'] += box_x
		xpbc_left.loc[:,'X coords'] -= box_x
		ypbc_up.loc[:,'Y coords'] += box_y
		ypbc_down.loc[:,'Y coords'] -= box_y
		
		concatenate.extend([xpbc_right, xpbc_left,ypbc_up,ypbc_down])
	
	padded = pd.concat(concatenate, ignore_index=True)
	return padded


def crop_grid(grid, pad_col, pad_row):
	
	grid = np.delete(grid, np.s_[:pad_row], axis=0)
	grid = np.delete(grid, np.s_[-pad_row:], axis=0)
	grid = np.delete(grid, np.s_[:pad_col], axis=1)
	grid = np.delete(grid, np.s_[-pad_col:], axis=1)
	
	return grid
"""
def pad_axes(x_values, y_values, box_x, box_y, pad_x, pad_y):
	x_values = np.array(x_values)
	y_values = np.array(y_values)
	
	padx_begin = x_values[-pad_x-1:-1]
	padx_end = x_values[1:pad_x+1]
	print(x_values)
	print(y_values)
	pady_begin = np.subtract(y_values[-pad_y-1:-1], box_y) 
	pady_end = np.add(y_values[1:pad_y+1], box_y)
	
	
	# padd the x axis
	x_values = np.insert(x_values, [0], padx_begin)
	xval_len = x_values.size
	x_values = np.insert(x_values, [xval_len], padx_end)
	
	# padd the y axis
	y_values = np.insert(y_values, [0], pady_begin)
	yval_len = y_values.size
	y_values = np.insert(y_values, [yval_len], pady_end)
	print(x_values)
	print()
	print(y_values)
	print("\n\n\n")
	return x_values, y_values
"""
	
def interpolate_grid(dataframe, bdim, name, pad_x, pad_y):
	time = int(bdim[0])
	box_x, box_y, box_z = [float(val) for val in bdim[1:]]
	
	## Building data from plotting
	grid_x, grid_y = np.mgrid[0:box_x:50j, 0:box_y:50j]
	
	grids_per_leaflet = {'lower leaflet': {'grid':None, 'plot':0},
						'upper leaflet': {'grid':None, 'plot':1}}
	
	property_name = dataframe.columns[-1]
	
	np.set_printoptions(threshold=np.inf, linewidth=np.inf)
	
	#padding array for periodic boundary conditions
	grid_x, grid_y, x_right_min, x_left_max, y_up_min, y_down_max = padding_grid(grid_x, grid_y, box_x, box_y, pad_x, pad_y)
	
	header = "#TIME (ps) {0:d}\n#\n".format(time)
	
	dataframe = padding_values(dataframe, pad_x, pad_y, box_x, box_y, x_right_min, x_left_max, y_up_min, y_down_max)
	
	for leaflet in grids_per_leaflet:
		x_values = dataframe[ dataframe['leaflet'] == leaflet ]['X coords']
		y_values = dataframe[ dataframe['leaflet'] == leaflet ]['Y coords']
		prop_values = dataframe.loc[ dataframe['leaflet'] == leaflet ][property_name]
		
		##Getting the values per leaflet
		points = np.stack((np.array(x_values).T, np.array(y_values).T), axis=-1)
		values = np.array(prop_values)
		
		# Interpolation of property values on grid of box_x and box_y dimensions
		grid = griddata(points, values, (grid_x, grid_y))
		
		#add the grid to the dict for later
		grids_per_leaflet[leaflet]['grid'] = grid
		
		
		grid_to_write = header + np.array2string(grid, separator=',')
		with open(name + leaflet.replace(' ','_') + '.grid','w') as grid_file:
			grid_file.write(grid_to_write)
	
	gridx_to_write = header + np.array2string(grid_x, separator=',')
	with open(name[:-1] + '.gridx','w') as gridx_file:
			gridx_file.write(gridx_to_write)
	
	gridy_to_write = header + np.array2string(grid_y, separator=',')
	with open(name[:-1] + '.gridy','w') as gridy_file:
			gridy_file.write(gridy_to_write)
	
	return [grid_x, grid_y , grids_per_leaflet, time]
	
	
def plot_grid(grid, prop, name, bdim):
	grid_x = grid[0]
	grid_y = grid[1]
	grid_dict = grid[2]
	time = grid[3]
	
	box_x, box_y, box_z = [float(val) for val in bdim[1:]]
	
	#should limit the range due to 'periodic' grid for interpolation
	bounds = None
	levels = None
	if prop == 'ORDER':
		bounds = np.linspace(-0.5, 1, 15)
		levels = np.linspace(-0.5, 1, 15)
		colormap = mcm.viridis
	if prop == 'APL':
		bounds = np.linspace(0.55, 0.75, 20)
		levels = np.linspace(0.55, 0.75, 20)
		colormap = mcm.plasma
	if prop == 'THICKNESS':
		bounds = np.linspace(1.0, 5.0, 20)
		levels = np.linspace(1.0, 5.0, 20)
		colormap = mcm.Greys_r
	#figure for plots
	fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,6), dpi=72, gridspec_kw = {'width_ratios':[5,5,1]})
	
	cs = [None, None]
	for leaflet in grid_dict:
		# Setting the axes on which to plot
		Nbleaflet = grid_dict[leaflet]['plot']
		grid_dict[leaflet]['plot'] = axs[ grid_dict[leaflet]['plot'] ]
		# Plot map
		cs[ Nbleaflet ] = grid_dict[leaflet]['plot'].contourf(grid_x, grid_y, grid_dict[leaflet]['grid'], cmap = colormap, levels=levels)
		grid_dict[leaflet]['plot'].set_title(prop+' for '+leaflet)
		
		grid_dict[leaflet]['plot'].set_xlabel("Box-X (nm)")
		grid_dict[leaflet]['plot'].set_ylabel("Box-Y (nm)")
		
	
	CB1 = fig.colorbar(cs[0],cax=axs[2],orientation = 'vertical', boundaries=bounds, ticks=bounds, extend='both', extendfrac='auto')
	#Set the time as title
	currentTime = "t = {0:d}".format(time)
	plt.suptitle(currentTime)
	
	for ax in axs[:-1]:
		ax.autoscale(axis='x', tight=True)
		ax.set_xlim(0.0, box_x)
		ax.set_ylim(0.0, box_y)
	
	
	
	plt.savefig(name)
	plt.close()
	
def plot_mean_grid(mean_list, std_list, prop, graph_name, grid_name, beginning_time, ending_time):
	grid_x = mean_list[0]
	grid_y = mean_list[1]
	grid_dict = mean_list[2]
	bilayer_mean = mean_list[3]
	mean_box = mean_list[4]
	
	#grid_x_std = mean_list[0]
	#grid_y_std = mean_list[1]
	#grid_dict_std = mean_list[2]
	#bilayer_std = mean_list[3]
	
	#should also limit the range
	bounds = None
	levels = None
	if prop == 'ORDER':
		plot_title = r"$\mathsf{P_2 for "
		bounds = [-0.5, 0.0, 0.5, 1.0]
		levels = np.linspace(-0.5, 1, 15)
		colormap = mcm.viridis
	if prop == 'APL':
		plot_title = r"$\mathsf{\frac{Area}{lipid} for "
		bounds = np.linspace(0.55, 0.75, 5)
		levels = np.linspace(0.55, 0.75, 20)
		colormap = mcm.plasma
	if prop == 'THICKNESS':
		plot_title = r"$\mathsf{Thickness for "
		bounds = np.linspace(2.0, 6.0, 5)
		levels = np.linspace(2.0, 6.0, 15)
		colormap = mcm.BuPu_r
	#figure for plots
	fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15,4.9), dpi=72, gridspec_kw = {'width_ratios':[5,5,5,1]})
	
	cs = [None, None, None]
	#sp = [None, None, None]
	
	header = "#Averaged on {0:d} ps from {1:d} to {2:d}\n#".format(ending_time - beginning_time,
																beginning_time, ending_time)
	L_y_set = False
	for leaflet in grid_dict:
		# Setting the axes on which to plot
		Nbleaflet = grid_dict[leaflet]['plot']
		grid_dict[leaflet]['plot'] = axs[ grid_dict[leaflet]['plot'] ]
		# Plot map
		cs[ Nbleaflet ] = grid_dict[leaflet]['plot'].contourf(grid_x, grid_y, grid_dict[leaflet]['grid'], cmap = colormap, levels=levels)
		#sp [ Nbleaflet ] = grid_dict_std[leaflet]['plot'].contourf(grid_x, grid_y, grid_dict[leaflet]['grid'], cmap = colormap, levels=levels)
		name = None
		
		if leaflet == 'upper leaflet':
			name = "upper\ leaflet"
		elif leaflet == "lower leaflet":
			name = "lower\ leaflet"
		
		grid_dict[leaflet]['plot'].set_title( r"$\mathsf{"+ name +r" }$") #prop+' for '+leaflet)
		grid_dict[leaflet]['plot'].set_xlabel(r"$\mathsf{L_x\ (\SI{}{\nano\metre})}$")
		if leaflet == 'lower leaflet':
			grid_dict[leaflet]['plot'].set_ylabel(r"$\mathsf{L_y\ (\SI{}{\nano\metre})}$")
			L_y_set = True
		
		grid_dict[leaflet]['plot'].set_xlim(0.0, mean_box[0])
		grid_dict[leaflet]['plot'].set_ylim(0.0, mean_box[1])
		
		grid_to_write = header + np.array2string(grid_dict[leaflet]['grid'], separator=',')
		with open(grid_name + leaflet.replace(' ','_') + '_mean.grid','w') as grid_file:
			grid_file.write(grid_to_write)
	
	gridx_to_write = header + np.array2string(grid_x, separator=',')
	with open(grid_name + 'mean.gridx','w') as gridx_file:
			gridx_file.write(gridx_to_write)
	
	gridy_to_write = header + np.array2string(grid_y, separator=',')
	with open(grid_name + 'mean.gridy','w') as gridy_file:
			gridy_file.write(gridy_to_write)
	
	gridbi_to_write = header + np.array2string(bilayer_mean, separator=',')
	with open(grid_name + 'bilayer_mean.grid','w') as gridbi_file:
			gridbi_file.write(gridbi_to_write)
	
	# Plot map
	cs[2] = axs[2].contourf(grid_x, grid_y, bilayer_mean, cmap = colormap, levels=levels)
	
	axs[2].set_title( r"$\mathsf{ bilayer }$")
	axs[2].set_xlabel(r"$\mathsf{L_x\ (\SI{}{\nano\metre})}$")
	#axs[2].set_ylabel(r"$\mathsf{L_y\ (\SI{}{\nano\metre}})$")
	
	ticks = [r"$\mathsf{"+str(round(b,2))+"}$" for b in bounds]
	CB1 = fig.colorbar(cs[0],cax=axs[3],orientation = 'vertical', boundaries=bounds, ticks=bounds, extend='both', extendfrac='auto')
	CB1.set_ticklabels(ticks)
	
	
		
	#Set the time as title
	title = "Averaged on {0:d} ps from {1:d} to {2:d}".format(ending_time - beginning_time,
																beginning_time, ending_time)
	plt.suptitle(title)
	
	for ax in axs[:-1]:
		ax.autoscale(axis='x', tight=True)
		ax.set_xlim(0.0, mean_box[0])
		ax.set_ylim(0.0, mean_box[1])
		
		x_ticks = [0.0,5.0,10.0,15.0]
		y_ticks = [0.0,5.0,10.0,15.0]
		
		ax.set_xticks(x_ticks)
		ax.set_yticks(y_ticks)

		y_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(y,1)) for y in y_ticks]
		ax.set_yticklabels(y_ticks_label)

		x_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(x,1)) for x in x_ticks]
		ax.set_xticklabels(x_ticks_label)

	
		
	plt.tight_layout()
	plt.savefig(graph_name)
	plt.close()

def compute_mean_radial_distribution2d(mean_list, radial_increment, prop, graph_name, file_name,
									   beginning_time, ending_time, pad_x, pad_y):
	
	
	
	grid_x = mean_list[0]
	grid_y = mean_list[1]
	grid_dict = mean_list[2]
	bilayer_mean = mean_list[3]
	box_mean = mean_list[4]
	box_x, box_y, box_z = [val for val in box_mean]
	
	#Getting the extrema for x axis
	#replace by box dimensions
	min_valx = 0.0
	max_valx = box_x
	
	#Getting the extrema for y axis
	min_valy = 0.0
	max_valy = box_y
	
	#Setting the center of the averaged box
	center_x = (max_valx - min_valx)/ 2.0
	center_y = (max_valy - min_valy)/ 2.0
	
	#Mean values
	low_leaflet_data = []
	upp_leaflet_data = []
	bil_leaflet_data = []
	
	#Errors
	low_leaflet_err = []
	upp_leaflet_err = []
	bil_leaflet_err = []
	
	radii = []
	
	
	prop_values_low = grid_dict['lower leaflet']['grid']
	prop_values_up = grid_dict['upper leaflet']['grid']
	prop_values_bilayer = bilayer_mean
	
	#unit vectors for box
	x_unit_vec = np.array([box_x, 0.0])
	y_unit_vec = np.array([0.0, box_y])
	#Need to remove the pad_x and pad_y first and last values
	grid_x_cell = crop_grid(grid_x, pad_x, pad_y)
	grid_y_cell = crop_grid(grid_y, pad_x, pad_y)
	#print(grid_x_cell)
	#print(grid_x_cell.shape)
	
	prop_values_low_cell = crop_grid(prop_values_low, pad_x, pad_y)
	prop_values_up_cell = crop_grid(prop_values_up, pad_x, pad_y)
	bilayer_mean_cell = crop_grid(bilayer_mean, pad_x, pad_y)
	
	gridx_from_center = np.subtract(grid_x_cell, center_x)
	gridy_from_center = np.subtract(grid_y_cell, center_y)
	
	vectors_center = np.concatenate( np.stack( (gridx_from_center, gridy_from_center) , axis=-1) )
	
	
	"""
	UL | UU | UR
	LL | CC | RR
	DL | DD | DR
	"""
	
	#along axis
	vectors_up = np.add(vectors_center, y_unit_vec)
	vectors_right = np.add(vectors_center, x_unit_vec)
	vectors_down = np.subtract(vectors_center, y_unit_vec)
	vectors_left = np.subtract(vectors_center, x_unit_vec)
	
	#diagonals
	vectors_up_right = np.add(vectors_up, x_unit_vec)
	vectors_down_right = np.subtract(vectors_right, y_unit_vec)
	vectors_down_left = np.subtract(vectors_left, y_unit_vec)
	vectors_up_left = np.subtract(vectors_up, x_unit_vec)
	
	## DISTANCES
	distances_center = np.sqrt(np.sum(vectors_center**2, axis=1))
	#print(vectors_center)
	#print(distances_center)
	
	"""
	#along axis
	distances_up = np.sqrt(np.sum(vectors_up**2, axis=1))
	distances_right = np.sqrt(np.sum(vectors_right**2, axis=1))
	distances_down = np.sqrt(np.sum(vectors_down**2, axis=1))
	distances_left = np.sqrt(np.sum(vectors_left**2, axis=1))
	
	#diagonals
	distances_up_right = np.sqrt(np.sum(vectors_up_right**2, axis=1))
	distances_down_right = np.sqrt(np.sum(vectors_down_right**2, axis=1))
	distances_down_left = np.sqrt(np.sum(vectors_down_left**2, axis=1))
	distances_up_left = np.sqrt(np.sum(vectors_up_left**2, axis=1))
	"""
	
	#flattened array of values
	prop_values_low = np.array(prop_values_low_cell).flatten()
	prop_values_up = np.array(prop_values_up_cell).flatten()
	prop_values_bil = np.array(bilayer_mean_cell).flatten()
	
	#print(prop_values_bil)
	
	radius = 0.0
	dr = radial_increment
	diagonal = math.sqrt((box_x/2.)**2 + (box_y/2.)**2)
	print(diagonal)
	while( radius <= diagonal ):
		condition = (distances_center >= radius) * (distances_center < radius + dr)
		
		indexes = np.where(condition)
		
		radii.append(radius)
		
		low_leaflet_data.append(np.mean(prop_values_low[indexes]))
		upp_leaflet_data.append(np.mean(prop_values_up[indexes]))
		bil_leaflet_data.append(np.mean(prop_values_bil[indexes]))
		
		low_leaflet_err.append(np.std(prop_values_low[indexes]))	#/ math.sqrt(len(prop_values_low[indexes]))	)
		upp_leaflet_err.append(np.std(prop_values_up[indexes]))		#/ math.sqrt(len(prop_values_up[indexes]))	)
		bil_leaflet_err.append(np.std(prop_values_bil[indexes]))	#/ math.sqrt(len(prop_values_bil[indexes]))	)
		
		radius += dr
	
	"""
	#print(radii)
	#print(low_leaflet_data)
	#print(center_x, center_y)
	#print(box_x, box_y)
	#zipped = zip(radii, low_leaflet_data, upp_leaflet_data, bil_leaflet_data)
	#zipped.sort()
	#radii, low_leaflet_data, upp_leaflet_data, bil_leaflet_data = zip(*sipped)
	#np.sqrt(    np.add(  np.power(np.subtract(grid_x, center_x), 2), np.power(np.subtract(grid_y, center_y), 2)  )    )
	
	#for i in range(0, 50):
		#for j in range(0, 50):
			#R = sqrt( (grid_x[i][j]-center_x)**2 + (grid_y[i][j]-center_y)**2)
			#points.append([])
	
	
	grid_x = mean_list[0]
	grid_y = mean_list[1]
	grid_dict = mean_list[2]
	bilayer_mean = mean_list[3]
	box_mean = mean_list[4]
	
	#Getting the extrema for x axis
	#replace by box dimensions
	min_valx = box_mean[]
	max_valx = np.amax(grid_x)
	
	#Getting the extrema for y axis
	min_valy = np.amin(grid_y)
	max_valy = np.amax(grid_y)
	
	#Setting the center of the averaged box
	center_x = (max_valx - min_valx)/ 2.0
	center_y = (max_valy - min_valy)/ 2.0
	
	radius_x = 0.0
	radius_y = 0.0
	
	low_leaflet_data = []
	upp_leaflet_data = []
	bil_leaflet_data = []
	
	radii = []
	
	# We collect data till we hit the box border
	#divided by 2.0 before
	while radius_x <= max_valx/2.0 and radius_y <= max_valy/2.0:
		value_lower_leaflet = 0.0
		value_upper_leaflet = 0.0
		value_bilayer = 0.0
		
		count = 0
		for i in range(0, 50):
			for j in range(0, 50):
				if(  (grid_x[i][j] >= center_x + radius_x and grid_x[i][j] < center_x + radius_x + radial_increment) and 
				(grid_y[i][j] >= center_y + radius_y and grid_y[i][j] < center_y + radius_y + radial_increment)  ):
					#Checking for nan values
					if math.isnan(grid_dict['lower leaflet']['grid'][i][j]):
						pass
					else:
						value_lower_leaflet += grid_dict['lower leaflet']['grid'][i][j]
						
					if math.isnan(grid_dict['upper leaflet']['grid'][i][j]):
						pass
					else:
						value_upper_leaflet += grid_dict['upper leaflet']['grid'][i][j]
						
					if math.isnan(bilayer_mean[i][j]):
						pass
					else:
						value_bilayer += bilayer_mean[i][j]

					count += 1
					
		if count != 0.0:
			low_leaflet_data.append(value_lower_leaflet/count)
			upp_leaflet_data.append(value_upper_leaflet/count)
			bil_leaflet_data.append(value_bilayer/count)
		else:
			low_leaflet_data.append(0.0)
			upp_leaflet_data.append(0.0)
			bil_leaflet_data.append(0.0)
		
		# storing the radius
		radius2d = math.sqrt(radius_x*radius_x + radius_y*radius_y)
		radii.append(radius2d)
		
		# Increment the radius
		radius_x += radial_increment
		radius_y += radial_increment
		#if <= max_valx and radius_y <= max_valy
	
	"""
	
	low_leaflet_data = np.array(low_leaflet_data)
	upp_leaflet_data = np.array(upp_leaflet_data)
	bil_leaflet_data = np.array(bil_leaflet_data)
	
	low_leaflet_err = np.array(low_leaflet_err)
	upp_leaflet_err = np.array(upp_leaflet_err)
	bil_leaflet_err = np.array(bil_leaflet_err)
		
	#plotting the function for each leaflet and total bilayer
	fig, ax = plt.subplots(nrows=1, ncols=1, dpi=72)
	
	#plot the curves
	ax.plot(radii, low_leaflet_data, alpha=1.0, color='blue', label=r"$\mathsf{lower\ leaflet}$",linewidth=3.0)
	ax.plot(radii, upp_leaflet_data, alpha=1.0, color='red', label=r'$\mathsf{upper\ leaflet}$',linewidth=3.0)
	ax.plot(radii, bil_leaflet_data, alpha=1.0, color='purple', label=r'$\mathsf{bilayer}$', linewidth=3.0)
	
	#plot the error
	ax.fill_between(radii, low_leaflet_data-low_leaflet_err, low_leaflet_data+low_leaflet_err, 
						alpha=0.2, facecolor='blue', linewidth=0)
	
	ax.fill_between(radii, upp_leaflet_data-upp_leaflet_err, upp_leaflet_data+upp_leaflet_err,
						alpha=0.2, facecolor='red', linewidth=0)
	
	ax.fill_between(radii, bil_leaflet_data-bil_leaflet_err, bil_leaflet_data+bil_leaflet_err,
						alpha=0.2, facecolor='purple', linewidth=0)
	
	ax.set_xlabel(r"$\mathsf{radius\ (\SI{}{\nano\metre}})$")
	
	y_name = None
	if prop == "ORDER":
		y_name = r"$\mathsf{P_2}$"
		ax.set_ylabel(y_name,rotation=0)
	elif prop == "APL":
		y_name = r"$\mathsf{ \frac{\displaystyle Area }{\displaystyle lipid }\ (\SI{}{\nano\metre^2}) } $"
		ax.set_ylabel(y_name)
	elif prop == "THICKNESS":
		y_name = r"$\mathsf{Thickness\ (\SI{}{\nano\meter})}$"
		ax.set_ylabel(y_name)
		
	ax.grid('on')
	ax.legend()
	
	#title = "Radial distribution of {0} \n averaged on {1:d} ps from {2:d} to {3:d} ps".format(prop, ending_time - beginning_time, beginning_time, ending_time)
	#plt.suptitle(title)
	
	x_ticks = ax.get_xticks()
	y_ticks = ax.get_yticks()

	
	y_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(y,1)) for y in y_ticks]
	ax.set_yticklabels(y_ticks_label)

	x_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(x,1)) for x in x_ticks]
	ax.set_xticklabels(x_ticks_label)
	
	plt.tight_layout()
	#plt.show()
	plt.savefig(graph_name)
	plt.close()
	
	header="radius,lower leaflet,upper leaflet,bilayer"
	data = np.array([radii, low_leaflet_data, upp_leaflet_data, bil_leaflet_data])
	data = data.T
	np.savetxt(file_name, data, header=header ,delimiter=",")
	
	
	
def histogram(dataframe, prop, csv_name, svg_name, index, last_histo, time, axs):
	hist_range = None
	if prop == 'ORDER':
		hist_range = (-0.5, 1.0)
	if prop == 'APL':
		hist_range = (0.4, 0.8)
	if prop == 'THICKNESS':
		hist_range = (0.0, 6.0)
		
	fig_hist = plt.figure(figsize=(15,10))
	
	upper_histo = plt.subplot(221)
	lower_histo = plt.subplot(223)
	bilayer_histo = plt.subplot(122)
	
	property_name = dataframe.columns[-1]
	
	low_prop = dataframe[dataframe['leaflet'] == 'lower leaflet'][property_name]
	upp_prop = dataframe[dataframe['leaflet'] == 'upper leaflet'][property_name]
	bil_prop = dataframe[property_name]
	
	hist_prop_low, low_bins, patches = lower_histo.hist(low_prop, bins='auto', range=hist_range,
														facecolor='blue')
	lower_histo.set_title(prop+' for '+'lower leaflet')
	lower_histo.set_xlabel(prop)
	lower_histo.set_xlabel('# lipids')
	lower_histo.grid('on')
	
	hist_prop_upp, upp_bins, patches = upper_histo.hist(upp_prop, bins='auto', range=hist_range,
														facecolor='red')
	upper_histo.set_title(prop+' for '+'upper leaflet')
	upper_histo.set_xlabel(prop)
	upper_histo.set_xlabel('# lipids')
	upper_histo.grid('on')
	
	hist_prop_bil, bil_bins, patches = bilayer_histo.hist(bil_prop, bins='auto', 
														range=hist_range, facecolor='purple')
	bilayer_histo.set_title(prop+' in bilayer')
	bilayer_histo.set_xlabel(prop)
	bilayer_histo.set_ylabel('# lipids')
	bilayer_histo.grid('on')
	
	title = "Histogram of {0} \n at t= {1:d} ps".format(prop, time)
	plt.suptitle(title)

	lower_histo.autoscale(axis='x', tight=True)
	upper_histo.autoscale(axis='x', tight=True)
	bilayer_histo.autoscale(axis='x', tight=True)
	
	plt.savefig(svg_name)
	plt.close(fig_hist)
	
	data = np.asarray([low_bins, hist_prop_low, upp_bins, hist_prop_upp, bil_bins, hist_prop_bil])
	data = data.T
	
	header = str("lower leaflet bins, lower leaflet nb,upper leaflet bins, upper leaflet nb,"
					" bilayer bins, bilayer nb, #AT TIME {0:d} ps").format(time)
	#with open(csv_name,'wb') as csv_file:
	np.savetxt(csv_name, data, header=header, delimiter=",", fmt="%s")
	
	if index == 0 or index == last_histo:
		axs[0].hist(low_prop, bins='auto',range=hist_range, facecolor='blue', alpha=0.75)
		axs[1].hist(upp_prop, bins='auto',range=hist_range, facecolor='red', alpha=0.75)
		axs[2].hist(bil_prop, bins='auto',range=hist_range, facecolor='purple', alpha=0.75)
	else:
		axs[0].hist(low_prop, bins='auto',range=hist_range, facecolor='blue', alpha=0.1)
		axs[1].hist(upp_prop, bins='auto',range=hist_range, facecolor='red', alpha=0.1)
		axs[2].hist(bil_prop, bins='auto',range=hist_range, facecolor='purple', alpha=0.1)




def plotting_histo_trace(fig, axs, name, beginning_time, ending_time):
	axs[0].set_title(prop+' for '+'lower leaflet')
	axs[0].set_xlabel(prop)
	axs[0].set_xlabel('# lipids')
	axs[0].grid('on')
	axs[0].autoscale(axis='x', tight=True)
	
	axs[1].set_title(prop+' for '+'upper leaflet')
	axs[1].set_xlabel(prop)
	axs[1].set_xlabel('# lipids')
	axs[1].grid('on')
	axs[1].autoscale(axis='x', tight=True)
	
	axs[2].set_title(prop+' in bilayer')
	axs[2].set_xlabel(prop)
	axs[2].set_ylabel('# lipids')
	axs[2].grid('on')
	axs[2].autoscale(axis='x', tight=True)
	
	title = "histograme trace on {0:d} ps from {1:d} to {2:d}".format(ending_time - beginning_time,
															beginning_time, ending_time)
	plt.suptitle(title)

	plt.savefig(name)
	plt.close(fig)
		
		
		
		
		
		
		
		
			
def PDB_out(dataframe, name):
	pdb_file = name
	HG_PDB = 'HEAD'
	
	property_name = dataframe.columns[-1]
	
	with open(pdb_file,'a') as PDB:
		for resid, leaflet, x, y, z, val in zip(dataframe['resid'], dataframe['leaflet'],dataframe['X coords'], dataframe['Y coords'],
										dataframe['Z coords'], dataframe[property_name]):

			x *= 10.
			y *= 10.
			z *= 10.
			leaflet = leaflet.upper()[:3]
			line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(
			"ATOM", resid, HG_PDB, ' ', leaflet, ' ', 0, ' ', x, y, z, 1.0, val, ' ',' ')

			PDB.write(line)
		
def scatter_plot(dataframe_x, dataframe_y, name):
	"""
	Makes a scatter plot.
	
	dataframe_x : properties to plot on x axis
	dataframe y : properties to plot on y axis
	
	"""
	fig, ax = plt.subplots()
	
	property_x_name = dataframe_x.columns[-1]
	property_y_name = dataframe_y.columns[-1]
	
	property_x = dataframe_x[property_x_name]
	property_y = dataframe_y[property_y_name]
	
	if property_x.shape[0] < property_y.shape[0]:
		property_x = np.pad(property_x, (0, (property_y.shape[0] - property_x.shape[0])) ,mode='constant', constant_values=np.nan)
	if property_y.shape[0] < property_x.shape[0]:
		property_y = np.pad(property_y, (0, (property_x.shape[0] - property_y.shape[0])) ,mode='constant', constant_values=np.nan)

	pearson_coeff = Pearson(property_x, property_y)
	spearman_coeff = Spearman(property_x, property_y)
	
	coefficients = "p = {0:.1}\ns = {1:.1}".format(pearson_coeff,
															spearman_coeff)
	
	color = None
	colormap = None
	if 'order' in property_y_name.lower():
		color = 'blue'
		colormap = mcm.Blues
	if 'area' in property_y_name.lower():
		color = 'red'
		colormap = mcm.Reds
	if 'thickness' in property_y_name:
		color = 'green'
		colormap = mcm.Greens
		
	scatter = plt.scatter(property_x, property_y, alpha=0.3, label=coefficients)
	ax.set_xlabel(property_x_name)
	ax.set_ylabel(property_y_name)
	ax.grid('on')
	ax.legend()
	plt.savefig(name+'.svg')
	
	plt.close()
	
	fig, ax = plt.subplots()
	
	hexbin = plt.hexbin(property_x, property_y, cmap=colormap,label=coefficients)
	ax.set_xlabel(property_x_name)
	ax.set_ylabel(property_y_name)
	ax.grid('on')
	ax.legend()
	plt.savefig(name+'_hex.svg')
	
	plt.close()
	
	
def Var(xs, mu=None, ddof=0):
	xs = np.asarray(xs)
	
	if mu is None:
		mu = xs.mean()
	
	ds = xs - mu
	return np.dot(ds, ds) / (len(xs) - ddof)

def MeanVar(xs, ddof=0):
	"""
	Computes mean and variance
	"""
	xs = np.asarray(xs)
	mean = xs.mean()
	s2 = Var(xs, mean, ddof)
	return mean, s2

def Cov(xs,ys, meanx=None, meany=None):
	xs = np.asarray(xs)
	ys = np.asarray(ys)
	
	if meanx is None:
		meanx = np.mean(xs)
	if meany is None:
		meany = np.mean(ys)
	
	cov = np.dot(xs-meanx, ys-meany)/ len(xs)
	return cov

def Pearson(property_x, property_y):
	xs = np.asarray(property_x)
	ys = np.asarray(property_y)
	
	meanx, varx = MeanVar(xs)
	meany, vary = MeanVar(ys)
	
	corr = Cov(xs, ys, meanx, meany)/math.sqrt(varx * vary)
	return corr

def Spearman(xs, ys):
	xs = pd.Series(xs)
	ys = pd.Series(ys)
	
	xranks = xs.rank()
	yranks = ys.rank()
	return xs.corr(ys, method='spearman')

def transform_to_unit_box(grid, bdim):
	box_x, box_y, box_z = [float(val) for val in bdim[1:]]
	
	grid[0] = grid[0]/box_x
	grid[1] = grid[1]/box_y
	
#def norm_vec2d(vector):
	#return sqrt(vector[0]**2 + vector[1]**2)
	
	
	
def scatter_cmd(csv_found):
	data_scatter = {}
	#Scatter plot
	for job_name in data:
		start_time = time.time()
		data_scatter[job_name] = {}
		for prop in data[job_name]:
			frames = data[job_name][prop]
			data_scatter[job_name][prop] = pd.concat(frames)
		print( "--- gathering scatter data for {0} done in {1:f} seconds ---".format(job_name ,time.time() - start_time) )
			
	for job_name in data_scatter:
		start_time = time.time()
		couple = []
		
		for prop_x in data_scatter[job_name]:
			for prop_y in data_scatter[job_name]:
				if [prop_y, prop_x] in couple or prop_x == prop_y:
					continue
				name = "{0}/{1}/graphs/scatter_X-{2}_Y-{3}".format(analysis_folder, job_name, prop_x, prop_y)
				
				X = data_scatter[job_name][prop_x]
				Y = data_scatter[job_name][prop_y]
				
				scatter_plot(X, Y, name)
				couple.append([prop_x, prop_y])
		
		print( "--- Scatter for {0} done in {1:f} seconds ---".format(job_name ,time.time() - start_time) )
		gc.collect()
			
	#for job_name in csv_found:
		#for prop in csv_found[job_name]:
			#property_values = data[job_name][prop].columns[-1]
			#x_axes = [val for val in data[job_name][prop].columns[:-1] if val !='leaflet']
			#for x_axis in x_axes:
				#data[job_name][prop].plot.scatter(x=x_axis,y=property_values)
	#plt.show()

def pdb_cmd(data):
	#Creates pdb file with headgroups
	nb_frames = None
	
	for job_name in data:
		start_time = time.time()
		for prop in data[job_name]:
			pdb_folder = "{0}/{1}/PDB/{2}/".format(analysis_folder, job_name, prop)
			
			#list previous files and remove them
			files = glob.glob(pdb_folder+'*.pdb')
			for f in files:
				os.remove(f)
			
			nb_frames = len(str( len(data[job_name][prop]) ))

			for index, frame in enumerate(data[job_name][prop]):
				name = "{0}/{1}_frame_{2}.pdb".format(pdb_folder, prop, str(index).zfill(nb_frames))
				PDB_out(frame, name)
		print( "--- PDB for {0} done in {1:f} seconds ---".format(job_name ,time.time() - start_time) )
		gc.collect()





def map_cmd(data):
	#values for padding
	pad_x = pad_y = 10
	for job_name in data:
		start_time = time.time()
		box_file_path = "{0}/{1}/*BOX.xvg".format(analysis_folder, job_name)
		box_file = glob.glob(box_file_path)[0]
		
		box_dim = []
		nb_frames = None
		with open(box_file,'r') as boxF:
			for line in boxF:
				if line.startswith('@') or line.startswith('#'):
					continue
				else:
					box_dim.append(line.split()[:4])
		
		box_dim = box_dim[BEGIN:END+1:STRIDE]
		#convert to float
		for index, box in enumerate(box_dim):
			box_dim[index] = [float(i) for i in box]
		
		beginning_time = int(box_dim[0][0])
		ending_time = int(box_dim[-1][0])
		
		#compute the mean box
		box_dimension = np.array([i[1:] for i in box_dim])
		mean_box = np.mean(box_dimension, axis=0)
		
		for prop in data[job_name]:
			
			#list containing grid data
			grid_list = []
			grid_folder = "{0}/{1}/GRID_data/{2}".format(analysis_folder, job_name, prop)
			grid_graph_folder = "{0}/{1}/GRID_data/graphs/{2}".format(analysis_folder, job_name, prop)
			
			if not os.path.isdir(grid_folder):
				os.makedirs(grid_folder, exist_ok=True)
			else:
				files = glob.glob(grid_folder+'/*.grid*')
				for f in files:
					os.remove(f)
			
			if not os.path.isdir(grid_graph_folder):
				os.makedirs(grid_graph_folder, exist_ok=True)
			else:
				files = glob.glob(grid_graph_folder+'/*.svg')
				for f in files:
					os.remove(f)
			
			nb_frames = len(str( len(box_dim) ))
			box_and_frame = zip(box_dim, data[job_name][prop])
			
			for index, bandf in enumerate(box_and_frame):
				bdim = bandf[0]
				frame = bandf[1]
				
				name = "{0}/{1}_frame_{2}_".format(grid_folder, prop, str(index).zfill(nb_frames))
				
				#Output the interpolation
				grid_data = interpolate_grid(frame, bdim, name, pad_x, pad_y)
				
				#appending data to the grid data list to compute mean later
				grid_list.append(grid_data)
			
			#plot the interpolations every STRIDE
			#nb_plot_frames = len(box_dim)
			
			grid_box = zip(grid_list, box_dim)
			for index, g_b in enumerate(grid_box):
				grid = g_b[0]
				bdim = g_b[1]
				name = "{0}/{1}_frame_{2}.svg".format(grid_graph_folder, prop, str(index).zfill(nb_frames))
				#plot_grid(grid, prop, name, bdim)
			
			##computing the mean
			# mean grid_x
			grid_x_array = np.array([i[0] for i in grid_list])
			grid_x_mean = np.mean(grid_x_array, axis=0)
			#grid_x_std = np.sqrt(np.var(grid_x_array, axis=0))
			
			# mean grid_y
			grid_y_array = np.array([i[1] for i in grid_list])
			grid_y_mean = np.mean(grid_y_array, axis=0)
			#grid_y_std = np.sqrt(np.var(grid_y_array, axis=0))
			
			#mean property for lower leaflet
			grid_lower_leaflet_array = np.array([i[2]['lower leaflet']['grid'] for i in grid_list])
			grid_lower_leaflet_mean = np.mean(grid_lower_leaflet_array, axis=0)
			#grid_lower_leaflet_var = np.var(grid_lower_leaflet_array, axis=0)
			#grid_lower_leaflet_std = np.sqrt(grid_lower_leaflet_var)
			
			#mean property for upper leaflet
			grid_upper_leaflet_array = np.array([i[2]['upper leaflet']['grid'] for i in grid_list])
			grid_upper_leaflet_mean = np.mean(grid_upper_leaflet_array, axis=0)
			#grid_upper_leaflet_var = np.var(grid_upper_leaflet_array, axis=0)
			#grid_upper_leaflet_std = np.sqrt(grid_upper_leaflet_var)
			
			
			#mean property for bilayer
			grid_bilayer_mean = np.add(grid_lower_leaflet_mean, grid_upper_leaflet_mean)/2.0
			#grid_bilayer_std = np.sqrt( (grid_lower_leaflet_var + grid_upper_leaflet_var)/2.0 )
			
			
			graph_name = "{0}/{1}_mean.eps".format(grid_graph_folder, prop)
			grid_name = "{0}/{1}_mean_".format(grid_folder, prop)
			
			"""
			grid_x_mean = np.zeros(shape=(50,50), dtype=float)
			grid_y_mean = np.zeros(shape=(50,50), dtype=float)
			grid_lower_leaflet_mean = np.zeros(shape=(50,50), dtype=float)
			grid_upper_leaflet_mean = np.zeros(shape=(50,50), dtype=float)
			
			grid_x_array = np.array(grid[])
			grid_x_mean = np.mean()
			for i in range(0, 50):
				for j in range(0, 50):
					for grid in grid_list:
						grid_x_mean[i][j] += grid[0][i][j]
						grid_y_mean[i][j] += grid[1][i][j]
						grid_lower_leaflet_mean[i][j] += grid[2]['lower leaflet']['grid'][i][j]
						grid_upper_leaflet_mean[i][j] += grid[2]['upper leaflet']['grid'][i][j]
						
						
						
						
			nbgrid = float(len(grid_list))
			grid_x_mean = grid_x_mean/nbgrid
			grid_y_mean = grid_y_mean/nbgrid
			grid_lower_leaflet_mean = grid_lower_leaflet_mean/nbgrid
			grid_upper_leaflet_mean = grid_upper_leaflet_mean/nbgrid
			"""

			mean_list = [grid_x_mean, grid_y_mean, 
						{'lower leaflet' : {'grid':grid_lower_leaflet_mean, 'plot':0},
						'upper leaflet': {'grid':grid_upper_leaflet_mean, 'plot':1}}
						, grid_bilayer_mean, mean_box]
						
			std_list = 0
			"""
			[grid_x_std, grid_y_std, 
						{'lower leaflet' : {'grid':grid_lower_leaflet_std, 'plot':0},
						'upper leaflet': {'grid':grid_upper_leaflet_std, 'plot':1}}
						, grid_bilayer_std]
			"""
			plot_mean_grid(mean_list, std_list, prop, graph_name, grid_name, beginning_time, ending_time)
			
			if RADIAL is not None:
				graph_name = "{0}/{1}_radial.eps".format(grid_graph_folder, prop)
				file_name = "{0}/{1}.rad".format(grid_folder, prop)
				compute_mean_radial_distribution2d(mean_list, RADIAL, prop, graph_name, file_name,
													beginning_time, ending_time, pad_x, pad_y)
				
		print( "--- Maps for {0} done in {1:f} seconds ---".format(job_name ,time.time() - start_time) )
		gc.collect()
















def hist_cmd(data):
	for job_name in data:
		start_time = time.time()
		box_file_path = "{0}/{1}/*BOX.xvg".format(analysis_folder, job_name)
		box_file = glob.glob(box_file_path)[0]
		
		box_dim = []
		nb_frames = None
		with open(box_file,'r') as boxF:
			for line in boxF:
				if line.startswith('@') or line.startswith('#'):
					continue
				else:
					box_dim.append(line.split()[:4])
		
		box_dim = box_dim[BEGIN:END+1:STRIDE]
		
		beginning_time = int(box_dim[0][0])
		ending_time = int(box_dim[-1][0])
		
		for prop in data[job_name]:
			#list containing grid data
			hist_list = []
			hist_folder = "{0}/{1}/HISTO/{2}".format(analysis_folder, job_name, prop)
			hist_graph_folder = "{0}/{1}/HISTO/graphs/{2}".format(analysis_folder, job_name, prop)
			
			if not os.path.isdir(hist_folder):
				os.makedirs(hist_folder, exist_ok=True)
			else:
				files = glob.glob(hist_folder+'/*.csv*')
				for f in files:
					os.remove(f)
			
			if not os.path.isdir(hist_graph_folder):
				os.makedirs(hist_graph_folder, exist_ok=True)
			else:
				files = glob.glob(hist_graph_folder+'/*.svg')
				for f in files:
					os.remove(f)
			
			nb_frames = len(str( len(box_dim) ))
			
			#plot the interpolations every STRIDE
			fig_trace = plt.figure(figsize=(15,10))
			
			lower_histo_trace = plt.subplot(223)
			upper_histo_trace = plt.subplot(221)
			bilayer_histo_trace = plt.subplot(122)
			
			axs = [lower_histo_trace, upper_histo_trace, bilayer_histo_trace]
			
			last_histo = len(hist_list) - 1
			
			box_and_frame = zip(box_dim, data[job_name][prop])
			for index, b_a_f in enumerate(box_and_frame):
				bdim = b_a_f[0]
				frame = b_a_f[1]
				
				csv_name = "{0}/{1}_frame_{2}.csv".format(hist_folder, prop, str(index).zfill(nb_frames))
				
				svg_name = "{0}/{1}_frame_{2}.svg".format(hist_graph_folder, prop, str(index).zfill(nb_frames))
				histogram(frame, prop, csv_name, svg_name, index, last_histo, int(bdim[0]), axs)
				
			name = "{0}/{1}_trace.svg".format(hist_graph_folder, prop)
			plotting_histo_trace(fig_trace, axs, name, beginning_time, ending_time)
			
		print( "--- Histograms for {0} done in {1:f} seconds ---".format(job_name ,time.time() - start_time) )
		gc.collect()
		
















#=======================##=======================##=======================#
#=======================##THE SCRIPT STARTS HERE ##=======================#
#=======================##=======================##=======================#

default_properties = ['APL','ORDER','THICKNESS']

# Command parser #
parser = argparse.ArgumentParser(description='Compute the desired properties with fatslim '
								 'and plot them.\n'
								 'Use -r and -f for automatic file finding in the specified '
								 'folder.\n'
								 'Else if you want to compute properties on a peculiar file '
								 'use -c')


sub_parser = parser.add_subparsers()

# Adding parser options ##############################################################
compute = sub_parser.add_parser('compute')

genOpt = compute.add_argument_group("General options")
genOpt.add_argument('-g','--gmx', dest='gmx_path',
					type=str, default='gmx',
                    help='Set the gromacs path')

genOpt.add_argument('-i', '--input', dest='md_run_input',
					type=str, nargs='*', default=None,
                    help='Set the name of MD run to analyse')

genOpt.add_argument('-f', '--folder', dest='folder',
					type=str, default=None,
                    help='Set the folder where to search for data')

genOpt.add_argument('--fatslim', dest='fatslim',
					type=str, default='fatslim',
                    help='Set the path to FATSLIM')

# parameters for fatslim ##############################################################
fatslimOpt = compute.add_argument_group("Options for fatslim")

#option for headgroup
fatslimOpt.add_argument('--hg-group', dest='head_group',
					type=str, default="PO4",
                    help='Set the headgroup for fatslim')

#option for configuration file
fatslimOpt.add_argument('-c', '--conf', dest='conf',
					type=str, default=None,
                    help='Set the name of the configuration file to use')

#option for trajectory file
fatslimOpt.add_argument('-t', '--trajectory', dest='trajectory',
					type=str, nargs='*', default=None,
                    help='Set the name of MD run for trajectory')

#option for begin-frame
fatslimOpt.add_argument('-b', '--begin', dest='begin_frame',
					type=int, default=None,
                    help='Set the first frame for the trajectories')

#option for end-frame
fatslimOpt.add_argument('-e', '--end', dest='end_frame',
					type=int, default=None,
                    help='Set the last frame for the trajectories')

#option for index file
fatslimOpt.add_argument('-n', '--index', dest='index',
					type=str, default='index.ndx',
                    help='Set the name of MD run to analyse')

# cutoff options
fatslimOpt.add_argument('--cutoff', dest='cutoff',
					type=float, default=6.0,
                    help='Set the global cutoff for fatslim')

fatslimOpt.add_argument('--cutapl', dest='cutapl',
					type=float, default=3.0,
                    help='Set the apl cutoff for fatslim')

fatslimOpt.add_argument('--cutthick', dest='cutthick',
					type=float, default=6.0,
                    help='Set the thickness cutoff for fatslim')

# axis for order parameter
fatslimOpt.add_argument('--main-axis', dest='main_axis',
					type=int, nargs=3, default=None,
                    help='Set the membrane cutoff for fatslim')

# option for frequency
fatslimOpt.add_argument('--idfreq', dest='idfreq',
					type=int, default=None,
                    help='Set the frequency for fatslim')

## parameters for analysis ##############################################################
analyse = sub_parser.add_parser('analyse')

genOptAnalysis = analyse.add_argument_group("General options")

genOptAnalysis.add_argument('-f', '--folder', dest='folder',
					type=str, default=None,
                    help='Set the folder where to search for data')

genOptAnalysis.add_argument('-n', '--name', dest='name',
					type=str, nargs='*', default=None,
                    help='Name of the job to analyse (if specific)')

analysisOpt = analyse.add_argument_group("Options for analysis")

analysisOpt.add_argument('--hist', action='store_true',
                    help='Plot histograms from .csv files')

analysisOpt.add_argument('--map', action='store_true',
                    help='Plot 2D maps from .csv files')

analysisOpt.add_argument('--radial','-r', dest='radial',
					type=float, default= None,
                    help='Set the infinitesimal length dr [nm] for position-based histogram.')

analysisOpt.add_argument('--scatter', action='store_true',
                    help='Plot scatter graph with apl, order and thickness')

analysisOpt.add_argument('--prop', dest='properties',
					type=str, nargs='*', default= default_properties,
                    help='Set the name of the properties to analyse. (default = False)')

analysisOpt.add_argument('--pdb', action='store_true',
                    help='Write out pdb from .csv files. (default = False)')

analysisOpt.add_argument('--stride', dest='stride',
					type=int, default=None,
                    help='Set the stride for selecting .csv files. (default=1).')

# stride for each mode
analysisOpt.add_argument('-sh', dest='hist_stride',
					type=int, default=None,
                    help='Set the stride for histogram. (default= stride).')

analysisOpt.add_argument('-sm', dest='map_stride',
					type=int, default=None,
                    help='Set the stride for maps. (default= stride).')

analysisOpt.add_argument('-ss', dest='scatter_stride',
					type=int, default=None,
                    help='Set the stride for scatter. (default= stride).')

# analysis timespan
analysisOpt.add_argument('--begin', dest='begin',
					type=int, default=None,
                    help='Set the beginning for .csv files. (default= first).')

analysisOpt.add_argument('--end', dest='end',
					type=int, default=None,
                    help='Set the end for selecting .csv files. Applies to --map and '
                    '--hist (default= last).')

#begin and end for each mode
analysisOpt.add_argument('-bh', dest='hist_begin',
					type=int, default=None,
                    help='Set the beginning for histogram. (default= begin).')

analysisOpt.add_argument('-bm', dest='map_begin',
					type=int, default=None,
                    help='Set the beginning for maps. (default= begin).')

analysisOpt.add_argument('-bs', dest='scatter_begin',
					type=int, default=None,
                    help='Set the beginning for scatter. (default= begin).')

analysisOpt.add_argument('-eh', dest='hist_end',
					type=int, default=None,
                    help='Set the end for histogram. (default= end).')

analysisOpt.add_argument('-em', dest='map_end',
					type=int, default=None,
                    help='Set the end for maps. (default= end).')

analysisOpt.add_argument('-es', dest='scatter_end',
					type=int, default=None,
                    help='Set the end for scatter. (default= end).')


## parameters for animate ##############################################################
animate = sub_parser.add_parser('animate')

genOptAnimate = animate.add_argument_group("General options")

genOptAnimate.add_argument('-f', '--folder', dest='folder',
					type=str, required=True,
                    help='Set the folder where to search for data')

genOptAnimate.add_argument('-t', '--type', dest='type',
					type=str, default=None,
                    help='Set the image type that will be used')

## parameters for mayavi #
mayavi = sub_parser.add_parser('mayavi')

genOptMayavi = mayavi.add_argument_group("General options")

genOptMayavi.add_argument('-f', '--filename', dest='filename',
					type=str, required=True,
                    help='Set the filename to plot using mayavi')

genOptMayavi.add_argument('-t', '--type', dest='type',
					type=str, required=True,
                    help='Set the type to plot using mayavi: loctemp or fatslim')

## parameters for mayavi #
visual = sub_parser.add_parser('visual')

genOptVisual = visual.add_argument_group("General options")

genOptVisual.add_argument('-f', '--filename', dest='filename',
				type=str, required=True,
				help='Set the filename to visualise with vpython')


# Adding parser options to gromacs analysis##############################################################
gromacs = sub_parser.add_parser('gromacs')

gmxOpt = gromacs.add_argument_group("General options")
gmxOpt.add_argument('-g','--gmx', dest='gmx_path',
					type=str, default='gmx',
                    help='Set the gromacs path')

gmxOpt.add_argument('-i', '--input', dest='md_run_input',
					type=str, nargs='*', default=None,
                    help='Set the name of MD run to analyse')

gmxOpt.add_argument('--ext', dest='extension',
					type=str, default='xtc',
                    help='Set the name of the extension for trajecotry : gro , xtc (default = xtc)')

gmxOpt.add_argument('-f', '--folder', dest='folder',
					type=str, default=None,
                    help='Set the folder where to search for data')

#option for begin-frame
gmxOpt.add_argument('-b', '--begin', dest='begin_frame',
					type=int, default=None,
                    help='Set the first frame for the trajectories')

#option for end-frame
gmxOpt.add_argument('-e', '--end', dest='end_frame',
					type=int, default=None,
                    help='Set the last frame for the trajectories')

gmxOpt.add_argument('-dt', dest='dt',
					type=int, default=None,
                    help='Set the stride for density computation')

gmxOpt.add_argument('--radius', dest='radius',
					type=float, default=8.0,
                    help='Set the radius outside which the membrane will be analysed')

gmxOpt.add_argument('--density', action='store_true',
                    help='Compute density using gmx density')

gmxOpt.add_argument('--energy', action='store_true',
                    help='Compute energy using gmx energy')

gmxOpt.add_argument('--loctemp', action='store_true',
                    help='Compute local temperature in the box')

gmxOpt.add_argument('--locpres', action='store_true',
                    help='Compute local pressure in the box')


gmxOpt.add_argument('--select', action='store_true',
                    help='Compute solvent per lipid')

gmxOpt.add_argument('--headgroup', dest='headgroup',
					type=str, nargs='*', default=None,
                    help='Set the name of headgroup')

gmxOpt.add_argument('--linkgroup', dest='linkgroup',
					type=str, nargs='*', default=None,
                    help='Set the name of linkgroup')

gmxOpt.add_argument('--tailgroup', dest='tailgroup',
					type=str, nargs='*', default=None,
                    help='Set the name of tailgroup')


# Adding parser options to gromacs analysis##############################################################
mda = sub_parser.add_parser('mda')

mdaOpt = mda.add_argument_group("General options")

mdaOpt.add_argument('-i', '--input', dest='md_run_input',
					type=str, nargs='*', default=None,
                    help='Set the name of MD run to analyse')

mdaOpt.add_argument('--ext', dest='extension',
					type=str, default='xtc',
                    help='Set the name of the extension for trajecotry : gro , xtc (default = xtc)')

mdaOpt.add_argument('-f', '--folder', dest='folder',
					type=str, default=None,
                    help='Set the folder where to search for data')

#option for begin-frame
mdaOpt.add_argument('-b', '--begin', dest='begin_frame',
					type=int, default=None,
                    help='Set the first frame for the trajectories')

#option for end-frame
mdaOpt.add_argument('-e', '--end', dest='end_frame',
					type=int, default=None,
                    help='Set the last frame for the trajectories')

mdaOpt.add_argument('-dt', dest='stride',
					type=int, default=None,
                    help='Set the stride (in outputfrequency * timestep)')

mdaOpt.add_argument('--radius', dest='radius',
					type=float, default=None,
                    help='Set the radius outside which the membrane will be analysed')

mdaOpt.add_argument('--density', dest='density',
					type=float, default=None,
                    help='Set the bin size for computing density')

mdaOpt.add_argument('--grid', dest='grid',
					type=float, default=None,
                    help='Set the bin size for computing grid')

mdaOpt.add_argument('--local', dest='local',
					type=float, nargs=3, default=None,
                    help='Set the bin size for computing local properties')

mdaOpt.add_argument('--select', action='store_true',
                    help='Compute solvent per lipid')

mdaOpt.add_argument('--headgroup', dest='headgroup',
					type=str, nargs='*', default=None,
                    help='Set the name of headgroup')

mdaOpt.add_argument('--tailgroup', dest='tailgroup',
					type=str, nargs='*', default=None,
                    help='Set the name of tailgroup')

mdaOpt.add_argument('--sugroup', dest='sugroup',
					type=str, nargs='*', default=None,
                    help='Set the name of sugroup')


# Adding parser options to gromacs analysis##############################################################
xvgplot = sub_parser.add_parser('xvgplot')

xvgOpt = xvgplot.add_argument_group("General options")

xvgOpt.add_argument('-q', '--quantities', dest='quantities',
					type=str, nargs='*', default=None,
                    help='Set the name of the quantities to plot')

xvgOpt.add_argument('-f', '--file', dest='xvg_file',
					type=str, required=True,
                    help='Set the name of the file to plot')

xvgOpt.add_argument('--integrate', dest='limits',
					type=float, nargs=2, default=None,
                    help='Set the limits for integration')

xvgOpt.add_argument('--plot', action='store_true',
                    help='Plot the selected quantities')

xvgOpt.add_argument('--mean', action='store_true',
                    help='Return the mean, std and stderr values of the selected quantities')

# Adding parser options to reflectometry analysis##############################################################
reflectometry = sub_parser.add_parser('reflectometry')

reflectometryOpt = reflectometry.add_argument_group("General options")

reflectometryOpt.add_argument('-sl', dest='sl_file',
					type=str, required=True,
					help='Set the name of the file to use for scattering length data (.sl)')

reflectometryOpt.add_argument('-ref', dest='reflectivity',
					type=str, default=None,
					help='Set the quantity name to compute reflectivity')

reflectometryOpt.add_argument('-f', '--file', dest='xvg_file',
					type=str, required=True,
					help='Set the name of the file to use for sld (density.xvg)')

reflectometryOpt.add_argument('-m','--molecules', dest='molecules',
					type=str, nargs='*', required=True,
					help='Set the molecules used to compute sld')

reflectometryOpt.add_argument('--limits', dest='limits',
					type=float, nargs=2, default=None,
					help='Set the limits for the sld')

reflectometryOpt.add_argument('-sub','--sublayers', dest='sublayers',
					type=str, nargs='*', default=None,
                    help='Set the sub layers composition to compute the reflectivity [[mol density thickness rugosity], ...]')

reflectometryOpt.add_argument('-sup','--superlayers', dest='superlayers',
					type=str, nargs='*', default=None,
                    help='Set the super layers composition to compute the reflectivity [[mol density thickness rugosity], ...]')

reflectometryOpt.add_argument('-q', dest='qrange',
					type=float, nargs=3, default=None,
                    help='Set the Qmin, Qmax and Qgrid size to compute reflectometry curve')

reflectometryOpt.add_argument('--check', action='store_true',
                    help='Show the sld curve before computing reflectometry')

reflectometryOpt.add_argument('-res', dest='resolution',
					type=int, default=None,
                    help='Set the number of points for gaussian convolution')

reflectometryOpt.add_argument('-exp', dest='experimental',
					type=str, default=None,
                    help='File with data from experimental neutron reflectometry')

#dlambda/lambda = 0.1 for Koutsioubas

cmdParam = parser.parse_args(sys.argv[1:])
print(cmdParam)


if 'compute' in sys.argv:
	start_time = time.time()
	# Setting general options
	GMX = cmdParam.gmx_path
	MD_RUN_INPUT = cmdParam.md_run_input
	FOLDER = cmdParam.folder
	FATSLIM = cmdParam.fatslim

	# Setting parameters for fatslim
	HEAD_GROUP = cmdParam.head_group
	CONF_FILE = cmdParam.conf
	TRAJECTORY = cmdParam.trajectory
	BEGIN_FRAME = cmdParam.begin_frame
	END_FRAME = cmdParam.end_frame
	
	INDEX = cmdParam.index
	IDFREQ = cmdParam.idfreq
	CUTOFF = cmdParam.cutoff
	CUTOFF_THICK = cmdParam.cutthick
	CUTOFF_APL = cmdParam.cutapl
	MAIN_AXIS = cmdParam.main_axis

	filesFound = {}

	if MD_RUN_INPUT is not None and FOLDER is not None:
		"""
		These options enables membrane property computation for a lot of data
		If you want only a peculiar file set it using -c and not -r and -f
		"""
		if FOLDER[-1] == '/': FOLDER = FOLDER[:-1]
		#print(FOLDER)
		
		for Run in MD_RUN_INPUT:
			for path in glob.glob(FOLDER +'/**/*'+ Run +'*gro', recursive=True):
				if 'fatslimAnalysis' in path: continue
			
				path_test = path.replace(FOLDER,'').split('/')
				path = path.split('/')
				filename = '/'.join(path)
				
				project_name = None
				job_name = None
				
				if len(path_test) >= 3 : 
					project_name = path[-3]
					job_name = path[-2]
					"""
					#if project_name not in filesFound: 
						#filesFound[project_name] = { job_name: [filename] }
						
					#else:
						#if job_name not in filesFound[project_name]:
							#filesFound[project_name][job_name] = [filename]
							
						#else:
							#filesFound[project_name][job_name].append(filename)
					"""
				else:
					project_name = path[-3]
					job_name = path[-2]
					
				if 'JOB' in job_name and 'OUTPUT' in job_name:
					job_name = path[-3]
					project_name = path[-4]
				
				if project_name not in filesFound: 
					filesFound[project_name] = { job_name: [filename] }
					
				else:
					if job_name not in filesFound[project_name]:
						filesFound[project_name][job_name] = [filename]
						
					else:
						filesFound[project_name][job_name].append(filename)
						
		if TRAJECTORY is not None:
			for Traj in TRAJECTORY:
				for path in glob.glob(FOLDER +'/**/*'+ Traj +'*xtc', recursive=True):
					if 'fatslimAnalysis' in path: continue
				
					path_test = path.replace(FOLDER,'').split('/')
					path = path.split('/')
					filename = '/'.join(path)
					
					project_name = None
					job_name = None
					
					if len(path_test) >= 3 : 
						project_name = path[-3]
						job_name = path[-2]
					else:
						project_name = path[-3]
						job_name = path[-2]
					
					if 'JOB' in job_name and 'OUTPUT' in job_name:
						job_name = path[-3]
						project_name = path[-4]
					
					if project_name not in filesFound: 
						filesFound[project_name] = { job_name: [filename] }
						
					else:
						if job_name not in filesFound[project_name]:
							filesFound[project_name][job_name] = [filename]
							
						else:
							filesFound[project_name][job_name].append(filename)
		
		#print(filesFound)
		#Get the depth of the dictionnary
		DEPTH = depth(filesFound)
		print(filesFound)
		
		#Creates the base directory
		analysisFolder = FOLDER+'/fatslimAnalysis'
		#print(analysisFolder)
		if not os.path.isdir(analysisFolder):
			os.makedirs(analysisFolder, exist_ok=True)
		else:
			bakCount = 0
			backupExt = '-bak'
			
			for path in glob.glob(FOLDER +'/*fatslimAnalysis*'):
				print(path)
				tmp = path.split('/')[-1]
				
				if backupExt in tmp:
					bakCount += 1
					
			
			backupFile = analysisFolder + backupExt+str(bakCount)
			
			try:
				os.rename(analysisFolder, backupFile)
			except (OSError, IOError) as e:
				print(sys.stderr, 'Error moving {0} to {1}: {2}'.format(analysisFolder, backupFile, e))
				
			os.makedirs(analysisFolder, exist_ok=True)
			
			
		# Set the make index command
		makeIndexOpt = "EOF\na {0}\na DEF\nq\nEOF".format(HEAD_GROUP)
		
		if DEPTH >= 3:
			for project_name in filesFound:
				for job_name, file_names in filesFound[project_name].items():
					os.makedirs(analysisFolder+'/'+job_name, exist_ok=True)
					
					gro_file = None
					xtc_file = None
					ndx_file = None
					destination = None
					
					plotname = None
					csvname = None
					
					# Directory for NDX leaflets
					ndx_dir = analysisFolder+'/'+job_name+'/NDX'
					os.makedirs(ndx_dir)
					# Directory for XVG files
					xvg_dir = analysisFolder+'/'+job_name+'/XVG'
					os.makedirs(xvg_dir)
					
					#Directories for PDB files
					pdb_apl_dir = analysisFolder+'/'+job_name+'/PDB/APL'
					pdb_order_dir = analysisFolder+'/'+job_name+'/PDB/ORDER'
					pdb_thickness_dir = analysisFolder+'/'+job_name+'/PDB/THICKNESS'
					
					os.makedirs(pdb_apl_dir)
					os.makedirs(pdb_order_dir)
					os.makedirs(pdb_thickness_dir)
					
					#Directories for CSV files
					csv_apl_dir = analysisFolder+'/'+job_name+'/CSV/APL'
					csv_order_dir = analysisFolder+'/'+job_name+'/CSV/ORDER'
					csv_thickness_dir = analysisFolder+'/'+job_name+'/CSV/THICKNESS'
					
					os.makedirs(csv_apl_dir)
					os.makedirs(csv_order_dir)
					os.makedirs(csv_thickness_dir)
					
					print(file_names)
					for file_name in file_names:
						name = file_name.split('/')[-1]
						print(name)
						
						#Copying the files
						destination = analysisFolder+'/'+job_name+'/'+name
						#shutil.copyfile(file_name, destination)
						
						if TRAJECTORY is None:
							# Setting the file names for fatslim
							gro_file = file_name
							
							#File name for fatslim index
							hg_membrane_ndx = ndx_dir +'/'+ name.replace('.gro','_HG.ndx')
							
							ndx_file = destination.replace('.gro','.ndx')
							#put a temp name for plot file name
							plotname = xvg_dir +'/'+ name.replace('.gro','tmp')
							#put a temp name for csv file name
							csv_apl_name = csv_apl_dir +'/'+ name.replace('.gro','_APL.csv')
							csv_order_name = csv_order_dir +'/'+ name.replace('.gro','_ORDER.csv')
							csv_thickness_name = csv_thickness_dir +'/'+ name.replace('.gro','_THICKNESS.csv')
							
						else:
							if 'gro' in name:
								gro_file = file_name
								
								#File name for fatslim index
								hg_membrane_ndx = ndx_dir +'/'+ name.replace('.gro','_HG.ndx')
								
							if 'xtc' in name: 
								##Copying the .edr file for box dimensions
								#edr_file = file_name.replace('xtc','edr')
								#shutil.copyfile(edr_file, destination.replace('.xtc', '.edr'))
								
								xtc_file = file_name
								
								ndx_file = destination.replace('.xtc','.ndx')
								#put a temp name for plot file name
								plotname = destination.replace('.xtc','tmp')
								#put a temp name for csv file name
								csvname = destination.replace('.xtc','tmp')
								
								#put a temp name for plot file name
								plotname = xvg_dir +'/'+ name.replace('.xtc','tmp')
								#put a temp name for csv file name
								csv_apl_name = csv_apl_dir +'/'+ name.replace('.xtc','_APL.csv')
								csv_order_name = csv_order_dir +'/'+ name.replace('.xtc','_ORDER.csv')
								csv_thickness_name = csv_thickness_dir +'/'+ name.replace('.xtc','_THICKNESS.csv')
								
					
					#Making the index
					makeIndexCmd = "{0} make_ndx -f {1} -o {2} << {3}".format(GMX, gro_file, ndx_file,
																				makeIndexOpt)
					try:
						sub.call(makeIndexCmd, shell=True)
					except CalledProcessError as e:
						print("Error: {0}".format(e))
					
					
					
					# Getting the box dimensions for the plots
					box_file = ndx_file.replace('.ndx','_BOX.xvg')
					if TRAJECTORY is not None:
						#Need the box dimensions for the grid when plotting
						trajBoxCmd = "{0} traj -s {1} -f {2} -ob {3} ".format(GMX, gro_file, xtc_file, box_file)
						if BEGIN_FRAME is not None:
							trajBoxCmd += "-b {0} ".format(BEGIN_FRAME)
						if END_FRAME is not None:
							trajBoxCmd += "-e {0} ".format(END_FRAME)
						
						trajBoxCmd += "<< EOF\nSystem\nEOF"
						
						try:
							sub.call(trajBoxCmd, shell=True)
						except CalledProcessError as e:
							print("Error: {0}".format(e))
						
					else:
						with open(gro_file) as fp:
							for line in fp:
								line = line.strip()

								if len(line) == 0:
									continue

								last_line = line

							box_x, box_y, box_z = [float(val) for val in line.split()[:3]]
						with open(box_file,'w') as boxF:
							boxF.write(' 0     {0} {1} {2}'.format(box_x, box_y, box_z))
							
					
					
					#Look in the index to find DEF atoms
					interacting_group = None
					with open(ndx_file,'r') as index:
						indexStr = index.read()
						if 'DEF' in indexStr:
							interacting_group = 'DEF'
									
					
					
					fatslim_memb_cmd = "{0} membranes -c {1}  -n {2}  --output-index-hg {3} --cutoff {4} --hg-group {5} ".format(
																										FATSLIM,
																										gro_file, 
																										ndx_file, 
																										hg_membrane_ndx,
																										CUTOFF,
																										HEAD_GROUP)
					
					fatslim_apl_cmd = "{0} apl -c {1}  -n {2} --cutoff {3} --apl-cutoff {4} --hg-group {5} ".format(FATSLIM, gro_file, 
																													ndx_file,
																													CUTOFF,
																													CUTOFF_APL,
																													HEAD_GROUP)
					
					fatslim_thick_cmd = "{0} thickness -c {1}  -n {2} --cutoff {3} --thickness-cutoff {4} --hg-group {5} ".format(FATSLIM, 
																																gro_file,
																																ndx_file,
																																CUTOFF,
																																CUTOFF_THICK, 
																																HEAD_GROUP)
					
					fatslim_order_cmd = "{0} order -c {1}  -n {2} --cutoff {3} --hg-group {4} ".format(FATSLIM, gro_file, ndx_file, CUTOFF,
																										HEAD_GROUP)
					
					if MAIN_AXIS is not None:
						fatslim_order_cmd += "--main-axis {0} {1} {2} ".format(MAIN_AXIS[0], MAIN_AXIS[1], MAIN_AXIS[2])
						
					if TRAJECTORY is not None: 
						fatslim_memb_cmd += "-t {0} ".format(xtc_file)
						fatslim_apl_cmd += "-t {0} ".format(xtc_file)
						fatslim_thick_cmd += "-t {0} ".format(xtc_file)
						fatslim_order_cmd += "-t {0} ".format(xtc_file)
						
						if BEGIN_FRAME is not None: 
							fatslim_memb_cmd += "--begin-frame {0} ".format(BEGIN_FRAME)
							fatslim_apl_cmd += "--begin-frame {0} ".format(BEGIN_FRAME)
							fatslim_thick_cmd += "--begin-frame {0} ".format(BEGIN_FRAME)
							fatslim_order_cmd += "--begin-frame {0} ".format(BEGIN_FRAME)
							
						if END_FRAME is not None: 
							fatslim_memb_cmd += "--end-frame {0} ".format(END_FRAME)
							fatslim_apl_cmd += "--end-frame {0} ".format(END_FRAME)
							fatslim_thick_cmd += "--end-frame {0} ".format(END_FRAME)
							fatslim_order_cmd += "--end-frame {0} ".format(END_FRAME)
						
					#Add interacting group for def found previously
					if interacting_group is not None:
						fatslim_memb_cmd += "--interacting-group {0} ".format(interacting_group)
						fatslim_apl_cmd += "--interacting-group {0} ".format(interacting_group)
						fatslim_thick_cmd += "--interacting-group {0} ".format(interacting_group)
						fatslim_order_cmd += "--interacting-group {0} ".format(interacting_group)
						
					# Adding plot for properties
					fatslim_apl_cmd += "--plot-apl {0} --plot-area {1} ".format(plotname.replace('tmp','_APL.xvg'),
																				plotname.replace('tmp','_AREA.xvg'))
					fatslim_thick_cmd += "--plot-thickness {0} ".format(plotname.replace('tmp','_THICKNESS.xvg'))
					fatslim_order_cmd += "--plot-order {0} ".format(plotname.replace('tmp','_ORDER.xvg'))
					
					# Adding plot for properties
					fatslim_apl_cmd += "--export-apl-raw {0} ".format(csv_apl_name)
					fatslim_thick_cmd += "--export-thickness-raw {0} ".format(csv_thickness_name)
					fatslim_order_cmd += "--export-order-raw {0} ".format(csv_order_name)
					
					#cmd_for_all = "#! /bin/bash -x\n"
					cmd_for_all = "## Script generated from the following command:\n##"
					cmd_for_all += "## {0} \n##\n\n".format(' '.join(sys.argv))
					cmd_for_all += fatslim_memb_cmd+"\n\n"+fatslim_apl_cmd+"\n\n"
					cmd_for_all += fatslim_thick_cmd+"\n\n"
					cmd_for_all += fatslim_order_cmd+"\n\n"
					
					
					script_file = '/'.join(ndx_file.split('/')[:-1])+"/fatslim-script.sh"
					
					with open(script_file,'w') as script:
						script.write(cmd_for_all)
					
					try:
						sub.call(cmd_for_all, shell=True)
					except CalledProcessError as e:
						print("Error: {0}".format(e))
					
					
					#Creating directories for property under graphs
					apl_graph_dir = analysisFolder+'/'+job_name+'/graphs/APL'
					order_graph_dir = analysisFolder+'/'+job_name+'/graphs/ORDER'
					thickness_graph_dir = analysisFolder+'/'+job_name+'/graphs/THICKNESS'
					
					os.makedirs(apl_graph_dir, exist_ok=True)
					os.makedirs(order_graph_dir, exist_ok=True)
					os.makedirs(thickness_graph_dir, exist_ok=True)
					
					##list of csv_files
					csv_files = [csv_apl_name, csv_order_name, csv_thickness_name]
					graph_dirs = [apl_graph_dir, order_graph_dir, thickness_graph_dir]
					pdb_dirs = [pdb_apl_dir, pdb_order_dir, pdb_thickness_dir]
					
					#Create a list with box dimensions for each step
					
								
					#print(box_dim)
					
					#for csvf, graphd, pdbd in zip(csv_files, graph_dirs, pdb_dirs):
						#csv_frame_base = csvf.replace('.csv','') + '*'
						#print(csv_frame_base)
						#csv_frames = glob.glob(csv_frame_base)
						
						
						#for csv_f, bdim in zip(csv_frames, box_dim):
							#try:
									#plotCSV(csv_f, bdim, graphd)
									#PDB_out(csv_f, pdbd)
							#except:
								#exctype, value = sys.exc_info()[:2]
								#print( "Error: {0}\n".format(exctype) )
								#print("Value: {0}\n\n".format(value))
								#traceback.print_exc()
						##getMeanMap()Greys
						##
					
		else:
			for job_name, file_names in filesFound.items():
				
					os.makedirs(analysisFolder+'/'+job_name, exist_ok=True)
					
					for file_name in file_names:
						destination = analysisFolder+'/'+job_name+'/'+file_name.split('/')[-1]
						shutil.copyfile(file_name, destination)
						
	
	print( "--- Computing done in {0:f} seconds ---".format(time.time() - start_time) )
	
	
elif 'analyse' in sys.argv:
	
	#for multiprocessing
	import multiprocessing as multproc
	
	
	try:
		import numpy as np
		import numpy.ma as ma
	except ImportError as e:
		print("You need numpy to use this script")
		print("Error {0}".format(e))
		
	try:
		import matplotlib.pyplot as plt
		import matplotlib.cm as mcm
		from matplotlib import rcParams
		from matplotlib import rc
		import matplotlib.ticker as ticker #to change the ticks
		
		params = {
		'pgf.texsystem': 'xelatex',        # change this if using xetex or lautex
		'font.size': 16,
		'figure.figsize': [8,4],
		'axes.labelsize': 20,
		'legend.fontsize': 20,
		'xtick.labelsize': 18,
		'ytick.labelsize': 18,
		'text.usetex': True,
		'font.family': 'sans-serif',
		"text.latex.preamble": [
			r"\usepackage[utf8]{inputenc}",    # use utf8 fonts 
			r"\usepackage[T1]{fontenc}",        # plots will be generated

			r"\usepackage[detect-all]{siunitx}",         # load additional packages
			r"\usepackage{amsmath}",
			r"\usepackage{amssymb}",
			r"\sisetup{mode = math,math-rm=\mathsf}"
			]
		}
		rcParams.update(params)
		
	except ImportError as e:
		print("You need matplotlib.pyplot to use this script")
		print("Error {0}".format(e))
		
	try:
		from scipy.interpolate import griddata
	except ImportError as e:
		print("You need scipy.interpolate to use this script")
		print("Error {0}".format(e))
		
	try:
		# for data plotting
		import pandas as pd
	except ImportError as e:
		print("You need pandas to use this script")
		print("Error {0}".format(e))
	
	
	FOLDER = cmdParam.folder
	if FOLDER[-1] == '/':
		FOLDER = FOLDER[:-1]
	NAMES = cmdParam.name
	STRIDE = cmdParam.stride
	
	BEGIN = None
	if cmdParam.begin is not None: 
		BEGIN = cmdParam.begin
		
	END = None
	if cmdParam.end is not None:
		END = cmdParam.end + 1
	
	PDB_SWITCH = cmdParam.pdb
	MAP_SWITCH = cmdParam.map
	HIST_SWITCH = cmdParam.hist
	SCATTER_SWITCH = cmdParam.scatter
	RADIAL = cmdParam.radial
	
	PROPERTIES = [prop.upper() for prop in cmdParam.properties]
	
	"""
	Not used yet 
	For each mode
	HIST_STRIDE = cmdParam.hist_stride
	if HIST_STRIDE is None : 
		HIST_STRIDE = STRIDE
	
	HIST_BEGIN = cmdParam.hist_begin
	if HIST_BEGIN is None : HIST_BEGIN = BEGIN
	
	HIST_END = cmdParam.hist_end
	if HIST_END is None : HIST_END = END
	
	
	MAP_STRIDE = cmdParam.map_stride
	if MAP_STRIDE is None : 
		MAP_STRIDE = STRIDE
	
	MAP_BEGIN = cmdParam.map_begin
	if MAP_BEGIN is None : 
		MAP_BEGIN = BEGIN
	else:
		MAP_BEGIN = cmdParam.map_begin + 1
	
	MAP_END = cmdParam.map_end
	if MAP_END is None :
		MAP_END = END
	else:
		MAP_END = cmdParam.map_end + 1
	
	
	SCATTER_STRIDE = cmdParam.scatter_stride
	if SCATTER_STRIDE is None : 
		SCATTER_STRIDE = STRIDE
	
	SCATTER_BEGIN = cmdParam.scatter_begin
	if SCATTER_BEGIN is None : 
		SCATTER_BEGIN = BEGIN
	
	SCATTER_END = cmdParam.scatter_end
	if SCATTER_END is None : SCATTER_END = END
	"""
	
	
	
	csv_found = OrderedDict()
	data = OrderedDict()
	
	analysis_folder = FOLDER
	# Finding the .csv files
	print("Finding csv files from FatSlim ...")
	for prop in PROPERTIES:
		tmp_found = glob.glob(analysis_folder+'/**/*'+prop+'*csv', recursive=True)
		if NAMES is not None:
			for files in tmp_found:
				for name in NAMES:
					if name.lower() in files.lower():
						job_name = files.split('/')[-4]
						
						if job_name in csv_found:
							if prop in csv_found[job_name]:
								csv_found[job_name][prop].append(files)
								
							else:
								csv_found[job_name].update( { prop : [files] } )
								data[job_name][prop] = None
						else:
							csv_found[job_name] = { prop: [files] }
							data[job_name] ={prop : None}
		else:
			for files in tmp_found:
				job_name = files.split('/')[-4]
				if job_name in csv_found:
					if prop in csv_found[job_name]:
						csv_found[job_name][prop].append(files)
						
					else:
						csv_found[job_name].update( { prop : [files] } )
						data[job_name][prop] = None
				else:
					csv_found[job_name] = { prop: [files] }
					data[job_name] ={ prop : None }
		
	
	
	start_time = time.time()
	for job_name in csv_found:
		for prop in csv_found[job_name]:
			frames = [pd.read_csv(csv) for csv in csv_found[job_name][prop][BEGIN:END:STRIDE]]
			data[job_name][prop] = frames
	print( "--- gathering data done in {1:f} seconds ---".format(job_name ,time.time() - start_time) )
	
	if SCATTER_SWITCH:
		try:
			scatter = multproc.Process(target=scatter_cmd, args=(data,), name='Process-SCATTER')
			scatter.start()
		except multproc.ProcessError as error:
			print(error)
		
	if PDB_SWITCH:
		try:
			scatter = multproc.Process(target=pdb_cmd, args=(data,), name='Process-PDB')
			scatter.start()
		except multproc.ProcessError as error:
			print(error)
					
	if MAP_SWITCH:
		try:
			scatter = multproc.Process(target=map_cmd, args=(data,), name='Process-MAP')
			scatter.start()
		except multproc.ProcessError as error:
			print(error)

	if HIST_SWITCH:
		try:
			scatter = multproc.Process(target=hist_cmd, args=(data,), name='Process-HIST')
			scatter.start()
		except multproc.ProcessError as error:
			print(error)
	
	


elif 'animate' in sys.argv:
	
	try:
		import imageio
	except ImportError as e:
		print("You need imageio to use this command")
		print("Error {0}".format(e))
	try:
		import io
	except ImportError as e:
		print("You need io to use this command")
		print("Error {0}".format(e))
	try:
		import cairo
	except ImportError as e:
		print("You need cairo to use this command")
		print("Error {0}".format(e))
	try:
		from gi import require_version
		require_version('Rsvg', '2.0')
		from gi.repository import Rsvg as rsvg
	except ImportError as e:
		print("You need gir to use this command")
		print("Error {0}".format(e))
	
	start_time = time.time()

	FOLDER = cmdParam.folder
	if FOLDER[-1] == '/': FOLDER = FOLDER[:-1]
	
	filenames = glob.glob(FOLDER+'/*frame*.svg')
	filenames.sort()
	
	svg_files = [ rsvg.Handle().new_from_file(file_name) for file_name in filenames]
	
	with imageio.get_writer(FOLDER+'/movie.mp4', mode='I') as writer:
		print('Creating movie.mp4 ...')
		for svg in svg_files:
			width = svg.props.width
			height = svg.props.height
			img = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
			
			ctx = cairo.Context(img)
			svg.render_cairo(ctx)
			
			img.write_to_png('tmp.png')
			
			image = imageio.imread('tmp.png')
			writer.append_data(image)
		
	print('Done !')
	os.remove('tmp.png')
	
	print( "--- Movie made in {0:f} ---".format(time.time() - start_time) )
	
elif 'mayavi' in sys.argv:
	
	try:
		from mayavi import mlab
	except ImportError as e:
		print("If you want to see 3d representation of the membrane please install mayavi")
	else:
		
		FILENAME	= cmdParam.filename
		TYPE		= cmdParam.type
		filename = None
		
		try:
			filename = glob.glob(FILENAME)[0]
		except IOError as e:
			print(e)
			
		if TYPE == 'fatslim':
			Membrane = {'lower leaflet': {'x':[],'y':[],'z':[],'xd':[],'yd':[],'zd':[],'xn':[],'yn':[],'zn':[],'prop value':[]},
			'upper leaflet':{'x':[],'y':[],'z':[],'xd':[],'yd':[],'zd':[],'xn':[],'yn':[],'zn':[],'prop value':[]}}
			with open(filename) as fp:
				for lino, line in enumerate(fp):
					if lino == 0:
						membrane_property = line.split(",")[-1].strip()

					else:
						line = line.strip()

						if len(line) == 0:
							continue

						resid, leaflet, x, y, z, xd, yd, zd, xn, yn, zn, value = line.split(",")
						Membrane[leaflet]['x'].append(float(x))
						Membrane[leaflet]['y'].append(float(y))
						Membrane[leaflet]['z'].append(float(z))

						Membrane[leaflet]['xd'].append(- float(xd))
						Membrane[leaflet]['yd'].append(- float(yd))
						Membrane[leaflet]['zd'].append(- float(zd))

						Membrane[leaflet]['xn'].append(float(xn))
						Membrane[leaflet]['yn'].append(float(yn))
						Membrane[leaflet]['zn'].append(float(zn))
						Membrane[leaflet]['prop value'].append(float(value))
			
			#Lower leaflet
			lower_leaflet_pts = mlab.points3d(Membrane['lower leaflet']['x'], Membrane['lower leaflet']['y'], Membrane['lower leaflet']['z'], Membrane['lower leaflet']['prop value'])
			
			lower_leaflet_directions = mlab.quiver3d(Membrane['lower leaflet']['x'], Membrane['lower leaflet']['y'], Membrane['lower leaflet']['z'], 
							Membrane['lower leaflet']['xd'], Membrane['lower leaflet']['yd'], Membrane['lower leaflet']['zd'], mode='arrow', scale_factor=1.5, resolution=32, color=(1,0,0))
			
			lower_leaflet_normals = mlab.quiver3d(Membrane['lower leaflet']['x'], Membrane['lower leaflet']['y'], Membrane['lower leaflet']['z'], 
							Membrane['lower leaflet']['xn'], Membrane['lower leaflet']['yn'], Membrane['lower leaflet']['zn'], mode='arrow', scale_factor=1.5, resolution=32, color=(0,0,1))
			
			#Upper leaflet
			upper_leaflet_pts = mlab.points3d(Membrane['upper leaflet']['x'], Membrane['upper leaflet']['y'], Membrane['upper leaflet']['z'], Membrane['upper leaflet']['prop value'])
			
			upper_leaflet_directions = mlab.quiver3d(Membrane['upper leaflet']['x'], Membrane['upper leaflet']['y'], Membrane['upper leaflet']['z'], 
							Membrane['upper leaflet']['xd'], Membrane['upper leaflet']['yd'], Membrane['upper leaflet']['zd'], mode='arrow', scale_factor=1.5, resolution=32, color=(1,0,0))
			
			upper_leaflet_normals = mlab.quiver3d(Membrane['upper leaflet']['x'], Membrane['upper leaflet']['y'], Membrane['upper leaflet']['z'], 
							Membrane['upper leaflet']['xn'], Membrane['upper leaflet']['yn'], Membrane['upper leaflet']['zn'], mode='arrow', scale_factor=1.5, resolution=32, color=(0,0,1))

			
			lower_leaflet_pts.glyph.glyph.scale_mode = "data_scaling_off"
			lower_leaflet_pts.glyph.glyph.scale_factor = 1.0
			lower_leaflet_pts.glyph.glyph_source.glyph_source.phi_resolution = 32
			lower_leaflet_pts.glyph.glyph_source.glyph_source.theta_resolution = 32
			lower_leaflet_pts.glyph.glyph.color_mode = "color_by_scalar"
			
			upper_leaflet_pts.glyph.glyph.scale_mode = "data_scaling_off"
			upper_leaflet_pts.glyph.glyph.scale_factor = 1.0
			upper_leaflet_pts.glyph.glyph_source.glyph_source.phi_resolution = 32
			upper_leaflet_pts.glyph.glyph_source.glyph_source.theta_resolution = 32
			upper_leaflet_pts.glyph.glyph.color_mode = "color_by_scalar"
			
			lower_leaflet_directions.glyph.glyph_source.glyph_source.shaft_radius = 0.07
			lower_leaflet_directions.glyph.glyph_source.glyph_source.tip_radius = 0.2
			
			upper_leaflet_directions.glyph.glyph_source.glyph_source.shaft_radius = 0.07
			upper_leaflet_directions.glyph.glyph_source.glyph_source.tip_radius = 0.2
			
			lower_leaflet_normals.glyph.glyph_source.glyph_source.shaft_radius = 0.07
			lower_leaflet_normals.glyph.glyph_source.glyph_source.tip_radius = 0.2
			
			upper_leaflet_normals.glyph.glyph_source.glyph_source.shaft_radius = 0.07
			upper_leaflet_normals.glyph.glyph_source.glyph_source.tip_radius = 0.2
			
			mlab.colorbar()
			mlab.show()
			
		elif TYPE in ['loctemp', 'locpres']:
			
			try:
				import numpy as np
			except ImportError as e:
				print("You need numpy to use this command")
				print("Error {0}".format(e))
		
			loctemp_data = np.load(filename)
			
			grid_x = loctemp_data['mgx']
			grid_y = loctemp_data['mgy']
			grid_z = loctemp_data['mgz']
			grid = loctemp_data['mg']
			
			if TYPE == 'loctemp':
				mlab.figure(figure="local temperature", bgcolor=(1.0, 1.0, 1.0))
			elif TYPE == 'locpres':
				mlab.figure(figure="local pressure", bgcolor=(1.0, 1.0, 1.0))
				
			# Volume data
			local_property = mlab.pipeline.scalar_field(grid_x, grid_y, grid_z, grid)
			nb_points1D = len(grid)
			volume_data = mlab.pipeline.volume(local_property)
			
			# Planes
			mlab.pipeline.image_plane_widget(local_property,
											plane_orientation='x_axes',
											slice_index=nb_points1D/2.,
										)
			mlab.pipeline.image_plane_widget(local_property,
										plane_orientation='y_axes',
										slice_index=nb_points1D/2.,
									)
			volume_data.lut_manager.show_scalar_bar		= True
			if TYPE == 'loctemp':
				volume_data.lut_manager.scalar_bar.title	= 'T(K)'
			elif TYPE == 'locpres':
				volume_data.lut_manager.scalar_bar.title	= 'P(bar)'
				
			
			volume_data.lut_manager.title_text_property.bold	= True
			volume_data.lut_manager.label_text_property.bold	= True
			
			volume_data.lut_manager.label_text_property.color	= (0.0, 0.0, 0.0)
			volume_data.lut_manager.title_text_property.color	= (0.0, 0.0, 0.0)
			
			volume_data.volume_mapper.sample_distance		= 0.1
			volume_data.volume_property.diffuse				= 0.9
			volume_data.volume_property.specular = 1.0
			volume_data.volume_property.scalar_opacity_unit_distance	= 0.2
			
			mlab.show()
			
		elif TYPE in ['density']:
			
			try:
				import numpy as np
				#import matplotlib.pyplot as plt
				from scipy.interpolate import interp1d
			except ImportError as e:
				print("You need numpy to use this command")
				print("Error {0}".format(e))
			try:
				from tvtk.util.ctf import ColorTransferFunction
				from tvtk.util.ctf import PiecewiseFunction
				from vtk import vtkLookupTable
			except ImportError as e:
				print("Error {0}".format(e))
		
			density_data = np.load(filename)
			
			grid_x	= density_data['mgx']
			grid_y	= density_data['mgy']
			grid_z	= density_data['mgz']
			grid	= density_data['mg']
			
			if TYPE == 'density':
				mlab.figure(figure="Density", bgcolor=(1.0, 1.0, 1.0))
				
			# Volume data
			density = mlab.pipeline.scalar_field(grid_x, grid_y, grid_z, grid)
			nb_points1D = len(grid)
			volume_data = mlab.pipeline.volume(density)
			
			ctf	= ColorTransferFunction()
			ctf.add_rgb_point(0.0, 1.0,  1.0, 1.0)
			ctf.add_rgb_point(1.0, 0.0,  0.0, 1.0)
			ctf.add_rgb_point(2.0, 0.34, 1.0, 0.0)
			ctf.add_rgb_point(3.0, 1.0,  0.0, 0.0)
			ctf.add_rgb_point(4.0, 1.0,  1.0, 0.0)
			volume_data._volume_property.set_color(ctf)
			volume_data._ctf	= ctf
			volume_data.update_ctf = True
			
			otf = PiecewiseFunction()
			otf.add_point(0.0, 0.01)
			otf.add_point(1.0, 0.005)
			otf.add_point(2.0, 0.5)
			otf.add_point(3.0, 0.5)
			otf.add_point(4.0, 0.5)
			volume_data._otf	= otf
			volume_data._volume_property.set_scalar_opacity(otf)
			
			# Planes
			x_plane = mlab.pipeline.image_plane_widget(density,
											plane_orientation='x_axes',
											slice_index=nb_points1D/2.,
										)
			y_plane = mlab.pipeline.image_plane_widget(density,
										plane_orientation='y_axes',
										slice_index=nb_points1D/2.,
									)
			
			z_plane = mlab.pipeline.image_plane_widget(density,
										plane_orientation='z_axes',
										slice_index=nb_points1D/2.,
									)
			
			planes = [x_plane, y_plane, z_plane]
			
			R = [255,   0,   0, 255, 255]
			G = [255,   0, 255,   0, 255]
			B = [255, 255,   0,   0,   0]
			A = [255, 255, 255, 255, 255]
			X = [0.0, 1.0, 2.0, 3.0, 4.0]
			x_interp = np.linspace(0, 4.0, 256, endpoint=True)
			
			R_interp = interp1d(X, R, kind='linear')
			G_interp = interp1d(X, G, kind='linear')
			B_interp = interp1d(X, B, kind='linear')
			A_interp = interp1d(X, A, kind='linear')

			R_interp = np.clip(R_interp(x_interp),0,255)
			G_interp = np.clip(G_interp(x_interp),0,255)
			B_interp = np.clip(B_interp(x_interp),0,255)
			A_interp = A_interp(x_interp)
			
			#plt.plot(x_interp, R_interp,'r')
			#plt.plot(x_interp, G_interp,'g')
			#plt.plot(x_interp, B_interp,'b')
			#plt.plot(x_interp, A_interp,'k')
			#plt.show()
			RGBA = np.matrix([R_interp, G_interp, B_interp, A_interp])

			
			lut = RGBA.T
			for plane in planes:
				plane.module_manager.scalar_lut_manager.lut.table = lut
			
			x_plane.module_manager.scalar_lut_manager.use_default_range = False
			x_plane.module_manager.scalar_lut_manager.data_range = [0.0, 4.0]

			volume_data.volume_mapper.sample_distance		= 0.1
			volume_data.volume_property.diffuse				= 0.9
			volume_data.volume_property.specular = 1.0
			volume_data.volume_property.scalar_opacity_unit_distance	= 0.1
			volume_data.volume_property.shade = False
			volume_data.volume_property.interpolation_type = 'nearest'
			
			mlab.draw()
			mlab.show()
			
			
elif 'vpython' in sys.argv:
	try:
		import numpy as np
	except ImportError as e:
		print("You need numpy to use this command")
		print("Error {0}".format(e))
		
	try:
		import vpython as vp
	except ImportError as e:
		print("You need vpython to use this command")
		print("Error {0}".format(e))
		
	FILENAME = cmdParam.filename
	datatype = np.dtype({'names' : ('resid', 'leaflet','X','Y','Z','XD','YD','ZD','XN','YN','ZN', 'property'),
				'formats' : ('i', 'S13','f','f','f','f','f','f','f','f','f','f')})
	
	data = np.loadtxt(FILENAME, skiprows=1, dtype=datatype, delimiter=',')
	
	box_filename = "{0}/*BOX*".format( '/'.join(FILENAME.split('/')[:-3]) )
	box_filename = glob.glob(box_filename)[0]
	
	box_dim = []
	nb_frames = None
	with open(box_filename,'r') as boxF:
		for line in boxF:
			if line.startswith('@') or line.startswith('#'):
				continue
			else:
				box = (int(line.split()[0]),)
				box = box + tuple( [float(DIM) for DIM in line.split()[1:]] )
				box_dim.append(box)
	
	boxdatatype = np.dtype({'names' : ('time', 'XX','YY','ZZ','YX','ZX','ZY'),
				'formats' : ('i', 'f','f','f','f','f','f')})
	box_dim = np.array(box_dim, dtype=boxdatatype)
	
	print(len(box_dim))
	lower_leaflet_atoms = []
	lower_leaflet_normals = []
	lower_leaflet_directions = []

	upper_leaflet_atoms = []
	upper_leaflet_normals = []
	upper_leaflet_directions = []

	scene = vp.canvas(width=1280, height=720)
	
	index = 20000
	vector_centering = vp.vector(0.5*box_dim['XX'][index], 0.5*box_dim['ZZ'][index], 0.5*box_dim['YY'][index])
	box_size = vp.vector(box_dim['XX'][index], box_dim['ZZ'][index], box_dim['YY'][index])
	MDbox = vp.box(pos=vp.vector(0,0,0), size=box_size, opacity=0.1)
	for dat in data:
		position = vp.vector(dat['X'], dat['Z'], dat['Y']) - vector_centering
		normal_axis = vp.vector(dat['XN'], dat['ZN'], dat['YN'])
		direction_axis = vp.vector(-2*dat['XD'], -2*dat['ZD'], -2*dat['YD'])
		atom = vp.sphere( pos=position, radius=0.47, color=vp.color.green )
		normal = vp.arrow(pos=position, axis=normal_axis, color=vp.color.blue)
		direction = vp.arrow(pos=position, axis=direction_axis, color=vp.color.red)
		
		if 'lower leaflet' in dat:
			lower_leaflet_atoms.append(atom)
			lower_leaflet_normals.append(normal)
			lower_leaflet_directions.append(direction)
		
		if 'upper leaflet' in dat:
			upper_leaflet_atoms.append(atom)
			upper_leaflet_normals.append(normal)
			upper_leaflet_directions.append(direction)
			
			

elif 'gromacs' in sys.argv:
		
	start_time = time.time()
	
	# Setting general options
	GMX				=	cmdParam.gmx_path
	MD_RUN_INPUT	=	cmdParam.md_run_input
	EXTENSION		=	cmdParam.extension
	FOLDER			=	cmdParam.folder
	
	# Setting parameters for gmx energy and density
	BEGIN_FRAME	=	cmdParam.begin_frame
	END_FRAME	=	cmdParam.end_frame
	DT			=	cmdParam.dt
	RADIUS		=	cmdParam.radius
	
	#Options to compute
	DENSITY		=	cmdParam.density
	ENERGY		=	cmdParam.energy
	SELECT		=	cmdParam.select
	LOCAL_TEMP	=	cmdParam.loctemp
	LOCAL_PRESS	=	cmdParam.locpres
	
	#OPTIONS FOR LIPIDS
	HEADGROUP	= cmdParam.headgroup
	LINKGROUP	= cmdParam.linkgroup
	TAILGROUP	= cmdParam.tailgroup
	
	
	

	filesFound = {}

	if MD_RUN_INPUT is not None and FOLDER is not None:
		"""
		These options enables membrane property computation for a lot of data
		If you want only a peculiar file set it using -c and not -r and -f
		"""
		
		if FOLDER[-1] == '/': FOLDER = FOLDER[:-1]
		
		
		for Run in MD_RUN_INPUT:
			for path in glob.glob(FOLDER +'/**/*'+ Run +'*'+EXTENSION, recursive=True):
				if 'fatslimAnalysis' in path: continue
			
				path_test = path.replace(FOLDER,'').split('/')
				path = path.split('/')
				filename = '/'.join(path)
				
				project_name = None
				job_name = None
				
				if len(path_test) >= 3 : 
					project_name = path[-3]
					job_name = path[-2]
				else:
					project_name = path[-3]
					job_name = path[-2]
				
				if 'JOB' in job_name and 'OUTPUT' in job_name:
					job_name = path[-3]
					project_name = path[-4]
					
				if project_name not in filesFound: 
					filesFound[project_name] = { job_name: filename }
					
				else:
					if job_name not in filesFound[project_name]:
						filesFound[project_name][job_name] = filename
						
					#else:
						#filesFound[project_name][job_name].append(filename)
		
		DEPTH = depth(filesFound)
		print(filesFound)
		
		#Creates the base directory
		analysisFolder = FOLDER+'/GMX_Analysis'
		#print(analysisFolder)
		if not os.path.isdir(analysisFolder):
			os.makedirs(analysisFolder, exist_ok=True)
		else:
			bakCount = 0
			backupExt = '-bak'
			
			for path in glob.glob(FOLDER +'/*GMX_Analysis*'):
				print(path)
				tmp = path.split('/')[-1]
				
				if backupExt in tmp:
					bakCount += 1
					
			
			backupFile = analysisFolder + backupExt+str(bakCount)
			
			try:
				os.rename(analysisFolder, backupFile)
			except (OSError, IOError) as e:
				print(sys.stderr, 'Error moving {0} to {1}: {2}'.format(analysisFolder, backupFile, e))
				
			os.makedirs(analysisFolder, exist_ok=True)
		
		
		# Directory for gmx energy
		if ENERGY:
			energy_dir = analysisFolder+'/GMX_ENERGY'
			os.makedirs(energy_dir)
		
		# Directory for gmx density
		if DENSITY:
			density_dir = analysisFolder+'/GMX_DENSITY'
			os.makedirs(density_dir)
		
		# Directory for gmx density
		if SELECT:
			select_dir = analysisFolder+'/GMX_SELECT'
			os.makedirs(select_dir)
		
		if HEADGROUP is not None or LINKGROUP is not None or TAILGROUP is not None:
			ndx_dir = analysisFolder +'/GMX_NDX'
			os.makedirs(ndx_dir)
		
		if DEPTH >= 3:
			for project_name in filesFound:
				for job_name, file_name in filesFound[project_name].items():
					

					xtc_file = None
					edr_file = None
					tpr_file = None
					ndx_file = None
					destination = None
					
					plotname = None
					csvname = None
					
					name = file_name.split('/')[-1]
					path = '/'.join(file_name.split('/')[:-1])
					
					xtc_file = file_name
					edr_file = file_name.replace(EXTENSION, 'edr')
					tpr_file = file_name.replace(EXTENSION, 'tpr')
					
					
					ndx_file = path +'/' + '_'.join(name.split('_')[:-1]) + ".ndx"
					
					if '_WALL' in name:
						ndx_file = '_'.join(ndx_file.split('_')[:-1]) + ".ndx"
					
					

					# Get the number of group to output
					read_index_cmd = """cat {0} | grep "\[" """.format(ndx_file)
					read_index_proc = sub.Popen(read_index_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
					read_index_out = read_index_proc.stdout.read()
					nb_index = len(read_index_out.splitlines())
					
					#Making new index including heads, tails and linkage
					create_ndx_head_link_tail = ""
					ndx_file_out = ndx_file
					if HEADGROUP is not None:
						create_ndx_head_link_tail += repr("a {0}\nname {1} head\n".format(' | a '.join(HEADGROUP), nb_index))
						nb_index += 1
						
					if LINKGROUP is not None:
						create_ndx_head_link_tail += repr("a {0}\nname {1} link\n".format(' | a '.join(LINKGROUP), nb_index))
						nb_index += 1
						
					if TAILGROUP is not None:
						create_ndx_head_link_tail += repr("a {0}\nname {1} tail\n".format(' | a '.join(TAILGROUP), nb_index))
						nb_index += 1
					create_ndx_head_link_tail += repr("q\n")
					
					if HEADGROUP is not None or LINKGROUP is not None or TAILGROUP is not None:
						job_ndx_dir = analysisFolder +'/GMX_NDX/'
						
						ndx_file_out = job_ndx_dir + job_name + '.ndx'
						
						input_file = file_name
						if EXTENSION != 'gro':
							input_file = file_name.replace(EXTENSION, 'gro')
						
						try:
							make_ndx_cmd = """ printf {0} | {1} make_ndx -f {2} -n {3} -o {4} \n\n""".format(create_ndx_head_link_tail,
																											GMX, input_file,
																											ndx_file , ndx_file_out)
							sub.call(make_ndx_cmd, shell=True)
						
						except OSError as e:
							print(e)
						
					# Setting the begining and end of analysis
					begin_end = ""
					
					if BEGIN_FRAME is not None:
						begin_end +=" -b {0}".format(BEGIN_FRAME)
					if END_FRAME is not None:
						begin_end +=" -e {0}".format(END_FRAME)
						
					if ENERGY:
						job_energy_dir = analysisFolder+'/GMX_ENERGY/'
						
						energy_xvg_file = job_name + '-' + xtc_file.split('-')[-1].replace('.'+EXTENSION, '_EN.xvg')
						energy_avg_file = energy_xvg_file.replace('xvg','avg')
						
						# To output all quantities
						variables_for_output = "echo "
						for i in range(1, 1000):
							variables_for_output += repr('{0}\n'.format(i))
						variables_for_output += repr('\n\n')
						
						# Call gmx energy
						try:
							energy_cmd = "{0} | {1} energy {2} -f {3} -s {4} -o {5} > {6}".format(variables_for_output, GMX, begin_end, edr_file, tpr_file,
																							job_energy_dir+energy_xvg_file, job_energy_dir+energy_avg_file)
							sub.call(energy_cmd, shell=True)
						
						except OSError as e:
							print(e)
					
					
					
					# Calling gmx density
					if DENSITY:
						job_density_dir = analysisFolder+'/GMX_DENSITY/'
						
						density_xvg_file = job_name + '-' + xtc_file.split('-')[-1].replace('.'+EXTENSION, '_DENS.xvg')
						
						
						
						groups_for_output = "echo "
						for i in range(0, nb_index):
							groups_for_output += repr('{0}\n'.format(i))
						groups_for_output += repr('\n\n')
						
						try:
							density_cmd = "{0} | {1} density -dens number {2} -f {3} -s {4} -o {5} -ng {6} -n {7}".format(groups_for_output, GMX, begin_end, xtc_file, tpr_file,
																											job_density_dir + density_xvg_file, 
																											nb_index, ndx_file_out)
							sub.call(density_cmd, shell=True)
						
						except OSError as e:
							print(e)
					
					
					if SELECT:
						job_select_dir = analysisFolder+'/GMX_SELECT/' + job_name
						os.makedirs(job_select_dir, exist_ok=True)
						
						remove_pertubed_membrane = ""
						above_su = ""
						solvent = "W"
						read_index_out = str(read_index_out)
						
						flipflop_xvg_file = job_name + '-' + xtc_file.split('-')[-1].replace('.'+EXTENSION, '_FLIPFLOP.xvg')
						spl_xvg_file = job_name + '-' + xtc_file.split('-')[-1].replace('.'+EXTENSION, '_SOLPERLEAF.xvg')
						
						if 'defo' in  read_index_out:
							remove_pertubed_membrane = "and (not within {0} of (name DEF))".format(RADIUS)
						if 'su' in read_index_out:
							above_su = " and (z > z of com of group su)"
						if 'PW' in read_index_out:
							solvent = "PW"
							
						mono = ""
						if 'mono' in read_index_out:
							mono = str("and (z < z of com of (resname {2} and (z > z of com of group bilayer)))"
										"").format(remove_pertubed_membrane, solvent)
							
						lower_leaflet_hg_bil = "(name PO4 {0}) and (z < z of com of group bilayer)".format(remove_pertubed_membrane)
						upper_leaflet_hg_bil = "(name PO4 {0}) and (z > z of com of group bilayer) {2}".format(remove_pertubed_membrane, solvent, mono)
						
						lower_leaflet_solvent_bil = "(resname {0} {1}) and (z < z of com of group bilayer) {2}".format(solvent, above_su, remove_pertubed_membrane)
						upper_leaflet_solvent_bil = "(resname {0} {1}) and (z > z of com of group bilayer) {2}".format(solvent, mono, remove_pertubed_membrane)
						
						leaflet_hg_mono = "name PO4 {0} and (z > z of com of (resname {1} and (z > z of com of group bilayer)))".format(remove_pertubed_membrane, solvent)
						leaflet_solvent_mono = str("resname {1} {0} and (z > z of com of group bilayer)"
													" and (z > z of com of (resname {1} and (z > z of com of group bilayer)))").format(remove_pertubed_membrane, solvent)
						
						## Call gmx select
						try:
							select_cmd = """{0} select {1} -f {2} -s {3} -n {4} -oi {5} -select "{6}" \n""".format(GMX, begin_end, xtc_file, tpr_file,
																											ndx_file,
																											job_select_dir+'/'+ 'lower_leaflet_hg_bil.dat',
																											lower_leaflet_hg_bil)
							
							select_cmd += """{0} select {1} -f {2} -s {3} -n {4} -oi {5} -select "{6}" \n""".format(GMX, begin_end, xtc_file, tpr_file,
																											ndx_file,
																											job_select_dir+'/'+ 'upper_leaflet_hg_bil.dat',
																											upper_leaflet_hg_bil)
							
							select_cmd += """{0} select {1} -f {2} -s {3} -n {4} -oi {5} -select "{6}" \n""".format(GMX, begin_end, xtc_file, tpr_file,
																											ndx_file,
																											job_select_dir+'/'+ 'lower_leaflet_solvent_bil.dat',
																											lower_leaflet_solvent_bil)
							
							select_cmd += """{0} select {1} -f {2} -s {3} -n {4} -oi {5} -select "{6}" \n""".format(GMX, begin_end, xtc_file, tpr_file,
																											ndx_file,
																											job_select_dir+'/'+ 'upper_leaflet_solvent_bil.dat',
																											upper_leaflet_solvent_bil)
							
							if 'mono' in read_index_out:
								select_cmd += """{0} select {1} -f {2} -s {3} -n {4} -oi {5} -select "{6}" \n""".format(GMX, begin_end, xtc_file, tpr_file,
																											ndx_file,
																											job_select_dir+'/'+ 'leaflet_hg_mono.dat',
																											leaflet_hg_mono)
								
								select_cmd += """{0} select {1} -f {2} -s {3} -n {4} -oi {5} -select "{6}" \n""".format(GMX, begin_end, xtc_file, tpr_file,
																											ndx_file,
																											job_select_dir+'/'+ 'leaflet_solvent_mono.dat',
																											leaflet_solvent_mono)
							
							
							print(select_cmd)
							sub.call(select_cmd, shell=True)
						
						except OSError as e:
							print(e)
						
						#Need to add mono later
						
						llhb_file = open(job_select_dir+'/'+ 'lower_leaflet_hg_bil.dat','r')
						ulhb_file = open(job_select_dir+'/'+ 'upper_leaflet_hg_bil.dat','r') 
						llsb_file = open(job_select_dir+'/'+ 'lower_leaflet_solvent_bil.dat','r')
						ulsb_file = open(job_select_dir+'/'+ 'upper_leaflet_solvent_bil.dat','r')
						
						if 'mono' in read_index_out:
							lhm_file = open(job_select_dir+'/'+ 'leaflet_hg_mono.dat','r') 
							lsm_file = open(job_select_dir+'/'+ 'leaflet_solvent_mono.dat','r')
						
						#arrays for data
						time_arr = []
						llhb_nb = []
						ulhb_nb = []
						llsb_nb = []
						ulsb_nb = []
						
						for line in llhb_file:
							timestep = float(line.strip().split()[0])
							number = float(line.strip().split()[1])
							
							time_arr.append(timestep)
							llhb_nb.append(number)
							
						file_dict = {'ulhb': {'file': ulhb_file, 'array': ulhb_nb},
										'llsb': {'file': llsb_file, 'array': llsb_nb},
										'ulsb': {'file': ulsb_file, 'array': ulsb_nb},
										}
						
						for key in file_dict:
							dat_file = file_dict[key]['file']
							
							for line in dat_file:
								number = float(line.strip().split()[1])
								file_dict[key]['array'].append(number)
						
						ulhb_nb = file_dict['ulhb']['array']
						llsb_nb = file_dict['llsb']['array']
						ulsb_nb = file_dict['ulsb']['array']
						
						with open(job_select_dir+'/{0}'.format(spl_xvg_file),'a+') as xvg_out:
							header = """
									@    title ""
									@    xaxis  label "Time (ps)"
									@    yaxis  label ""
									@TYPE xy
									@ view 0.15, 0.15, 0.75, 0.85
									@ legend on
									@ legend box on
									@ legend loctype view
									@ legend 0.78, 0.8
									@ legend length 2
									@ s0 legend "lower_leaflet_hg"
									@ s1 legend "upper_leaflet_hg"
									@ s2 legend "lower_leaflet_sol"
									@ s3 legend "upper_leaflet_sol"
									@ s4 legend "sol_per_lower_hg"
									@ s5 legend "sol_per_upper_hg"
									"""
							
							if 'mono' in read_index_out:
								header += """
									@ s6 legend "mono_hg"
									@ s7 legend "mono_sol"
									@ s8 legend "sol_per_mono_hg"
									"""
							
							xvg_out.write(ut.RemoveUnwantedIndent(header)+"\n")
							
							for T, llhb, ulhb, llsb, ulsb in zip(time_arr, llhb_nb, ulhb_nb, llsb_nb, ulsb_nb):
								row = "    {0}    {1}    {2}    {3}    {4}    {5}    {6}\n".format(T, llhb, ulhb, llsb, ulsb, llsb/llhb, ulsb/ulhb)
								xvg_out.write(row)
						
						flip_time = []
						llhb_flip = []
						ulhb_flip = []
						llsb_flip = []
						ulsb_flip = []
						
						llhb_flip_speed = []
						ulhb_flip_speed = []
						llsb_flip_speed = []
						ulsb_flip_speed = []
						
						for i in range(0, len(time_arr), 2):
							if i == (len(time_arr) - 1): break
					
							flip_time.append( (time_arr[i+1] - time_arr[i])/2 + time_arr[i] )
							llhb_flip.append( llhb_nb[i+1] - llhb_nb[i])
							ulhb_flip.append( ulhb_nb[i+1] - ulhb_nb[i])
							llsb_flip.append( llsb_nb[i+1] - llsb_nb[i])
							ulsb_flip.append( ulsb_nb[i+1] - ulsb_nb[i])
							
							llhb_flip_speed.append( (llhb_nb[i+1] - llhb_nb[i])/(time_arr[i+1] - time_arr[i]) )
							ulhb_flip_speed.append( (ulhb_nb[i+1] - ulhb_nb[i])/(time_arr[i+1] - time_arr[i]) )
							llsb_flip_speed.append( (llsb_nb[i+1] - llsb_nb[i])/(time_arr[i+1] - time_arr[i]) )
							ulsb_flip_speed.append( (ulsb_nb[i+1] - ulsb_nb[i])/(time_arr[i+1] - time_arr[i]) )
							
						with open(job_select_dir+'/{0}'.format(flipflop_xvg_file),'a+') as xvg_out:
							header = """
									@    title ""
									@    xaxis  label "Time (ps)"
									@    yaxis  label ""
									@TYPE xy
									@ view 0.15, 0.15, 0.75, 0.85
									@ legend on
									@ legend box on
									@ legend loctype view
									@ legend 0.78, 0.8
									@ legend length 2
									@ s0 legend "lower_leaflet_hg_flip"
									@ s1 legend "upper_leaflet_hg_flip"
									@ s2 legend "lower_leaflet_sol_flip"
									@ s3 legend "upper_leaflet_sol_flip"
									@ s4 legend "lower_leaflet_hg_flip_speed"
									@ s5 legend "upper_leaflet_hg_flip_speed"
									@ s6 legend "lower_leaflet_sol_flip_speed"
									@ s7 legend "upper_leaflet_sol_flip_speed"
									"""
							
							xvg_out.write(ut.RemoveUnwantedIndent(header)+"\n")
							
							for T, llhb, ulhb, llsb, ulsb, llhbspd, ulhbspd, llsbspd, ulsbspd in zip(flip_time, llhb_flip, ulhb_flip, llsb_flip, ulsb_flip, 
																llhb_flip_speed, ulhb_flip_speed, llsb_flip_speed,ulsb_flip_speed):
								row = "    {0}    {1}    {2}    {3}    {4}    {5}    {6}    {7}    {8}\n".format(T, llhb, ulhb, llsb, ulsb,
																												llhbspd, ulhbspd, llsbspd, ulsbspd)
								xvg_out.write(row)
							
						llhb_file.close()
						ulhb_file.close()
						llsb_file.close()
						ulsb_file.close()
					
					
					
		else:
			for job_name, file_names in filesFound.items():
				
					os.makedirs(analysisFolder+'/'+job_name, exist_ok=True)
					
					for file_name in file_names:
						destination = analysisFolder+'/'+job_name+'/'+file_name.split('/')[-1]
						shutil.copyfile(file_name, destination)
						
	
	print( "--- gromacs done in {0:f} seconds ---".format(time.time() - start_time) )

elif 'mda' in sys.argv:
	""" Compute local temperature and temperature, dynamic selected density using MDAnalysis"""
	start_time = time.time()
	try:
		import numpy as np
	except ImportError as e:
		print("You need numpy to use this command")
		print("Error {0}".format(e))
	try:
		import pandas as pd
	except ImportError as e:
		print("If you want to see 3d representation of the membrane please install mayavi")

	try:
		import MDAnalysis as mda
	except ImportError as e:
		print("You need MDAanalysis to use this command")
		print("Error {0}".format(e))
		
	try:
		from scipy.interpolate import griddata
	except ImportError as e:
		print("You need scipy.interpolate to use this script")
		print("Error {0}".format(e))
	
	#
	MD_RUN_INPUT	=	cmdParam.md_run_input
	EXTENSION		=	cmdParam.extension
	FOLDER			=	cmdParam.folder
	
	# Setting parameters for gmx energy and density
	BEGIN_FRAME	=	cmdParam.begin_frame
	END_FRAME	=	cmdParam.end_frame
	STRIDE		=	cmdParam.stride
	RADIUS		=	cmdParam.radius
	
	#Options to compute
	LOCAL		=	cmdParam.local
	DENSITY		=	cmdParam.density
	GRID		=	cmdParam.grid
	SELECT		=	cmdParam.select
	
	#OPTIONS FOR LIPIDS
	HEADGROUP	= cmdParam.headgroup
	TAILGROUP	= cmdParam.tailgroup
	SUGROUP		= cmdParam.sugroup
	
	filesFound = {}
	pad_x = pad_y = pad_z = 10
	if MD_RUN_INPUT is not None and FOLDER is not None:
		"""
		These options enables membrane property computation for a lot of data
		If you want only a peculiar file set it using -c and not -r and -f
		"""
		
		if FOLDER[-1] == '/': FOLDER = FOLDER[:-1]
		
		
		for Run in MD_RUN_INPUT:
			for path in glob.glob(FOLDER +'/**/*'+ Run +'*'+EXTENSION, recursive=True):
				if 'fatslimAnalysis' in path: continue
			
				path_test = path.replace(FOLDER,'').split('/')
				path = path.split('/')
				filename = '/'.join(path)
				
				project_name = None
				job_name = None
				
				if len(path_test) >= 3 : 
					project_name = path[-3]
					job_name = path[-2]
				else:
					project_name = path[-3]
					job_name = path[-2]
				
				if 'JOB' in job_name and 'OUTPUT' in job_name:
					job_name = path[-3]
					project_name = path[-4]
					
				if project_name not in filesFound: 
					filesFound[project_name] = { job_name: filename }
					
				else:
					if job_name not in filesFound[project_name]:
						filesFound[project_name][job_name] = filename
		
		DEPTH = depth(filesFound)
		print(filesFound)
		
		#Creates the base directory
		analysisFolder = FOLDER+'/MD_Analysis'
		#print(analysisFolder)
		if not os.path.isdir(analysisFolder):
			os.makedirs(analysisFolder, exist_ok=True)
		else:
			bakCount = 0
			backupExt = '-bak'
			
			for path in glob.glob(FOLDER +'/*MD_Analysis*'):
				print(path)
				tmp = path.split('/')[-1]
				
				if backupExt in tmp:
					bakCount += 1
					
			
			backupFile = analysisFolder + backupExt+str(bakCount)
			
			try:
				os.rename(analysisFolder, backupFile)
			except (OSError, IOError) as e:
				print(sys.stderr, 'Error moving {0} to {1}: {2}'.format(analysisFolder, backupFile, e))
				
			os.makedirs(analysisFolder, exist_ok=True)
		
		
		# Directory for local temperature
		if LOCAL is not None:
			loctemp_dir = analysisFolder+'/MD_LOC_TEMP'
			os.makedirs(loctemp_dir)
		
		## Directory for local pressure
		#if LOCAL_PRES:
			locpres_dir = analysisFolder+'/MD_LOC_PRES'
			os.makedirs(locpres_dir)
			
		# Directory for density
		if DENSITY is not None:
			density_dir = analysisFolder+'/MD_DENSITY'
			os.makedirs(density_dir)
		
		# Directory for gmx density
		if SELECT:
			select_dir = analysisFolder+'/MD_SELECT'
			os.makedirs(select_dir)
		
		
		if DEPTH >= 3:
			for project_name in filesFound:
				for job_name, file_name in filesFound[project_name].items():
					

					xtc_file = None
					edr_file = None
					tpr_file = None
					ndx_file = None
					destination = None
					
					plotname = None
					csvname = None
					
					name = file_name.split('/')[-1]
					path = '/'.join(file_name.split('/')[:-1])
					
					xtc_file = file_name
					edr_file = file_name.replace(EXTENSION, 'edr')
					tpr_file = file_name.replace(EXTENSION, 'tpr')
					
					begin_end = ""
					
					trr_file = file_name.replace('xtc','trr')
					is_trr_file = os.path.isfile(trr_file)
					
					coord		= mda.Universe(tpr_file, xtc_file, convert_units=False)
					if is_trr_file:
						vel_force	= mda.Universe(tpr_file, trr_file, convert_units=False)
					
					# Array to get the types of atoms
					atom_types = {}
					multiple_masses = False
					atom_name_field = {'W': 1, 'NC3': 2, 'PO4':3, 'GL1':4, 'GL2':5, 'C1A':6, 'C2A':7, 'C3A':8, 'C4A':9, 'C5A':10,
																					'C1B':11, 'C2B':12, 'C3B':13, 'C4B':14, 'C5B':15,
																					'DEF':0, 'SUN': 50, 'SUP':51}
					
					lipid_groups_field = None
					if HEADGROUP is not None and TAILGROUP is not None:
						lipid_groups_field = atom_name_field
						
						for name in lipid_groups_field:
							if name in HEADGROUP:
								lipid_groups_field[name] = 2
								
							if name in TAILGROUP:
								lipid_groups_field[name] = 3
					else:
						lipid_groups_field = atom_name_field
						
					if SUGROUP is not None:
						for name in lipid_groups_field:
							if name in SUGROUP:
								lipid_groups_field[name] = 4
					
					print(lipid_groups_field)
					for atom in coord.atoms:
						if atom.name not in atom_types:
							atom_types.update( { atom.name : { 'mass' : atom.mass , 'field': float(lipid_groups_field[atom.name]) } } )
							
					
					
					nb_frames	= len(coord.trajectory)
					nb_atoms	= len(coord.atoms)
					
					print(job_name)
				
					if LOCAL is not None:
						
						grid_x_array		= []
						grid_y_array		= []
						grid_z_array		= []
						
						grid_array_temp		= []
						grid_array_pres		= []
						
						
						vel_for_pos = zip(vel_force.trajectory[BEGIN_FRAME:END_FRAME:STRIDE], coord.trajectory[BEGIN_FRAME:END_FRAME:STRIDE])
						for index, vfp in enumerate(vel_for_pos):
							
							print("Computing for timestep = {0}".format(vel_force.atoms.ts.time))
							#Creating the box and grid
							box_xx, box_yy, box_zz, box_xy, box_xz, box_yz = vel_force.atoms.ts.dimensions
							
							dx = box_xx/LOCAL[0]
							dy = box_yy/LOCAL[1]
							dz = box_zz/LOCAL[2]
							
							dx = complex(0,dx)
							dy = complex(0,dy)
							dz = complex(0,dz)
							
							grid_x, grid_y, grid_z = np.mgrid[0:box_xx:dx, 0:box_yy:dy, 0:box_zz:dz]
							#print(np.shape(grid_x),np.shape(grid_y),np.shape(grid_z))
							
							# Computing the local temperature using velocities
							sqrd = np.power(vel_force.atoms.velocities,2)
							summed = np.sum(sqrd, axis=1)
							
							#If atoms have different masses
							if multiple_masses:
								masses = [atom.mass for atom in vel_force.atoms]
								summed = np.multiply(summed, masses)
							else:
								summed *= 72.0
								
							# second term on right hand side to convert to Kelvin/mass
							in_kelvin = summed * 40.09074226350679338
							
							
							# Computing the local pressure using velocities, forces and coordinates.
							part_volume		= box_xx * box_yy * box_zz / nb_atoms
							pos_dot_force	= np.sum( np.multiply(vel_force.atoms.atoms.forces, coord.atoms.positions), axis=1 )
							pressure		= summed/part_volume + pos_dot_force/part_volume
							pressure		= pressure / 16.6053886
							
							points		= coord.atoms.positions
							temp_values	= np.array(in_kelvin)
							pres_values	= np.array(pressure)

							grid_x_array.append(grid_x)
							grid_y_array.append(grid_y)
							grid_z_array.append(grid_z)
							
							grid_array_temp.append(griddata(points, temp_values, (grid_x, grid_y, grid_z) , method='nearest'))
							grid_array_pres.append(griddata(points, pres_values, (grid_x, grid_y, grid_z) , method='nearest'))
						
						
						#print(np.shape(grid_x_array))
						#print(np.shape(grid_y_array))
						#print(np.shape(grid_z_array))
						mean_grid_x		= np.mean(grid_x_array,axis=0)
						mean_grid_y		= np.mean(grid_y_array,axis=0)
						mean_grid_z		= np.mean(grid_z_array,axis=0)
						mean_grid_temp	= np.mean(grid_array_temp,axis=0)
						mean_grid_pres	= np.mean(grid_array_pres,axis=0)
						
						loctemp_npz_file = job_name + '-' + xtc_file.split('-')[-1].replace('.xtc', '_LOCTMP.npz')
						locpres_npz_file = job_name + '-' + xtc_file.split('-')[-1].replace('.xtc', '_LOCPRES.npz')
						np.savez_compressed(loctemp_dir+'/{0}'.format(loctemp_npz_file), mgx=mean_grid_x, mgy=mean_grid_y, mgz=mean_grid_z, mg=mean_grid_temp)
						np.savez_compressed(locpres_dir+'/{0}'.format(locpres_npz_file), mgx=mean_grid_x, mgy=mean_grid_y, mgz=mean_grid_z, mg=mean_grid_pres)
						
						
					if DENSITY is not None:
						dens_start_time = time.time()
						
						frame_density		= []
						grid_x_array		= []
						grid_y_array		= []
						grid_z_array		= []
						
						grid_array_density	= []
						
						#For dataframe
						atom_names = [atom_type for atom_type in atom_types]
						
						for index, p in enumerate(coord.trajectory[BEGIN_FRAME:END_FRAME:STRIDE]):
							
							print("Computing for timestep = {0}".format(coord.atoms.ts.time))
							
							#Creating the box and grid
							box_xx, box_yy, box_zz, box_xy, box_xz, box_yz = coord.atoms.ts.dimensions
							#print(vel_force.atoms.ts)
							
							dx = dy = dz = None
							if GRID is not None:
								dx = dy = dz = GRID#math.floor(box_zz/DENSITY[2])
								dx = complex(0,dx)
								dy = complex(0,dy)
								dz = complex(0,dz)
								grid_x, grid_y, grid_z = np.mgrid[0:box_xx:dx, 0:box_yy:dy, 0:box_zz:dz]
							
							DZ = box_zz/DENSITY
							
							
							
							##Array keeping all the data per frame
							points = np.empty(shape=(1,3))
							values = []
							
							if RADIUS is not None:
								slice_volume = (box_xx * box_yy - math.pi*RADIUS**2) * DZ
								#DEF = coord.select_atoms("name DEF")
								#DEF_COG = DEF.center_of_geometry()
								#z_DEF = DEF_COG[2]
							else:
								slice_volume = (box_xx * box_yy) * DZ
							
								
							density_per_name = {}
							
							for atom_type in atom_types:
								atom_names.append(atom_type)
								selection = None
								if RADIUS is not None and atom_type != 'DEF':
									selection	= coord.select_atoms("name {0} and not (cyzone {1} {2} -{2} (name DEF))".format(atom_type, RADIUS, box_zz), update=True)
									#(box_zz-0.5)/2.0, (box_zz)/2.0))
									
								else:
									selection	= coord.select_atoms("name {0}".format(atom_type))
								
								z = 0.0
								z_bin = []
								nb_atoms = []
								while( z < box_zz ):
									condition = (selection.positions[:,2] >= z) * (selection.positions[:,2] < z + DZ)
									
									indexes = np.where(condition)
									z_bin.append(z+DZ/2.)
									
									nb_atoms.append(len(indexes[0]) / slice_volume)
									
									z += DZ
								
								if 'zbin' not in density_per_name:
									density_per_name.update({'zbin':z_bin})
								
								density_per_name.update({atom_type : nb_atoms})
								points = np.append(points, selection.atoms.positions, axis=0)
								
								if GRID is not None:
									values.extend([atom_types[atom_type]['field']]*len(selection.atoms.positions))
							
							# Creating a dataframe dor density and adding it to the frames
							density_per_name = pd.DataFrame(density_per_name)
							frame_density.append(density_per_name)
							
							if GRID is not None:
								if RADIUS is not None:
									selection	= coord.select_atoms("cyzone {1} {2} -{2} (name DEF)".format(atom_type, RADIUS, box_zz), update=True)
									#(box_zz-0.5)/2.0, (box_zz)/2.0))
									points		= np.append(points, selection.atoms.positions, axis=0)
									values.extend([0.0]*len(selection.atoms.positions))
									
								points = np.delete(points,0,0)
								values = np.array(values)

								grid = griddata(points, values, (grid_x, grid_y, grid_z), method='nearest')
								grid_x_array.append(grid_x)
								grid_y_array.append(grid_y)
								grid_z_array.append(grid_z)
								grid_array_density.append(grid)
						
						if GRID is not None:
							mean_grid_x			= np.mean(grid_x_array,axis=0)
							mean_grid_y			= np.mean(grid_y_array,axis=0)
							mean_grid_z			= np.mean(grid_z_array,axis=0)
							mean_grid_density	= np.mean(grid_array_density,axis=0)
							
							density_npz_file = job_name + '-' + xtc_file.split('-')[-1].replace('.xtc', '_DENS.npz')
							np.savez_compressed(density_dir+'/{0}'.format(density_npz_file), mgx=mean_grid_x, mgy=mean_grid_y, mgz=mean_grid_z, mg=mean_grid_density)
						
						df_concat = pd.concat(frame_density)
						mean = df_concat.groupby(level=0).mean()
						
						if HEADGROUP is not None:
							heads = intersect(HEADGROUP, list(mean))
							mean['head'] = mean[heads].sum(axis=1)
						if TAILGROUP is not None:
							tails = intersect(TAILGROUP, list(mean))
							mean['tail'] = mean[tails].sum(axis=1)
						if SUGROUP is not None:
							su = intersect(SUGROUP, list(mean))
							mean['su'] = mean[su].sum(axis=1)
							
						# get a list of columns
						cols = list(mean)
						# move the column to head of list using index, pop and insert
						cols.insert(0, cols.pop(cols.index('zbin')))
						# use ix to reorder
						mean = mean.ix[:, cols]
						
						# FAIRE LA MOYENNE AVEC PANDA SUR CHAQUE BIN ! plus incertitude !!!!
						header = """
								@    title ""
								@    xaxis  label "z-coordinate (nm)"
								@    yaxis  label ""
								@TYPE xy
								@ view 0.15, 0.15, 0.75, 0.85
								@ legend on
								@ legend box on
								@ legend loctype view
								@ legend 0.78, 0.8
								@ legend length 2\n
								"""
						header = ut.RemoveUnwantedIndent(header)
						with open(density_dir+'/'+job_name + '-' + xtc_file.split('-')[-1].replace('.xtc', '_MDADENS.xvg'),'a+') as xvg_out:
							for index, atom_name in enumerate(mean):
								if atom_name != "zbin":
									header += """@ s{0} legend "{1}"\n""".format(index-1, atom_name)
								
							for row in mean.itertuples():
								header += '   '.join(map(str, row[1:]))+'\n'
							
							xvg_out.write(header+"\n")
							
						print( "--- density for {0} done in {1:f} seconds ---".format(job_name, time.time() - dens_start_time) )
	
	print( "--- mda done in {0:f} seconds ---".format(time.time() - start_time) )
				

elif 'reflectometry' in sys.argv:
	"""
		Function to computer the scattering length density (SLD) of the molecules defined in .sl file
		and the compute the the reflectometry curve.
	"""
	try:
		import numpy as np
		import pandas as pd
		import matplotlib as mpl
		from matplotlib import rcParams
		from matplotlib import rc
		import matplotlib.ticker as ticker #to change the ticks
		
		#print(rcParams.keys())
		#mpl.use('pgf')
		params = {
		'pgf.texsystem': 'xelatex',        # change this if using xetex or lautex
		'font.size': 16,
		'figure.figsize': [8,4],
		'axes.labelsize': 20,
		'legend.fontsize': 20,
		'xtick.labelsize': 16,
		'ytick.labelsize': 16,
		'text.usetex': True,
		'font.family': 'sans-serif',
		"text.latex.preamble": [
			r"\usepackage[utf8]{inputenc}",    # use utf8 fonts 
			r"\usepackage[T1]{fontenc}",        # plots will be generated
			#r"\usepackage{anyfontsize}",
			#r"\fontsize{16}{19.20}",
			#r"\DeclareMathSizes{10}{10}{10}{10}",
			r"\usepackage[detect-all]{siunitx}",         # load additional packages
			r"\usepackage{amsmath}",
			r"\usepackage{amssymb}",
			r"\sisetup{mode = math,math-rm=\mathsf}"
			#r"\usepackage{unicode-math}",  # unicode math setup
			#r"\setmathfont{xits-math.otf}",
			]
}
			
		rcParams.update(params)
		import matplotlib.pyplot as plt
		from scipy.interpolate import CubicSpline
	except ImportError as e:
		print("Error {0}".format(e))
	
	SL_FILE			= cmdParam.sl_file
	XVG_FILE		= cmdParam.xvg_file
	
	MOLECULES		= cmdParam.molecules
	SUBLAYERS		= cmdParam.sublayers
	SUPLAYERS		= cmdParam.superlayers
	LIMITS			= cmdParam.limits
	
	REFLECTOMETRY	= cmdParam.reflectivity
	Q_RANGE			= cmdParam.qrange
	EXPERIMENTAL	= cmdParam.experimental
	RESOLUTION		= cmdParam.resolution
	
	CHECK			= cmdParam.check
	sl_data = {}
	
	#Reading data from .sl file
	with open(SL_FILE,'r') as sl:
		molecule = None
		for line in sl:
			sl_atom = line.strip().split()
			if '[' in line and ']' in line:
				molecule = line.strip()[1:-2].strip()
				sl_data.update({molecule : {}})
			elif len(sl_atom) == 2:
				atom_name		= sl_atom[0]
				sl_atom_name	= sl_atom[1]
				sl_data[molecule][atom_name] = float(sl_atom_name)
			elif len(sl_atom) == 1:
				sl_atom_name	= sl_atom[0]
				sl_data.update({molecule : float(sl_atom_name)})
				
	#print(sl_data)
	
	#Reading data from the .xvg file
	density_data = {}
	density_data.update( { 'z' : 0} )
	col_ndx = 1
	with open(XVG_FILE,'r') as data:
		for line in data:
			if line.startswith('@'):
				#Extract the legend to get the Column/quantity match
				if 's{0}'.format(col_ndx-1) in line:
					atom_name = line[  line.find('"')+1 : line.rfind('"')  ]
					density_data.update( { atom_name : col_ndx} )
					col_ndx += 1
				
				continue
			
			elif line.startswith('#'):
				pass
			
			else:
				break
	
	output_str = ""
	with open(XVG_FILE,'r') as inputFile:
		for line in inputFile:
			if '#' not in line:
				if '@' not in line:
					output_str += line
	
	readFrom = StringIO(output_str)
	readFrom.seek(0)
	
	for atom_name in density_data:
		density_data[atom_name] = np.loadtxt(readFrom, dtype='float', usecols = (int(density_data[atom_name]),) )
		readFrom.seek(0)
	
	density_data = pd.DataFrame(density_data)
	
	for mol in MOLECULES:
		for atom_name in  sl_data[mol]:
			# converting to (10^-6 A^-2) with 1e-8
			density_data[atom_name] = density_data[atom_name]* sl_data[mol][atom_name]*1e-8 
	
	all_atoms = []
	for mol in MOLECULES:
		all_atoms.extend(sl_data[mol])
	
	total = intersect(all_atoms, list(density_data))
	density_data['total'] = density_data[total].sum(axis=1)
	
	# get a list of columns
	cols = list(density_data)
	# move the column to head of list using index, pop and insert
	cols.insert(0, cols.pop(cols.index('z')))
	# use ix to reorder
	density_data = density_data.ix[:, cols]
	
	if LIMITS is not None:
		lower_limit = density_data['z'] >= LIMITS[0]
		upper_limit = density_data['z'] <= LIMITS[1]
	
	density_data = density_data[lower_limit & upper_limit]
	
	dz = None
	if SUBLAYERS is not None or SUPLAYERS is not None:
		dz = density_data['z'].diff().mean()
		
	sub_layers_dict = {}
	sup_layers_dict = {}
	
	print(density_data)
	
	if SUBLAYERS is not None:
		
		sub_list = []
		with_density = False
		if len(SUBLAYERS) % 3 == 0:
			for i in range(0, len(SUBLAYERS), 3):
				sub_list.append(SUBLAYERS[i:i+3])

			for index, sub in enumerate(sub_list):
				mol, thick, rug = sub
				#assert(mol not in sl_data), "The molecule you choosed to add below the sample is not in your .sl data file"
				sub_layers_dict.update({ index : {'molecule': mol,'thickness': float(thick),
									  								'rugosity': float(rug)}})
		elif len(SUBLAYERS) % 4 == 0:
			for i in range(0, len(SUBLAYERS), 4):
				sub_list.append(SUBLAYERS[i:i+4])

			for index, sub in enumerate(sub_list):
				mol, thick, rug, dens = sub
				molecule, atom = mol.split('_')
				#assert(mol not in sl_data), "The molecule you choosed to add below the sample is not in your .sl data file"
				sub_layers_dict.update({ index : {'molecule': molecule, 'atom': atom,
									  								'thickness': float(thick),
									  								'rugosity': float(rug),
									  								'density': float(dens)
									  								}})
			with_density = True
		
		min_z	= density_data['z'][0]
		
		total_value_at_sup	= density_data['total'][0]
		total_new_z_bin		= []
		total_new_sld		= []
		
		for index in range(0, len(sub_layers_dict)):
			molecule	= sub_layers_dict[index]['molecule']
			thickness	= sub_layers_dict[index]['thickness']
			rugosity	= sub_layers_dict[index]['rugosity']
			
			sld = None
			if not with_density:
				sld		= sl_data[molecule]*1e-6 # in (10^-6 A^-2)
			else:
				density	= sup_layers_dict[index]['density']
				atom	= sup_layers_dict[index]['atom']
				sld		= sl_data[molecule][atom]*density*1e-8
			
			number_bin		= None
			rough_bin		= None
			
			new_z_bin		= None
			sld_values		= None
			rough_interface = None
			
			if rugosity == 0.0:
				number_bin	= 2#int(thickness / dz)
				new_z_bin	= np.linspace(min_z, -thickness + min_z, number_bin)
				min_z		= -thickness + min_z
				
				sld_values = np.empty(number_bin)
				sld_values.fill(sld)
				
				# Updating the value at sup layer
				total_value_at_sup	= sld
			
			else:
				# values to interpolate for layers roughness
				x	= [min_z-rugosity, min_z-(rugosity/2.), min_z]
				y	= [sld, (sld+total_value_at_sup)/2., total_value_at_sup]
				
				
				# number of bins in the rough part and the layer
				rough_interp	= CubicSpline(x,y, extrapolate=True, bc_type='clamped')
				number_bin		= 2#int(thickness / dz)
				rough_bin		= int(rugosity / dz)
				
				#Creating the bins for nex layers with roughness
				rough_interface	= np.linspace(min_z, min_z-rugosity, 10)[1:-1]
				new_z_bin		= np.linspace(min_z-rugosity,-thickness+min_z, number_bin)
				
				#Interpolation on the bins
				rough_sld		= rough_interp(rough_interface)
				#Changing the minimum z values for next loop
				min_z			= -thickness + min_z
			
				sld_values	= np.empty(number_bin)
				sld_values.fill(sld)
				sld_values	= np.append(rough_sld,sld_values)
				
				new_z_bin	= np.append(rough_interface, new_z_bin)
				
				# Updating the value at sup layer
				total_value_at_sup	= sld
				
			total_new_z_bin	= np.append(total_new_z_bin, new_z_bin)
			total_new_sld	= np.append(total_new_sld, sld_values)
		
		# Removes the first value to avoid repetition
		total_new_z_bin	= total_new_z_bin[1:]
		total_new_sld	= total_new_sld[1:]
		# Reverse the curves to have the right z-direction
		total_new_z_bin	= total_new_z_bin[::-1]
		total_new_sld	= total_new_sld[::-1]
		
		#Create a new total with the added layers
		extract_total = density_data[['z','total']]
		new_total = pd.DataFrame({'z': total_new_z_bin, 'total': total_new_sld})
		new_total = pd.concat([new_total, extract_total], ignore_index=True)
		new_total = new_total.rename(columns={'total' : 'fittotal'})
		
		val_to_add = new_total.shape[0] - density_data.shape[0]
		density_data = density_data.append(pd.DataFrame([0.0]*val_to_add), ignore_index=True)
		density_data = density_data.shift(periods=val_to_add)
		density_data = density_data.assign(fittotal=new_total.fittotal, z=new_total.z)
	
	print(density_data)
	#density_data.plot(x='z', y='fittotal')
	#plt.show()
	
	if SUPLAYERS is not None:
		
		sup_list = []
		with_density = False
		if len(SUPLAYERS) % 3 == 0:
			for i in range(0, len(SUPLAYERS), 3):
				sup_list.append(SUPLAYERS[i:i+3])

			for index, sub in enumerate(sup_list):
				mol, thick, rug = sub
				#assert(mol not in sl_data), "The molecule you choosed to add below the sample is not in your .sl data file"
				sup_layers_dict.update({ index : {'molecule': mol,'thickness': float(thick),
									  								'rugosity': float(rug)}})
		elif len(SUPLAYERS) % 4 == 0:
			for i in range(0, len(SUPLAYERS), 4):
				sup_list.append(SUPLAYERS[i:i+4])

			for index, sub in enumerate(sup_list):
				mol, thick, rug, dens = sub
				molecule, atom = mol.split('_')
				#assert(mol not in sl_data), "The molecule you choosed to add below the sample is not in your .sl data file"
				sup_layers_dict.update({ index : {'molecule': molecule, 'atom': atom,
									  								'thickness': float(thick),
									  								'rugosity': float(rug),
									  								'density': float(dens)
									  								}})
			with_density = True
		
		min_z	= density_data['z'][density_data.index[-1]]
		
		total_value_at_sub	= density_data['total'][density_data.index[-1]]
		total_new_z_bin		= []
		total_new_sld		= []
		
		for index in range(0, len(sup_layers_dict)):
			molecule	= sup_layers_dict[index]['molecule']
			thickness	= sup_layers_dict[index]['thickness']
			rugosity	= sup_layers_dict[index]['rugosity']
			
			sld = None
			if not with_density:
				sld		= sl_data[molecule]*1e-6 # in (10^-6 A^-2)
			else:
				density	= sup_layers_dict[index]['density']
				atom	= sup_layers_dict[index]['atom']
				sld		= sl_data[molecule][atom]*density*1e-8
			
			number_bin		= None
			rough_bin		= None
			
			new_z_bin		= None
			sld_values		= None
			rough_interface = None
			
			if rugosity == 0.0:
				number_bin	= 2 #int(thickness / dz)
				new_z_bin	= np.linspace(min_z, thickness+min_z, number_bin)
				min_z		= thickness + min_z
				
				sld_values = np.empty(number_bin)
				sld_values.fill(sld)
				
				# Updating the value at sup layer
				total_value_at_sub	= sld
			
			else:
				# values to interpolate for layers roughness
				x	= [min_z, min_z+(rugosity/2.), min_z+rugosity]
				y	= [total_value_at_sub, (sld+total_value_at_sub)/2., sld]
				
				
				# number of bins in the rough part and the layer
				rough_interp	= CubicSpline(x,y, extrapolate=True, bc_type='clamped')
				number_bin		= 2 #int(thickness / dz)
				rough_bin		= int(rugosity / dz)
				
				#Creating the bins for nex layers with roughness
				rough_interface	= np.linspace(min_z, min_z+rugosity, 10)[-1:1]
				new_z_bin		= np.linspace(min_z+rugosity,thickness+min_z, number_bin)
				
				#Interpolation on the bins
				rough_sld		= rough_interp(rough_interface)
				
				#Changing the minimum z values for next loop
				min_z			= thickness + min_z
			
				sld_values	= np.empty(number_bin)
				sld_values.fill(sld)
				sld_values	= np.append(rough_sld,sld_values)
				
				new_z_bin	= np.append(rough_interface, new_z_bin)
				
				# Updating the value at sup layer
				total_value_at_sub	= sld
				
			total_new_z_bin	= np.append(total_new_z_bin, new_z_bin)
			total_new_sld	= np.append(total_new_sld, sld_values)
		
		total_new_z_bin	= total_new_z_bin[1:]
		total_new_sld	= total_new_sld[1:]
		
		#Create a new total with the added layers
		new_total = None
		if SUBLAYERS is not None:
			extract_total	= density_data[['z','fittotal']]
			new_total		= pd.DataFrame({'z': total_new_z_bin, 'fittotal': total_new_sld})
			new_total		= pd.concat([extract_total, new_total], ignore_index=True)
			
			val_to_add = new_total.shape[0] - density_data.shape[0]
			density_data = density_data.append(pd.DataFrame([0.0]*val_to_add), ignore_index=True)
			
		else:
			extract_total	= density_data[['z','total']]
			new_total		= pd.DataFrame({'z': total_new_z_bin, 'total': total_new_sld})
			new_total		= pd.concat([extract_total, new_total], ignore_index=True)
			new_total		= new_total.rename(columns={'total' : 'fittotal'})
			
			val_to_add = new_total.shape[0] - density_data.shape[0]
			density_data = density_data.append(pd.DataFrame([0.0]*val_to_add), ignore_index=True)
		
		
		density_data = density_data.assign(fittotal=new_total.fittotal, z=new_total.z).drop(0,1)
	
	with pd.option_context('display.max_rows', None):
		print(density_data)
	#print(density_data)
	#density_data.plot(x='z', y='fittotal')
	#density_data.plot(x='z', y='total')
	#plt.show()
	
	# get a list of columns
	cols = list(density_data)
	# move the column to head of list using index, pop and insert
	cols.insert(0, cols.pop(cols.index('z')))
	# use ix to reorder
	density_data = density_data.ix[:, cols]
	
	
	#density_data.fillna(0.0, inplace=True)
	header = """
			@    title ""
			@    xaxis  label "z (nm)"
			@    yaxis  label "sld (A-1)"
			@TYPE xy
			@ view 0.15, 0.15, 0.75, 0.85
			@ legend on
			@ legend box on
			@ legend loctype view
			@ legend 0.78, 0.8
			@ legend length 2\n
			"""
	
	if SUPLAYERS is not None:
		adding_to_header = """
			#SLD fittotal was made using the following parameters for sub layers:
			# {0}
			#""".format(str(SUPLAYERS))
		
		header = adding_to_header + header
	
	if SUBLAYERS is not None:
		adding_to_header = """
			#SLD fittotal was made using the following parameters for super layers:
			# {0}
			#""".format(str(SUBLAYERS))
		
		header = adding_to_header + header
		
	
	header = ut.RemoveUnwantedIndent(header)
	
	xvg_out = None
	if XVG_FILE.endswith('MDADENS.xvg'):
		xvg_out = open(XVG_FILE.replace('MDADENS.xvg', 'MDASLD.xvg'),'w')
		
	if XVG_FILE.endswith('DENS.xvg'):
		xvg_out = open(XVG_FILE.replace('DENS.xvg', 'SLD.xvg'),'w')

	for index, atom_name in enumerate(density_data):
		if atom_name != "z":
			header += """@ s{0} legend "{1}"\n""".format(index-1, atom_name)
		
	for row in density_data.itertuples():
		header += '   '.join(map(str, row[1:]))+'\n'
	
	xvg_out.write(header+"\n")
	xvg_out.close()
	
	if CHECK:
		if 'fittotal' in density_data:
			density_data.plot(x='z', y='fittotal', style=['b'], linewidth=3.0, label='Sim. SLD')
		else:
			density_data.plot(x='z', y='total')
		
		y_ticks	= plt.yticks()[0]
		y_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(y*1e6,2)) for y in y_ticks]
		plt.yticks(y_ticks, y_ticks_label)
		#plt.ylim(0.01e-8, 5.0e-8)
		
		x_ticks	= plt.xticks()[0]
		x_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(x,2)) for x in x_ticks]
		plt.xticks(x_ticks, x_ticks_label)
		#plt.xlim(0.0, 0.25)
		
		
		plt.ylabel(r"$ \mathsf{SLD\ ( \times \SI{E-06}{\angstrom^{-2} } ) } $")
		plt.xlabel(r"$\mathsf{ z\ (\SI{}{\angstrom}) }$")
		#legend = plt.legend()
		#legend.get_texts()[0].set_text('Sim. SLD')
		#legend.draggable(True)
		plt.tight_layout()
		plt.show()
			
	
	if REFLECTOMETRY is not None:
		try:
			import numpy.matlib as npm
			import cmath as cm
		except ImportError as e:
			print("Error {0}".format(e))
			
		q_min	= None
		q_max	= None
		q_grid	= None
		q_range	= None
		
		PI = math.pi
		if EXPERIMENTAL is not None:
			if EXPERIMENTAL.endswith('.csv'):
				q_range	= np.loadtxt(EXPERIMENTAL, dtype='float', usecols = (0,) )
				Rq4_exp	= np.loadtxt(EXPERIMENTAL, dtype='float', usecols = (1,) )
				Rq4_std			= None
				
				if 'data' in EXPERIMENTAL:
					Rq4_std = np.loadtxt(EXPERIMENTAL, dtype='float', usecols = (2) )
					
		elif Q_RANGE is not None:
			q_min, q_max, q_grid = Q_RANGE
			q_range	= np.linspace(q_min, q_max, q_grid)
		
		sld_0	= density_data[REFLECTOMETRY].iloc[0]
		
		# Array for the reflectometry data
		R	= []
		for q in q_range:
			#Initial wave vector and phase factor 
			k_0	= complex(q/2., 0)
			b_0	= complex(0,0)
			
			# Initial matrix before multiplication (diagonal)
			M = npm.matrix([[complex(1,0),complex(0,0)],[complex(0,0), complex(1,0)]], 
								dtype=complex)
			
			# Wave vector and phase factor for the previous layer (start)
			# Changed in loop
			k_j				= k_0
			phase_factor_j	= b_0
			
			for i in range(1, len(density_data.index)):
				# prec layer
				j = i-1
				# Computing the wave vector of the ith layer using SLD value
				sld_i	= density_data[REFLECTOMETRY].iloc[i]
				#sld_j	= density_data[REFLECTOMETRY].iloc[j]
				k	= (k_0.real)**2 - 4*PI*(sld_i - sld_0)
				
				k_i	= None
				if k >= 0:
					k_i = complex(math.sqrt(k), 0)
				else:
					k_i = complex(0, math.sqrt(-k))
				
				# Computing the phase factor for the current layer
				layer_thick		= abs(density_data['z'].iloc[i] - density_data['z'].iloc[j]) * 10.
				#print("index: ",i)
				#print("z_i:",density_data['z'].iloc[i])
				#print("z_j:",density_data['z'].iloc[j])
				#print("thick: ", layer_thick)
				phase_factor_i	= k_i * complex(layer_thick,0)
				
				#print(layer_thick)
				# Computing the fresnel coefficent
				# No roughtness as already included in the sld (might need to add them though)
				k_sub		= k_j - k_i
				k_sum		= k_j + k_i
				k_prod		= k_j * k_i
				
				#taking into account roughness
				pre_exp		= complex(-2,0)
				exp_k_prod	= cmath.exp(pre_exp * k_prod * 1)
				
				R_fresnel	= (k_sub / k_sum) #* exp_k_prod
				
				# Matrix for reflection between the 2 layers
				ib			= complex(0,1) * phase_factor_j
				ibneg		= complex(-1,0) * ib
				exp_ib		= cmath.exp(ib)
				exp_ibneg	= cmath.exp(ibneg)
				
				C_mat	= npm.matrix([[exp_ib, R_fresnel * exp_ib],[R_fresnel * exp_ibneg, exp_ibneg]],
										dtype=complex)
				
				M = M @ C_mat
				
				# Putting the top k and phase factor for the next bottom one 
				k_j				= k_i
				phase_factor_j	= phase_factor_i
				
			
			M11		= M[0,0]
			M11conj	= M11.conjugate()
			M21		= M[1,0]
			M21conj	= M21.conjugate()
			
			Rm_q	= (M21 * M21conj)/(M11 * M11conj) 
			Rm_q	= Rm_q.real
			
			R.append(Rm_q)
			
		
		q4	= np.power(q_range,4)
		Rq4	= np.multiply(R,q4)
		reflectivity_data	= pd.DataFrame({'R': R, 'Rq4': Rq4,'q': q_range})
		
		Chi = None
		if EXPERIMENTAL is not None:
			reflectivity_data	= reflectivity_data.assign(Rq4_exp = Rq4_exp)
			if Rq4_std is not None:
				reflectivity_data	= reflectivity_data.assign(Rq4_std = Rq4_std)
			
			#Computing Chi
			diff		= Rq4_exp - reflectivity_data['Rq4']
			diff_sqrd	= np.power(diff,2)
			sig_sqrd	= np.power(Rq4_std,2)
			num			= np.sum(np.divide(diff_sqrd,sig_sqrd))
			denum		= len(Rq4_exp) -1
			
			Chi_sqrd	= num/denum
			Chi			= np.sqrt(Chi_sqrd)
			
		
		if RESOLUTION is not None:
			# Formula from DOI: 10.1021/acs.jpcb.6b05433#
			p = RESOLUTION
			
			nb_points			= np.arange(-p,p+1,1)
			term_in_exp			= -2 * np.power(nb_points/p ,2 )
			weight_terms		= np.exp(term_in_exp)
			normalising_term	= np.sum(weight_terms)
			
			
			norm_weights		= weight_terms / normalising_term
			
			
			R_prime	= []
			
			R		= reflectivity_data['R']
			q_range	= reflectivity_data['q']
			#q_delta	= (q_max - q_min)/q_grid
			
			R_interp	= CubicSpline(q_range,R, extrapolate=True, bc_type='clamped')
			
			values_interp	= np.multiply(R_interp(q_range),q4)
			reflectivity_data	= reflectivity_data.assign(interp=values_interp)
			
			
			q_prec	= 0.0
			for q in q_range:
				Rpq	= 0.0
				dq	= 0.02 #(q - q_prec)/ q#q_delta/q
				#print(dq)
				for a, w in zip(nb_points, norm_weights):
					#print(a)
					Rpq += w * R_interp(q + a/p * dq)
				
				R_prime.append(Rpq)
				q_prec	= q
				
			
			Rpq4	= np.multiply(R_prime,q4)
			reflectivity_data	= reflectivity_data.assign(Rpq4=Rpq4)
			
			
			reflectivity_data.plot(x='q', y=['Rq4','Rpq4', 'interp'])
			plt.show()
			
		
		if CHECK:
			with pd.option_context('display.max_rows', None):
				print(reflectivity_data)
			
			if EXPERIMENTAL is not None:
				
				
				reflectivity_data.plot(x='q', y=['Rq4'], style=['b'], logy=True, linewidth=3.0, label='Sim.')
				if Rq4_std is not None:
					plt.errorbar('q', 'Rq4_exp', yerr='Rq4_std', data=reflectivity_data, fmt='ro', label='Exp.')
					
					anot = plt.annotate(r'$\chi = {0}$'.format(round(Chi,2)), xy=(0.10, 2e-8), xytext=(0.10, 3e-8))
			else:
				reflectivity_data.plot(x='q', y='Rq4')
			
			y_ticks	= plt.yticks()[0]
			y_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(y*1e8,2)) for y in y_ticks]
			plt.yticks(y_ticks, y_ticks_label)
			maximum = np.amax(y_ticks) #5.0e-8 
			plt.ylim(0.01e-8, maximum)
			
			x_ticks	= plt.xticks()[0]
			x_ticks_label	= [r'$\mathsf{{ {0} }}$'.format(round(x,2)) for x in x_ticks]
			plt.xticks(x_ticks, x_ticks_label)
			plt.xlim(0.0, 0.25)
			
			
			plt.ylabel(r"$ \mathsf{Rq^4 \times \SI{E-08}{(\angstrom^{-4}) } } $")
			plt.xlabel(r"$\mathsf{ q\ (\SI{}{\angstrom^{-1} }) }$")
			legend = plt.legend()
			legend.get_texts()[0].set_text('Sim.')
			legend.draggable(True)
			anot.draggable(True)
			plt.tight_layout()
			plt.show()
			
		
		# get a list of columns
		cols = list(reflectivity_data)
		# move the column to head of list using index, pop and insert
		cols.insert(0, cols.pop(cols.index('q')))
		# use ix to reorder
		reflectivity_data = reflectivity_data.ix[:, cols]
		
		
		#density_data.fillna(0.0, inplace=True)
		header = """
				@    title ""
				@    xaxis  label "z (nm)"
				@    yaxis  label " Rq4 (A-4)"
				@TYPE xy
				@ view 0.15, 0.15, 0.75, 0.85
				@ legend on
				@ legend box on
				@ legend loctype view
				@ legend 0.78, 0.8
				@ legend length 2\n
				"""
		
		if SUPLAYERS is not None:
			adding_to_header = """
				#Reflectivity curve computed using {0} of SLD file with:
				# {1}
				#""".format(REFLECTOMETRY, str(SUPLAYERS))
			
			header = adding_to_header + header
		
		if SUBLAYERS is not None:
			adding_to_header = """
				#Reflectivity curve computed using {0} of SLD file with:
				# {1}
				#""".format(REFLECTOMETRY, str(SUBLAYERS))
			
			header = adding_to_header + header
			
		
		header = ut.RemoveUnwantedIndent(header)
		
		xvg_out = None
		if XVG_FILE.endswith('MDADENS.xvg'):
			xvg_out = open(XVG_FILE.replace('MDADENS.xvg', 'MDAREF.xvg'),'w')
			
		if XVG_FILE.endswith('DENS.xvg'):
			xvg_out = open(XVG_FILE.replace('DENS.xvg', 'REF.xvg'),'w')

		for index, qt in enumerate(reflectivity_data):
			if qt != "q":
				header += """@ s{0} legend "{1}"\n""".format(index-1, qt)
			
		for row in reflectivity_data.itertuples():
			header += '   '.join(map(str, row[1:]))+'\n'
		
		xvg_out.write(header+"\n")
		xvg_out.close()
	
	

elif 'xvgplot' in sys.argv:
	
	try:
		import numpy as np
	except ImportError as e:
		print("You need numpy to use this command")
		print("Error {0}".format(e))
	
	try:
		from io import StringIO
	except ImportError as e:
		print("You need numpy to use this command")
		print("Error {0}".format(e))
		
	try:
		import matplotlib.pyplot as plt
		import matplotlib.cm as mcm
		from matplotlib import rc
	except ImportError as e:
		print("You need matplotlib.pyplot to use this script")
		print("Error {0}".format(e))
	
	QUANTITIES = cmdParam.quantities
	XVG_FILE = cmdParam.xvg_file
	INTEGRATE_LIMITS = cmdParam.limits
	PLOT = cmdParam.plot
	MEAN = cmdParam.mean
	
	quantities_in_file = {}
	if XVG_FILE.endswith('EN.xvg'):
		quantities_in_file.update({ 'X': {'values': 0, 'unit': 'time (ps)'} })
		
	elif XVG_FILE.endswith('DENS.xvg'):
		quantities_in_file.update({ 'X': {'values': 0, 'unit': 'position (nm)'} })
		
	else:
		quantities_in_file.update({ 'X': {'values': 0, 'unit': ''} })
		
	col_ndx = 1
									
	with open(XVG_FILE,'r') as data:
		for line in data:
			if line.startswith('@'):
				#Extract the legend to get the Column/quantity match
				if 's{0}'.format(col_ndx-1) in line:
					Qtty = line[  line.find('"')+1 : line.rfind('"')  ]
					quantities_in_file.update( { Qtty : {'values': col_ndx, 'unit': ""}} )
					col_ndx += 1
				
				continue
			
			elif line.startswith('#'):
				pass
			
			else:
				break

	avg_file = XVG_FILE.replace('xvg','avg')
	if os.path.isfile(avg_file):
		for qtty in quantities_in_file:
			linenumber = quantities_in_file[qtty]['values']
			line = linecache.getline(avg_file, linenumber+6)
			unit = '.'.join(line.split(qtty)[1].split()[4:])
			
			quantities_in_file[qtty]['unit'] = unit
			
	#Remove the '@'s and '#'s
	output_str = ""
	
	with open(XVG_FILE, 'r') as input_file:
		for line in input_file:
			if '#' not in line:
				if '@' not in line:
					output_str += line
		#Modified the # position recently
				
	output_file = StringIO(output_str)
	
	quantities_in_file['X']['values'] = np.loadtxt(output_file, dtype='float', usecols=(0,) )
	output_file.seek(0)
	
	
	if QUANTITIES is not None:
		for qtty in QUANTITIES:
			if qtty in quantities_in_file:
				quantities_in_file[qtty]['values'] = np.loadtxt(output_file, dtype='float', usecols=(quantities_in_file[qtty]['values'],) )
				output_file.seek(0)
			else:
				print("{0} was not in {1}".format(qtty, XVG_FILE))
			
	else:
		for qtty in quantities_in_file:
			if qtty != 'X':
				quantities_in_file[qtty]['values'] = np.loadtxt(output_file, dtype='float', usecols =(quantities_in_file[qtty]['values'],) ) 
				output_file.seek(0)
	
	#plotting the function for each leaflet and total bilayer
	if PLOT:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,6), dpi=72)
		ax.set_xlabel(quantities_in_file['X']['unit'])
		if QUANTITIES is not None:
			for qtty in QUANTITIES:
				if qtty in quantities_in_file:
					ax.plot(quantities_in_file['X']['values'], quantities_in_file[qtty]['values'], alpha=1.0,
														label=qtty+' '+quantities_in_file[qtty]['unit'])
		else:
			for qtty in quantities_in_file:
				if qtty != 'X':
					ax.plot(quantities_in_file['X']['values'], quantities_in_file[qtty]['values'], alpha=1.0,
															label=qtty+' '+quantities_in_file[qtty]['unit'])
		
		if XVG_FILE.endswith('DENS.xvg'):
			ax.set_ylabel('Density (CG/nm^3)')
			ax.set_xlabel('Z-coord (nm)')
			
		elif XVG_FILE.endswith('SLD.xvg'):
			#rc('text', usetex=True)
			#ax.yaxis.major.formatter._useMathText = True
			ax.set_ylabel('SLD (10^{-6} A^{-1})')
			ax.set_xlabel('Z-coord (nm)')
			ax.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
		elif XVG_FILE == 'flipflop.xvg' or XVG_FILE == 'solvent_per_leaflet.xvg':
			ax.set_xlabel('Time (ps)')
			
		
		ax.grid('on')
		ax.legend()
		ax.legend().draggable()
		
		#plt.show()
		plt.savefig(XVG_FILE.replace('xvg','svg'))
		plt.close()
	
	
	if INTEGRATE_LIMITS is not None:
		print("##################################################")
		print("##################   INTEGRATE   #################")
		print("##################################################")
		density_slice = quantities_in_file['X']['values'][1] - quantities_in_file['X']['values'][0]
		if QUANTITIES is not None:
			for qtty in QUANTITIES:
				if qtty in quantities_in_file:
					quantities_in_file[qtty]['integrate'] = 0.0
					
					for X, qt in zip(quantities_in_file['X']['values'], quantities_in_file[qtty]['values']):
						
						if X > INTEGRATE_LIMITS[1]:
							break
						elif X >= INTEGRATE_LIMITS[0] and X <= INTEGRATE_LIMITS[1]:
							quantities_in_file[qtty]['integrate'] += qt*density_slice
					
					
					print("# Integration	of	{0}	between	{1}	and	{2}	:	{3}".format(qtty, 
																			INTEGRATE_LIMITS[0],
																			INTEGRATE_LIMITS[1],
																			quantities_in_file[qtty]['integrate']))
				
			
		else:
			for qtty in quantities_in_file:
				if qtty != 'X':
					quantities_in_file[qtty]['integrate'] = 0.0
					
					for X, qt in zip(quantities_in_file['X']['values'], quantities_in_file[qtty]['values']):
						
						if X > INTEGRATE_LIMITS[1]:
							break
						elif X >= INTEGRATE_LIMITS[0] and X <= INTEGRATE_LIMITS[1]:
							quantities_in_file[qtty]['integrate'] += qt*density_slice
					
					
					print("# Integration	of	{0}	between	{1}	and	{2}	:	{3}".format(qtty, 
																			INTEGRATE_LIMITS[0],
																			INTEGRATE_LIMITS[1],
																			quantities_in_file[qtty]['integrate']))
		
		print("##################################################")
		print("##################   INTEGRATE   #################")
		print("##################################################")
		
	if MEAN:
		print("##################################################")
		print("##############   MEAN, STD, STDERR   #############")
		print("##################################################")
		density_slice = quantities_in_file['X']['values'][1] - quantities_in_file['X']['values'][0]
		if QUANTITIES is not None:
			for qtty in QUANTITIES:
				if qtty in quantities_in_file:
					quantities_in_file[qtty]['mean'] = 0.0
					quantities_in_file[qtty]['std'] = 0.0
					quantities_in_file[qtty]['stderr'] = 0.0
					
					quantities_in_file[qtty]['mean'] = np.mean(quantities_in_file[qtty]['values'], dtype=np.float64)
					quantities_in_file[qtty]['std'] = np.std(quantities_in_file[qtty]['values'], dtype=np.float64)
					quantities_in_file[qtty]['stderr'] = quantities_in_file[qtty]['std']/ math.sqrt(len(quantities_in_file[qtty]['values']))
					
					
					print("# For	{0}	:	{1}		{2}		{3}".format(qtty, quantities_in_file[qtty]['mean'],
																quantities_in_file[qtty]['std'],
																quantities_in_file[qtty]['stderr']))
				
			
		else:
			for qtty in quantities_in_file:
				if qtty != 'X':
					quantities_in_file[qtty]['mean'] = 0.0
					quantities_in_file[qtty]['std'] = 0.0
					quantities_in_file[qtty]['stderr'] = 0.0
					
					quantities_in_file[qtty]['mean'] = np.mean(quantities_in_file[qtty]['values'], dtype=np.float64)
					quantities_in_file[qtty]['std'] = np.std(quantities_in_file[qtty]['values'], dtype=np.float64)
					quantities_in_file[qtty]['stderr'] = quantities_in_file[qtty]['std']/ math.sqrt(len(quantities_in_file[qtty]['values']))
					
					
					print("# For	{0}	:	{1}		{2}		{3}".format(qtty, quantities_in_file[qtty]['mean'],
																quantities_in_file[qtty]['std'],
																quantities_in_file[qtty]['stderr']))
		
		print("##################################################")
		print("##############   MEAN, STD, STDERR   #############")
		print("##################################################")
	
	
	
else:
	print("You need to provide one of the following command:\n")
	print("	compute - command to compute properties using fatslim")
	print("	analyse - command to analyse the properties")
	print("	animate - command to create a movie from analyse output")
	print("	mayavi  - command to see a representation of the membrane with mayavi")
	print("	vpython - command to see a representation of the membrane with vpython")
	print("	gromacs - command to compute gmx energy and density")



