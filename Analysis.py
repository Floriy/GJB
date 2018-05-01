#!/usr/bin/python3
# -*- coding: utf-8 -*-


import argparse #Argument parser module !
import sys
import glob
import os
import shutil
import subprocess as sub
import traceback
import math
from collections import OrderedDict
import gc #garbage collector
import time


#global variables
extensions = ['ndx','gro','xtc']

#def check_existing_dir
def depth(d, level=1):
	if not isinstance(d, dict) or not d:
		return level
	
	return max(depth(d[k], level + 1) for k in d)

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
		bounds = np.linspace(0.0, 6.0, 30)
		levels = np.linspace(0.0, 6.0, 30)
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
		bounds = np.linspace(-0.5, 1, 15)
		levels = np.linspace(-0.5, 1, 15)
		colormap = mcm.viridis
	if prop == 'APL':
		bounds = np.linspace(0.55, 0.75, 20)
		levels = np.linspace(0.55, 0.75, 20)
		colormap = mcm.plasma
	if prop == 'THICKNESS':
		bounds = np.linspace(0.0, 6.0, 30)
		levels = np.linspace(0.0, 6.0, 30)
		colormap = mcm.BuPu_r
	#figure for plots
	fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15,6), dpi=72, gridspec_kw = {'width_ratios':[5,5,5,1]})
	
	cs = [None, None, None]
	#sp = [None, None, None]
	
	header = "#Averaged on {0:d} ps from {1:d} to {2:d}\n#".format(ending_time - beginning_time,
																beginning_time, ending_time)
	for leaflet in grid_dict:
		# Setting the axes on which to plot
		Nbleaflet = grid_dict[leaflet]['plot']
		grid_dict[leaflet]['plot'] = axs[ grid_dict[leaflet]['plot'] ]
		# Plot map
		cs[ Nbleaflet ] = grid_dict[leaflet]['plot'].contourf(grid_x, grid_y, grid_dict[leaflet]['grid'], cmap = colormap, levels=levels)
		#sp [ Nbleaflet ] = grid_dict_std[leaflet]['plot'].contourf(grid_x, grid_y, grid_dict[leaflet]['grid'], cmap = colormap, levels=levels)
		grid_dict[leaflet]['plot'].set_title(prop+' for '+leaflet)
		grid_dict[leaflet]['plot'].set_xlabel("Box-X (nm)")
		grid_dict[leaflet]['plot'].set_ylabel("Box-Y (nm)")
		
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
	axs[2].set_title(prop+' for bilayer')
	axs[2].set_xlabel("Box-X (nm)")
	axs[2].set_ylabel("Box-Y (nm)")
	
	CB1 = fig.colorbar(cs[0],cax=axs[3],orientation = 'vertical', boundaries=bounds, ticks=bounds, extend='both', extendfrac='auto')
	#Set the time as title
	title = "Averaged on {0:d} ps from {1:d} to {2:d}".format(ending_time - beginning_time,
																beginning_time, ending_time)
	plt.suptitle(title)
	
	for ax in axs[:-1]:
		ax.autoscale(axis='x', tight=True)
		ax.set_xlim(0.0, mean_box[0])
		ax.set_ylim(0.0, mean_box[1])
	
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
	
	low_leaflet_data = []
	upp_leaflet_data = []
	bil_leaflet_data = []
	
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
	
	#plotting the function for each leaflet and total bilayer
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,6), dpi=72)
	
	ax.fill(radii, low_leaflet_data, alpha=0.4, color='blue', label='lower leaflet')
	ax.fill(radii, upp_leaflet_data, alpha=0.4, color='red', label='upper leaflet')
	ax.plot(radii, bil_leaflet_data, alpha=1.0, color='purple', label='bilayer')
	ax.set_xlabel("r (nm)")
	ax.set_ylabel(prop)
	ax.grid('on')
	ax.legend()
	
	title = "Radial distribution of {0} \n averaged on {1:d} ps from {2:d} to {3:d} ps".format(prop, ending_time - beginning_time, beginning_time, ending_time)
	plt.suptitle(title)
	
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
				plot_grid(grid, prop, name, bdim)
			
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
			
			
			graph_name = "{0}/{1}_mean.svg".format(grid_graph_folder, prop)
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
				graph_name = "{0}/{1}_radial.svg".format(grid_graph_folder, prop)
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

gmxOpt.add_argument('-r','--radius', dest='radius',
					type=int, default=None,
                    help='Set the radius outside which the membrane will be analysed')

# Adding parser options to gromacs analysis##############################################################
xvgplot = sub_parser.add_parser('xvgplot')

xvgOpt = xvgplot.add_argument_group("General options")

xvgOpt.add_argument('-q', '--quantities', dest='quantities',
					type=str, nargs='*', default=None,
                    help='Set the name of the quantities to plot')

xvgOpt.add_argument('-f', '--file', dest='xvg_file',
					type=str, default=None,
                    help='Set the name of the file to plot')


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
							hg_membrane_ndx = ndx_dir +'/'+ name.replace('_out.gro','_HG.ndx')
							
							ndx_file = destination.replace('_out.gro','.ndx')
							#put a temp name for plot file name
							plotname = xvg_dir +'/'+ name.replace('_out.gro','tmp')
							#put a temp name for csv file name
							csv_apl_name = csv_apl_dir +'/'+ name.replace('_out.gro','_APL.csv')
							csv_order_name = csv_order_dir +'/'+ name.replace('_out.gro','_ORDER.csv')
							csv_thickness_name = csv_thickness_dir +'/'+ name.replace('_out.gro','_THICKNESS.csv')
							
						else:
							if 'gro' in name:
								gro_file = file_name
								
								#File name for fatslim index
								hg_membrane_ndx = ndx_dir +'/'+ name.replace('_out.gro','_HG.ndx')
								
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
		
		FILENAME = cmdParam.filename
		filename = None
		try:
			filename = glob.glob(FILENAME)[0]
		except Error as e:
			print(e)
			
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
		
elif 'vpython' in sys.argv:
	try:
		import numpy as np
	except ImportError as e:
		print("You need vpython to use this command")
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
	GMX = cmdParam.gmx_path
	MD_RUN_INPUT = cmdParam.md_run_input
	FOLDER = cmdParam.folder
	BEGIN_FRAME = cmdParam.begin_frame
	END_FRAME = cmdParam.end_frame
	DT = cmdParam.dt
	RADIUS =cmdParam.radius

	# Setting parameters for gmx energy and density
	
	

	filesFound = {}

	if MD_RUN_INPUT is not None and FOLDER is not None:
		"""
		These options enables membrane property computation for a lot of data
		If you want only a peculiar file set it using -c and not -r and -f
		"""
		if FOLDER[-1] == '/': FOLDER = FOLDER[:-1]
		#print(FOLDER)
		
		for Run in MD_RUN_INPUT:
			for path in glob.glob(FOLDER +'/**/*'+ Run +'*xtc', recursive=True):
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
		energy_dir = analysisFolder+'/GMX_ENERGY'
		os.makedirs(energy_dir)
		# Directory for gmx density
		density_dir = analysisFolder+'/GMX_DENSITY'
		os.makedirs(density_dir)
		
		if DEPTH >= 3:
			for project_name in filesFound:
				for job_name, file_name in filesFound[project_name].items():
					
					job_energy_dir = analysisFolder+'/GMX_ENERGY/'+job_name
					job_density_dir = analysisFolder+'/GMX_DENSITY/'+job_name
					
					os.makedirs(job_energy_dir, exist_ok=True)
					os.makedirs(job_density_dir, exist_ok=True)
					
					xtc_file = None
					edr_file = None
					tpr_file = None
					ndx_file = None
					destination = None
					
					plotname = None
					csvname = None
					
					
					
					name = file_name.split('/')[-1]
					print(name)
					
					xtc_file = file_name
					edr_file = file_name.replace('xtc','edr')
					tpr_file = file_name.replace('xtc','tpr')
					energy_xvg_file = xtc_file.split('/')[-1].replace('.xtc','_EN.xvg')
					energy_avg_file = energy_xvg_file.replace('xvg','avg')
					density_xvg_file = xtc_file.split('/')[-1].replace('.xtc','_DENS.xvg')
					ndx_file = '_'.join(file_name.split('_')[:-1]) + ".ndx"
					
					
					#Look in the index to find DEF atoms
					is_def_in_sample = False
					#with open(ndx_file,'r') as index:
						#indexStr = index.read()
						#if 'DEF' in indexStr:
							#is_def_in_sample = True
							
					#if is_def_in_sample:
						#select_cmd = "gmx2016 select -f BILAYER_740DSPC_14800W_21DEFv1.0C_1-NPT1.xtc -n BILAYER_740DSPC_14800W_21DEFv1.0C.ndx -s BILAYER_740DSPC_14800W_21DEFv1.0C_1-NPT1.tpr -on index.ndx"
									
					begin_end = ""
					if BEGIN_FRAME is not None:
						begin_end +=" -b {0}".format(BEGIN_FRAME)
					if END_FRAME is not None:
						begin_end +=" -e {0}".format(END_FRAME)
					
					#Computing gmx energy
					# To output all quantities
					variables_for_output = "echo "
					for i in range(1,100):
						variables_for_output += repr('{0}\n'.format(i))
					variables_for_output += repr('\n\n')
					
					try:
						energy_cmd = "{0} | {1} energy {2} -f {3} -s {4} -o {5} > {6}".format(variables_for_output, GMX, begin_end, edr_file, tpr_file,
																						job_energy_dir+'/'+ energy_xvg_file, job_energy_dir + '/' + energy_avg_file)
						sub.call(energy_cmd, shell=True)
					
					except OSError as e:
						print(e)
					
					
					
					# Computing gmx density
					# Get the number of group to output
					read_index_cmd = """cat {0} | grep "\[" """.format(ndx_file)
					read_index_proc = sub.Popen(read_index_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
					read_index_out = read_index_proc.stdout.read()
					
					nb_index = len(read_index_out.splitlines())
					
					groups_for_output = "echo "
					for i in range(0, nb_index):
						groups_for_output += repr('{0}\n'.format(i))
					groups_for_output += repr('\n\n')
					
					# Call gmx density
					try:
						density_cmd = "{0} | {1} density {2} -f {3} -s {4} -o {5} -ng {6} -n {7}".format(groups_for_output, GMX, begin_end, xtc_file, tpr_file,
																										job_density_dir+'/'+ density_xvg_file, 
																										nb_index, ndx_file)
						sub.call(density_cmd, shell=True)
					
					except OSError as e:
						print(e)
					
					
		else:
			for job_name, file_names in filesFound.items():
				
					os.makedirs(analysisFolder+'/'+job_name, exist_ok=True)
					
					for file_name in file_names:
						destination = analysisFolder+'/'+job_name+'/'+file_name.split('/')[-1]
						shutil.copyfile(file_name, destination)
						
	
	print( "--- gromacs done in {0:f} seconds ---".format(time.time() - start_time) )

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
	except ImportError as e:
		print("You need matplotlib.pyplot to use this script")
		print("Error {0}".format(e))
	
	QUANTITIES = cmdParam.quantities
	XVG_FILE = cmdParam.xvg_file
	
	quantities_in_file = {}
	if XVG_FILE.endswith('EN.xvg'):
		quantities_in_file.update({ 'X': {'values': 0, 'unit': 'time (ps)'} })
		
	elif XVG_FILE.endswith('DENS.xvg'):
		quantities_in_file.update({ 'X': {'values': 0, 'unit': 'position (nm)'} })
		
	else:
		quantities_in_file.update({ 'X': {'values': 0, 'unit': ''} })
		
	ColumnIndex = 1
									
	with open(XVG_FILE,'r') as data:
		for line in data:
			if line.startswith('@'):
				#Extract the legend to get the Column/quantity match
				if 's{0}'.format(ColumnIndex-1) in line:
					Qtty = line[  line.find('"')+1 : line.rfind('"')  ]
					quantities_in_file.update( { Qtty : {'values': ColumnIndex, 'unit': ""}} )
					ColumnIndex += 1
				
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
			if not '#' in line:
				if not '@' in line:
					output_str += line
				
	output_file = StringIO(output_str)
	
	quantities_in_file['X']['values'] = np.loadtxt(output_file, dtype='float', usecols=(0,) )
	output_file.seek(0)
	
	print(quantities_in_file)
	
	if QUANTITIES is not None:
		for qtty in QUANTITIES:
			if qtty in quantities_in_file:
				quantities_in_file[qtty]['values'] = np.loadtxt(output_file, dtype='float', usecols=(quantities_in_file[qtty]['values'],) )
				output_file.seek(0)
			else:
				print("{0} was not in {1}".format(qtty, xvg_file))
	else:
		for qtty in quantities_in_file:
			if qtty != 'X':
				quantities_in_file[qtty]['values'] = np.loadtxt(output_file, dtype='float', usecols =(quantities_in_file[qtty]['values'],) ) 
				output_file.seek(0)
	
	#plotting the function for each leaflet and total bilayer
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
		ax.set_ylabel('Density (kg/m^3)')
		
	ax.grid('on')
	ax.legend()
	ax.legend().draggable()
	
	plt.show()
	plt.close()
	
	
	
else:
	print("You need to provide one of the following command:\n")
	print("	compute - command to compute properties using fatslim")
	print("	analyse - command to analyse the properties")
	print("	animate - command to create a movie from analyse output")
	print("	mayavi  - command to see a representation of the membrane with mayavi")
	print("	vpython - command to see a representation of the membrane with vpython")
	print("	gromacs - command to compute gmx energy and density")



