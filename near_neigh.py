#!/usr/bin/env python
import numpy as np
from numpy.linalg import norm
from numpy import *
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
import matplotlib
import matplotlib.pyplot as plt
import time
import xlsxwriter

#Calculates the nearest neighbour between oxygen in Prodan and Phosphate by finding the distance between O and each P and then finding the minimum
def nns(u, Prod):
	Phos = u.select_atoms("resid 12490 and name P8")
	Dist = distance_array(Phos.positions,Prod.positions,box = np.array([130.2080,130.2080,82.9370, 90., 90., 90.], dtype=float32))
	min_dis = np.amin(Dist)
	return min_dis

if __name__ == "__main__":
	
	u=MDAnalysis.Universe("npt_gs.gro", "md_dopc_gs.xtc")
	domains ={
		'Prodan1': u.select_atoms("resid 24738 and type O"),
		'Prodan2': u.select_atoms("resid 18553 and type O")
		}
	Prod=dict((name,g) for name,g in domains.iteritems())
	
	#Creates an empty dictionary into which results would be inserted
	nframes=len(u.trajectory)
	results = dict((name, np.zeros((nframes, 2), dtype=np.float32)) for name in domains)	
	
	#Iteration over all timesteps in trajector, finds the nearest neighbour for each frame
	for iframe,ts in enumerate(u.trajectory):
		for name, g in domains.iteritems():
			results[name][iframe, :] = u.trajectory.time, nns(u,Prod[name])
	#plotting
	i=0 
	fig = plt.figure(figsize=(24,8))	
	for name in "Prodan1","Prodan2":
		i=i+1	#done to create the appropriate number of subplots
		data = results[name]
		ax = fig.add_subplot(1,2,i)
		ax.plot(data[:,0],data[:,1], 'b-', lw=1)
		ax.set_xlabel(r"time (ns)")
		ax.set_ylabel(r"distance")
		ax.set_title("Nearest Neighbours between C terminus and Oxygen in %(x)s " % {'x':name})
	fig.savefig("Dis_P2.png")

	#saving into file
	row = 0
	col = 0

	workbook = xlsxwriter.Workbook('Dis_P2.xlsx')
	worksheet = workbook.add_worksheet()
	columns = 2
	for name in (results):
		data=results[name]
		worksheet.write(row, col, name)
		for i in range(columns):
			worksheet.write_column(row+1,col,data[:,i])
			col+=1
		col +=1
	workbook.close()
