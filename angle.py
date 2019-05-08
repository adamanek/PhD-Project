#!/usr/bin/env python
# coding=utf-8
import numpy as np
from numpy.linalg import norm
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis.leaflet import LeafletFinder
from MDAnalysis.analysis.distances import dist, distance_array 
import xlsxwriter
import time
start_time=time.time()

def theta(x,y,z):
	V1 = x-y
	V2 = y-z
	theta1 = np.arccos(np.dot(V1,V2)/(norm(V1)*norm(V2)))
	return np.rad2deg(theta1)
	
if __name__ == "__main__":
	u = MDAnalysis.Universe("npt_gs.gro","md_dopc_gs.xtc")

	distances = np.array([])
	angles = np.array([])
	for ts in u.trajectory[30000::5]:
		prodan1 = u.select_atoms("resid 24738 and name O5").positions
		phos = u.select_atoms("resname DOPC and name P8 and around 15 (resid 24738 and name O5)")
		dis=distance_array(prodan1,phos.positions)
		distances = np.append(distances,dis)
		for i in phos:
			selection = "resname DOPC and name P8 and around 15 (resid 24738 and name O5)"
			oxy = u.select_atoms("resname DOPC and name O11 and byres resid %(x)s" %{'x':i.resid})
			nit = u.select_atoms("resname DOPC and name N4 and byres resid %(x)s" %{'x':i.resid})
			l = np.squeeze(np.asarray(oxy.positions))
			k = np.squeeze(np.asarray(nit.positions))
			ang = theta(l,i.position,k)
			angles = np.append(angles,ang)

	

	workbook = xlsxwriter.Workbook('headgroup_angles_40ns.xlsx')
	worksheet = workbook.add_worksheet()
	columns = len(distances)
	k=0
	for i in range(columns):
		worksheet.write(i,0,distances[i])
		worksheet.write(i,1,angles[i])
	workbook.close()

	elapsed_time = time.time() - start_time
	print(elapsed_time)
