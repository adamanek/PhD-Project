#!/usr/bin/env python
# coding=utf-8
import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis.leaflet import LeafletFinder
from MDAnalysis.analysis.distances import dist, distance_array 
import xlsxwriter
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

"""def nh20(u,p,m):
	NW=np.array(len(p),2)
	for i in p:
		NW[i][0]=distance_array(m,p)
		NW[i][1] = len(u.select_atoms("resname SOL and around 3.2 %(x)s" % {'x':i})
	return NW"""


if __name__ == "__main__":
	u = MDAnalysis.Universe("pull_dian2.gro","after_pull_dian2.trr")
	"""distances = distance_array(prodan1,phos)
	length=len(phos)
	NW=np.zeros((1,length),dtype=np.int8)
	k=0
	for i in phos:
		#NW[i,0]=distance_array(prodan1,phos)
		j = np.array_str(i).strip("[]")
		NW[0,k]=len(u.select_atoms("resname SOL and point %(x)s 3.2" % {'x':j}))
		k=k+1"""

	distances = np.array([])
	NW = np.array([],dtype=np.int8)
	snap = np.array([])
	select1=u.select_atoms("resname DPPC and name N and (around 25 (resname DIAN and name N1))", updating=True)
	"""exl_sel=u.select_atoms("resname DOPC and name P and (around 15 (resid 25482 and name N))", updating=True)"""
	for ts in u.trajectory[1500::5]:
		prodan1 = u.select_atoms("resname DIAN and name N1").positions
		atoms=select1
		phos1 = atoms.positions
		dis1=distance_array(prodan1,phos1)		
		distances = np.append(distances,dis1)
		snap = np.append(snap,len(phos1))
		for i in phos1:	
			j = np.array_str(i).strip("[]")
			NW1=len(u.select_atoms("name OH2 and point %(x)s 6.4" % {'x':j}))
			NW = np.append(NW,NW1)
	
	workbook = xlsxwriter.Workbook('hydration_dppc_dian_N_2.xlsx')
	worksheet = workbook.add_worksheet()
	columns = len(distances)
	print(columns)
	k=0
	for i in range(columns):
		worksheet.write(i,0,distances[i])
	"""	worksheet.write(i,1,NW[i])"""
	worksheet.write_column('B1',NW)	
	worksheet.write_column('C1',snap)
	workbook.close()

