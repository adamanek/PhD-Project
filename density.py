#!/usr/bin/env python
# coding=utf-8
import numpy as np
from numpy.linalg import norm
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis.leaflet import LeafletFinder
from MDAnalysis.analysis.distances import dist, distance_array 
from MDAnalysis.analysis.density import density_from_Universe, Density
import xlsxwriter

if __name__ == "__main__":
	u = MDAnalysis.Universe("npt_gs.gro","md_dopc_gs.xtc")

	distances = np.array([])
	x_coord = np.array([])
	y_coord = np.array([])
	for ts in u.trajectory[30000::10]:
		prodan1 = u.select_atoms("resid 24738 and name O5").positions
		phos = u.select_atoms("resname DOPC and around 20 (resid 24738 and name O5)").positions
		dis=distance_array(prodan1,phos)
		distances = np.append(distances,dis)
		for i in range(len(phos)):		
			x_coord = np.append(x_coord,phos[i,0])
			y_coord = np.append(y_coord,phos[i,1])

	workbook = xlsxwriter.Workbook('DOPC_P_dens_whole_lipid.xlsx')
	worksheet = workbook.add_worksheet()
	columns = 4
	k=0
	worksheet.write_column(0,0,distances)
	for i in range(len(x_coord)):	
		worksheet.write(i,1,x_coord[i])
		worksheet.write(i,2,y_coord[i])
	workbook.close()
