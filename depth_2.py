#!/usr/bin/env python
# coding=utf-8

import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis.leaflet import LeafletFinder
import xlsxwriter

def dpth(u,pro_z):
	L = LeafletFinder(u, 'name P')
	L0 = L.group(0)
	L1 = L.group(1)
	"""bi_cntr = (np.mean(L0.positions,axis=0)[-1]+np.mean(L1.positions,axis=0)[-1])/2"""
	int_cntr=L1.centroid()[-1]
	bi_cntr_dpth = (pro_z-int_cntr)
	return bi_cntr_dpth


if __name__ == "__main__":
	u = MDAnalysis.Universe("pull_dian2.gro","after_pull_dian2.trr")
	pro1=u.select_atoms("resid 25387 and (name C20 or name C19)")
	data=np.array([])
	for ts in u.trajectory[0::10]:
		pro_z=pro1.centroid()[-1]
		dis=dpth(u,pro_z)
		data=np.append(dis,data)

	#plotting
	i=0
	fig = plt.figure(figsize=(24,8))	
	for name in "Prodan1","Prodan2":
		i=i+1
		ax = fig.add_subplot(1,2,i)
		ax.plot(data, 'b-', lw=1)
		ax.set_xlabel(r"time (ns)")
		ax.set_ylabel(r"distance (Angstrom)")
		ax.set_title("Distance between the bilayer interface and %(x)s " % {'x':name})
	fig.savefig("Laurdan_depth_DOPC.png")
	#saving into file
	row = 0
	col = 0

	workbook = xlsxwriter.Workbook('Laurdan_depth_DOPC.xlsx')
	worksheet = workbook.add_worksheet()
	columns = 2
	worksheet.write_column('A1', data)
	print(data)
	workbook.close()
