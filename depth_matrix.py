#!/usr/bin/env python
# coding=utf-8

import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis.leaflet import LeafletFinder
import xlsxwriter
import seaborn as sns
import tqdm

def dpth(u):
	L = LeafletFinder(u, 'name P')
	L0 = L.group(0)
	L1 = L.group(1)
	"""bi_cntr = (np.mean(L0.positions,axis=0)[-1]+np.mean(L1.positions,axis=0)[-1])/2"""
	int_cntr=L0.centroid()[-1]
	return int_cntr


if __name__ == "__main__":
    u = MDAnalysis.Universe("gs_dopc_laur.gro","gs_dopc_laur_small.trr")
    #dye=u.select_atoms("resid 28662 and (type N or type O or name C24)", updating = True)
    dye=u.select_atoms("resid 28663 and (not type H)", updating = True)
    leaflet_P = u.select_atoms("name P", updating = True)
    L = LeafletFinder(u, 'name P')
    L0 = L.group(0)
    L1 = L.group(1)
    data=np.empty([len(dye)+1,len(u.trajectory)])  
    for i, ts in enumerate(tqdm.tqdm(u.trajectory[0::1])):
        leaflet_pos_z = L1.centroid()[-1]        
        for j, atom in enumerate(dye):
            data[0,i] = u.trajectory.time
            atom_z = atom.position[-1]
            distance_atom_bound = - atom_z + leaflet_pos_z
            data[j+1,i] = distance_atom_bound
    ax = sns.heatmap(data[1:,:], xticklabels = 1000)
    
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
