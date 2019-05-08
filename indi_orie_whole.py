#!/usr/bin/env python
import numpy as np
from numpy.linalg import norm
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
import time
import xlsxwriter

start_time=time.time()

Origin = np.array([0,0,0])

def theta1(u,ax_vec):
	"""Calculates the angle between axis and vector going from Oxygen to Nitrogen in Prodan 1"""

	Opro1 = u.select_atoms("resid 24738 and type O").center_of_mass()
	Npro1 = u.select_atoms("resid 24738 and type N").center_of_mass()
	V1 = Opro1-Npro1
	theta1 = np.arccos(np.dot(V1,ax_vec)/(norm(V1)*norm(ax_vec)))
	return np.rad2deg(theta1)
	
def theta2(u,ax_vec):
	Opro2 = u.select_atoms("resid 18553 and type O").center_of_geometry()
	Npro2 = u.select_atoms("resid 18553 and type N").center_of_geometry()
	V2 = Opro2-Npro2
	theta2 = np.arccos(np.dot(V2,ax_vec)/(norm(V2)*norm(ax_vec)))
	return np.rad2deg(theta2)

u = MDAnalysis.Universe("npt_gs.gro","md_dopc_gs_fixed.xtc")

domains ={
	'X': np.array([1,0,0]),
	'Y': np.array([0,1,0]),
	'Z': np.array([0,0,1]),
	}
#Calculating axis vector
ax_vec = dict((name, g) for name,g in domains.iteritems())
nframes=len(u.trajectory)
results = dict((name, np.zeros((nframes, 3), dtype=np.float32)) for name in domains)

#plotting
fig = plt.figure(figsize=(24,16))

for iframe,ts in enumerate(u.trajectory):
	for name, g in domains.iteritems():
		results[name][iframe, :] = u.trajectory.time,theta1(u,ax_vec[name]), theta2(u,ax_vec[name])

i=0	
for name in "X","Y","Z":
	i=i+1
	data = results[name]
	ax = fig.add_subplot(2,3,i)
	ax.plot(data[:,0],data[:,1], 'b-', lw=1)
	ax.set_xlabel(r"time (ns)")
	ax.set_ylabel(r"angle $\theta$")
	ax.set_title("Orientation of Prodan 1 wrt %(x)s axis" % {'x':name})

for name in "X","Y","Z":
	i=i+1
	data = results[name]
	ax = fig.add_subplot(2,3,i)
	ax.plot(data[:,0],data[:,2], 'r-', lw=1)
	ax.set_xlabel(r"time (ns)")
	ax.set_ylabel(r"angle $\theta$")
	ax.set_title("Orientation of Prodan 2 wrt %(x)s axis" % {'x':name})

fig.savefig("Orientations_DEG_fixed.png")

#saving into file
row = 0
col = 0

workbook = xlsxwriter.Workbook('Orientations_deg.xlsx')
worksheet = workbook.add_worksheet()
columns = 3
for name,g in results.iteritems():
	data=results[name]
	worksheet.write(row, col, name)
	for i in range(columns):
		worksheet.write_column(row+1,col,data[:,i])
		col+=1
	col +=1

workbook.close()
elapsed_time = time.time() - start_time
print(elapsed_time)

