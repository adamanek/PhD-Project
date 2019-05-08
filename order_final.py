#!/usr/bin/env python
# coding=utf-8
import numpy as np
from numpy import absolute
import MDAnalysis
from MDAnalysis import *
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
from numpy.linalg import norm
import xlsxwriter

u = Universe("pull_dian2.gro","after_pull_dian2.trr")
mem_normal = np.array([0,0,1])
tail_carbons = np.arange(2, 17)
diam = ([0,7,12,16,20]) #sets the ranges of distances from dye
orders = np.zeros((len(diam)-1,len(tail_carbons)))
orders_stderr = np.zeros((len(diam)-1,len(tail_carbons)))

#Loop that goes over all the concentric spheres
for k in range(len(diam)-1): 
    #Initiating variables
    x_axis = []
    y_axis = []
    Szz = np.array([])
    Syy = np.array([])
    Syz = np.array([])
    order_param = np.zeros(len(tail_carbons))
    order_param_stderr = np.zeros(len(tail_carbons))
    #Loop that goes over each carbon in the phospholipid tail
    for i, carbon in enumerate(tail_carbons):
        S_cd = np.array([])

        #Loop that goes over the trajectory file
        for ts in u.trajectory[1000::100]:

            #Selection of the two C-H bonds around the dye
            selection = "(resname DPPC and ( name C2{0:d} or name H{1:d}R or name H{2:d}S or name C3{0:d} or name H{4:d}X or name H{5:d}Y) and byres ((resname DPPC and name P) and (sphlayer %(x)s %(y)s (name N1 and resname DIAN))))".format(* \
                        ((carbon,) * 6)) % {'x':diam[k], 'y':diam[k+1]}
            #This selection is used for the C=C bond as they only have one H atom each
            selZ9 = "(resname DPPC and (name C28 or name C38 or name C29 or name C39 or name C210 or name C310) and byres ((resname DPPC and name P) and (sphlayer %(x)s %(y)s (name N1 and resname DIAN))))" % {'x':diam[k], 'y':diam[k+1]}
            selZ = u.select_atoms(selZ9, format="afc").positions
            data = u.select_atoms(selection, format="afc").positions

            #Following if conditions check if the 'tail_carbons' loops has reached the C=C bond, if yes uses the C=C selection and calculates order parameters for the two carbons, else uses normal C-H selection
            if len(data) % 3 != 0:
                if carbon == 9:
                    selZ9 = "(resname DPPC and (name C29 or name C39 or name H9R or name H9X) and byres ((resname DPPC and name P) and (sphlayer %(x)s %(y)s (name N1 and resname DIAN))))" % {'x':diam[k], 'y':diam[k+1]}
                    data = u.select_atoms(selZ9, format="afc").positions
                    cd = np.concatenate((data[1::2] - data[0::2]), axis=0)
                    for l in range(len(cd)):
                        cos_theta = np.dot(cd[l],mem_normal)/(norm(cd[l])*norm(mem_normal))
                        S_cd = np.append(S_cd, -0.5 * (3 * np.square(cos_theta) - 1))
                if carbon == 10:
                    selZ9 = "(resname DPPC and (name C210 or name C310 or name H10R or name H10X) and byres ((resname DPPC and name P) and (sphlayer %(x)s %(y)s (name N1 and resname DIAN))))" % {'x':diam[k], 'y':diam[k+1]}
                    data = u.select_atoms(selZ9, format="afc").positions
                    cd = np.concatenate((data[1::2] - data[0::2]), axis=0)
                    for l in range(len(cd)):
                        cos_theta = np.dot(cd[l],mem_normal)/(norm(cd[l])*norm(mem_normal))
                        S_cd = np.append(S_cd, -0.5 * (3 * np.square(cos_theta) - 1))
            else: 
                cd = np.concatenate((data[1::3] - data[0::3], data[2::3] - data[0::3]), axis=0) 
                cd_r = np.sqrt(np.sum(np.power(cd, 2), axis=-1))

            # Dot product with the z axis
                cos_theta = cd[..., 2] / cd_r
                S_cd1 = -0.5 * (3. * np.square(cos_theta) - 1)
                S_cd = np.append(S_cd, S_cd1)

        del data
        del cd
        order_param[i] = np.average(S_cd)
        order_param_stderr[i] = np.std(S_cd) / np.sqrt(len(S_cd))
    orders[k] = order_param
    orders_stderr[k] = order_param_stderr
    print(orders)

workbook = xlsxwriter.Workbook('order_dianDPPC.xlsx')
worksheet = workbook.add_worksheet()
for i in range(len(diam)-1):
    worksheet.write_column(0, i+1, orders[i])
    worksheet.write_column(0, i+5, orders_stderr[i])   
workbook.close()	
