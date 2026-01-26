#!/usr/bin/env python

# B. Dutagaci, 2025

import os,sys
import numpy as np
import math
import random

import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--input_path',type=str, default='input.pdb',
    help='Input structure.')
parser.add_argument(
    '--output_path',type=str, default='output.pdb',
    help='Output structure.')
parser.add_argument(
    '--cutoff',type=float, default=5,
    help='Distance cutoff between carbon atoms.')
parser.add_argument(
    '--box',type=float, default=100,
    help='Box size.')
parser.add_argument(
    '--particle_number',type=int, default=4,
    help='Number of particles.')
arg = parser.parse_args()

outputfile = open(arg.output_path,"w")

inputfile = open(arg.input_path,"r")
lines = inputfile.readlines()
inputfile.close()

atomlist = ["CA","CB","CG","CD","CZ","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12"]

xca = []
yca = []
zca = []
cog = [0,0,0]
for line in lines:
    if "ATOM" in line:
        if line.split()[2] in atomlist:
            xca.append(float(line.split()[6]))
            yca.append(float(line.split()[7]))
            zca.append(float(line.split()[8]))
            cog[0]+=float(line.split()[6])*12.011
            cog[1]+=float(line.split()[7])*12.011
            cog[2]+=float(line.split()[8])*12.011

cog[0]=cog[0]/(len(xca)*12.011)
cog[1]=cog[1]/(len(xca)*12.011)
cog[2]=cog[2]/(len(xca)*12.011)

r = arg.cutoff
box = arg.box
np = arg.particle_number
maxh = 10000
new_cog_list = []
new_angle_list = []

newx_list = []
newy_list = []
newz_list = []

for i in range(1,np+1):
    print("np index = ",i)
    h=0
    while h<maxh:
        newx = []
        newy = []
        newz = []

        x = random.uniform(0.0,box)
        y = random.uniform(0.0,box)
        z = random.uniform(0.0,box)

        tethax = random.uniform(0.0,360.0)
        tethay = random.uniform(0.0,360.0)
        tethaz = random.uniform(0.0,360.0)

        tethax_rad = math.radians(tethax)
        tethay_rad = math.radians(tethay)
        tethaz_rad = math.radians(tethaz)

        for j in range(len(xca)):
            x_center = xca[j]-cog[0]
            y_center = yca[j]-cog[1]
            z_center = zca[j]-cog[2]

            xi = x_center
            yi = y_center*math.cos(tethax_rad) - z_center*math.sin(tethax_rad)
            zi = y_center*math.sin(tethax_rad) + z_center*math.cos(tethax_rad)

            xii = xi*math.cos(tethay_rad) + zi*math.sin(tethay_rad)
            yii = yi
            zii = -1*xi*math.sin(tethay_rad) + zi*math.cos(tethay_rad)

            xiii = xii*math.cos(tethaz_rad) - yii*math.sin(tethaz_rad)
            yiii = xii*math.sin(tethaz_rad) + yii*math.cos(tethaz_rad)
            ziii = zii


            newx.append(xiii+x)
            newy.append(yiii+y)
            newz.append(ziii+z)
        if i==1:
            new_angle_list.append([tethax_rad,tethay_rad,tethaz_rad])
            new_cog_list.append([x,y,z])
            newx_list.append(newx)
            newy_list.append(newy)
            newz_list.append(newz)
            h=maxh
        else:
          break_con = False
          for j in range(len(newx)):
            for k in range(len(newx_list)):
              for t in range(len(newx_list[k])):
                d = ((newx[j]-newx_list[k][t])**2 + (newy[j]-newy_list[k][t])**2 + (newz[j]-newz_list[k][t])**2)**0.5
                if d < r:
                    h=h+1
                    break_con = True
                    break
              if break_con:
                break
            if break_con:
              break
          if not break_con:
            new_angle_list.append([tethax_rad,tethay_rad,tethaz_rad])
            new_cog_list.append([x,y,z])
            newx_list.append(newx)
            newy_list.append(newy)
            newz_list.append(newz)
            print("h = ",h)
            h=maxh

j=1
print(cog)
for i in range(len(new_cog_list)):
    for line in lines:
        if "ATOM" in line:
            x_center = float(line.split()[6])-cog[0]
            y_center = float(line.split()[7])-cog[1]
            z_center = float(line.split()[8])-cog[2]

            xi = x_center
            yi = y_center*math.cos(new_angle_list[i][0]) - z_center*math.sin(new_angle_list[i][0])
            zi = y_center*math.sin(new_angle_list[i][0]) + z_center*math.cos(new_angle_list[i][0])

            xii = xi*math.cos(new_angle_list[i][1]) + zi*math.sin(new_angle_list[i][1])
            yii = yi
            zii = -1*xi*math.sin(new_angle_list[i][1]) + zi*math.cos(new_angle_list[i][1])

            xiii = xii*math.cos(new_angle_list[i][2]) - yii*math.sin(new_angle_list[i][2])
            yiii = xii*math.sin(new_angle_list[i][2]) + yii*math.cos(new_angle_list[i][2])
            ziii = zii


            newx = xiii+new_cog_list[i][0]
            newy = yiii+new_cog_list[i][1]
            newz = ziii+new_cog_list[i][2]
            print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%11.3f"%(newx),"%7.3f"%(newy),"%7.3f"%(newz), " 1.00  0.00 ",file=outputfile)
            j=j+1

outputfile.close()
