#!/usr/bin/env python

# B. Dutagaci, 2025

import os,sys
import numpy as np
import math
import random


#DEFINING FUNCTIONS

#Defining the function of center of mass calculation.
def center(x,y,z,xcom,ycom,zcom):
    xi = x - xcom
    yi = y - ycom
    zi = z - zcom

    return xi,yi,zi

#Defining the function of rotation around x axis.
def rotatex(x,y,z,td):

    t = math.radians(td)

    xi = x
    yi = y*math.cos(t) - z*math.sin(t)
    zi = y*math.sin(t) + z*math.cos(t)

    return xi,yi,zi

#Defining the function of rotation around y axis.
def rotatey(x,y,z,td):

    t = math.radians(td)

    xi = x*math.cos(t) + z*math.sin(t)
    yi = y
    zi = -1*x*math.sin(t) + z*math.cos(t)

    return xi,yi,zi

#Defining the function of rotation around z axis.
def rotatez(x,y,z,td):

    t = math.radians(td)

    xi = x*math.cos(t) - y*math.sin(t)
    yi = x*math.sin(t) + y*math.cos(t)
    zi = z

    return xi,yi,zi

#Defining the function of translation.
def translate(x,y,z,xcom,ycom,zcom,xt,yt,zt):
    xi = x+xcom+xt
    yi = y+ycom+yt
    zi = z+zcom+zt

    return xi,yi,zi

#GENERATING STRUCTURES

#Generating the dimer structure of triple helices from the initial triple helix structure.

outputfile = open("assembly_double.pdb","w")

inputfile = open("model_ccbuilder_prl.pdb","r")
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
cog[1]=cog[1]/(len(yca)*12.011)
cog[2]=cog[2]/(len(zca)*12.011)

j=1

for line in lines:
  if "ATOM" in line:
      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%11.3f"%(float(line.split()[6])),"%7.3f"%(float(line.split()[7])),"%7.3f"%(float(line.split()[8])), " 1.00  0.00 ",file=outputfile)
      j=j+1

angle_y = 180
angle_z = -35
trans = [0,18,15]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line.split()[6]),float(line.split()[7]),float(line.split()[8]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_rot2,y_rot2,z_rot2 = rotatez(x_rot,y_rot,z_rot,angle_z)
      x_new,y_new,z_new = translate(x_rot2,y_rot2,z_rot2,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%11.3f"%(x_new),"%7.3f"%(y_new),"%7.3f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

outputfile.close()

#Generating the hexamer structure of triple helices from the dimer structure.

outputfile = open("assembly_hexa.pdb","w")

inputfile = open("assembly_double.pdb","r")
lines = inputfile.readlines()
inputfile.close()

xca = []
yca = []
zca = []
cog = [0,0,0]
for line in lines:
    if "ATOM" in line:
        if line.split()[2] in atomlist:
            xca.append(float(line.split()[5]))
            yca.append(float(line.split()[6]))
            zca.append(float(line.split()[7]))
            cog[0]+=float(line.split()[5])*12.011
            cog[1]+=float(line.split()[6])*12.011
            cog[2]+=float(line.split()[7])*12.011

cog[0]=cog[0]/(len(xca)*12.011)
cog[1]=cog[1]/(len(yca)*12.011)
cog[2]=cog[2]/(len(zca)*12.011)

j=1

for line in lines:
  if "ATOM" in line:
      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[4])),"%11.3f"%(float(line.split()[5])),"%7.3f"%(float(line.split()[6])),"%7.3f"%(float(line.split()[7])), " 1.00  0.00 ",file=outputfile)
      j=j+1

angle_x = -60
trans = [0,55,-30]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatex(x_center,y_center,z_center,angle_x)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[4])),"%11.3f"%(x_new),"%7.3f"%(y_new),"%7.3f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

angle_x = 60
trans = [0,-55,-30]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatex(x_center,y_center,z_center,angle_x)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[4])),"%11.3f"%(x_new),"%7.3f"%(y_new),"%7.3f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

outputfile.close()

#Generating the dodecamer structure of triple helices from the hexamer structure.

outputfile = open("assembly_dodeca.pdb","w")

inputfile = open("assembly_hexa.pdb","r")
lines = inputfile.readlines()
inputfile.close()

xca = []
yca = []
zca = []
cog = [0,0,0]
for line in lines:
    if "ATOM" in line:
        if line.split()[2] in atomlist:
            xca.append(float(line.split()[5]))
            yca.append(float(line.split()[6]))
            zca.append(float(line.split()[7]))
            cog[0]+=float(line.split()[5])*12.011
            cog[1]+=float(line.split()[6])*12.011
            cog[2]+=float(line.split()[7])*12.011

cog[0]=cog[0]/(len(xca)*12.011)
cog[1]=cog[1]/(len(yca)*12.011)
cog[2]=cog[2]/(len(zca)*12.011)

j=1

for line in lines:
  if "ATOM" in line:
      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[4])),"%11.3f"%(float(line.split()[5])),"%7.3f"%(float(line.split()[6])),"%7.3f"%(float(line.split()[7])), " 1.00  0.00 ",file=outputfile)
      j=j+1

angle_x = 180
trans = [0,0,-80]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatex(x_center,y_center,z_center,angle_x)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[4])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

outputfile.close()
sys.exit()

