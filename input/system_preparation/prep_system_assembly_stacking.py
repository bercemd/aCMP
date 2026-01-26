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

#Generating the assembly of 12 hexagonal units from a single hexagonal assembly unit

outputfile = open("assembly_dodeca_stacked.pdb","w")

inputfile = open("assembly_dodeca_min.pdb","r")
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
            xca.append(float(line[30:38]))
            yca.append(float(line[38:46]))
            zca.append(float(line[46:54]))
            cog[0]+=float(line[30:38])*12.011
            cog[1]+=float(line[38:46])*12.011
            cog[2]+=float(line[46:54])*12.011

cog[0]=cog[0]/(len(xca)*12.011)
cog[1]=cog[1]/(len(yca)*12.011)
cog[2]=cog[2]/(len(zca)*12.011)

print(cog)
print(len(xca))
print(len(yca))
print(len(zca))


spacing = 16.5
j=1

#Unit 1

for line in lines:
  if "ATOM" in line:
      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(float(line[30:38])),"%7.2f"%(float(line[38:46])),"%7.2f"%(float(line[46:54])), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 2

angle_y = 180
d = 1*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 3

angle_y = 180
d = -1*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 4

d = 2*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_new,y_new,z_new = translate(x_center,y_center,z_center,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 5

d = -2*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_new,y_new,z_new = translate(x_center,y_center,z_center,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 6

angle_y = 180
d = 3*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 7

angle_y = 180
d = -3*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 8

d = 4*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_new,y_new,z_new = translate(x_center,y_center,z_center,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 9

d = -4*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_new,y_new,z_new = translate(x_center,y_center,z_center,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 10

angle_y = 180
d = 5*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 11

angle_y = 180
d = -5*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_rot,y_rot,z_rot = rotatey(x_center,y_center,z_center,angle_y)
      x_new,y_new,z_new = translate(x_rot,y_rot,z_rot,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1

#Unit 12

d = 6*spacing
trans = [d,0,0]

for line in lines:
  if "ATOM" in line:
      x_center,y_center,z_center = center(float(line[30:38]),float(line[38:46]),float(line[46:54]),cog[0],cog[1],cog[2])
      x_new,y_new,z_new = translate(x_center,y_center,z_center,cog[0],cog[1],cog[2],trans[0],trans[1],trans[2])

      print("ATOM","%6i"%j,line[12:21],"%4i"%(int(line.split()[5])),"%10.2f"%(x_new),"%7.2f"%(y_new),"%7.2f"%(z_new), " 1.00  0.00 ",file=outputfile)
      j=j+1



outputfile.close()
sys.exit()


