# import system modules 
import arcpy
from arcpy import env
import os
from os.path import basename, splitext


#first split aphid traps into individual shapefiles, then create the various buffer distances needed with this script
# Set environment settings
env.workspace = "C:/testGISout/buffers/buffers1km/"
inDir = "C:/testGIS/points/"
outDir = "C:/testGISout/buffers/buffers1km/"

dirList = []

# fill dirList with tif files in 'inDir'
for files in os.listdir(inDir): 
    if files.endswith(".shp"): 
        dirList.append(files)
        print(files)   

# iterate over dirList
for i in dirList:

	SuctionTrap = inDir + i
	fileName, fileExtension = os.path.splitext(i)
	PointBuffer = fileName +"_"+ "1kmbuff.shp"
	bufferDistance = "1 Kilometer"
	
	arcpy.Buffer_analysis(SuctionTrap, PointBuffer, bufferDistance)