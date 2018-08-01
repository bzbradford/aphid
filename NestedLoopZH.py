import arcpy
from arcpy import env
from arcpy.sa import *
import os
from os.path import basename, splitext
# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
arcpy.ListFeatureClasses()

#allows to overwrite files if rerunning script 
arcpy.env.overwriteOutput = True

# Set environment settings
env.workspace = "C:/testGIS/"
buffers = "C:/testGISout/buffers/"
outDir = "C:/testGISout/ZH/"
CDLs = "C:/testGIS/CDLRaster/"


RasterList = []

for files in os.listdir(CDLs): 
    if files.endswith(".tif"): 
        RasterList.append(files)
        print(files)

for j in RasterList:


# initialize empty dirList
	dirList = []

	# fill dirList with shp files in 'inDir'
	for files in os.listdir(buffers): 
	    if files.endswith(".shp"): 
	        dirList.append(files)
	        print(files)   

	# iterate over dirList
	for i in dirList:
	   	print("Calculating Zonal Histogram for " + buffers + i) 

	 	# Set local variables
	   	inZoneData = buffers + i
	   	zoneField = "Site2"
		inValueRaster = CDLs + j
	   	fileName, fileExtension = os.path.splitext(i)
	   	RasterYear, fileExtension = os.path.splitext(j[:8]) # this takes the filename (i) and splits it at the '.', assigning the first part to an object named 'fileName' and the second part (the extension) to an object named 'fileExtension'. You could call these objects 'mork' and 'mindy' - just thought these names would be more intuitive
	   	outTable = fileName + RasterYear + "_" + "ZH" + ".dbf"
	  
	  # Execute ZonalHistogram
	   	ZonalHistogram (inZoneData, zoneField, inValueRaster, outTable)

	   	print("Zonal histogram calculated for " + i + j + ", output table written to " + outTable)
	  