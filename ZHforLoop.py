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
inDir = "C:/testGIS/"
outDir = "C:/testGISout/"

# initialize empty dirList
dirList = []

# fill dirList with tif files in 'inDir'
for files in os.listdir(inDir): 
    if files.endswith(".shp"): 
        dirList.append(files)
        print(files)   

# iterate over dirList
for i in dirList:
    print("Calculating Zonal Histogram for " + inDir + i) 

  # Set local variables
    inZoneData = inDir + i
    zoneField = "Location"
    inValueRaster = "C:/testGIS/CDL_2014_30m_F1.tif"
    fileName, fileExtension = os.path.splitext(i) # this takes the filename (i) and splits it at the '.', assigning the first part to an object named 'fileName' and the second part (the extension) to an object named 'fileExtension'. You could call these objects 'mork' and 'mindy' - just thought these names would be more intuitive
    outTable = fileName + "_" + "2014_20km_ZH" + ".dbf"
  
  # Execute ZonalHistogram
    ZonalHistogram (inZoneData, zoneField, inValueRaster, outTable)

    print("Zonal histogram calculated for " + i + ", output table written to " + outTable)
  
print "Process Complete!"