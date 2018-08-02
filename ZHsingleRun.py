import arcpy
from arcpy import env
from arcpy.sa import *
import os
arcpy.CheckOutExtension("Spatial") # Check out the ArcGIS Spatial Analyst extension license

#allows to overwrite files if rerunning script 
arcpy.env.overwriteOutput = True

# Set environment settings

inDir = "C:/testGIS/"
outDir = "C:/testGISout/"
env.workspace = "C:/testGIS/"


            # Set local variables
inZoneData = "C:/testGIS/Antigo.shp" 
zoneField = "Location"
inValueRaster = "C:/testGIS/CDL_2014_30m_F1.tif" #inDir+i   I need this to say 'for all files in the directory ending in .tif'
        #inRaster = inDir + i  
outTable = outDir + "antigo" + "2014_20km_ZH" + ".dbf"
print("OutTable = " + outTable) 
               

            # Execute ZonalHistogram
print ("inValueRaster " + inValueRaster)
ZonalHistogram (inZoneData, zoneField, inValueRaster, outTable)