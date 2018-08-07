### setup ###
import arcpy
import os

# set environment variables
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False
arcpy.env.workspace = "C:/ArcGIS/Default.gdb"
pts = "L:/aphids/gis/aphid_pts.shp"
sites = "L:/aphids/gis/sites"
buffers = "L:/aphids/gis/buffers"
tables = "L:/aphids/gis/tables"
rasters = "D:/gis_data/cropscape/national/filtered"


### split trap shapefile ###

inFC = pts
field = "SiteID"
outNamePrefix = "Site_"
outDir = sites
arcpy.ImportToolbox("L:/arc/tools/SplitLayerByAttributes.tbx")
arcpy.SplitLayerByAttributes(inFC, field, outNamePrefix, outDir)


### Generate buffers ###

# set variables
dists = [500, 1000, 2000, 5000, 10000]
dists = [20000]
distUnit = " Meters"
bufferDistances = [i + distUnit for i in map(str, dists)]

# create buffers
for shp in os.listdir(sites):
	if shp.endswith(".shp"):
		for dist in bufferDistances:
			inFile = os.path.join(sites, shp)
			outFileName = shp.replace("Site_", "").replace(".shp", "_") + dist.replace(distUnit, "m") + ".shp"
			outFile =  os.path.join(buffers, outFileName)
			print("Creating buffer for " + inFile)
			if arcpy.Exists(outFile):
				print(outFile + " exists!")
			else:
				arcpy.Buffer_analysis(inFile, outFile, dist)
				print(outFile + " created!")


### Run Zonal Histogram on buffers with rasters

# generate buffer file list
buffer_list = []
for file in os.listdir(buffers):
	if file.endswith(".shp"):
		buffer_list.append(os.path.join(buffers, file))

# generate raster file list
raster_list = []
for file in os.listdir(rasters):
	if file.endswith(".tif"):
		raster_list.append(os.path.join(rasters, file))

# run zonal histograms
for raster in raster_list:
	rasterYear = os.path.basename(raster).split("_")[1]
	inValueRaster = raster
	for shp in buffer_list:
		inZoneData = shp
		ZoneField = "SiteID"
		outTable = os.path.join(tables, os.path.basename(shp).replace(".shp", "") + "_" + rasterYear + ".dbf")
		if arcpy.Exists(outTable):
			print(outTable + " exists!")
		else:
			arcpy.sa.ZonalHistogram(inZoneData, ZoneField, inValueRaster, outTable)
			print(outTable + " created!")

print("Processing complete.")
