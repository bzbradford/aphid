### setup ###
import arcpy
import os

# set environment variables
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False
arcpy.env.workspace = "c:/arc/default.gdb"
pts = "c:/arc/aphids/buffers/aphid_pts.shp"
sites = "c:/arc/aphids/buffers/sites"
buffers = "c:/arc/aphids/buffers/buffers"
tables = "c:/arc/aphids/buffers/tables"
rasters = "d:/gis_data/cropscape/national/filtered"


### section 1 ###

# split trap shapefile
inFC = pts
field = "SiteID"
outNamePrefix = "Site_"
outDir = sites
arcpy.ImportToolbox("D:/arc/tools/SplitLayerByAttributes.tbx")
arcpy.SplitLayerByAttributes(inFC, field, outNamePrefix, outDir)

# iterate trap shps and buffer them
dists = [500, 1000, 2000, 5000, 10000]
distUnit = " Meters"
bufferDistances = [i + distUnit for i in map(str, dists)]

for shp in os.listdir(sites):
	if shp.endswith(".shp"):
		for dist in bufferDistances:
			inFile = os.path.join(sites, shp)
			outFileName = shp.replace("Site_", "").replace(".shp", "_") + dist.replace(distUnit, "m") + ".shp"
			outFile =  os.path.join(buffers, outFileName)
			if arcpy.Exists(outFile):
				print(outFile + " exists!")
			else:
				arcpy.Buffer_analysis(inFile, outFile, dist)
				print(outFile + " created!")

# iterate rasters and ZH them
# generate file lists
buffer_list = []
for file in os.listdir(buffers):
	if file.endswith(".shp"):
		buffer_list.append(os.path.join(buffers, file))

raster_list = []
for file in os.listdir(rasters):
	if file.endswith(".tif"):
		raster_list.append(os.path.join(rasters, file))

# run zh
for raster in raster_list:
	rasterYear = os.path.basename(raster).split("_")[1]
	inValueRaster = raster
	for shp in buffer_list:
		inZoneData = shp
		ZoneField = "SiteID"
		outTable = os.path.join(tables, os.path.basename(shp).replace(".shp", "") + "_" + rasterYear + ".dbf")
		# print(inZoneData)
		# print(ZoneField)
		# print(inValueRaster)
		# print(outTable)
		# print("\n")
		if arcpy.Exists(outTable):
			print(outTable + " exists!")
		else:
			arcpy.sa.ZonalHistogram(inZoneData, ZoneField, inValueRaster, outTable)
			print(outTable + " created!")
print("Processing complete.")

