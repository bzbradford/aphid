# import system modules 
import arcpy
from arcpy import env


#first split aphid traps into individual shapefiles, then create the various buffer distances needed with this script
# Set environment settings
env.workspace = "C:/testGIS/"
inDir = "C:/testGIS/points/"
outDir = "C:/testGISout/buffers/buffers1km/"



SuctionTrap = "C:/testGIS/points/IA_Ames.shp"
PointBuffer = "C:/testGISout/buffers/buffer1km/buffer1km.gdb/buffer_output4"
bufferDistance = "1 Kilometer"
line_side = "FULL"
endType = "ROUND"
dissolveOption = "NONE"
dissolveField = "Distance"

#arcpy.Buffer_analysis(SuctionTrap, PointBuffer, bufferDistance, line_side, endType, dissolveOption, dissolveField)

arcpy.Buffer_analysis("C:/testGIS/points/IA_Ames.shp", "C:/testGISout/buffers/buffer1km/buffer1km.gdb/buffer_output3", "1 Kilometer")
arcpy.Buffer_analysis(SuctionTrap, PointBuffer, bufferDistance)