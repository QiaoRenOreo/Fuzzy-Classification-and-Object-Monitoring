import gdal
from osgeo import gdal
from osgeo import ogr
from osgeo import gdal_array as gdarr
import matplotlib.pyplot as plt
import os
import numpy

# step 1. read transectline.shp
dataDirectory=r'D:\Study\Module8\IndividualAssignment\ObjectMonit_step1_output' #set data directory
os.chdir(dataDirectory)

datasource_TransectLine = ogr.Open("transectline.shp") # open raster file, set datasource
layerCount = datasource_TransectLine.GetLayerCount()
layer_TransectLine = datasource_TransectLine.GetLayer(0) # transectline only has one layer

# the buffer should have the same spatial reference as the transect line.
# so we need to fetch the spatial reference of transect line
srs=layer_TransectLine.GetSpatialRef()

# step2. get feature from transect line
for feature in layer_TransectLine:
    transcetLine_feature=feature
geometryRef_TransectLine = transcetLine_feature.GetGeometryRef()

# step3. create buffer
bufferDistance = 60 # buffer distance is set to be 60 meters
buffer = geometryRef_TransectLine.Buffer(bufferDistance) #create the buffer
areaBuffer=buffer.Area()

#create a file to save buffer
driver = ogr.GetDriverByName("ESRI Shapefile")
data_source = driver.CreateDataSource("transectLine_buffer.shp")
layer = data_source.CreateLayer("transectLine_buffer", srs, ogr.wkbPolygon)
field_name = ogr.FieldDefn("Name", ogr.OFTString)
field_name.SetWidth(24)
layer.CreateField(field_name)
field_area = ogr.FieldDefn("Area", ogr.OFTReal)
field_area.SetWidth(32)
field_area.SetPrecision(2)
layer.CreateField(field_area)
feature = ogr.Feature(layer.GetLayerDefn())
feature.SetField("Name", 'name1') # go to the attribute field "Name", set its value as 'name1'
feature.SetField("Area", areaBuffer) # go to the attribute field "Area", set its value as areaBuffer
feature.SetGeometry(buffer)
layer.CreateFeature(feature)
feature = None
data_source = None
##########################

# step5. subset 2013 and 2017 image from the buffer
# firstly read file 2013 and 2017
# secondly, use buffer as a mask. mask the pixels inside buffer as 1. outside buffer as 0.
# thirdly, rasterize the mask
# fourthly, do multiplication to execute subset

nodatavalue = 0
tempDataset13 = gdal.Open('FCM_lc820130827_s5_run_1.tif')
tempDataset17 = gdal.Open('FCM_lc820170619_s5_run_1.tif')

driver = gdal.GetDriverByName('MEM') # memory raster

# get raster size
x = tempDataset13.RasterXSize
y = tempDataset13.RasterYSize

# define the new output raster (I want an integer raster with 1 bands)
maskRaster = driver.Create('Memory_Image', x, y, 1, gdal.GDT_Int16)

maskRaster.SetProjection(tempDataset13.GetProjection())
maskRaster.SetGeoTransform(tempDataset13.GetGeoTransform())

#Get mask band and set the nodata value
maskBand = maskRaster.GetRasterBand(1)
maskBand.SetNoDataValue(nodatavalue)

# open the vector file and get the first layer
vector_ds = ogr.Open("transectLine_buffer.shp")
vectorLayer = vector_ds.GetLayer()
vectorLayer.SetAttributeFilter("Name = 'name1'")


# rasterize the mask layer. Because tif is raster, buffer is vector.
# we have to convert vector to raster.
# Otherwise we can not do the further computation.

gdal.RasterizeLayer(maskRaster, [1], vectorLayer, burn_values=[1])

maskPX = gdarr.DatasetReadAsArray(maskRaster, 0, 0, x, y)
tempPx13 = gdarr.DatasetReadAsArray(tempDataset13, 0, 0, x, y)
tempPx17 = gdarr.DatasetReadAsArray(tempDataset17, 0, 0, x, y)

# do the multiplication:
# membership value* mask = either store as 0 or store the membership value.
# output is an array

pxResult13=numpy.multiply(maskPX,tempPx13)
pxResult17=numpy.multiply(maskPX,tempPx17)
print("hey!")


# remove zeros from the array. convert array to a one dimentional list
pxResult13=numpy.extract(pxResult13>0,pxResult13)
pxResult17=numpy.extract(pxResult17>0,pxResult17)


# step6. plot and add labels and legends
plt.xlabel("pixel number")
plt.ylabel("membership value")
plt.grid()
plt.plot(pxResult13,":",color='red',label="2013")
plt.plot(pxResult17,":",color='green',label="2017")
plt.title("Plot Of Membership Value Along The Transect Line In 2013 And 2017")
plt.legend(loc="lower left")
plt.show()


# step7. flush the cache and clean memory
layer = None
maskBand.FlushCache()
maskBand = None
dataset = None
outraster = None
source_ds = None

