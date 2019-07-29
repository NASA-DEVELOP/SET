from osgeo import gdal
#from osgeo import gdal_array
import numpy as np
from PIL import Image
#from matplotlib import pyplot as plt

fp000 = "kernel_30.23151_1.0_80.0_0.0.tif"
fp180 = "kernel_30.23151_1.0_80.0_-180.0.tif"

ds_k000 = gdal.Open(fp000)
ds_k180 = gdal.Open(fp180)

print("GUIS Kernel 000 Metadata")
print(ds_k000.GetMetadata())
print("GUIS Kernel 180 Metadata")
print(ds_k180.GetMetadata())

print("[ RASTER BAND COUNT ]: ", ds_k000.RasterCount)
bnd000 = ds_k000.GetRasterBand(1)
nodata000 = bnd000.GetNoDataValue()
print(type(bnd000))
print("[ RASTER BAND COUNT ]: ", ds_k180.RasterCount)
bnd180 = ds_k180.GetRasterBand(1)
nodata180 = bnd180.GetNoDataValue()
print(type(bnd180))

arr000 = bnd000.ReadAsArray()
print(type(arr000))
# arr000 = np.ma.masked_equal(arr000, nodata000)
arr000[np.isnan(arr000)] = 9999
print(arr000.min())
arr180 = bnd180.ReadAsArray()
print(type(arr180))
# arr180 = np.ma.masked_equal(arr180, nodata180)
arr180[np.isnan(arr180)] = 999
print(arr180.min())

arr000flip = np.flipud(arr000)
arr_symchk = np.divide(arr180, arr000flip)
print(arr_symchk.min())
print(arr_symchk.max())

# img = Image.fromarray(arr_symchk)
# #img.save('symmetry_check.png')
# img.show()
