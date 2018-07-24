from osgeo import gdal

def main():
  lat, lon = 32, 24
  filein = 'GroundTruth/anthlightmags_wgs.tif'
  hem_raster = gdal.Open(filein)
  transform = hem_raster.GetGeoTransform()
  data = hem_raster.ReadAsArray()
  x_origin, y_origin = transform[0], transform[3]
  px_width, px_height = transform[1], transform[5]
  print(x_origin, y_origin, px_width, px_height)


  row = int((lat - y_origin)/px_height)
  col = int((lon - x_origin)/px_width)
  val = data[row][col]
  print(val)


main()
