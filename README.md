# map2oruxmaps
Convert raster maps to OruxMaps format. Supports assignment of overview maps.

Input maps are read by GDAL (version >= 2) and converted to JPEG tiles. The maps are sorted by scale and an appropriate zoom level is assigned.
Datum and projection strings are translated to OruxMap's format. I have included UTM and Gauss-Kruger. Feel free to add more.
The code is a fork of map2rmap from QLandkarteGT by Oliver Eichler.

```
usage: map2oruxmaps -r <algo> -q <1..100> -s <411|422|444> -m <memory> -c -n <mapname> <file1> <file2> ... <fileN> <output directory>

  -r    Resampling algorithm (nearest, bilinear, cubic, cubicspline, lanczos, average, mode, gauss) 
  -q    The JPEG quality from 1 to 100. Default is 75
  -s    The chroma subsampling. Default is 411
  -m    Sets GDAL_CACHEMAX memory (for big maps).
  -c    Crop all zoom levels to the base map.
```
