# GeoRasterGrid

Many geospatial datasets are very large and are often broken down into tiles 
so that users may load only the geographic region relevant to their work.

A few data products that implement this are:

* [NASA's Blue Marble DEM](https://visibleearth.nasa.gov/images/73934/topography) (included; see `./data`)
* [GMTED2010](https://topotools.cr.usgs.gov/gmted_viewer/viewer.htm)
* [Global Surface Water](https://global-surface-water.appspot.com/download)
* and many more

This tiled approach is more memory-efficient, but it can be annoying to work with 
for users who only need to look up data by lat-lon.  This class provides a 
facility to access tiled geospatial data products by lat-lon, and internally 
it handles tile lookup, loading, caching, and interpolation.

## Quickstart

The easiest way to start is with the Blue Marble dataset included in the `./data` folder.
if no arguments are provided to `GeoRasterGrid`, it will look for data in that folder. So, 
these two commands will both produce the same result:

```
blue_marble = fileparts(which('gebco_08_rev_elev_A1_grey_geo.tif'))
assert(~isempty(blue_marble), 'Expected ./data/blue-marble to be on the MATLAB path');

% load the map by providing the folder containing the GeoTIFF files
map = GeoRasterGrid(blue_marble)

% with no input args, it defaults to the included dataset (same as above)
map = GeoRasterGrid()
```

With the map object created, you can call `map.show` to see the tile boundaries 
overlaid on an RGB map of the Earth:

