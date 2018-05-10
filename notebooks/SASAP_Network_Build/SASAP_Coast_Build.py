import pandas as pd
import geopandas as gpd
import os
from pyriv.coastal import Land
import networkx as nx

## Setup: ######
data_dir = '/home/shares/scientist/pyriv/shapefiles/'
fullpath = lambda s: os.path.join(data_dir, s)
full_land_fn = fullpath('SASAP_Coastline_Polygon.shp')
simp_num = 1000
n_processors = 38

land_out_fn = fullpath('SASAP_coast_poly{}m.shp'.format(simp_num))
gpic_out_fn = fullpath('SASAP_coast{}m.gpickle'.format(simp_num))

# read in the full resolution land shapefile
full_land = gpd.read_file(full_land_fn)
# dump small polygons and simplify the remains
lowres_land = full_land.set_geometry(full_land.simplify(simp_num).buffer(0))
lowres_land.to_file(land_out_fn)
# generate the graph
lnd = Land(land_out_fn)
ecgraph = lnd.graph(n_jobs=n_processors)
nx.write_gpickle(ecgraph, gpic_out_fn)