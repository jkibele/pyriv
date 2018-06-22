import geopandas as gpd
import pandas as pd
from river_graph import point_to_tuple
from shapely.geometry import Point

def river_distances(pnts_gdf, riv_graph, node_distance=False):
    """
    Assumes projection units are meters for the moment.
    """
    if pnts_gdf.__class__.__name__ != 'GeoDataFrame':
        pnts = gpd.read_file(shp_file)
    else:
        pnts = pnts_gdf
    # function to find closest node for a point geometry
    cl_nd = lambda p: riv_graph.closest_node(point_to_tuple(p))
    if riv_graph.coastline is None:
        path_find = lambda p: riv_graph.shortest_path_to_deadend(cl_nd(p))
    else:    
        path_find = lambda p: riv_graph.shortest_path_to_coast(cl_nd(p))
    pnts['path'] = pnts.geometry.apply(path_find)
    pth_dist = lambda p: p.length * 0.001
    pnts['rivdist_km'] = pnts.path.apply(pth_dist)
    if node_distance:
        pnts['nearest_node'] = pnts.geometry.apply(lambda p: Point(cl_nd(p)))
        d = {}
        for i, row in pnts.iterrows():
            d[i] = row['geometry'].distance(row['nearest_node'])
        dser = pd.Series(d, name='node_dist')
        pnts = pnts.join(dser)
    return pnts

def save_points(riv_dist_gdf, outfile):
    obj_cols = riv_dist_gdf.select_dtypes(include=['object']).columns
    cols_to_drop = [c for c in ['path', 'nearest_node'] if c in obj_cols]
    riv_dist_gdf.drop(cols_to_drop, axis=1).to_file(outfile)

def save_paths(riv_dist_gdf, outfile):
    outdf = riv_dist_gdf.set_geometry('path', drop=True)
    if 'nearest_node' in riv_dist_gdf.columns:
        outdf.drop(['nearest_node'], axis=1, inplace=True)
    zero_ind = outdf.index[outdf.rivdist_km==0.0]
    outdf.drop(zero_ind).to_file(outfile)