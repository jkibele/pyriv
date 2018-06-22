import geopandas as gpd
import pandas as pd
from river_graph import point_to_tuple
from shapely.geometry import Point

class RiverDist(gpd.GeoDataFrame):
    def __init__(self, *args, **kwargs):
        self = super(RiverDist, self).__init__(*args, **kwargs)

    @property
    def points_df(self):
        obj_cols = self.select_dtypes(include=['object']).columns
        cols_to_drop = [c for c in ['path', 'nearest_node'] if c in obj_cols]
        return self.drop(cols_to_drop, axis=1)
    
    def save_points(self, outfile):
        self.points_df.to_file(outfile)

    @property
    def paths_df(self):
        outdf = self.set_geometry('path', drop=True)
        if 'nearest_node' in self.columns:
            outdf.drop(['nearest_node'], axis=1, inplace=True)
        zero_ind = outdf.index[outdf.rivdist_km==0.0]
        outdf.drop(zero_ind)
        return outdf  

    def save_paths(self, outfile):
        self.paths_df.to_file(outfile)

    def summary_plot(self, **kwargs):
        ax = self.paths_df.plot(color='g', zorder=9, label='Paths', **kwargs)
        ax = self.points_df.plot(ax=ax, marker='*', markersize=150, c='r', 
                                label='Start Point', linewidth=0, zorder=8)
        return ax