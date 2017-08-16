import sys
sys.path.append("/Users/haleighwright/Desktop/NCEAS/pyriv")
sys.path.append("/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages")
import networkx as nx
import numpy as np
import sys
import geopandas as gpd
import seaborn # makes matplotlib graphs look prettier
from shapely.geometry import Point
from multiprocessing import Pool

from pyriv import graph_prep as GraphBuilder
from pyriv import river_graph as RiverGraph
from pyriv import snapping_processor as SnapTool

%matplotlib inline
%config InlineBackend.figure_format = 'retina' # makes matplotlib graphs 2x resolution for retina display