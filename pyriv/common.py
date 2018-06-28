import networkx as nx
import numpy as np
import json
from shapely import ops
from shapely.geometry import LineString, MultiLineString, point, Point, LinearRing
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
import geopandas as gpd
import pandas as pd

def node_rounding(graph, decimal_places=3):
    myround = lambda f: round(f, decimal_places)
    tround = lambda t: tuple([myround(f) for f in t])
    graph = nx.relabel_nodes(graph, tround)
    return graph

def add_geom_edge(graph, path, reverse_too=True):
    """
    Add an edge to graph with the geometry of the path going into
    the edge attributes. The end and/or start nodes should be
    rounded to the same precision as the network you want to attach
    it to.
    
    Parameters
    ----------
      graph : NetworkX.Graph or NetworkX.DiGraph
        A NetworkX Graph that you want to add an edge to using a Shapely geometry class.
        Define reverse_too=True to add reverse edge (reverse_too=False establishes directionality).
      path : Shapely.LineString or Shapely.MultiLineString
        A shapely linestring that gives attributes to the new edge.
    
    Returns
    ----------
      graph : NetworkX.Graph or NetworkX.DiGraph
        The graph with new edge and attributes from line geometry.
    """
    vertlist = list([list(t) for t in path.coords])
    startnode = tuple(vertlist[0])
    endnode = tuple(vertlist[-1])
    pjson = {
        'type': 'LineString',
        'coordinates': vertlist
    }
    edict = {
        'distance': path.length * 0.001,
        'Json': json.dumps(pjson)
    }
    graph.add_edge(startnode, endnode, attr_dict=edict)
    if reverse_too:
        # now add reverse edge
        edict['Json'] = json_linestring_reverse(edict['Json'])
        graph.add_edge(endnode, startnode, attr_dict=edict)
    return graph

def json_linestring_reverse(ls_json):
    """
    Reverses vertices within a JSON-encoded linestring. Used to edit attribute
    dictionary (redefine 'Json' attribute) in reverse edge of NetworkX.Graph
    when you are reversing an edge that's been added by a Shapely Linestring geometry.
    
    Parameters
    ----------
      ls_json : GeoJSON-encoded linestring
        Attribute 'Json' of NetworkX.Graph edge, containing coordinates of a Shapely LineString.
    
    Returns
    ----------
      rev_json : GeoJSON-encoded linestring
        Reverses 'coordinates' attribute of 'Json' within input GeoJSON object
        Defines attributes for a new Networkx.Graph edge when edge has been reversed.
    """
    ed_attr_dict = json.loads(ls_json)
    ed_attr_dict["coordinates"] = ed_attr_dict["coordinates"][::-1]
    rev_json = json.dumps(ed_attr_dict)
    return rev_json

def full_reverse(G):
    """
    This will reverse all of the linestring path between nodes as well as the
    order of all nodes.
    
    Parameters
    ----------
      G : NetworkX.DiGraph
        A directed NetworkX graph.
    
    Returns
    ----------
      G: NetworkX.DiGraph
        A directed NetworkX graph, the reverse of all edges and
        order of nodes in the original directed graph G.
    """
    G = G.reverse()
    for n0,n1 in G.edges_iter():
        try:
            ls_json = G[n0][n1]['Json']
            G[n0][n1]['Json'] = json_linestring_reverse(ls_json)
        except KeyError:
            # This means there's no json path so we don't need to 
            # reverse it.
            pass
    return G

def add_reverse(G):
    """
    Make a fully reversed directional copy and add it to the original. 
    This will let us find the proper path up and down the same river.
    
    Parameters
    ----------
      G: NetworkX.DiGraph
        A directed NetworkX graph. This function reverses the direction of all graph
        edges in a directed graph.
    
    Returns
    ----------
      nx.compose(G, rG) : NetworkX.Graph
        An undirected NetworkX.Graph. The returned graph is undirected because it is
        the union of two directed graphs (G and the reverse of G).
    """
    rG = full_reverse(G)
    return nx.compose(G, rG)

def explode(indf):
    """
    Change a multipolygon geodataframe into a single polygon geodataframe.
    
    Code borrowed from: https://gist.github.com/mhweber/cf36bb4e09df9deee5eb54dc6be74d26
    """
    #indf = gpd.GeoDataFrame.from_file(indata)
    outdf = gpd.GeoDataFrame(columns=indf.columns)
    for idx, row in indf.iterrows():
        if type(row.geometry) == Polygon:
            outdf = outdf.append(row,ignore_index=True)
        if type(row.geometry) == MultiPolygon:
            multdf = gpd.GeoDataFrame(columns=indf.columns)
            recs = len(row.geometry)
            multdf = multdf.append([row]*recs,ignore_index=True)
            for geom in range(recs):
                multdf.loc[geom,'geometry'] = row.geometry[geom]
            outdf = outdf.append(multdf,ignore_index=True)
    return outdf

def path_completion(land_df, path_geom):
    """
    Don't pass in a point that's already outside the land polygon.
    """
    if type(path_geom).__name__ == 'LineString':
        riv_mouth_pnt = Point(path_geom.coords[-1])
        ncp = nearest_coast_pnt(land_df, riv_mouth_pnt)
        ret_path = ops.linemerge((path_geom, LineString((riv_mouth_pnt, ncp))))
    else: # sometimes I want to just pass in a point
        ncp = nearest_coast_pnt(land_df, path_geom)
        ret_path = LineString((path_geom, ncp))
    return ret_path
    
    

def nearest_coast_pnt(land_df, riv_mouth_pnt, return_dist=False):
    """
    in pseudo-code:
        if riv mouth point is outside land poly, just return it.
        if in land poly find nearest point on outer ring of that poly
    """
    in_poly = land_df.intersects(riv_mouth_pnt)
    if not in_poly.any():
        c_pnt = riv_mouth_pnt
        c_dist = 0.0
    else:
        lglr = LinearRing(land_df[in_poly].geometry.iloc[0].exterior)
        c_dist = lglr.project(riv_mouth_pnt)
        c_pnt = lglr.interpolate(c_dist)
    
    if return_dist:
        return riv_mouth_pnt.distance(c_pnt)
    else:
        return c_pnt

def nearest_river_pnt(riv_df, pnt, return_dist=False, threshold=None):
    """
    Find nearest point on a line in riv_df.
    """
    nearest_riv = riv_df.loc[riv_df.distance(pnt).argmin()]
    rivdist = nearest_riv.geometry.distance(pnt)
    
    r_geom = nearest_riv.geometry
    r_dist = r_geom.project(pnt)
    r_pnt = r_geom.interpolate(r_dist)
    
    if return_dist:
        return pnt.distance(r_pnt)
    elif threshold:
        if pnt.distance(r_pnt) > threshold:
            return None
        else:
            return r_pnt
    else:
        return r_pnt

def deadend_distances(shp_file, riv_graph, node_distance=False):
    """
    I think I've superseded this with RiverGraph.river_distances. I'll probably
    delete this. 
    
    Assumes projection units are meters for the moment.
    """
    if type(shp_file).__name__ == 'GeoDataFrame' :
         pnts = shp_file
    else:
         pnts = gpd.read_file(shp_file)
    # function to find closest node for a point geometry
    cl_nd = lambda p: riv_graph.closest_node(point_to_tuple(p))
    path_find = lambda p: riv_graph.shortest_path_to_deadend(cl_nd(p))
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

def point_to_tuple(g):
    """
    Convert `shapely.geometry.point.Point` or geopandas geoseries of length 1
    to x,y coordinate tuple. If given a tuple, list, or numpy array a tuple will
    be returned.
    """
    # this will handle geoseries
    if g.__class__.__name__ == 'GeoSeries':
        if len(g) == 1:
            g = g.iloc[0]
        else:
            raise Exception("Can not convert {} length geoseries to tuple. Only length 1.".format(len(g)))
    elif g.__class__.__name__ == 'GeoDataFrame':
        if len(g) == 1:
            g = g.geometry.iloc[0]
        else:
            raise Exception("Can not convert {} length geodataframe to tuple. Only length 1.".format(len(g)))

    if g.__class__.__name__ in ['tuple', 'list', 'ndarray']:
        result = tuple(g)
    else:
        # this will handle shapely `Point` geometries
        result = tuple(np.array(g))
    return result

def extract_poly_exterior_lines(geom):
    """Extract the exterior lines from a polygon or multipolygon.
    
    Returns
    -------
    LineString or MultiLineString
        All the exterior rings of geom.
    
    """
    if geom.type == 'Polygon':
        exterior_lines = np.array(geom.exterior.coords)
        return LineString(exterior_lines)
    elif geom.type == 'MultiPolygon':
        exterior_lines = []
        for part in geom:
            epc = extract_poly_exterior_lines(part)  # Recursive call
            exterior_lines.append(epc)
        return MultiLineString(exterior_lines)
    else:
        raise ValueError('Unhandled geometry type: ' + repr(geom.type))

def get_coastline_geom(shape):
    """
    Read a shapefile, convert polygons to lines if necessary, run unary_union 
    on the geometries and return the resulting geometry. In the case of a 
    coastline, this will be a multilinestring of the coast.

    Parameters
    ----------
      shape : geopandas.GeoDataFrame or string
        A geodataframe or the filepath to a shapefile. Can be a polygon or a 
        line shapefile.

    Returns
    -------
      shapely.geometry.MultiLinestring
        Just the geometry. Ready to use for distance calculations.
      crs
        The crs of the geometry.
    """
    if shape is not None:
        # handle filepath strings or geodataframes
        if type(shape) == str:
            cldf = gpd.read_file(shape_fn)
        else:
            # assume it's a GeoDataFrame
            cldf = shape
        
        # handle lines or polygons
        if cldf.geom_type.apply(lambda s: s.find("Polygon") >= 0).all():
            cldf = cldf.set_geometry(cldf.geometry.apply(extract_poly_exterior_lines))
        return cldf.unary_union, cldf.crs
    else:
        return None, None

