import geopandas as gpd

from BZ import regions, bbGIS

from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

"""
Clean the bridges inventory:
    https://hifld-geoplatform.opendata.arcgis.com/datasets/a9b05a595ff94f3fa3888d1240545740

Column names:
    https://www.fhwa.dot.gov/bridge/nbi/format.cfm
"""


if __name__== '__main__':
    reg = 'EC'

    ## this takes awhile its big
    src = Path(os.getenv('dataroot')) / 'GIS' / 'Bridges' / 'National_Bridge_Inventory.GeoJSON'
    gdf = gpd.read_file(src)

    # keep these columns
    to_keep = []
    for col in [col.lower() for col in gdf.columns]:
        for k in 'bridge clr type kind feat deck lat lon geometry'.split():
            if k in col:
                to_keep.append(col)
    to_keep  = pd.Series(to_keep).drop_duplicates().to_list()

    gdf_keep = gdf[to_keep]

    gdf_crop = bbGIS.in_dfll(gdf, SNWE=regions[reg])

    dst = src.parent / f'National_Bridge_Inventory_{reg}.GeoJSON'
    gdf_crop.to_file(dst)
    print ('Wrote:', dst)