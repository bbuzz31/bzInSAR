#!/usr/bin/env python3

"""
Use the asf_search module to find and download sentinel data
Essentially a clone of ariaDownload
"""
import os, os.path as op
import argparse
import math
import re
from datetime import datetime
import asf_search as asf

def createParser():
    """Download Sentinel-1 products using asf_search

    see: https://github.com/asfadmin/Discovery-asf_search
    """
    parser = argparse.ArgumentParser(description=
        'Command line interface to download GUNW products from the ASF DAAC. '
        'GUNW products are hosted at the NASA ASF DAAC.\nDownloading them '
        'requires a NASA Earthdata URS user login and requires users to add '
        '"GRFN Door (PROD)" and "ASF Datapool Products" to their URS '
        'approved applications.',
        epilog='Examples of use:\n\t ASF_DL-search.py --track 004 '
                '--output count'
                '\n\t ASF_DL-search.py --bbox "36.75 37.225 -76.655 -75.928"'
                '\n\t ASF_DL-search.py -t 004,077 --start 20190101 -o count',
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', dest='output', default='Download', \
                        type=str,
        help='Output type, default is "Download". "Download", "Count", and "Url"'
             '"Kmz" are currently supported. Use "Url" for ingestion to '
             'aria*.py')
    parser.add_argument('-t', '--track', dest='track', default=None, type=str,
        help='track to download; single number (including leading zeros) or '
             'comma separated')
    parser.add_argument('-b', '--bbox', dest='bbox',  default=None, type=str,
        help='Lat/Lon Bounding SNWE, or GDAL-readable file containing '
             'POLYGON geometry.')
    parser.add_argument('-f', '--frame', dest='frame', default=None, type=str,
        help='frame(s); can be a comma separated string')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products',
                        type=str, help='Specify directory to deposit all outputs. Default '\
                     'is "products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default=None, type=str,
        help='Start date as YYYYMMDD; If none provided, starts at beginning '
             'of Sentinel record (2014).')
    parser.add_argument('-e', '--end', dest='end', default=None, type=str,
        help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument('-u', '--user', dest='user', default=None, type=str,
        help='NASA Earthdata URS user login. Users must add "GRFN Door '
             '(PROD)" and "ASF Datapool Products" to their URS approved '
             'applications.')
    parser.add_argument('-p', '--pass', dest='passw', default=None, type=str,
        help='NASA Earthdata URS user password. Users must add "GRFN Door '
             '(PROD)" and "ASF Datapool Products" to their URS approved '
             'applications.')
    parser.add_argument('-nt', '--num_threads', dest='num_threads',
                        default='8', type=str,
        help='Specify number of threads for multiprocessing '
             'download. By default "1". Can also specify "All" to use all '
             'available threads.')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None,
                        type=str, help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument('-i', '--ignore-dupes', action='store_false',
                    help='Download all, regardless of duplicates')
    parser.add_argument('--version', dest='version',  default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods; default: '
             'newest')
    parser.add_argument('-v', '--verbose', dest='v', action='store_true',
        help='Print products to be downloaded to stdout')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(2)

    inps = parser.parse_args(args=iargs)

    if not inps.track and not inps.bbox:
        raise Exception('Must specify either a bbox or track')

    if not inps.output.lower() in ['count', 'kmz', 'kml', 'url', 'download']:
        raise Exception ('Incorrect output keyword. Choose "count", "kmz", '
                         '"url", or "download"')

    inps.output = 'Kml' if inps.output.lower() == 'kmz' or \
            inps.output.lower() == 'kml' else inps.output.title()
    return inps


def make_bbox(inp_bbox):
    """Make a WKT from SNWE or a shapefile"""
    if inp_bbox is None:
        return None
    from shapely.geometry import Polygon
    if op.exists(op.abspath(inps.bbox)):
        from ARIAtools.shapefile_util import open_shapefile
        ring = open_shapefile(inps.bbox, 0, 0).exterior
        poly = Polygon(ring)
    else:
        try:
            S, N, W, E = [float(i) for i in inps.bbox.split()]
            ## adjust for degrees easting / northing (0 - 360 / 0:180)
            if W > 180: W -= 360; print('AdjustedW')
            if E > 180: E -= 360; print('AdjustedE')
            if N > 90: N-=90; S-=90; print('Adjusted N/S')
        except:
            raise Exception('Cannot understand the --bbox argument. '
            'Input string was entered incorrectly or path does not '
            'exist.')

        poly = Polygon([(W, N), (W,S), (E,S), (E, N)])
    return poly


def fmt_dst(inps):
    """Format the save name"""
    ext  = '.kmz' if inps.output == 'Kml' else '.txt'

    if inps.track is not None:
        fn_track = f'track{inps.track}'.replace(',', '-')
    else:
        fn_track = ''

    if inps.bbox is not None:
        WSEN = make_bbox(inps.bbox).bounds
        WSEN_fmt = []
        for i, coord in enumerate(WSEN):
            if i < 2:
                WSEN_fmt.append(math.floor(float(coord)))
            else:
                WSEN_fmt.append(math.ceil(float(coord)))
        fn_bbox = f'_bbox{WSEN_fmt[0]}W{WSEN_fmt[1]}S{WSEN_fmt[2]}E{WSEN_fmt[3]}N'
    else:
        fn_bbox = ''

    dst   = op.join(inps.wd, f'{fn_track}{fn_bbox}_0{ext}'.lstrip('_'))
    count = 1 # don't overwrite if already exists
    while op.exists(dst):
        basen  = f'{re.split(str(count-1)+ext, op.basename(dst))[0]}' \
                 f'{count}{ext}'
        dst    = op.join(op.dirname(dst), basen)
        count += 1
    return dst


class Downloader(object):
    def __init__(self, inps):
        self.inps            = inps
        self.inps.output     = self.inps.output.title()
        self.inps.wd         = op.abspath(self.inps.wd)
        os.makedirs(self.inps.wd, exist_ok=True)


    def __call__(self):
        scenes   = self.query_asf()

        if self.inps.output == 'Count':
            print(f'Found -- {len(scenes)} -- products')

        elif self.inps.output == 'Kml':
            dst    = fmt_dst(inps)
            print('Kml option is not yet supported.')

        elif self.inps.output == 'Url':
            urls =  [s.geojson()['properties']['url'] for s in scenes]
            dst  = fmt_dst(inps)
            with open(dst, 'w') as fh: [print(url, sep='\n', file=fh) \
                                        for url in urls]
            print(f'Wrote -- {len(urls)} -- product urls to: {dst}')

        elif self.inps.output == 'Download':
            ## turn the list back into an ASF object
            fns    =  [s.geojson()['properties']['sceneName'] for s in scenes]
            if self.inps.ignore_dupes:
                new    = []
                for i, f in enumerate(fns):
                    dst = op.join(self.inps.wd, f'{f}.zip')
                    if op.exists(dst):
                        print (f, 'already exists, skipping...')
                        continue
                    else:
                        new.append(scenes[i])

            scenes = asf.ASFSearchResults(new)
            nt     = int(self.inps.num_threads) # so legacy works
            ## allow a user to specify username / password
            if self.inps.user is not None:
                session = asf.ASFSession()
                session.auth_with_creds(self.inps.user, self.inps.passw)
                scenes.download(self.inps.wd, processes=nt, session=session)

            else:
                scenes.download(self.inps.wd, processes=nt)
            print(f'Wrote -- {len(scenes)} -- products to: {self.inps.wd}')

        return


    def query_asf(self):
        """Get the scenes from ASF"""
        st     = datetime.strptime(self.inps.start, '%Y%m%d') \
                    if self.inps.start is not None else None
        en     = datetime.strptime(self.inps.end, '%Y%m%d') \
                    if self.inps.end is not None else None
        bbox   = make_bbox(self.inps.bbox)

        if self.inps.track is not None:
            tracks = self.inps.track.split(',')
            tracks = [int(track) for track in tracks]
        else:
            tracks = self.inps.track

        if self.inps.frame is not None:
            frames = self.inps.frame.split(',')
            frames = [int(frame) for frame in frames]
        else:
            frames = self.inps.frame

        dct_kw = dict(platform=asf.constants.SENTINEL1,
                    processingLevel=asf.PRODUCT_TYPE.SLC,
                    relativeOrbit=tracks, frame=frames,
                    start=st, end=en,
                    lookDirection=self.inps.flightdir,
                    intersectsWith=bbox)
        # print (dct_kw)
        scenes = asf.geo_search(**dct_kw)
        return scenes


if __name__ == '__main__':
    inps = cmdLineParse()
    Downloader(inps)()
