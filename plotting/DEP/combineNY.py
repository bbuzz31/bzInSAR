"""
Combine the VLM from NYC and NNJ to get Sandy Hook
"""

from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE.plotting.plotExp import PlotExp
import xarray as xr
from BZ import bbPlot


NY  = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
NYP = PlotExp(NY, kind='unc', show_gps=True)

NNJ = ExpBase(NNJ_SR, 'ERA5_SET_PM_ex_Fast', 'NJHT', neofs=15)
NJP = PlotExp(NNJ, kind='rate', show_gps=True)

# make grids that spans both datasets
da_ny0 = NYP.da.copy()
da_ny = NYP.da
da_nj = NJP.da

coords = {'lat': da_ny.lat.astype(np.float32), 'lon': da_ny.lon.astype(np.float32)}
da_ny  = da_ny.assign_coords(coords)

coords2 = {'lat': da_nj.lat.astype(np.float32), 'lon': da_nj.lon.astype(np.float32)}
da_nj  = da_nj.assign_coords(coords2)

xs = (da_ny.lon[1] - da_ny.lon[0]).item()
ys = (da_ny.lat[1] - da_ny.lat[0]).item()

W1, E1 = np.min([da_ny.lon.min(), da_nj.lon.min()]), np.max([da_ny.lon.max(), da_nj.lon.max()])
S1, N1 = np.min([da_ny.lat.min(), da_nj.lat.min()]), np.max([da_ny.lat.max(), da_nj.lat.max()])

lons = np.arange(W1, E1+xs, xs)
lats = np.arange(S1, N1+ys, ys)
da_grid = xr.DataArray(np.meshgrid(lons, lats)[0], coords={'lat': lats, 'lon':lons},
            dims='lat lon'.split())

da_nyi  = da_ny.interp_like(da_grid, method='nearest')
da_nji  = da_nj.interp_like(da_grid, method='nearest')

da_nji  = da_nji.where(da_nji.lat<40.55)
# da_nji  = da_nji.where(da_nji.lon>-74.0254)
da_nji  = da_nji.where(da_nji.lon>-74.145226)


## check; should be the same but i'll let them be sorta close i guess
lat, lon = 40.74829, -73.901642
print (da_ny.sel(lat=lat, lon=lon, method='nearest'))
print (da_nyi.sel(lat=lat, lon=lon, method='nearest'))

da_m = da_nyi.where(~da_nyi.isnull(), da_nji)
# breakpoint()

# fig, axes = NYP.plot_basic(da=da_m)
fig, axes = NYP.plot_basic(da=da_m)
S, N, W, E = NYP.SNWE
axes.set_extent([W, E, 40.4370, N])
# axes.set_extent([W, E, S, N])
fig.set_label(f'{fig.get_label()}+nnj')
bbPlot.savefigs(PATH_RES, True, True)
plt.show()
