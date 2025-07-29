from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
import xarray as xr
import pygmt

gdal.UseExceptions()


class PlotVLM(object):
    def __init__(self, region, path_vup, ref_sta=[]):
        self.region         = region
        self.path_vup       = path_vup
        self.path_transects = op.join(self.path_vup, 'Transects')
        self.ref_stas       = DCT_REG[region][3] if not ref_sta else ref_sta

        self.SNWE = DCT_REG[region][0]
        self.WESN = np.array([*self.SNWE[2:], *self.SNWE[:2]])

        if region == 'HR':
            buff = 0.0025
            xy   = 'n0.75/0.5'
            self.cbar_pos  = f'{xy}+w2.0c/0.20c'

            self.clims  = [(-5, 3, 2), (0, 3, 1), (0, 200, 50)] # rate, unc, residue
            self.ti_pos = [(-76.20, 37.2), (-76.20, 37.18)]
            ng   = 0
            ngU  = 0

            # gps markers
            self.just_ref  = 'BL'
            self.just_non  = 'BR'
            self.ref_stas.extend('LS03 SPVA LOYZ'.split())
            parms  = {'L': True}
            parmsU = {'L': True}
            parmsR = {'L': True}

        elif region in 'Charleston SC'.split():
            self.WESN = np.array([-80.047996, -79.732874, 32.641118, 32.892924])
            buff = 0.0075 # buffer around map
            xy   = 'n0.475/0.05'                                 # loc of colorbar
            ng   = 1
            ngU  = 1
            self.cbar_pos  = f'{xy}+w3.0c/0.20c+h-Np'                 # for multipanel

            offsets   = [['-0.1/-0', '0.025/0.025', '0.025/0.025'], # transect label st
                        ['0.15/-0.1', '0.2/0.05', '0.025/0.025']]   # transect label en

            self.clims = ([-3, 3, 2], [0, 3, 1], [0, 200, 50])
            self.ti_pos = [(-79.85, 32.66), (-79.9, 32.89)] # Vertical rate mm/yr

            self.just_ref  = 'TL' # justification of GPS reference labels
            self.just_non  = 'BR' # justification of GPS non-reference labels
            parms  = {'L': True} #if not  dict(frame=f'xa{ng}')
            parmsU  = {'L': True} #if not  dict(frame=f'xa{ng}')
            parmsR = {'L': True}


        elif region == 'Savannah':
            buff      = 0.0025
            self.WESN = np.array([-81.196693, -80.98298,32.006252, 32.151274])
            # xy        = 'n0.075/-0.125'
            xy        = 'n0.250/0.955'
            self.cbar_pos  = f'{xy}+w2.50c/0.15c+h -Np'
            offsets     = [['-0.015/-0.04', '-0.02/0', '0.03/-0.035'],
                        ['0.05/0.1', '0.05/0', '0.02/-0.010']]
            self.clims  = ([-3, 3, 1], [0, 3, 1], [0, 200, 50])
            self.ti_pos = [(-81.130, 32.133), (-81.100, 32.133)] # Vertical rate mm/yr


            self.just_ref  = 'BL' # justification of GPS reference labels
            self.just_non  = 'BR' # justification of GPS non-reference labels
            parms  = {'L': True} #if not  dict(frame=f'xa{ng}')
            parmsU  = {'L': True} #if not  dict(frame=f'xa{ng}')
            parmsR = {'L': True}


        elif region in 'NYC'.split():
            buff = 0.01 # buffer around map
            xy   = 'n0.475/0.05'                                 # loc of colorbar
            self.cbar_pos  = f'{xy}+w3.0c/0.20c+h-Np'            # for multipanel

            offsets   = [['-0.1/-0', '0.025/0.025', '0.025/0.025'], # transect label st
                        ['0.15/-0.1', '0.2/0.05', '0.025/0.025']]   # transect label en

            # self.clims = ([-5, 5, 2], [0, 2, 0.5], [0, 200, 50])
            self.clims = ([-4, 4, 1], [0, 1.5, 0.25], [0, 200, 50])
            self.ti_pos = [(-79.85, 32.66), (-79.9, 32.89)] # Vertical rate mm/yr
            self.just_ref  = 'TL' # justification of GPS reference labels
            self.just_non  = 'BR' # justification of GPS non-reference labels
            parms  = {'equalsize': True}
            parmsU = {'equalsize': True}
            parmsR = {'equalsize': True}



        #### SHARED
        self.size, self.hwid = 6, 0.5
        self.WESN[[0, 2]]   -= buff
        self.WESN[[1, 3]]   += buff
        self.fc    = 'white' # fontcolor
        self.cmap  = 'vik'
        self.cmapU = 'lajolla'
        self.cmapR = 'imola'
        self.rev   = False # reverse colormap
        self.revU  = False # reverse colormap
        self.revR  = False # reverse colormap
        self.lcol  = 'gray30'
        self.wcol  = 'black'
        self.proj  = f'M{self.size}c'
        self.frame  = ['ws', 'y0.15','x0.2','g0']
        self.frameL = self.frame.copy(); self.frameL[0] = 'SnWe'
        self.frameR = self.frame.copy(); self.frameR[0] = 'Snwe'
        self.framec = 'black' # color of divisions between colorbar
        self.parms  = parms
        self.parmsU = parmsU
        self.parmsR = parmsR

        #  gps markers
        self.gms       = [0.30] # gps marker size
        self.mbor      = '0.50p,white' # marker border
        self.tis    = ['Vertical Rate (mm/yr)', 'Uncertainty (mm/yr)', 'Residual (mm)']

        self.ill    = 0.35 # illumination; brighter for presentation

        return


    def plot_basic(self, da, kind='rate', w_mask=None, continuous=False, df_lalo=None):
        """ Basic rate/uncertainty plot.

        Optional w_mask for background. Optional continuous colorbar. Optional custom markers
        """
        da = da.astype(float)
        if da.units == 'm/y':
            da = da.copy() * 1000

        fig  = pygmt.Figure()
        ## shrink for jupyter, lalo ticks (an colorbar), lalo labels, black frame
        pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5, MAP_TICK_LENGTH=0,
                    MAP_FRAME_TYPE='plain', FONT_ANNOT_PRIMARY='6p')

        fig.coast(region=self.WESN, frame=self.frameL, projection=self.proj,
                            resolution='f', land=self.lcol, water=self.wcol)

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameL,
                projection=self.proj, cmap=f'{self.wcol},{self.wcol}', nan_transparent=True)

        ## this was for when water is wrong?
        # wcol = self.wcol if w_mask is not None else None # so coast doesnt overwrite
        # wcol = self.wcol

        ## colormap and rate
        if kind.lower() == 'rate':
            idx = 0
            cmap = self.cmap
        elif kind.lower() == 'unc':
            idx  = 1
            cmap = self.cmapU

        pygmt.makecpt(cmap=cmap, series=self.clims[idx], reverse=self.rev,
                                        continuous=continuous)
        fig.grdimage(grid=da, cmap=True, shading=self.ill, nan_transparent=True)

        parms = dict(frame='af') if continuous else dict(equalsize=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**parms, position=self.cbar_pos+'+ml', N='p', shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[idx],
                    font=f'5p,{self.fc}', justify='BL')

        if df_lalo is not None:
            ## add custom lalo markers
            fig.plot(x=df_lalo.lon[:1], y=df_lalo.lat[:1],
                style=f'c{self.gms[0]}c', pen='1.0p,black') # color=non_sta.u_sig

            # in case only one point
            try:
                fig.plot(x=df_lalo.lon[1:], y=df_lalo.lat[1:], size=self.gms*(df_lalo.shape[0]-1),
                    style='tc', pen='1.0p,black') # fill=non_sta.u_sig
            except:
                pass

        return fig


    def plot_together(self, da_rate, da_unc, w_mask=None, continuous=False):
        if da_rate.units == 'm/y':
            da_rate *= 1000

        if da_unc.units == 'm/y':
            da_unc *= 1000


        fig  = pygmt.Figure()
        ## shrink for jupyter, lalo ticks (an colorbar), lalo labels, black frame
        pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5, MAP_TICK_LENGTH=0,
                    MAP_FRAME_TYPE='plain', FONT_ANNOT_PRIMARY='6p')

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameL,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')

        wcol = self.wcol if w_mask is not None else None # so coast doesnt overwrite
        wcol = self.wcol

        fig.coast(region=self.WESN, frame=self.frameL, projection=self.proj,
                            resolution='f', land=self.lcol, water=wcol)

        ## colormap and rate
        pygmt.makecpt(cmap=self.cmap, series=self.clims[0], reverse=self.rev, continuous=continuous)

        fig.grdimage(grid=da_rate, cmap=True, shading=self.ill, nan_transparent=True)

        parms = dict(frame='af') if continuous else dict(equa_size=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**parms, position=self.cbar_pos+'+ml', N='p', shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[0],
                    font=f'5p,{self.fc}', justify='BL')


        ## now the uncertainty
        fig.shift_origin(xshift=self.size+self.hwid, yshift=f'a') # shift right

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameR,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')

        fig.coast(region=self.WESN, frame=self.frameR, projection=self.proj,
                            resolution='f', land=self.lcol, water=wcol)

        ## colormap and rate
        pygmt.makecpt(cmap=self.cmapU, series=self.clims[1], reverse=self.revU,
                                        continuous=continuous)

        fig.grdimage(grid=da_unc, cmap=True, shading=self.ill, nan_transparent=True)

        parms = dict(frame='af') if continuous else dict(equalsize=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**parms, position=self.cbar_pos+'+ml', N='p', shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[1], font=f'5p,{self.fc}', justify='BL')

        return fig


    def plot_together_gps(self, da_rate, da_unc, df_gps, w_mask=None, continuous=False):
        """ Mask should be JUST water and should have 1 at water, nan elsewhere """
        if da_rate.units == 'm/y':
            da_rate *= 1000

        if da_unc.units == 'm/y':
            da_unc *= 1000

        da_rate_reg = da_rate#.sel(lat=slice(*self.SNWE[:2]), lon=slice(*self.SNWE[2:]))
        da_unc_reg  = da_unc#.sel(lat=slice(*self.SNWE[:2]), lon=slice(*self.SNWE[2:]))
        print (f'Mean Rate: {da_rate_reg.mean():.2f} mm/yr')
        print (f'Mean Unc: {da_unc_reg.mean():.2f} mm/yr')

        if continuous:
            # cbar_parms = {'frame': 'g1f0.5'}
            cbar_parms = {'frame': 'fa1'}
            ## dont get the spacing
            # clims = self.clims[0][:2], self.clims[1][:2]
            clims = self.clims
            # clims[1][-1] = 0.1 # change interval
            # self.ill = None
        else:
            cbar_parms = {'equalsize': True}
            clims = self.clims

        fig  = pygmt.Figure()
        ## shrink for jupyter, lalo ticks (an colorbar), lalo labels, black frame
        pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5, MAP_TICK_LENGTH=0,
                    MAP_FRAME_TYPE='plain', FONT_ANNOT_PRIMARY='6p')

        fig.coast(region=self.WESN, frame=self.frameL, projection=self.proj,
                            resolution='f', water=self.wcol, land=self.lcol)

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameL,
                    projection=self.proj, cmap=f'{self.wcol},{self.wcol}', nan_transparent=True)


        ## colormap and rate
        pygmt.makecpt(cmap=self.cmap, series=clims[0], reverse=self.rev, continuous=continuous)

        fig.grdimage(grid=da_rate, cmap=True, shading=self.ill, nan_transparent=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**cbar_parms, position=self.cbar_pos, shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[0],
                    font=f'5p,{self.fc}', justify='BL')

        ## add the GPS
        ref_sta   = df_gps[df_gps.sta.isin(self.ref_stas)]
        non_sta   = df_gps[~df_gps.sta.isin(self.ref_stas + ['WLP2'])]

        ## convert m to mm
        units = 1000 if np.abs(ref_sta.u_vel.iloc[0]) < 1e-2 else 1


        # non reference stas
        fig.plot(x=non_sta.lon, y=non_sta.lat, size=self.gms*non_sta.shape[0],
                fill=non_sta.u_vel*units, style='cc', intensity=self.ill, cmap=True, pen=self.mbor)
        fig.text(x=non_sta.lon, y=non_sta.lat, text=non_sta.sta,
                         justify=self.just_non, font='5p,gray99=~1.0p,gray1')

        # reference stas
        fig.plot(x=ref_sta.lon, y=ref_sta.lat, size=self.gms*ref_sta.shape[0],
                fill=ref_sta.u_vel*units, style='sc', intensity=self.ill, cmap=True, pen=self.mbor)
        fig.text(x=ref_sta.lon, y=ref_sta.lat, text=ref_sta.sta,
                         justify=self.just_ref, font='5p,gray99=~1.0p,gray1')

        ## add different for LOYZ
        # fig.plot(x=[sta_LOYZ.lon], y=[sta_LOYZ.lat], size=gms,
        #         fill=[sta_LOYZ.u_vel], style='sc', intensity=ill, cmap=True, pen=mbor)
        # fig.text(x=sta_LOYZ.lon, y=sta_LOYZ.lat, text=sta_LOYZ.sta,
        #                  justify=just_ref, font='5p,gray99=~1.0p,gray1')


        ## now the uncertainty
        fig.shift_origin(xshift=self.size+self.hwid, yshift=f'a') # shift right

        fig.coast(region=self.WESN, frame=self.frameR, projection=self.proj,
                            resolution='f', land=self.lcol, water=self.wcol)

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameR,
                    projection=self.proj, cmap=f'{self.wcol},{self.wcol}', nan_transparent=True)

        pygmt.makecpt(cmap=self.cmapU, series=clims[1], reverse=self.revU, continuous=continuous)
        fig.grdimage(grid=da_unc, cmap=True, shading=self.ill, nan_transparent=True)


        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(position=self.cbar_pos, shading=self.ill, **cbar_parms)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[1], font=f'5p,{self.fc}', justify='BL')

        fig.plot(x=non_sta.lon, y=non_sta.lat, size=self.gms*non_sta.shape[0],
                fill=non_sta.u_sig*units, style='cc', intensity=self.ill, cmap=True, pen=self.mbor)
        # fig.text(x=non_sta.lon, y=non_sta.lat, text=non_sta.sta,
        #                  justify=self.just_non, font='5p,gray99=~1.0p,gray1')


        fig.plot(x=ref_sta.lon, y=ref_sta.lat, size=self.gms*ref_sta.shape[0],
                fill=ref_sta.u_sig*units, style='sc', intensity=self.ill, cmap=True, pen=self.mbor)
        # fig.text(x=ref_sta.lon, y=ref_sta.lat, text=ref_sta.sta,
        #                  justify=self.just_ref, font='5p,gray99=~1.0p,gray1')

        return fig


    def plot_together_residue(self, da_rate, da_res, w_mask=None, continuous=False):
        fig  = pygmt.Figure()
        ## shrink for jupyter, lalo ticks (an colorbar), lalo labels, black frame
        pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5, MAP_TICK_LENGTH=0,
                    MAP_FRAME_TYPE='plain', FONT_ANNOT_PRIMARY='6p')

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameL,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')

        ## not sure if this is necessary to so coast doesnt overwrite
        # wcol = self.wcol if w_mask is not None else None
        wcol = self.wcol

        fig.coast(region=self.WESN, frame=self.frameL, projection=self.proj,
                            resolution='f', land=self.lcol, water=wcol)

        ## colormap and rate
        pygmt.makecpt(cmap=self.cmap, series=self.clims[0], reverse=self.rev, continuous=continuous)

        fig.grdimage(grid=da_rate, cmap=True, shading=self.ill, nan_transparent=True)
        if continuous:
            frame = 'gx'

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**self.parms, position=self.cbar_pos, shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[0],
                    font=f'5p,{self.fc}', justify='BL')



        ## add different for LOYZ
        # fig.plot(x=[sta_LOYZ.lon], y=[sta_LOYZ.lat], size=gms,
        #         fill=[sta_LOYZ.u_vel], style='sc', intensity=ill, cmap=True, pen=mbor)
        # fig.text(x=sta_LOYZ.lon, y=sta_LOYZ.lat, text=sta_LOYZ.sta,
        #                  justify=just_ref, font='5p,gray99=~1.0p,gray1')


        ## now the uncertainty
        fig.shift_origin(xshift=self.size+self.hwid, yshift=f'a') # shift right

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameR,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')

        fig.coast(region=self.WESN, frame=self.frameR, projection=self.proj,
                            resolution='f', land=self.lcol, water=wcol)

        pygmt.makecpt(cmap=self.cmapR, series=self.clims[2], reverse=self.revR,
                                        continuous=continuous)

        fig.grdimage(grid=da_res, cmap=True, shading=self.ill, nan_transparent=True)


        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**self.parmsR, position=self.cbar_pos, shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[2], font=f'5p,{self.fc}', justify='BL')


        return fig


    def plot_together_custom(self, da_rate, da_unc, df_lalo, df_gps=None, w_mask=None, continuous=False):
        """ Add custom lalo """
        fig  = pygmt.Figure()
        ## shrink for jupyter, lalo ticks (an colorbar), lalo labels, black frame
        pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5, MAP_TICK_LENGTH=0,
                    MAP_FRAME_TYPE='plain', FONT_ANNOT_PRIMARY='6p')

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameL,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')

        wcol = self.wcol if w_mask is not None else None # so coast doesnt overwrite
        wcol = self.wcol

        fig.coast(region=self.WESN, frame=self.frameL, projection=self.proj,
                            resolution='f', land=self.lcol, water=wcol)

        ## colormap and rate
        pygmt.makecpt(cmap=self.cmap, series=self.clims[0], reverse=self.rev, continuous=continuous)

        fig.grdimage(grid=da_rate, cmap=True, shading=self.ill, nan_transparent=True)

        parms = dict(frame='af') if continuous else dict(equalsize=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**parms, position=self.cbar_pos+'+ml', N='p', shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[0],
                    font=f'5p,{self.fc}', justify='BL')

        if df_gps is not None:
            ## add the GPS
            ref_sta   = df_gps[df_gps.sta.isin(self.ref_stas)]
            non_sta   = df_gps[~df_gps.sta.isin(self.ref_stas + ['WLP2'])]

            fig.plot(x=ref_sta.lon, y=ref_sta.lat, size=self.gms*ref_sta.shape[0],
                    fill=ref_sta.u_vel*1000, style='sc', intensity=self.ill, cmap=True, pen=self.mbor)
            fig.text(x=ref_sta.lon, y=ref_sta.lat, text=ref_sta.sta,
                             justify=self.just_ref, font='5p,gray99=~1.0p,gray1')

            fig.plot(x=non_sta.lon, y=non_sta.lat, size=self.gms*non_sta.shape[0],
                    fill=non_sta.u_vel*1000, style='cc', intensity=self.ill, cmap=True, pen=self.mbor)
            fig.text(x=non_sta.lon, y=non_sta.lat, text=non_sta.sta,
                             justify=self.just_non, font='5p,gray99=~1.0p,gray1')


            ## add different for LOYZ
            # fig.plot(x=[sta_LOYZ.lon], y=[sta_LOYZ.lat], size=gms,
            #         fill=[sta_LOYZ.u_vel], style='sc', intensity=ill, cmap=True, pen=mbor)
            # fig.text(x=sta_LOYZ.lon, y=sta_LOYZ.lat, text=sta_LOYZ.sta,
            #                  justify=just_ref, font='5p,gray99=~1.0p,gray1')


        ## now the uncertainty
        fig.shift_origin(xshift=self.size+self.hwid, yshift=f'a') # shift right

        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameR,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')

        fig.coast(region=self.WESN, frame=self.frameR, projection=self.proj,
                            resolution='f', land=self.lcol, water=wcol)

        ## colormap and unc
        pygmt.makecpt(cmap=self.cmapU, series=self.clims[1], reverse=self.revU,
                                        continuous=continuous)

        fig.grdimage(grid=da_unc, cmap=True, shading=self.ill, nan_transparent=True)

        parms = dict(frame='af') if continuous else dict(equalsize=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**parms, position=self.cbar_pos+'+ml', N='p', shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[1], font=f'5p,{self.fc}', justify='BL')

        if df_gps is not None:
            fig.plot(x=ref_sta.lon, y=ref_sta.lat, size=self.gms*ref_sta.shape[0],
                    fill=ref_sta.u_sig*1000, style='sc', intensity=self.ill, cmap=True, pen=self.mbor)
            # fig.text(x=ref_sta.lon, y=ref_sta.lat, text=ref_sta.sta,
            #                  justify=self.just_ref, font='5p,gray99=~1.0p,gray1')

            fig.plot(x=non_sta.lon, y=non_sta.lat, size=self.gms*non_sta.shape[0],
                    fill=non_sta.u_sig*1000, style='cc', intensity=self.ill, cmap=True, pen=self.mbor)
            # fig.text(x=non_sta.lon, y=non_sta.lat, text=non_sta.sta,
            #                  justify=self.just_non, font='5p,gray99=~1.0p,gray1')

        ## add custom lalo markers
        fig.plot(x=df_lalo.lon, y=df_lalo.lat, size=self.gms*df_lalo.shape[0],
            style='tc', pen='1p,white') # fill=non_sta.u_sig
        # fig.savefig(dst, dpi=500, anti_alias=True, show=False)
        # print ('Wrote:', dst)
        return fig


    def plot_transects(self, da_rate, df_gps, w_mask=None, continuous=False):
        import string
        offsets   = [['-0.1/-0', '0.025/0.025', '0.025/0.025'],
                     ['0.15/-0.1', '0.2/0.05', '0.025/0.025']]


        fig  = pygmt.Figure()
        # pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5) # shrink for jupyter
        pygmt.config(PS_SCALE_X=0.5, PS_SCALE_Y=0.5, MAP_TICK_LENGTH=0,
                    MAP_FRAME_TYPE='plain', FONT_ANNOT_PRIMARY='6p')
        if w_mask is not None:
            fig.grdimage(grid=w_mask, region=self.WESN, frame=self.frameL,
                    projection=self.proj, cmap=f'{self.wcol},{self.lcol}')
            self.wcol=None # so that next coast doesnt overwrite

        wcol = self.wcol if w_mask is not None else None # so coast doesnt overwrite
        wcol = self.wcol

        fig.coast(region=self.WESN, frame=self.frameL, projection=self.proj,
                            resolution='f', land=self.lcol, water=self.wcol)

        pygmt.makecpt(cmap=self.cmap, series=self.clims[0], reverse=self.rev, continuous=continuous)

        fig.grdimage(grid=da_rate, cmap=True, shading=self.ill, nan_transparent=True)

        parms = dict(frame='af') if continuous else dict(equalsize=True)

        ## make sure the colorbar ticks are white
        with pygmt.config(FONT_ANNOT_PRIMARY=f'4p,{self.fc}', FONT_ANNOT_SECONDARY=f'4p,{self.fc}'):
            fig.colorbar(**parms, position=self.cbar_pos+'+ml', N='p', shading=self.ill)
            fig.text(x=self.ti_pos[0][0], y=self.ti_pos[0][1],
                    text=self.tis[0],
                    font=f'5p,{self.fc}', justify='BL')


        # ## Transects
        lss          = '..', '-', '' # pen styles for the line
        sty          = 'f0.75/1p'   ## plot as fault; prob unnecessary
        df_transects = make_transects(self.path_transects)

        for i, (id, df_t) in enumerate(df_transects.groupby('id')):
            path_arrows = op.join(self.path_transects, f'arrows{i}.xy')
            fig.plot(x=df_t.x, y=df_t.y, style=sty, pen=f'1.8p,white,{lss[i]}')
            fig.text(path_arrows, F='+f6.5p,white=black+a') # arrows

            # annotate start and end of track; outline is after font
            fig.text(x=df_t.x.iloc[0], y=df_t.y.iloc[0],
                        font='5p,Helvetica-Bold,white=~0.6,black',
                        text=string.ascii_uppercase[i],
                        offset=offsets[0][i])

            fig.text(x=df_t.x.iloc[-1], y=df_t.y.iloc[-1],
                        font='5p,Helvetica-Bold,white=~0.6,black',
                        text=f"{string.ascii_uppercase[i]}'",
                        offset=offsets[1][i])

        ## GPS velocity
        ref_sta   = df_gps[df_gps.sta.isin(self.ref_stas)]
        non_sta   = df_gps[~df_gps.sta.isin(self.ref_stas + ['WLP2'])]

        fig.plot(x=ref_sta.lon, y=ref_sta.lat, size=self.gms*ref_sta.shape[0],
                fill=ref_sta.u_vel*1000, style='sc', intensity=self.ill, cmap=True, pen=self.mbor)
        fig.text(x=ref_sta.lon, y=ref_sta.lat, text=ref_sta.sta,
                         justify=self.just_ref, font='5p,gray99=~1.0p,gray1')

        fig.plot(x=non_sta.lon, y=non_sta.lat, size=self.gms*non_sta.shape[0],
                fill=non_sta.u_vel*1000, style='cc', intensity=self.ill, cmap=True, pen=self.mbor)
        fig.text(x=non_sta.lon, y=non_sta.lat, text=non_sta.sta,
                         justify=self.just_non, font='5p,gray99=~1.0p,gray1')

        # dst = op.join(path_vup, 'Figures', f'HR-Transects-Map.png')
        # fig.savefig(dst, dpi=500, anti_alias=True, show=False)
        # print ('Wrote:', dst)
        return fig


def make_transects(path_transects):
    print (os.sys.path)
    from VLM.EC.Profiles.SE_Transects_Calcs import Midpoint, make_df_transects

    # create the arrows
    for i in range(3):
        Midpoint('HR', i, path_transects)()

    # organize the lines
    df_transects = make_df_transects('Charleston')
    return df_transects


if __name__ == '__main__':

    # Exp    = ExpBase(Charleston_SR, 'Base_Fast')
    Exp    = ExpBase(NYC_SR, 'Base', 'NJHT', neofs=15)
    rateM  = xr.open_dataset(Exp.path_rate_msk_nc)['Band1']*1000
    uncM   = xr.open_dataset(Exp.path_std_msk_nc)['Band1']*1000

    df_gps = prep_gps(Exp.path_gps, Exp.reg)

    Obj = PlotVLM(Exp.reg, Exp.path_mp_exp)
    fig = Obj.plot_together_gps(rateM, uncM, df_gps, continuous=False)
    ## its too small to see locally
    # fig.show(dpi=600, width=1000)
