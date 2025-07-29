%% Prep Transect, GNSS Data
close all
clear all
clc

% prep fringe
dct_expF = 'Charleston_Base';
dct_expA = 'Charleston_ARIA';
mp_exp   = 'Base';

%% ---------------------------------------------------------------  Prepare Data
% [dataA, data_stdA, lonlatA] = prep_transect('ARIA');
% [data, data_std, lonlat] = prep_transect('FRInGE');

% prep fringe
command = strcat(['python profile_meta.py ', dct_expF, ' ', mp_exp]);
system(command);
load profile_parms

DATA        = ncread(PATH_RATE, 'Band1') * 1000;
DATA_STD    = ncread(PATH_STD, 'Band1') * 1000;
lat_vector  = ncread(PATH_RATE, 'lat');
lon_vector  = ncread(PATH_RATE, 'lon');
[LON,LAT]   = meshgrid(lon_vector, lat_vector);
LON=LON'; LAT=LAT';

lonlat = [reshape(LON,[],1) reshape(LAT,[],1)];
data = reshape(DATA,[],1);
data_std = reshape(DATA_STD,[],1);
ix_drop = isnan(data);

% lonlat(ix_drop,:)=[];
% data_std(ix_drop)=[];
% data(ix_drop)=[];


command = strcat(['python profile_meta.py ', dct_expA, ' ', mp_exp]);
system(command);
load profile_parms

DATA        = ncread(PATH_RATE, 'Band1') * 1000;
DATA_STD    = ncread(PATH_STD, 'Band1') * 1000;
lat_vector  = ncread(PATH_RATE, 'lat');
lon_vector  = ncread(PATH_RATE, 'lon');
[LON,LAT]   = meshgrid(lon_vector, lat_vector);
LON=LON'; LAT=LAT';
lonlatA = [reshape(LON,[],1) reshape(LAT,[],1)];
dataA = reshape(DATA,[],1);
data_stdA = reshape(DATA_STD,[],1);
ix_drop = isnan(dataA);
% lonlat(ix_drop,:)=[];
% data_std(ix_drop)=[];
% dataA(ix_drop)=[];


disp(['Plotting transects in ', REG])
disp(['Experiment: ', mp_exp])


%% ----------------------------------------------------------  Prepare Transects
Points{1}  =[
  -75.9410   36.7870 % get rid of
  -76.2587   36.7603
  -76.3000   36.8040
  -76.3485   36.8816
  -76.3600   36.9044
  -76.4043   37.0586
  -76.5546   37.1605
  -76.5700   37.17];

Points{2} = [
  -75.9896   36.9374
  -76.1445   36.8501   % add the midpoint as it needs 3 for poly
  -76.2393   36.7628]; % remove outside of bounds
  % -76.4091   36.6633];

Points{3}=[
  -76.401    36.84
  -76.4043   37.0586
  -76.3897   37.0659
  -76.356    37.085
  -76.45     37.155];

Points = COORDS;

%% ----------------------------------------------------------------- Prepare GPS
GPS_tmp  = readtable(PATH_GPS);

% GPS_tmp.Properties.VariableNames = {'index', 'station_name', 'lat' ,'lon', 'tss', 'tse', 'v', 'vs'};
GPS_tmp.Properties.VariableNames{'u_vel'} = 'v';
GPS_tmp.Properties.VariableNames{'u_sig'} = 'vs';

GPS_tmp.station_lonlat = [GPS_tmp.lon, GPS_tmp.lat];

% these will be blue
GPS_ref = GPS_tmp(ismember(GPS_tmp.sta, REF_STAS), :);

% these will be gold
GPS_non = setdiff(GPS_tmp, GPS_ref);

% get rid of a few specific ones
GPS_non = GPS_non(~ismember(GPS_non.sta, EX_STAS), :);

% GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'WLP2'), :);
% GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'LNG4'), :);
%% ---------------------------------------------------------------- Prepare Figs

fontsize=22
radius = 750
nbins  = 15; %200

figprop.fontsize = 22;
figprop.ps_min=3;
figprop.units ='mm/yr';
figprop.n_crossplot =3;
figprop.ylimits =YLIMS;
figprop.xlimits =XLIMS;
figprop.binsize=RADIUS;               % in m
% figprop.poly = poly;
% figprop.poly_color = 'k';
% figprop.poly_style = '--';
% figprop.poly_linewidth=1;
figprop.data{1} = [GPS_ref.v GPS_ref.vs];            % stations for referencing
figprop.data{2} = [GPS_non.v GPS_non.vs];            % stations not used in refercing.
figprop.data_lonlat{1} = [GPS_ref.station_lonlat];
figprop.data_lonlat{2} = [GPS_non.station_lonlat];
figprop.data_color{1} = rgb('blue');
figprop.data_color{2} = rgb('Gold'); % deep pink
figprop.ext_col       = rgb('Violet');
figprop.ext_lw        = 2.25; % outline weight

figprop.data_draw_manual  ='y';
figprop.data_draw_manual_width = 1.15;
% figprop.data_draw_manual_linestyle
figprop.data_err_linewidth = 2.5;

figprop.names{1} = GPS_ref.sta;
figprop.names{2} = GPS_non.sta;
figprop.sourceData_marker_size = 3;
% figprop.sourceData_sigma_alpha = 0.1;

% only show 1 sigma
figprop.sourceData_tsigma_alpha= 0;
figprop.sourceData_tsigma_color='white';

clc
save=1; % save the data to be overlaid

data_profile_overlays(dataA,data_stdA,lonlatA,Points,10,radius,figprop,save);
% % close all


save=0;
data_profile_overlays(data,data_std,lonlat,Points,10,radius,figprop, save);

dst=fullfile(PATH_FIGS, 'Transect_Overlay.png');
% print(dst,'-depsc','-tiff')
saveas(gcf, dst)
disp (['Wrote: ', dst]); 
close all;
return

% %% plot a polygon
% load poly2;
% figure
% for i=1:size(polycoords,1)
%   x = polycoords(i, :);
%   idx = find(x,1,'first'):find(x,1,'last');
% 
%   X = x(idx);
%   Y = polyvals(i, idx);
% 
%   fill(X, Y,'red','EdgeColor','green','FaceAlpha',0.2);
%   hold on
% end