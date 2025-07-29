% function[] =  profile(path_vup)
% addpath(fullfile(getenv('HOME'), 'Google_Drive', 'BB_LIB', 'VLM', 'utils'));
reg='HR';
% exp='0.0-0.7-0.2-var-crop-GACOS-unw-derampL';
disp(['Plotting transects in ', reg])
disp(['Experiment: ', 'FRInGE'])

ref = 'LOY2 LS03 SPVA VAHP';
xax_lims = [0, 100]; % Km;
yax_lims = [-14, 14]; % mm/yr
radius     = 500; % m

% nodata     = 9.96920996838686905e+36;

root = strcat(getenv('dataroot'), '/VLM/Sentinel1/HR/MintPy_2alks_5rlks_33_15/LOY2_Base/geo')
% root = strcat(getenv('dataroot'), '/VLM/Sentinel1/track_004/Results/IGS14_GRL2020/Vup_MintPy/Vup')
% root = '/Users/buzzanga/Desktop';

Points{1}  =[
  -75.9410   36.7870
  -76.2587   36.7603
  -76.3000   36.8040
  -76.3485   36.8816
  -76.3600   36.9044
  -76.4043   37.0586
  -76.5546   37.1605
  -76.5700   37.17];

Points{2} = [
  -75.9896   36.9374
  -76.2393   36.7628
  -76.4091   36.6633];

Points{3}=[
  -76.401    36.84
  -76.4043   37.0586
  -76.3897   37.0659
  -76.356    37.085
  -76.45     37.155];

% create these with save_gdal (Jupyter Notebook)
network    = 'SR'
fnames     = [string(fullfile(root, strcat('geo_rate_',network,'_masked.nc')))...
            string(fullfile(root,  strcat('geo_std_',network,'_masked.nc')))];

DATA        = ncread(fnames(1), 'Band1')* 1000;
DATA_STD    = ncread(fnames(2), 'Band1')* 1000;
lat_vector  = ncread(fnames(1), 'lat');
lon_vector  = ncread(fnames(1), 'lon');
[LON,LAT] = meshgrid(lon_vector, lat_vector);

%% flip these because underlying data is flipped
LON=LON';
LAT=LAT';

figure;
  h = imagesc(lon_vector,lat_vector,DATA');
  set(h,'alphadata',~isnan(DATA'))
  colormap('turbo')
  caxis([-10 10])
  axis equal
  axis tight
  colorbar
  axis xy
  title('all rates')
  set(gca,'Color',[0.8 0.8 0.8]); % gray background on ratemap

  % plot the profile points
  hold on
  markers=['s', 'd', 'o'];
  for i=1:3
    t1 = cell2mat(Points(i));
    scatter(t1(:, 1), t1(:, 2), 'filled', markers(i))
  end


figure
  lon = reshape(LON,[],1);
  lat = reshape(LAT,[],1);
  data2 = reshape(DATA,[],1);
  lon(isnan(data2))=[];
  lat(isnan(data2))=[];
  data2(isnan(data2))=[];
  scatter(lon(1:100:end),lat(1:100:end),5,data2(1:100:end));
  view(0,90)
  axis equal

%% compute profiles
GPS_root ='~/data/VLM/GNSS/UNR';
GPS_tmp  = readtable(fullfile(GPS_root, strcat('Midas_GPS-HR_IGS14_crop.csv')));


% GPS_tmp.Properties.VariableNames = {'index', 'station_name', 'lat' ,'lon', 'tss', 'tse', 'v', 'vs'};
GPS_tmp.Properties.VariableNames{'u_vel'} = 'v';
GPS_tmp.Properties.VariableNames{'u_sig'} = 'vs';

%% get rid of these
GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'CHR2'), :);
GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'WLP2'), :);
% GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'LNG4'), :);

GPS_tmp.station_lonlat = [GPS_tmp.lon, GPS_tmp.lat];

%% these will be gold
GPS_non = [GPS_tmp(strcmp(GPS_tmp.sta, 'DRV5/6'), :); ...
           GPS_tmp(strcmp(GPS_tmp.sta,'LNG4'), :); ...
           GPS_tmp(strcmp(GPS_tmp.sta,'WLP2'), :); ...
           ];

%% these will be blue
GPS_ref = setdiff(GPS_tmp, GPS_non);

GPS_non
GPS_ref

% preping data for profiles
lonlat = [reshape(LON,[],1) reshape(LAT,[],1)];
data = reshape(DATA,[],1);
data_std = reshape(DATA_STD,[],1);
ix_drop = isnan(data);
% lonlat(ix_drop,:)=[];
% data_std(ix_drop)=[];
% data(ix_drop)=[];

fontsize=22;


clear figprop
figprop.fontsize = 22;
figprop.ps_min=3;
figprop.units ='mm/yr';
figprop.n_crossplot = 3;
figprop.ylimits = yax_lims;        % mm/yr
figprop.xlimits = xax_lims;        % Km
figprop.binsize = 500;             % m
% figprop.poly = poly;
% figprop.poly_color = 'k';
% figprop.poly_style = '--';
% figprop.poly_linewidth=1;
figprop.data{1} = [GPS_ref.v GPS_ref.vs];      % ref stations diff color
figprop.data{2} = [GPS_non.v GPS_non.vs];      % stations not used inref
figprop.data_lonlat{1} = [GPS_ref.station_lonlat];
figprop.data_lonlat{2} = [GPS_non.station_lonlat];
figprop.data_color{1} = rgb('blue'); % rgb('Red');  % ref
figprop.data_color{2} = rgb('Gold'); % rgb('Lime'); % nonref % deep pink
figprop.data_draw_manual  ='y';
figprop.data_draw_manual_width = 1.15;
% figprop.data_draw_manual_linestyle
figprop.data_err_linewidth = 2.5;

figprop.names{1} = GPS_ref.sta;
figprop.names{2} = GPS_non.sta;
% two sigma color; turn alpha to 0 (and color to white?) to turn off
figprop.sourceData_tsigma_alpha= 0.5;
figprop.sourceData_tsigma_color=[0.8 0.8 0.8];
data_profile_HR(data,data_std,lonlat,Points,200,radius,figprop);
%% temporarily take out of function for debugging
% data_or = data;
% data_or_std = data_std;
% data_or=data;data_std_or=data_std; lonlat_or=lonlat; points_list=Points;
% bins=200; R=radius; figprop=figprop;
% data_profile_HR;

% print(strcat('Transects_', reg),'-deps','-tiff')
dst = fullfile(root, 'Figures', [reg,'_', network, '_Transects']);
print(dst, '-dpng')
disp(['Wrote: ', dst, '.png'])

% copygraphics(gcf) % copy so you can paste fig into word/pages; slow
% close all
return
