% function[] =  profile(dct_exp mp_exp)
clear all; close all; clc;

currentFolder = pwd;
% run this from the bzFRInGE/profiles directory
[~, folderName] = fileparts(currentFolder);

if strcmp(folderName, 'profiles')
else
    disp('You are not in the ''profiles'' directory.');
end

dct_exp = 'Houston_SR'; % DC_SR NYC_SRc Charleston_SR
mp_exp  = 'ERA5_SET_PM_ex_Fast'; % ERA5_SET_PM_ex_Fast_2017 Base_Fast
ref_sta = 'NASA'; % 'USN8 NYBK NASA'
% transects = ['Shipping_Channel_South']; % Houston
transects = ['Shipping_Channel_South', ',NE_Transect']; % Houston

% transects = ['West_Subsidence', ',East_Uplift']; % DC
% transects = ['Transect_NJ_NYOB', ',Transect_Manhattan',
% ',Transect_Brooklyn']; % NYC

neofs   = 5;

command = strcat(['python profile_meta.py ', transects, ' ', dct_exp, ' ', mp_exp, ' ', ref_sta, ' ', num2str(neofs)]);

system(command);
load profile_parms

% addpath(fullfile(getenv('HOME'), 'Google_Drive', 'BB_LIB', 'VLM', 'utils'));
disp(['Plotting transects in ', REG])
disp(['Experiment: ', mp_exp])

DATA        = ncread(PATH_RATE, 'Band1') * 1000;
DATA_STD    = ncread(PATH_STD, 'Band1') * 1000;
lat_vector  = ncread(PATH_RATE, 'lat');
lon_vector  = ncread(PATH_RATE, 'lon');
[LON,LAT] = meshgrid(lon_vector, lat_vector);
LON=LON';
LAT=LAT';
if ~iscell(COORDS)
    COORDS = {squeeze(COORDS)};
end
%% Diagnostic Plots
% flip these because underlying data is flipped
figure;
  h = imagesc(lon_vector,lat_vector,DATA');
  set(h,'alphadata',~isnan(DATA'))
  colormap('jet')
  
  caxis([-3 3])
  axis equal
  axis tight
  colorbar
  axis xy
  title('all rates')
  set(gca,'Color',[0.8 0.8 0.8]); % gray background on ratemap

  % plot the profile points; different marker for each transect
  hold on
  if size(COORDS, 1) == 3
    marker=['s', 'd', 'o'];
  elseif size(COORDS, 1) == 2
    markers=['s', 'd'];
  elseif size(COORDS, 1) == 1
    markers = ['s'];
  end 
  cmap = flag(size(COORDS, 1));
  for i=1:size(COORDS, 1)
    t1 = cell2mat(COORDS(i));
    scatter(t1(:, 1), t1(:, 2), 'filled', markers(i), 'MarkerFaceColor', cmap(i,:))
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
  

  %% Prep GPS
GPS_tmp  = readtable(PATH_GPS);

% GPS_tmp.Properties.VariableNames = {'index', 'station_name', 'lat' ,'lon', 'tss', 'tse', 'v', 'vs'};
GPS_tmp.Properties.VariableNames{'u_vel'} = 'v';
GPS_tmp.Properties.VariableNames{'u_sig'} = 'vs';
GPS_tmp.v  = GPS_tmp.v*1000;
GPS_tmp.vs = GPS_tmp.vs*1000;

GPS_tmp.station_lonlat = [GPS_tmp.lon, GPS_tmp.lat];

% these will be blue
GPS_ref = GPS_tmp(ismember(GPS_tmp.sta, REF_STAS), :);

% these will be red
% GPS_non = setdiff(GPS_tmp.sta, GPS_ref.sta, 'rows'); % doesnt workanymore
GPS_non = setdiff(GPS_tmp.sta, GPS_ref.sta);
% get rid of a few specific ones
GPS_non = GPS_tmp(ismember(GPS_tmp.sta, GPS_non), :);
GPS_non = GPS_non(~ismember(GPS_non.sta, EX_STAS), :);

% GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'W
% LP2'), :);
% GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'LNG4'), :);
%% compute profiles

% preping data for profiles
lonlat = [reshape(LON,[],1) reshape(LAT,[],1)];
data = reshape(DATA,[],1);
data_std = reshape(DATA_STD,[],1);
ix_drop = isnan(data);
% lonlat(ix_drop,:)=[];
% data_std(ix_drop)=[];
% data(ix_drop)=[];

clear figprop
figprop.fontsize = 16;
figprop.ps_min=3;
figprop.units ='mm/yr';
figprop.n_crossplot = 3;
figprop.ylimits = YLIMS;        % mm/yr
figprop.xlimits = XLIMS;        % Km
figprop.binsize = RADIUS;       % m
% figprop.poly = poly;
% figprop.poly_color = 'k';
% figprop.poly_style = '--';
% figprop.poly_linewidt
figprop.data{1} = [GPS_ref.v GPS_ref.vs];      % ref stations diff color
figprop.data{2} = [GPS_non.v GPS_non.vs];      % stations not used inref
figprop.data_lonlat{1} = [GPS_ref.station_lonlat];
figprop.data_lonlat{2} = [GPS_non.station_lonlat];
figprop.data_color{1} = rgb('blue'); % rgb('Red');  % ref
figprop.data_color{2} = rgb('mediumblue'); % rgb('Lime'); % nonref % deep pink
figprop.data_draw_manual  ='y';
figprop.data_draw_manual_width = 1.15;
% figprop.data_draw_manual_linestyle
figprop.data_err_linewidth = 2.5;

figprop.names{1} = GPS_ref.sta;
figprop.name{2} = GPS_non.sta;
% two sigma color; turn alpha to 0 (and color to white?) to turn off
figprop.sourceData_tsigma_alpha= 0.0;
% figprop.sourceData_tsigma_color=[0.8 0.8 0.8];
figprop.sourceData_tsigma_color=[1.0 1.0 1.0];

% shortens graph
data_profile_Houston(data,data_std,lonlat,COORDS,15,RADIUS,figprop);

% data_profile_HR(data,data_std,lonlat,COORDS,200,RADIUS,figprop);
% data_profile_NYC(data,data_std,lonlat,COORDS,15,RADIUS,figprop);
% data_profile_SE(data,data_std,lonlat,COORDS,15,RADIUS,figprop);
% data_profile_DC(data,data_std,lonlat,COORDS,15,RADIUS,figprop);

%% temporarily take out of function for debugging
% data_or = data;
% data_or_std = data_std;
% data_or=data;data_std_or=data_std; lonlat_or=lonlat; points_list=Points;
% bins=200; R=radius; figprop=figprop;
% data_profile_HR;
% it's messed up on save; do it manually

% print(strcat('Transects_', reg),'-deps','-tiff')
% dst = fullfile(PATH_FIGS, [REG,'_Transects']);
% print(dst, ['-r' num2str(500)], '-dpng')
% disp(['Wrote: ', dst, '.png'])

% copygraphics(gcf) % copy so you can paste fig into word/pages; slow
% close all
% return
