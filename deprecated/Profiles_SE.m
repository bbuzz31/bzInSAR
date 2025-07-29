% function[] =  profile(reg, exp)
% reg = 'HR';
% exp = 'FRInGE';
addpath(fullfile(getenv('HOME'), 'Google_Drive', 'BB_LIB', 'VLM', 'utils'));

% reg='Charleston';
% exp='0.5-0.7-0.2-var-crop-GACOS';

reg='HR';
exp='0.0-0.45-0.0-var-GACOS';


disp(['Plotting transects in ', reg])
disp(['Experiment: ', exp])

if strcmp(reg, 'Charleston')
  ref = {'SCHA'};
  xax_lims = [0, 30]; % Km;
  yax_lims = [-10, 5]; % mm/yr
  radius     = 250; % m
  root       = sprintf('%s/VLM/Sentinel1/track_150/Results/%s/Vup', getenv('dataroot'), reg);
  fnames     = [string(fullfile(root, sprintf('%s-rate_%s_130.nc', reg, exp))), ...
            string(fullfile(root,  sprintf('%s-unc_%s_130.nc', reg, exp)))];
    Points{1} = [-80.0473513336293 32.8586780035838
        -79.9651212479214 32.8584202290205
        -79.9380549187699 32.7831500565229
        -79.9249084160391 32.7795412126361
        -79.8421627812045 32.7573726001882];

    Points{2} = [-79.8992706969624 32.892916893999995
        -80.0052770268589 32.7957671017285
        -79.9712749965147 32.7329633515634
        -79.89991 32.75112];

    Points{3} = [-79.7836637937923 32.892916893999995
                -79.9992766685629 32.7401637815186];

elseif strcmp(reg, 'Savannah')
  ref = {'SAVA'};
  xax_lims = [0, 18]; % Km;
  yax_lims = [-10, 5]; % mm/yr
  radius     = 250; % m
  root       = sprintf('%s/VLM/Sentinel1/track_150/Results/%s/Vup', getenv('dataroot'), reg);
  fnames     = [string(fullfile(root, sprintf('%s-rate_%s_130.nc', reg, exp))), ...
              string(fullfile(root,  sprintf('%s-unc_%s_130.nc', reg, exp)))];
  Points{1} = [-81.1165005262072 32.0067339047308
               -81.0831712834927 32.095904247081];

  Points{2} = [-81.173803434734 32.024860334979
               -81.053935105673 32.0236908878662
               -81.0077419447178 32.0497110861258];

  Points{3} = [-81.1007129901846 32.0073186282872
               -81.1776041378505 32.1505758996039];


elseif strcmp(reg, 'HR')
  ref = {'LOY2', 'LOYZ', 'LS03', 'SPVA', 'VAHP'};
  xax_lims = [0, 100]; % Km;
  yax_lims = [-14, 14]; % mm/yr
  radius   = 500; % m
  % root     = sprintf('%s/VLM/Sentinel1/HR/MintPy_2alks_5rlks_33_15/Vup', getenv('dataroot'));
  % fnames   = [string(fullfile(root, 'geo_rate_masked.nc')), ...
  %           string(fullfile(root,  'geo_rate_masked.nc'))];
  root       = sprintf('%s/VLM/Sentinel1/track_004/Results/%s/Vup', getenv('dataroot'), reg);
  fnames     = [string(fullfile(root, sprintf('%s-rate_%s_125.nc', reg, exp))), ...
              string(fullfile(root,  sprintf('%s-unc_%s_125.nc', reg, exp)))];
  Points{1}  =[    -75.9410   36.7870
    -76.2587   36.7603
    -76.3000   36.8040
    -76.3485   36.8816
    -76.3600   36.9044
    -76.4043   37.0586
    -76.5546   37.1605
    -76.5700   37.17];

  Points{2} = [  -75.9896   36.9374
    -76.2393   36.7628
    -76.4091   36.6633];

  Points{3}=[  -76.575   36.851
    -76.5546   36.9689
    -76.3897   37.0659
    -76.356    37.085
    -76.45     37.155];
else
  disp(['No instructions for ', reg])
  disp('Exiting')
  return
end

nodata     = 0;

% create these with save_gdal (Jupyter Notebook)


DATA        = ncread(fnames(1), 'Band1') * 1000;
DATA_STD    = ncread(fnames(2), 'Band1') * 1000;
lat_vector  = ncread(fnames(1), 'lat');
lon_vector  = ncread(fnames(1), 'lon');

[LON,LAT] = meshgrid(lon_vector, lat_vector);

figure;
h = imagesc(lon_vector,lat_vector,DATA');
set(h,'alphadata',~isnan(DATA'))
colormap('jet')
caxis([-15 15])
axis equal
axis tight
colorbar
axis xy
title('all rates')
set(gca,'Color',[0.8 0.8 0.8]); % gray background on ratemap
% plot the points
hold on;
markers=['s', 'd', 'o'];
for i=1:3
  t1 = cell2mat(Points(i));
  scatter(t1(:, 1), t1(:, 2), 'filled', markers(i))
end


%keyboard
% figure
% lon = reshape(LON,[],1);
% lat = reshape(LAT,[],1);
% data2 = reshape(DATA,[],1);
% lon(isnan(data2))=[];
% lat(isnan(data2))=[];
% data2(isnan(data2))=[];
%
% scatter(lon(1:10:end),lat(1:10:end),5,data2(1:10:end));
% view(0,90)
% axis equal
%  close all

%% compute profiles
% GPS_root ='/Volumes/BB_4TB/data/VLM/GNSS/UNR';
GPS_root ='~/data/VLM/GNSS/UNR';
GPS_tmp  = readtable(fullfile(GPS_root, strcat('Midas_GPS-', reg, '_IGS14_crop.csv')));


% GPS_tmp.Properties.VariableNames = {'index', 'station_name', 'lat' ,'lon', 'tss', 'tse', 'v', 'vs'};
GPS_tmp.Properties.VariableNames{'u_vel'} = 'v';
GPS_tmp.Properties.VariableNames{'u_sig'} = 'vs';
% GPS_tmp = GPS_tmp(strcmp(GPS_tmp.sta, 'SCHA'), :);

% get the gps stations rows that match the references
idx = zeros(size(GPS_tmp, 1), 1);
for i=1:numel(GPS_tmp.sta)
  if ismember(string(GPS_tmp.sta(i)), ref)
    idx(i) = i;
  end
end
% get just the good rows
idx = idx(idx>0);

% remove stations?
% GPS_tmp = GPS_tmp(~strcmp(GPS_tmp.sta, 'CHR2'), :);
GPS_tmp.station_lonlat = [GPS_tmp.lon, GPS_tmp.lat];

% set the ref / non ref stations
GPS_ref = GPS_tmp(idx, :);
% GPS_ref = GPS_tmp(strcmp(GPS_tmp.sta, ref), :);

GPS_non = setdiff(GPS_tmp, GPS_ref);

% preping data for profiles
lonlat = [reshape(LON,[],1) reshape(LAT,[],1)];
data = reshape(DATA,[],1);
data_std = reshape(DATA_STD,[],1);
ix_drop = isnan(data);
lonlat(ix_drop,:)=[];
data_std(ix_drop)=[];
data(ix_drop)=[];

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
figprop.data_color{1} = rgb('blue'); % rgb('Red');
figprop.data_color{2} = rgb('Gold'); % rgb('Lime'); % deep pink
figprop.data_draw_manual  ='y';
figprop.data_draw_manual_width = 1.15;
% figprop.data_draw_manual_linestyle
figprop.data_err_linewidth = 2.5;

figprop.names{1} = GPS_ref.sta;
figprop.names{2} = GPS_non.sta;
% two sigma color; turn alpha to 0 (and color to white?) to turn off
figprop.sourceData_tsigma_alpha= 0.5;
figprop.sourceData_tsigma_color=[0.8 0.8 0.8];


if strcmp(reg, 'HR')
  data_profile_HR(data,data_std,lonlat,Points,200,radius,figprop);
else
  data_profile_SE(data,data_std,lonlat,Points,200,radius,figprop);
end


% print(strcat('Transects_', reg),'-deps','-tiff')
dst = fullfile(root, 'Figures', [reg, '-Transects_', exp]);
print(dst, '-dpng')
disp(['Wrote: ', dst, '.png'])

% copygraphics(gcf) % copy so you can paste fig into word/pages; slow
% close all
return
