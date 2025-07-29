function [] = data_profile(data_or,data_std_or,lonlat_or,points_list,bins,R,figprop, save_flag)

% By David Bekaert - December 2010
% give as input the data, the lonlat (2 colums), the profile points [lon lat], the number of
% bins and the radius [m] of the tube used to project points
% by default a 15 bins and a radius of 100m is used.
% for debug figures set debugflag to 1
%
% options for figprop include:
% SPECIFIC PROPERTIES FOR FIGURE AXES
% figprop.ps_min                    Minimum number of points in a bin to draw an thicker line and uncertainty bounds
% figprop.fontsize                  Fontsize of the labels
% figprop.units                     Units of the y-axis
% figprop.n_crossplot               Number of pannels (vertical) to plot in one figure
% figprop.ylimits                   Y-axis limits
% figprop.yticks                    Y-axis ticks
% figprop.xlimits                   X-axis limits
% figprop.xticks                    X-axis ticks
% figprop.binsize                   Binsize in m units
% PROPERTIES OF THE SOURCE DATA USED TO CALCULATE THE PROFILES
% figprop.sourceData_marker_color   Marker color in the profiles of the source data
% figprop.sourceData_marker_size    Marker size in the profiles of the source data
% figprop.sourceData_marker_alpha   Marker transparancy of the source data in the profiles
% figprop.sourceData_tsigma_color   Two sigma color of the source data uncertainty in the profile plot
% figprop.sourceData_sigma_color    One sigma color of the source data uncertainty in the profile plot
% figprop.sourceData_tsigma_alpha   Transparancy for the two sigma uncertainty in the profile plot.
% figprop.sourceData_sigma_alpha    Transparancy for the one sigma uncertainty in the profile plot
%                                   Use negative number not to plot and 0 to only draw the line
% figprop.sourceData_mean_linewidth Linewdith of the binned average, use 0 when you do not want to plot it
% figprop.sourceData_mean_linecolor Linecolor of the binned average
% PROPERTIES OF POLYGONS THAT NEED TO BE DRAWN IF PROFILE INTERSECTS
% figprop.poly                      Polygon coordinates (2column). Specify as a cell when more polygons.
% figprop.poly_color                Color of the polygon. Specify as cell when to change individually
% figprop.poly_style                Style of the polygon. Specify as cell when to change individually
% figprop.poly_linewidth            Linewidth of the polygon
% figprop.poly_close_flag           Set to 'y' when the polygons needs to be closed
% PROPERTIES OF ADDITIONAL DATA TO BE PLOTTED
% figprop.data                      Data to be plotted (2-column data and uncertainty optional).
%                                   Specify as a cell when to plot more datasets with e.g. different colors.
% figprop.data_lonlat               Data coordinates (2column). Specify as a cell when more polygons.
% figprop.data_color                Dataset color. Can be a cell to specify color for different datasets
% figprop.data_draw_manual          plot the errobars manual not using matlab errorbar function
% figprop.data_draw_manual_width    width in km of the manual error bars to be
%                                   drawn, default is 2 times binwidth of profile
% figprop.data_draw_manual_linestyle Linestyle to be used when drawing the error bar manual.
% figprop.data_err_linewidth        Linewidth to be used when drawing error bar manual


% modifications
% 23 May 2017       remove nan in the data first
% 1 August 2017     manual implementation of error bar drawing
%                   change order of theobservatiosn versus the error bars.

debugflag=0;
% threshold for minimum PS inside a bin
ps_min = 3;

if nargin<8
 save_flag=0;
end

if nargin<7
    figprop = [];
end
if nargin<6
    R = 200;            % default value of R in m
end
if nargin<5
    bins = 25;          % default number of bins
end
if nargin<4
    error('Specify more input arguments \n');
end


% remove any nan vlaues
ix_drop = sum(isnan([data_or data_std_or]),2)>0;
data_or(ix_drop)=[];
data_std_or(ix_drop)=[];
lonlat_or(ix_drop,:)=[];


% defaults
fontsize = 15;          % label fontsize
n_crossplot=1;
units = 'rad';
tsigma_color = [0.8 0.8 0.8];
sigma_color =[0.6 0.6 0.6];
marker_color = [0.4 0.4 0.4];
sigma_alpha = 0.5;
tsigma_alpha = 0.5;
marker_alpha = 0.3;
marker_size = 10;
binned_linewidth =1;
binned_linecolor = 'k';
poly_flag = 'n';
poly_close_flag = 'n';
extra_data_flag = 'n';
n_profiles = length(points_list);

% units
if isfield(figprop,'sourceData_units')
    units=figprop.units;
else
    fprintf(['Will assume units are in ' units '\n'])
end
% maker colors, size, and transparancy
if isfield(figprop,'sourceData_marker_color')
    marker_color = figprop.sourceData_marker_color;
end
if isfield(figprop,'sourceData_marker_size')
    marker_size = figprop.sourceData_marker_size;
end
if isfield(figprop,'sourceData_marker_alpha')
    marker_alpha = figprop.sourceData_marker_alpha;
end
% mean of the bins properties
if isfield(figprop,'sourceData_mean_linewidth')
    binned_linewidth = figprop.sourceData_mean_linewidth;
end
if isfield(figprop,'sourceData_mean_linecolor')
    binned_linecolor = figprop.sourceData_mean_linecolor;
end


% error colors and transparancy
if isfield(figprop,'sourceData_tsigma_color')
    tsigma_color = figprop.sourceData_tsigma_color;
end
if isfield(figprop,'sourceData_tsigma_alpha')
    tsigma_alpha = figprop.sourceData_tsigma_alpha;
end
if isfield(figprop,'sourceData_sigma_color')
    sigma_color = figprop.sourceData_sigma_color;
end
if isfield(figprop,'sourceData_sigma_alpha')
    sigma_alpha = figprop.sourceData_sigma_alpha;
end

% figure properties
if isfield(figprop,'n_crossplot')
    n_crossplot=figprop.n_crossplot;
end
if isfield(figprop,'fontsize')
    fontsize=figprop.fontsize;
end
if isfield(figprop,'units')
    units=figprop.units;
end
if isfield(figprop,'poly')
    poly=figprop.poly;
    poly_flag = 'y';
end
if isfield(figprop,'poly_close_flag');
    poly_close_flag = figprop.poly_close_flag;
end
if isfield(figprop,'ps_min');
    ps_min = figprop.ps_min;
end
if isfield(figprop,'data')
    extra_data_flag = 'y';
end




% loop over the number of profiles to be made
counter = 1;
counter2=1;
% n_profiles=1;
for k_profile=1:n_profiles
    points = points_list{k_profile};

    % limit the observations to that in the bounding box of the line/ polygon
    ll_bbox = min(points,[],1)- R/1000*2/110;
    ur_bbox = max(points,[],1)+ R/1000*2/110;
    ix_bbox = find(lonlat_or(:,1)>ll_bbox(1,1) & lonlat_or(:,1)<ur_bbox(1,1) & lonlat_or(:,2)>ll_bbox(1,2) & lonlat_or(:,2)<ur_bbox(1,2));
    data = data_or(ix_bbox);
    data_std = data_std_or(ix_bbox);
    lonlat = lonlat_or(ix_bbox,:);

    X_offset = 0;
    XY_new_keep = [];
    data_new_keep = [];
    data_std_new_keep = [];

    XY_AUXdata_keep =  [];
    AUXdata_keep =[];
    counter_list = 1;
    for k_segments = 1:size(points,1)-1
        point1=points(k_segments,:);
        point2=points(k_segments+1,:);

        % crossection computation
        % transform to a local reference frame, unit becomes [km]
        lonlatXY = llh2local(lonlat',point1);	% transformation of the PS pixels
        P1 = llh2local(point1',point1);         % transformation of the line points
        P2 = llh2local(point2',point1);


        % rotate the local from such that line becomes horizontal
        alpha = atan((P2(2,1)-P1(2,1))/(P2(1,1)-P1(1,1)));		% angle [rad]
        rot = [cos(alpha)    -sin(alpha)
               sin(alpha)    cos(alpha)];
        XY_new = rot\lonlatXY;                  % X is the first row, Y the second row
        P1_new = rot\P1;
        P2_new = rot\P2;
        clear alpha  P1 P2 lonlatXY

        if debugflag==1
            figure('name',['Debug_plot [1]: P' num2str(k_profile)]);
            scatter(XY_new(1,1:10:end),XY_new(2,1:10:end),3,data(1:10:end,1),'filled')
            hold on
            scatter([P1_new(1) P2_new(1)],[P1_new(2) P2_new(2)],'ro')
        end


        % swap the direction in case the second point is negative
        if P2_new(1)<0
           P2_new(1)= -1.*P2_new(1);
           XY_new(1,:)= -1.*XY_new(1,:);
           flip_xflag = 'y';
        else
           flip_xflag = 'n';
        end
        % exclude the data outside the profile
        ix_data = XY_new(1,:)>=0 & XY_new(1,:)<=P2_new(1);

        XY_new = XY_new(:,ix_data);
        data_new = data(ix_data,:);
        data_std_new = data_std(ix_data,:);


        % Search for all the PS within distance R
        ix_data = find((abs(XY_new(2,:))<=R/1000));

        % updating vectors
        XY_new = XY_new(:,ix_data);
        data_new = data_new(ix_data,:);
        data_std_new = data_std_new(ix_data,:);

        % sort all the data
        [temp,ix]=sort(XY_new(1,:));
        XY_new = XY_new(:,ix);
        data_new = data_new(ix);
        data_std_new = data_std_new(ix);


        if debugflag==1
            figure('name',['Debug_plot [2]: P' num2str(k_profile)]);
            scatter(XY_new(1,:),XY_new(2,:),3,data_new(:,1),'filled')
            hold on
            scatter([P1_new(1) P2_new(1)],[P1_new(2) P2_new(2)],'ro')
        end

        % appending the data in case its a polygon and not a stright line
        XY_new(1,:) = XY_new(1,:) + X_offset(1);
        XY_new_keep = [XY_new_keep XY_new];
        data_new_keep = [data_new_keep ; data_new];
        data_std_new_keep = [data_std_new_keep; data_std_new];
        clear XY_new

        % check if there is additional data that needs to be plotted
        if strcmpi(extra_data_flag,'y')
            % check how many datasets there are
            data_cell_flag ='n';
            if iscell(figprop.data)
                data_cell_flag = 'y';
                n_data = length(figprop.data);
            else
                n_data = 1;
            end
            for k_data =1:n_data
                % load the data temporaly
                if strcmpi(data_cell_flag,'y')
                    AUXdata_temp = figprop.data{k_data};
                    AUXdata_lonlat_temp = figprop.data_lonlat{k_data};
                else
                    AUXdata_temp = figprop.data;
                    AUXdata_lonlat_temp = figprop.data_lonlat;
                end

                AUXdata_xy_temp = llh2local(AUXdata_lonlat_temp',point1)';
                AUXdata_xy_new = (rot\AUXdata_xy_temp');
                clear AUXdata_xy_temp AUXdata_lonlat_temp

                % flip the data if needed
                if strcmpi(flip_xflag,'y')
                   AUXdata_xy_new(1,:)= -1.* AUXdata_xy_new(1,:);
                end

                % exclude the data outside the profile
                ix_data_temp = AUXdata_xy_new(1,:)>=0 & AUXdata_xy_new(1,:)<=P2_new(1);


                AUXdata_xy_new = AUXdata_xy_new(:,ix_data_temp);
                AUXdata_temp = AUXdata_temp(ix_data_temp,:);
                clear ix_data_temp

                % Search for all the PS within distance R
                ix_data_temp = find((abs(AUXdata_xy_new(2,:))<=R/1000));

                % updating vectors
                AUXdata_xy_new = AUXdata_xy_new(:,ix_data_temp);
                AUXdata_temp = AUXdata_temp(ix_data_temp,:);
%                 labels = figprop.names(ix_data_temp)

                clear ix_data_temp

                % populate the colors
                if  isfield(figprop,'data_color')
                    if iscell(figprop.data_color)
                        data_color = figprop.data_color{k_data};
                    else
                        data_color = figprop.data_color;
                    end
                else
                    data_color = 'r';
                end

                for k_errorbars=1:size(AUXdata_xy_new,2)
                    AUXdata_color{counter_list}=data_color;
                    counter_list=counter_list+1;
                end

                % generate the data vector
                AUXdata_xy_new(1,:) = AUXdata_xy_new(1,:) + X_offset(1);
                XY_AUXdata_keep = [XY_AUXdata_keep AUXdata_xy_new];
                AUXdata_keep = [AUXdata_keep ; AUXdata_temp];
                clear AUXdata_xy_new AUXdata_temp



            end
        end

        % increment the offset in case this is a polygon profile not a straight line
        X_offset = X_offset+P2_new(1);


    end

    % output to the screen
    PS_used = size(data_new_keep,1);
    if (isempty(PS_used)==1 || PS_used==0)
        error('No PS are found within the tube crossection. \n')
    else
        fprintf([num2str(PS_used),' PS used to compute projection on crossection \n'])
    end

    clear data_new XY_new
    data_new =data_new_keep;
    data_std_new = data_std_new_keep;
    XY_new = XY_new_keep;


    % binning of the results
    if isfield(figprop,'binsize')
        binsize=figprop.binsize/1000;
        bins = ceil((max(XY_new(1,:))-min(XY_new(1,:)))/binsize);
    else
        binsize = (max(XY_new(1,:))-min(XY_new(1,:)))/bins;
    end

    data_cross_binned = zeros([1 bins]);
    bins_xy = zeros([1 bins]);
    clear ix

    bins_xy = zeros([1 bins]);
    data_binned_cov = zeros([size(data,2) bins]);
    data_binned_var = zeros([size(data,2) bins]);
    data_binned_var_std = zeros([size(data,2) bins]);
    for k=1:bins
        bl = (k-1)*binsize+min(XY_new(1,:));                        % lower bound
        bu = (k)*binsize+min(XY_new(1,:));                          % upper bound
        ix = find(bl<=XY_new(1,:)  & XY_new(1,:)<bu);
        % Only showing bins with a minimum of PS contained
        if size(ix,2)>=ps_min
            bins_xy(1,k) = bl+binsize/2;                            % position center of bin
            A = ones([size(ix,2) 1]);

            %data_var = diag(data_std(ix).^2);
            data_var = diag(data_std_new(ix).^2);

            data_binned_var(:,k) = (inv(A'*inv(data_var)*A)*A'*inv(data_var)*data_new(ix,:));
            data_binned_var_std(1,k) = sqrt(inv(A'*inv(data_var)*A));
            data_binned_var_std(1,k) = sqrt((mean(diag(data_var))));

        else
            bins_xy(:,k) = NaN;
            data_binned_var(:,k) = NaN;
            data_binned_var_std (1,k) = NaN;

        end
        clear ix bl bu
    end


    %%% PLOTTING
    if counter==1
        profile_fig = figure('position',[           1          79        1424         726]);
    else
        try
            figure(profile_fig)
        catch
            profile_fig  =  figure('position',[           1          79        1424         726]);
        end
    end
%     hsubplot(counter2) = subplot(n_crossplot,1,counter);

% not correct, need to pudate
% Add the labels ----------------------------------------- BZ
%     disp(counter2)
%     text(0.0, 1, char(counter2 + 64), 'Units', 'normalized', ...
%                     'FontSize', 20, 'Color', 'w', 'BackgroundColor', 'k') 
%     text(0.99, 1, strcat(char(counter2 + 64),"'"), 'Units', 'normalized', ... 
%                     'FontSize',20, 'Color', 'w', 'BackgroundColor', 'k')
% 
%     hold on

%     hsubplot(counter2) = subplot('position',[           1          79        1424         726]);

%%% turned this off to do each separately
%     hsubplot(counter2) = subplot(n_crossplot,1,counter);
%     if k_profile > 1
%     hsubplot(counter2) = subplot('Position',[1          79+k_profile*1424        1424         726])
%     end


    spwidth  = 0.8;
    spheight = 0.225;
    left     = 0.115;
    bot      = 0.1;
    sp  = 0.30;
    xlimits21=figprop.xlimits(1)./2;
    xlimits22=figprop.xlimits(2)./2;
    
    % override temporarily BZ
    xlimits21=figprop.xlimits(1);
    xlimits22=figprop.xlimits(2);

    if k_profile==1
        hsubplot(counter2) = subplot('Position',[left bot+2*sp spwidth spheight]);
        set(hsubplot(counter2),'XLim',figprop.xlimits)
    elseif k_profile == 2

        hsubplot(counter2) = subplot('Position',[left bot+1*sp spwidth/2 spheight]);
        set(hsubplot(counter2),'XLim', [xlimits21, xlimits22]);
        set(hsubplot(counter2),'xtick', 0:10:xlimits22);

    elseif k_profile == 3
        hsubplot(counter2) = subplot('Position',[left bot spwidth/2 spheight]);
        set(hsubplot(counter2),'XLim', [xlimits21, xlimits22]);
        set(hsubplot(counter2),'xtick', 0:10:xlimits22);

    end
    hold on


% % % % % % %     % plot the scatter points
% % % % % % %     %scatter(XY_new(1,:),data_new(:,1),15,'MarkerEdgeColor',marker_color)
% % % % % % %     htemp = plot(XY_new(1,:),data_new(:,1),'.','Color',marker_color);
% % % % % % %
% % % % % % %     % plot the mean
% % % % % % %     hold on
% % % % % % %     plot(bins_xy,data_binned_var,'.-','color','k','linewidth',2)
% % % % % % %     hold on
% % % % % % %     ylabel(units,'fontsize',fontsize)
% % % % % % %     set(gca,'fontsize',fontsize)
% % % % % % %     box on
% % % % % % %     grid on



    % the polygon plots requires special approach to split it into pieces
    clear ix
    ix = find(isnan(bins_xy)==1);

    % keep backup
    data_binned_var_or = data_binned_var;
    data_binned_var_std_or = data_binned_var_std;
    bins_xy_or = bins_xy;

    ti = strcat('poly', num2str(k_profile));
    polycoordsS = zeros(length(ix)+1, 5000);
    polyvalsS   = zeros(length(ix)+1, 5000);
    % they're created on the first round of this, loaded on second
    if isfile(strcat(ti, '.mat'))
      load (strcat(ti, '.mat'))
    end

    if ~isempty(ix)
        for k=1:length(ix)
            if ix(k)>1
                ix_current = ix(k)-1;

                if ix_current>2
                    % plotting 2 sigma
                    poly_val = [data_binned_var_or(1:ix_current)-2*data_binned_var_std_or(1:ix_current) fliplr(data_binned_var_or(1:ix_current)+2*data_binned_var_std_or(1:ix_current))];
                    poly_coord = [bins_xy_or(1:ix_current) fliplr(bins_xy_or(1:ix_current))];
                    h2= fill(poly_coord,poly_val,tsigma_color,'EdgeColor',tsigma_color,'FaceAlpha',tsigma_alpha);
                    % plotting 1 sigma
                    if sigma_alpha>=0
                        poly_val = [data_binned_var_or(1:ix_current)-data_binned_var_std_or(1:ix_current) fliplr(data_binned_var_or(1:ix_current)+data_binned_var_std_or(1:ix_current))];
                        poly_coord = [bins_xy_or(1:ix_current) fliplr(bins_xy_or(1:ix_current))];

                        h1 = fill(poly_coord,poly_val,sigma_color,'EdgeColor',sigma_color,'FaceAlpha',sigma_alpha);

                        polycoordsS(k, 1:size(poly_coord, 2)) = poly_coord;
                        polyvalsS(k, 1:size(poly_val, 2)) = poly_val;

                    end
                end

                if ix(k)>length(data_binned_var_std_or)
                    data_binned_var_std_or(1:end)=[];
                    data_binned_var_or(1:end)=[];
                    bins_xy_or(1:end)=[];
                else
                    data_binned_var_std_or(1:ix(k))=[];
                    data_binned_var_or(1:ix(k))=[];
                    bins_xy_or(1:ix(k))=[];
                end
                ix = ix-ix(k);
            else

                data_binned_var_std_or(1:ix(k))=[];
                data_binned_var_or(1:ix(k))=[];
                bins_xy_or(1:ix(k))=[];
                ix = ix-ix(k);
            end
            % make sure to plot the remainder
            if k==length(ix)
                % plotting 2 sigma
                try
                    poly_val = [data_binned_var_or-2*data_binned_var_std_or fliplr(data_binned_var_or+2*data_binned_var_std_or)];
                catch
                    keyboard
                end
                poly_coord = [bins_xy_or fliplr(bins_xy_or)];
                h2= fill(poly_coord,poly_val,tsigma_color,'EdgeColor',tsigma_color,'FaceAlpha',tsigma_alpha);
                % plotting 1 sigma
                if sigma_alpha>=0
                    poly_val = [data_binned_var_or-data_binned_var_std_or fliplr(data_binned_var_or+data_binned_var_std_or)];
                    poly_coord = [bins_xy_or fliplr(bins_xy_or)];

                    h1 = fill(poly_coord,poly_val,sigma_color,'EdgeColor',sigma_color,'FaceAlpha',sigma_alpha);
                    polycoordsS(k+1, 1:size(poly_coord, 2)) = poly_coord;
                    polyvalsS(k+1, 1:size(poly_val, 2)) = poly_val;
                end
            end
        end

    else % BB - doesn't get run
        % plotting 2 sigma
        poly_val = [data_binned_var_or-2*data_binned_var_std_or fliplr(data_binned_var_or+2*data_binned_var_std_or)];
        poly_coord = [bins_xy_or fliplr(bins_xy_or)];
        h2= fill(poly_coord,poly_val,tsigma_color,'EdgeColor',tsigma_color,'FaceAlpha',tsigma_alpha);

        if sigma_alpha>=0
            % plotting 1 sigma
            poly_val = [data_binned_var_or-data_binned_var_std_or fliplr(data_binned_var_or+data_binned_var_std_or)];
            poly_coord = [bins_xy_or fliplr(bins_xy_or)];
            h1 = fill(poly_coord,poly_val,sigma_color,'EdgeColor',sigma_color,'FaceAlpha',sigma_alpha);
        end
    end
    hold on
    if (save_flag)
      polycoords = polycoordsS; polyvals = polyvalsS;
      save(ti, 'polycoords', 'polyvals')
    end
    % ------------------------- plot the external data; BZ
    for k=1:size(polycoords, 1)
     % if (k<size(polycoords, 1)); continue; end
     % if (k<5); continue; end
      % extreme x values
      polycoord_ext = polycoords(k, :);
      polycoord_ext = polycoord_ext(1:find(polycoord_ext,1,'last'));
      polyval_ext   = polyvals(k, 1:size(polycoord_ext, 2)); % may eventually be a problem

      if (isequal(sum(polycoord_ext), 0)) continue; end
      h1 = fill(polycoord_ext,polyval_ext,sigma_color,'EdgeColor',figprop.ext_col, 'LineWidth', figprop.ext_lw, 'LineStyle', '-.','FaceAlpha',0);


    end


    % plot the scatter points
    scatter(XY_new(1,:),data_new(:,1),marker_size,'MarkerEdgeColor',marker_color,'MarkerfaceColor',marker_color,'MarkerEdgealpha',marker_alpha,'Markerfacealpha',marker_alpha)

    % plot the mean
    hold on
    if binned_linewidth~=0
        plot(bins_xy,data_binned_var,'.-','color',binned_linecolor,'linewidth',binned_linewidth)
    end
    hold on
    ylabel(units,'fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    box on
    grid on
    if isfield(figprop,'ylimits');
        ylim(figprop.ylimits)
    end
    ylimits = get(gca,'ylim');

    if counter==n_crossplot
        counter=1;
        xlh = xlabel('Distance along crossection (km)','fontsize',fontsize);

    else
        if counter==2
          set(hsubplot(counter),'Xticklabel',[]) % turn off ticklabels
        end
        xlabel('','fontsize',fontsize)

        counter= counter+1;
    end
    counter2 = counter2+1;

   % plotting axilary data
   if ~isempty(AUXdata_keep)
       for k_errorbars = 1:size(XY_AUXdata_keep,2)

           hold on
           % see if the error needs to be drawn manually or using erorobar
           % command of matlab
           if ~isfield(figprop,'data_draw_manual')
               figprop.data_draw_manual='n'
           end
           if strcmpi(figprop.data_draw_manual,'n')
               hold on
               herror2 = errorbar(XY_AUXdata_keep(1,k_errorbars),AUXdata_keep(k_errorbars,1),AUXdata_keep(k_errorbars,2)); % BZ ONE SIGMA
               % herror2 = errorbar(XY_AUXdata_keep(1,k_errorbars),AUXdata_keep(k_errorbars,1),2*AUXdata_keep(k_errorbars,2));
               set(herror2,'color',AUXdata_color{k_errorbars},'linewidth',4.2); % BZ increase GPS/GNSS line width



           else
               % horizontal lines draw them with the width based on the binsize
               if isfield(figprop,'data_draw_manual_width');
                  error_width = figprop.data_draw_manual_width;
               else
                   error_width = binsize;
               end
               if isfield(figprop,'data_draw_manual_linestyle');
                   error_linestyle = figprop.data_draw_manual_linestyle;
               else
                   error_linestyle='-';
               end
               if isfield(figprop,'data_err_linewidth');
                   data_err_linewidth=figprop.data_err_linewidth;
               else
                   data_err_linewidth=4;
               end
               x_coordinates = [XY_AUXdata_keep(1,k_errorbars)-error_width/2 XY_AUXdata_keep(1,k_errorbars)+error_width/2];
               y_coordinates = [AUXdata_keep(k_errorbars,1) - 2*AUXdata_keep(k_errorbars,2)  AUXdata_keep(k_errorbars,1) + 2*AUXdata_keep(k_errorbars,2)];

               % manual drawing of the error bar
               herror2 = plot(x_coordinates,[y_coordinates(1) y_coordinates(1)]);

               set(herror2,'color',AUXdata_color{k_errorbars},'linewidth',data_err_linewidth,'linestyle',error_linestyle);
               herror2 = plot(x_coordinates,[y_coordinates(2) y_coordinates(2)]);
               set(herror2,'color',AUXdata_color{k_errorbars},'linewidth',data_err_linewidth,'linestyle',error_linestyle);
               herror2 = plot([XY_AUXdata_keep(1,k_errorbars) XY_AUXdata_keep(1,k_errorbars)],[y_coordinates]);
               set(herror2,'color',AUXdata_color{k_errorbars},'linewidth',data_err_linewidth,'linestyle',error_linestyle);
           end
       end
   end
    %% computing polygons if requested
    if strcmpi(poly_flag,'y')
        X_offset = 0;
        for k_segments_profile = 1:size(points,1)-1
            point1=points(k_segments_profile,:);
            point2=points(k_segments_profile+1,:);

            % crossection computation
            % transform to a local reference frame, unit becomes [km]
            P1 = llh2local(point1',point1);                 % transformation of the line points
            P2 = llh2local(point2',point1);

            % rotate the local from such that line becomes horizontal
            alpha = atan((P2(2,1)-P1(2,1))/(P2(1,1)-P1(1,1)));		% angle [rad]
            rot = [cos(alpha)    -sin(alpha)
                   sin(alpha)    cos(alpha)];
            P1_new = rot\P1;
            P2_new = rot\P2;
            clear alpha  P1 P2

            % swap the direction in case the second point is negative
            if P2_new(1)<0
               P2_new(1)= -1.*P2_new(1);
               flip_xflag = 'y';
            else
               flip_xflag = 'n';
            end

            for k_poly=1:length(poly)
                polyXY = llh2local(poly{k_poly}',point1);	% transformation of the PS pixels

                % rotate the local from such that line becomes horizontal
                poly_new = (rot\polyXY)';                  % X is the first column, Y the second column
                clear polyXY

                % close the polygon if requested
                if strcmpi(poly_close_flag,'y')
                    poly_new = [poly_new; poly_new(1,:)];
                end

                % see if the x coordinates where flipped for the other data, if
                % so do it here too.
                if strcmpi(flip_xflag,'y')
                    poly_new(:,1)=-1.*poly_new(:,1);
                end

                % loop over each segement of the poly and see if its crossed at
                % the x-axis within the range of P2 x-coordinate

                for k_segments=1:size(poly_new,1)-1
                    % define a line
                    poly_P1 = poly_new(k_segments,:);
                    poly_P2 = poly_new(k_segments+1,:);
                    % compute the gradient
                    gradient = (poly_P1(2)-poly_P2(2))/(poly_P1(1)-poly_P2(1));
                    offset = poly_P1(2)-gradient*poly_P1(1);
                    % evaluate the line at zero y
                    x_coord = -1*offset/gradient;
                    % check if the coordinate is within the profile
                    if x_coord>0 && x_coord<P2_new(1) && (min(poly_new(k_segments:k_segments+1,2))<=0 & max(poly_new(k_segments:k_segments+1,2))>=0)

                        % check if the crossing is within the line segment
                        if ~( min([poly_P1(1) poly_P2(1)])>x_coord || max([poly_P1(1) poly_P2(1)])<x_coord || min([poly_P1(2) poly_P2(2)])>0  || max([poly_P1(2) poly_P2(2)])<0 )
                            % appending the data in case its a polygon and not a stright line
                            x_coord = x_coord +X_offset;



                            h_line = plot([x_coord x_coord],[ylimits],'--','linewidth',2);
                            if isfield(figprop,'poly_color');
                                if iscell(figprop.poly_color)
                                    if length(figprop.poly_color)>1
                                        set(h_line,'color',figprop.poly_color{k_poly})
                                    else
                                        set(h_line,'color',figprop.poly_color{1})
                                    end
                                else
                                    set(h_line,'color',figprop.poly_color)
                                end
                            else
                                set(h_line,'color','r')
                            end
                            if isfield(figprop,'poly_style');
                                 if iscell(figprop.poly_style)
                                     if length(figprop.poly_style)>1
                                        set(h_line,'linestyle',figprop.poly_style{k_poly})
                                    else
                                        set(h_line,'linestyle',figprop.poly_style{1})
                                     end
                                 else
                                     set(h_line,'linestyle',figprop.poly_style)
                                 end
                            else
                                set(h_line,'linestyle','--')
                            end
                            if isfield(figprop,'poly_linewidth');
                                 if iscell(figprop.poly_linewidth)
                                     if length(figprop.poly_linewidth)>1
                                        set(h_line,'linewidth',figprop.poly_linewidth{k_poly})
                                    else
                                        set(h_line,'linewidth',figprop.poly_linewidth{1})
                                     end
                                 else
                                     set(h_line,'linewidth',figprop.poly_linewidth)
                                 end
                            else
                                set(h_line,'linewidth',2)
                            end


                        end
                    end
                end

            end
            clear rot
            X_offset = X_offset+P2_new(1);
        end
    end
end

% return % BZ
% move xlabel to center and pad a bit
xlh_pos=get(xlh,'position');
xlh_pos(1)= xlh_pos(1) + 30;
xlh_pos(2)= xlh_pos(2) - 2;
set(xlh,'position',xlh_pos);
for k=1:length(hsubplot)

    if isfield(figprop,'xticks')
         set(hsubplot(k),'xtick',figprop.xticks);
    end
    if isfield(figprop,'yticks')
         set(hsubplot(k),'ytick',figprop.yticks);
    end

end
