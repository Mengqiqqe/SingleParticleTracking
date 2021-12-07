clear; clc; clearvars -global; close all; 

type = '1mM'; % clarify [urea]
concentration = '1mM'; % file name format
video = [18:41];
Dmax = 5; % um^2/s
gap = 2;
FittingType = 0; % unbrownian: FittingType=1; brownian: FittingType=0
Mini_Trajlength = 25; % define the mini trajectory length to rule out those short trajectories
R_threshold = 0.9; % R^2>0.9 for good trajectory
frame_interval = 0.06; % unit: s
clip_factor_shortTraj= 0.2;  % usually set as 0.25, only fit the first one quarter of the MSD plot
clip_factor_longTraj = 0.1; % clip_factor for very longTraj due to heterogeneity
R_thresholdS = sprintf('%.1f', R_threshold);

T = [];
for k = 1:numel(video)
    
    load(sprintf('%s_%03d_Tracked.mat',concentration, video(k)));

    n_tracks = size(trackedPar,2);
    indices = [1:n_tracks]';
    tracks={};

    for i = 1 : n_tracks
        
        index = indices(i);
  
        if numel(trackedPar(index).TimeStamp) > Mini_Trajlength
            tracks{index} = [trackedPar(index).TimeStamp trackedPar(index).xy];
        end      
        
    end
        
    tracks = tracks(~cellfun('isempty',tracks));

    ma = msdanalyzer(2, 'µm', 's');
    ma = ma.addAll(tracks);
    ma = ma.computeMSD;

    %fit MSD
    
    n_spots = numel(ma.msd);

    fprintf('Fitting %d curves of MSD = f(t), taking only the first %d%% - %d%% of each curve... ',...
    n_spots, ceil(100 * clip_factor_shortTraj), ceil(100 * clip_factor_longTraj) )
    
    %if clip_factor < 1
        %fprintf('Fitting %d curves of MSD = f(t), taking only the first %d%% of each curve... ',...
        %n_spots, ceil(100 * clip_factor) )
    %else
        %fprintf('Fitting %d curves of MSD = f(t), taking only the first %d points of each curve... ',...
        %n_spots, round(clip_factor) )
    %end

    a = NaN(n_spots, 1);
    b = NaN(n_spots, 1);
    r2fit = NaN(n_spots, 1);
    c = NaN(n_spots, 1); % Mengqi added
    ft2_r2fit = NaN(n_spots, 1); % Mengqi added
    d = NaN(n_spots, 1); % Mengqi added
    e = NaN(n_spots, 1); % Mengqi added
    ft3_r2fit = NaN(n_spots, 1); % Mengqi added
    Tralength = NaN(n_spots, 1); % Mengqi added
    ft1 = fittype('poly1');
    ft2 = fittype({'x'});
    
    fprintf('%4d/%4d', 0, n_spots);
    for i_spot = 1 : n_spots
    
        fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i_spot, n_spots);
    
        msd_spot = ma.msd{i_spot};
    
        time = msd_spot(:,1);
        yaix = msd_spot(:,2);
        stderaix = msd_spot(:,3);
        waix = msd_spot(:,4);
    
        % Thrash bad data
        nonnan = ~isnan(yaix);
        t = time(nonnan);
        y = yaix(nonnan);
        stder = stderaix(nonnan);
        w = waix(nonnan);
        Tralength(i_spot) = size(t, 1); % Mengqi added
        
        
        % Clip data, never take the first one dt = 0
        
        if Tralength(i_spot) > 300
            clip_factor = clip_factor_longTraj;
        else 
            clip_factor = clip_factor_shortTraj;
        end
       
        if clip_factor < 1
            t_limit = 2 : round(numel(t) * clip_factor);
        else
            t_limit = 2 : min(1+round(clip_factor), numel(t));
        end
        x = t(t_limit);
        Y = y(t_limit);
        W = w(t_limit);
        
        %added
%         timelag = 2:20;
%         logx = t(timelag);
%         logy = y(timelag);
%         logw = w(timelag);
        %added
        
        if numel(Y) < 2
            continue
        end
    
        [fo, gof] = fit(x, Y, ft1, 'Weights', W);
        [fo2, gof2] = fit(x, Y, ft2, 'Weights', W);
        [fo3, gof3] = fit(log(x), log(Y), ft1, 'Weights', W);
        % fit1 parameter y=ax+b
        a(i_spot) = fo.p1;
        b(i_spot) = fo.p2;
        r2fit(i_spot) = gof.adjrsquare;
        % fit2 parameter y=ax
        c(i_spot) = fo2.a;
        ft2_r2fit(i_spot) = gof2.adjrsquare;
        % fit3 parameter log-log
        d(i_spot) = fo3.p1;
        e(i_spot) = fo3.p2;
        ft3_r2fit(i_spot) = gof3.adjrsquare;
 
        
        % Mengqi addded for plotting and saving data
        if FittingType == 1
   
           if ft3_r2fit(i_spot) > R_threshold
            %save MSD data
            MSD_ind_data = [t,y,stder,w];
            T_MSD_ind_data = array2table(MSD_ind_data, 'VariableNames',{'t_s', 'msd_um2', 'std','N'});
            writetable(T_MSD_ind_data, sprintf('%s_%d_MSDdata_Traj%d.xlsx', concentration, video(k), i_spot)); 
        
            % save MSD plot
%             h = figure;
%             plot(fo,'r',x,Y);
%             hold on 
%             plot(fo2,'g',x,Y);
%             plot (t,y);
%             hold off
%             saveas(h, sprintf('%s_%d_MSD_Traj%d.png', concentration, video(k), i_spot));
            % save logMSD plot
%             h = figure;
%             plot(fo3,'r',log(x),log(Y));
%             hold on
%             plot(log(t),log(y));
%             hold off
%             saveas(h, sprintf('%s_%d_LOGMSD_Traj%d.png', concentration, video(k), i_spot));
             
              loglog(t,y);
              hold on 
            end
           
        else
        
         if r2fit(i_spot) > R_threshold
            if ft2_r2fit(i_spot) > R_threshold
   
            %save MSD data
            MSD_ind_data = [t,y,stder,w];
            T_MSD_ind_data = array2table(MSD_ind_data, 'VariableNames',{'t_s', 'msd_um2', 'std','N'});
            writetable(T_MSD_ind_data, sprintf('%s_%d_MSDdata_Traj%d.xlsx', concentration, video(k), i_spot)); 
        
            % save MSD plot
%             h = figure;
%             plot(fo,'r',x,Y);
%             hold on 
%             plot(fo2,'g',x,Y);
%             plot (t,y);
%             hold off
%             saveas(h, sprintf('%s_%d_MSD_Traj%d.png', concentration, video(k), i_spot));
            % save logMSD plot
%             h = figure;
%             plot(fo3,'r',log(x),log(Y));
%             hold on
%             plot(log(t),log(y));
%             hold off
%             saveas(h, sprintf('%s_%d_LOGMSD_Traj%d.png', concentration, video(k), i_spot));
              
              % plot in log-log scale
              %loglog(t,y);
              
              % plot in normal scale
              plot(t,y);
              
              hold on 
            end
         end
       end
    end
    
    fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
  
    %find good trajectory
    if FittingType == 1
      %find good trajectory unbrownian power law fit
      FittingTypeS = 'unbrownian'; % define fitting type string
      good_enough_fit = ft3_r2fit > R_threshold;
      ind_good = find(good_enough_fit == 1); % vector
    else
      %find good trajectory brownian fit
      FittingTypeS = 'brownian'; % define fitting type string
      good_enough_fit1 = r2fit > R_threshold;
      good_enough_fit2 = ft2_r2fit > R_threshold;
      good_enough_fit = good_enough_fit1 & good_enough_fit2;
      ind_good = find(good_enough_fit == 1); % vector
    end

  if size(ind_good, 1) > 0
    % results for fit1 (y = ax + b)
    All_D = a(good_enough_fit)/ 2 / ma.n_dim;
    All_D_r2fit = r2fit(good_enough_fit);
    All_LocError = b(good_enough_fit); % localization error
    % results for fit2 (y = ax)
    All_D_2 = c(good_enough_fit)/ 2 / ma.n_dim;
    All_D_r2fit_2 = ft2_r2fit(good_enough_fit);
    % results for alpha fit
    All_D_alpha = d(good_enough_fit);
    All_D_3= exp(e(good_enough_fit))/2 / ma.n_dim;
    All_D_r2fit_3 = ft3_r2fit(good_enough_fit);
    
    All_D_timelength = Tralength(good_enough_fit);

    AllInfo = [ind_good, All_D, All_LocError, All_D_r2fit, All_D_2, All_D_r2fit_2, All_D_timelength,All_D_alpha,All_D_3, All_D_r2fit_3];

    %save individual vedio MSD
    FinalTable = array2table(AllInfo, 'VariableNames',{'Traj','D_um2_per_s','LocError','R2','Dtwo_um2_per_s','R2two','Trajlength_frame','alpha','Dgeneral_um2_per_s','R2_alpha'})
    writetable(FinalTable, sprintf('%s_All_D_%d.xlsx', concentration, video(k))); 

    %store all D in a new folder
%     folder = sprintf('/Users/SuMaggie/Desktop/AllD_%s', concentration);
%     if ~exist(folder, 'dir')
%         mkdir(folder);
%     end
%     baseFileName = sprintf('%s_All_D_%d.csv',concentration, j);
%     fullFileName = fullfile(folder, baseFileName);
%     xlswrite(fullFileName,AllInfo);
 
      Video_ID= zeros(size(AllInfo,1),1);
      Video_ID = Video_ID + video(k);
      T = [T; Video_ID AllInfo];
  end
    
end

%added
hold off
saveas(gcf, sprintf('%s_logScale_MSD ALL.png', type)); 
%added

AllD = array2table(T, 'VariableNames',{'Video','Traj','D','LocError','R','Dtwo','Rtwo','length','alpha','Dgeneral','R2_alpha'})

writetable(AllD, sprintf('AllDcollections_Dmax=%d_Gap=%d_%s_%s_Traj>%d_R>%s.xlsx', Dmax, gap, type, FittingTypeS, Mini_Trajlength, R_thresholdS)); %concentration, FittingTypeS, Mini_Trajlength, R_thresholdS

% unit of D is um^2/s  

