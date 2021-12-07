close all

type = 'aldolase'; % clarify [urea]
concentration = 'aldolase_30ms'; % file name format
video = [107:159];
Dmax = 20; % um^2/s
gap = 3;
FittingType = 1; % unbrownian: FittingType=1; brownian: FittingType=0
Mini_Trajlength = 10; % define the mini trajectory length to rule out those short trajectories
R_threshold = 0.9; % R^2>0.9 for good trajectory

frame_interval = 0.06; % unit: s
clip_factor =0.25;  % usually set as 0.25, only fit the first one quarter of the MSD plot
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
%         trackName = sprintf('Track %d', index );
%         
%         x = tracks{index}(:,2);
%         y = tracks{index}(:,3);
%         t = tracks{index}(:,1);
%         
%         hps(i) =  plot(ha, x, y, ...
%             'Color', colors(i,:), ...
%             'DisplayName', trackName );
        
    end
         tracks = tracks(~cellfun('isempty',tracks));
         
         ma = msdanalyzer(2, 'µm', 's');
         ma = ma.addAll(tracks);
         ma = ma.computeMSD;
         
         n_spots = numel(ma.msd);

         if clip_factor < 1
            fprintf('Fitting %d curves of MSD = f(t), taking only the first %d%% of each curve... ',...
            n_spots, ceil(100 * clip_factor) )
         else
            fprintf('Fitting %d curves of MSD = f(t), taking only the first %d points of each curve... ',...
            n_spots, round(clip_factor) )
         end
         
         % define variables
         a = NaN(n_spots, 1);
         b = NaN(n_spots, 1);
         r2fit = NaN(n_spots, 1);
         c = NaN(n_spots, 1); 
         ft2_r2fit = NaN(n_spots, 1); 
         d = NaN(n_spots, 1); 
         e = NaN(n_spots, 1); 
         ft3_r2fit = NaN(n_spots, 1); 
         Tralength = NaN(n_spots, 1); 
         ft1 = fittype('poly1'); % y=ax+b
         ft2 = fittype({'x'}); % y=ax
    
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
    
            % Clip data, never take the first one dt = 0
            if clip_factor < 1
                t_limit = 2 : round(numel(t) * clip_factor);
            else
                t_limit = 2 : min(1+round(clip_factor), numel(t));
            end
            x = t(t_limit); % x => t
            Y = y(t_limit); % Y => MSD
            W = w(t_limit);
        
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
            % fit3 parameter power law log-log fit
            d(i_spot) = fo3.p1;
            e(i_spot) = fo3.p2;
            ft3_r2fit(i_spot) = gof3.adjrsquare;
    
            Tralength(i_spot) = size(t, 1); 
            
            %Mengqi addded for plotting and saving data
            
            if FittingType == 1
                
                %find good trajectory unbrownian power law fit
                FittingTypeS = 'unbrownian'; % define fitting type string
                if gof3.adjrsquare > R_threshold;
                    %save MSD data
                    MSD_ind_data = [t,y,stder,w];
                    T_MSD_ind_data = array2table(MSD_ind_data, 'VariableNames',{'t_s', 'msd_um2', 'std','N'});
                    writetable(T_MSD_ind_data, sprintf('%s_%s_%d_MSDdata_Traj%d.xlsx', concentration, FittingTypeS, video(k), i_spot));
                end
                
            else
                
                %find good trajectory brownian fit
                FittingTypeS = 'brownian'; % define fitting type string
                if gof.adjrsquare > 0.9
                    if gof2.adjrsquare > 0.9
                        %save MSD data
                        MSD_ind_data = [t,y,stder,w];
                        T_MSD_ind_data = array2table(MSD_ind_data, 'VariableNames',{'t_s', 'msd_um2', 'std','N'});
                        writetable(T_MSD_ind_data, sprintf('%s_%s_%d_MSDdata_Traj%d.xlsx', concentration, FittingTypeS, video(k), i_spot)); 
                    end
                end
            end
         end
         
         fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
  
         if FittingType == 1  
             %find good trajectory unbrownian power law fit
             FittingTypeS = 'unbrownian'; % define fitting type string
             good_enough_fit = ft3_r2fit > R_threshold;
             ind_good = find(good_enough_fit == 1); % vector
         else
             %find good trajectory brownian fit
             FittingTypeS = 'brownian'; % define fitting type string
             good_enough_fit1 = r2fit > 0.9;
             good_enough_fit2 = ft2_r2fit > 0.9;
             good_enough_fit = good_enough_fit1 & good_enough_fit2;
             ind_good = find(good_enough_fit == 1); % vector
         end
 
         n_goodTracks = numel(ind_good);
    
         if n_goodTracks > 0
             
             % fitting results...................................
             % results for fit1 (y = ax + b)
             All_D = a(good_enough_fit)/ 2 / ma.n_dim;
             All_D_r2fit = r2fit(good_enough_fit);
             % results for fit2 (y = ax)
             All_D_2 = c(good_enough_fit)/ 2 / ma.n_dim;
             All_D_r2fit_2 = ft2_r2fit(good_enough_fit);
             % results for alpha fit
             All_D_alpha = d(good_enough_fit);
             All_D_3= exp(e(good_enough_fit))/2 / ma.n_dim;
             All_D_r2fit_3 = ft3_r2fit(good_enough_fit);
    
             All_D_timelength = Tralength(good_enough_fit);

             AllInfo = [ind_good, All_D,All_D_2, All_D_r2fit, All_D_r2fit_2, All_D_timelength,All_D_alpha,All_D_3, All_D_r2fit_3];

             %save individual video fitting results
             FinalTable = array2table(AllInfo, 'VariableNames',{'Traj','D_um2_per_s','Dtwo_um2_per_s','R2','R2two','Trajlength_frame','alpha','Dgeneral_um2_per_s','R2_alpha'})
             writetable(FinalTable, sprintf('%s_%s_All_D_%d.xlsx', concentration, FittingTypeS, video(k))); 

             Video_ID= zeros(size(AllInfo,1),1);
             Video_ID = Video_ID + video(k);
             T = [T; Video_ID AllInfo];
             
             % plot good trajectories............................
             colors = jet(n_goodTracks);
             hps = NaN(n_goodTracks, 1);
        
             for i = 1 : n_goodTracks
                track = tracks{ind_good(i)};
        
                trackName = sprintf('Track %d', ind_good(i));
        
                xi = track(:,2);
                yi = track(:,3);
    
                hps(i) =  plot(xi, yi, ...
                    'Color', colors(i,:), ...
                    'DisplayName', trackName );
                hold on
             end
            
         end
         clear tracks
         
end

% save all Traj plots
title(sprintf('MTT %s: Dmax=%d um^2/s, gap=%d, %s, Traj>%d, R^2>%s', type, Dmax, gap, FittingTypeS, Mini_Trajlength, R_thresholdS));
xlabel('x (um)') 
ylabel('y (um)')
saveas(gcf, sprintf('A_PlotTrajFromMTT_%s_Dmax=%d_gap=%d_%s_Traj>%d_R>%s_All.fig', concentration, Dmax, gap, FittingTypeS, Mini_Trajlength, R_thresholdS));

% save all fitting results
AllD = array2table(T, 'VariableNames',{'Video','Traj','D','Dtwo','R','Rtwo','length','alpha','Dgeneral','R2_alpha'})

writetable(AllD, sprintf('A_AllDcollections_%s_Dmax=%d_gap=%d_%s_Traj>%d_R>%s.xlsx', concentration, Dmax, gap, FittingTypeS, Mini_Trajlength, R_thresholdS));

% unit of D is um^2/s
