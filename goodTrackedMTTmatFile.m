clear; clc; clearvars -global; close all; 

type = 'buffer'; % clarify [urea]
concentration = 'ImagesD2_9.5ms_stack'; % file name format
video = [1];
Dmax = 2.5; % um^2/s
gap = 2;
FittingType = 1; % unbrownian: FittingType=1; brownian: FittingType=0
Mini_Trajlength = 0; % define the mini trajectory length to rule out those short trajectories
R_threshold = 0; % R^2>0.9 for good trajectory

frame_interval = 0.0095; % unit: s
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
%          a = NaN(n_spots, 1);
%          b = NaN(n_spots, 1);
         r2fit = NaN(n_spots, 1);
%          c = NaN(n_spots, 1); 
         ft2_r2fit = NaN(n_spots, 1); 
%          d = NaN(n_spots, 1); 
%          e = NaN(n_spots, 1); 
         ft3_r2fit = NaN(n_spots, 1); 
%          Tralength = NaN(n_spots, 1); 
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
%             a(i_spot) = fo.p1;
%             b(i_spot) = fo.p2;
            r2fit(i_spot) = gof.adjrsquare;
            % fit2 parameter y=ax
%             c(i_spot) = fo2.a;
            ft2_r2fit(i_spot) = gof2.adjrsquare;
            % fit3 parameter power law log-log fit
%             d(i_spot) = fo3.p1;
%             e(i_spot) = fo3.p2;
            ft3_r2fit(i_spot) = gof3.adjrsquare;
    
            Tralength(i_spot) = size(t, 1); 
          
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
            trackedPar = struct;
         
            for i = 1:n_goodTracks 
                trackedPar(1,i).xy =  tracks{ind_good(i)}(:, 2:3);
                Frame = tracks{ind_good(i)}(:, 1)./frame_interval;
                trackedPar(i).Frame = fix(Frame); % round to integer
                trackedPar(i).TimeStamp = tracks{ind_good(i)}(:, 1);
            end 
            load(sprintf('%s_%03d_Tracked.mat',concentration, video(k)), 'settings');
            save( sprintf('good_%s_%03d_Tracked.mat',concentration, video(k)), 'trackedPar','settings');
         end
         
         %clear tracks
end

