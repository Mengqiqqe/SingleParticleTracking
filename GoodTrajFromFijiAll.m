clear; clc; clearvars -global; close all; 


type = '1mM'; % clarify [urea]
concentration = '1mM'; % file name format
video = [22:43];
remove = [];
FittingType = 1; % unbrownian: FittingType=1; brownian: FittingType=0
Mini_Trajlength = 15; % define the mini trajectory length to rule out those short trajectories
R_threshold = 0.9; % R^2>0.9 for good trajectory
frame_interval = 0.08; % unit: s
clip_factor =0.25;  % usually set as 0.25, only fit the first one quarter of the MSD plot
R_thresholdS = sprintf('%.1f', R_threshold);

for j = 1:numel(remove)
    video(video == remove(j)) = [];
end

for k = 1:numel(video)
    
    datasheet = csvread(sprintf('%s_%d.csv',concentration, video(k)),1,0); % File name:'%s_30ms_%d.csv';
    Results = array2table(datasheet, 'VariableNames',{'Trajectory','Frame','x','y','z','mm','mmm','mmmm','nnn','mmmmmmm','NPscore'}); % Transfer .csv array to table;
    Traj_id = unique(Results.Trajectory(:)); % Results is your table name. finds the Traj'ID = 1,2,3...
    num_Traj = numel(Traj_id); %this gives you the total number of Trajectories;
    tracks = {};

    for i=1:num_Traj

        ID = Traj_id(i); % 'ID' is the ID for each Traj;
        eachTraj_data = find(Results.Trajectory == ID); % 'eachTraj_data': collect data points with the same ID indices; 
    
        if numel(eachTraj_data) > Mini_Trajlength % 'numel()': counts # data points in each Traj, which shows the length of each Traj;
        %if numel(ind_par) < 150
        tracks{ID} = [Results.Frame(eachTraj_data).*frame_interval Results.x(eachTraj_data).*(1/9.375) Results.y(eachTraj_data).*(1/9.375)]; % 'tracks{ID}': data for each Traj after unit conversion;
        %end
        end
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
    
        % Clip data, never take the first one dt = 0
        if clip_factor < 1
            t_limit = 2 : round(numel(t) * clip_factor);
        else
            t_limit = 2 : min(1+round(clip_factor), numel(t));
        end
        x = t(t_limit);
        Y = y(t_limit);
        W = w(t_limit);
        
        if numel(Y) < 2
            continue
        end
    
        [fo, gof] = fit(x, Y, ft1, 'Weights', W);
        [fo2, gof2] = fit(x, Y, ft2, 'Weights', W);
        [fo3, gof3] = fit(log(x), log(Y), ft1, 'Weights', W);
        % fit1 parameter
        a(i_spot) = fo.p1;
        b(i_spot) = fo.p2;
        r2fit(i_spot) = gof.adjrsquare;
        % fit2 parameter
        c(i_spot) = fo2.a;
        ft2_r2fit(i_spot) = gof2.adjrsquare;
        % fit3 parameter
        d(i_spot) = fo3.p1;
        e(i_spot) = fo3.p2;
        ft3_r2fit(i_spot) = gof3.adjrsquare;
    
        Tralength(i_spot) = size(t, 1); % Mengqi added
    
        % Mengqi addded for plotting and saving data
%         if gof.adjrsquare > 0.9
%             if gof2.adjrsquare > 0.9
%             %save MSD data
%             MSD_ind_data = [t,y,stder,w];
%             T_MSD_ind_data = array2table(MSD_ind_data, 'VariableNames',{'t_s', 'msd_um2', 'std','N'});
%             writetable(T_MSD_ind_data, sprintf('%s_%d_MSDdata_Traj%d.xlsx', concentration, video(k), i_spot)); 
   
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
    
    n_tracks = numel(ind_good);
    
    if n_tracks > 0
        
        for i = 1:n_tracks
           goodtrack{i} = tracks{ind_good(i)};
        end
        
        save(sprintf('%s_%03d_goodTracked_Traj>%d_R2>%s.mat', concentration, video(k), Mini_Trajlength, R_thresholdS), 'tracks', 'ind_good', 'goodtrack'); % Mengqi
        
    end
    
    clearvars goodtrack tracks
end

