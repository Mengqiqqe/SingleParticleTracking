function obj = fitMSD(obj, clip_factor)
%%FITMSD Fit all MSD curves by a linear function.
%
% obj = obj.fitMSD fits all MSD curves by a straight line
%                      y = a * x + b.
% The fit is therefore rigorously valid only for purely
% diffusive behavior.
%
% Results are stored in the 'fit' field of the returned
% object. It is a structure with 2 fields:
% - a: all the values of the slope of the linear fit.
% - b: all the values for the intersect of the linear fit.
% - r2fit: the adjusted R2 value as a indicator of the goodness
% of the fit.
%
% obj = obj.fitMSD(clip_factor) does the fit, taking into
% account only the first potion of the average MSD curve
% specified by 'clip_factor' (a double between 0 and 1). If the
% value exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.


if nargin < 2
    clip_factor = 0.25;
end

if ~obj.msd_valid
    obj = obj.computeMSD;
end
n_spots = numel(obj.msd);

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
c = NaN(n_spots, 1);
ft2_r2fit = NaN(n_spots, 1);
Tralength = NaN(n_spots, 1); % Mengqi added
ft1 = fittype('poly1');
ft2 = fittype({'x'});

fprintf('%4d/%4d', 0, n_spots);
for i_spot = 1 : n_spots
    
    fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i_spot, n_spots);
    
    msd_spot = obj.msd{i_spot};
    
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
    
%     % Thrash bad data %% Mengqi shift above!!!
%     nonnan = ~isnan(y);
%     x = t(nonnan);
%     y = y(nonnan);
%     w = w(nonnan);
    
    if numel(Y) < 2
        continue
    end
    
    [fo, gof] = fit(x, Y, ft1, 'Weights', W);
    [fo2, gof2] = fit(x, Y, ft2, 'Weights', W);
    % fit1 parameter
    a(i_spot) = fo.p1;
    b(i_spot) = fo.p2;
    r2fit(i_spot) = gof.adjrsquare;
    % fit2 parameter
    c(i_spot) = fo2.a;
    ft2_r2fit(i_spot) = gof2.adjrsquare;
    
    Tralength(i_spot) = size(x, 1); % Mengqi added
    
    % Mengqi addded for plotting and saving data
    if gof.adjrsquare > 0.9
         if gof2.adjrsquare > 0.9
         %save MSD data
         MSD_ind_data = [t,y,stder,w];
         T_MSD_ind_data = array2table(MSD_ind_data, 'VariableNames',{'t_s', 'msd_um2', 'std','N'});
         writetable( T_MSD_ind_data, sprintf('MSDdata_%d.xlsx', i_spot)); 
        
         % save MSD plot
         h = figure;
         plot(fo,'r',x,Y);
         hold on 
         plot(fo2,'g',x,Y);
         plot (t,y)
         hold off
         saveas(h, sprintf('FIG_%d.png', i_spot));
         end
    end
end
fprintf('\b\b\b\b\b\b\b\b\bDone.\n')

obj.lfit = struct(...
    'a', a, ...
    'b', b, ...
    'r2fit', r2fit,...
    'c', c,...
    'ft2_r2fit', ft2_r2fit,...
    'Tralength', Tralength); % Mengqi added

end