clear; clc; clearvars -global; close all; 

type = 'No Threshold'; % clarify [urea]
concentration = 'buffer'; % file name format
video = [133];
Dmax = 5; % um^2/s
gap = 2;

ha = gca;
hold(ha, 'on');

for k = 1:numel(video)
    
load(sprintf('%s_%03d_Tracked.mat',concentration, video(k)));

%ha = gca;
n_tracks = size(trackedPar,2);
colors = jet(n_tracks);
indices = [1:n_tracks]';
%hold(ha, 'on');
hps = NaN(n_tracks, 1);


for i = 1 : n_tracks
        
        index = indices(i);
        track = trackedPar(index).xy;
        trackName = sprintf('Track %d', index );
        
        x = track(:,1);
        y = track(:,2);
        
        hps(i) =  plot(ha, x, y, ...
            'Color', colors(i,:), ...
            'DisplayName', trackName );
        
end

end
title(sprintf('MTT %s: Dmax=%d um^2/s, gap=%d', type, Dmax, gap));
xlabel('x (um)') 
ylabel('y (um)')
saveas(ha, sprintf('PlotTrajFromMTT_%s_Dmax=%d_gap=%d_All', type, Dmax, gap));
saveas(ha, sprintf('PlotTrajFromMTT_%s_Dmax=%d_gap=%d_All.jpg', type, Dmax, gap));
