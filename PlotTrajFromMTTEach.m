concentration = 'NOgly_no ND_h76_1to1000_4ms_buffer';
video = [34];

for k = 1:numel(video)
    
load(sprintf('%s_%03d_Tracked.mat',concentration, video(k)));

n_tracks = size(trackedPar,2);
colors = jet(n_tracks);
indices = [1:n_tracks]';
hps = NaN(n_tracks, 1);


for i = 1 : n_tracks
        
        index = indices(i);
        track = trackedPar(index).xy;
        trackName = sprintf('Track %d', index );
        
        x = track(:,1);
        y = track(:,2);
        
        hps(i) =  plot(x, y, ...
            'Color', colors(i,:), ...
            'DisplayName', trackName );
        hold on
        
end

hold off
saveas(gcf, sprintf('PlotTrajFromMTT_EACH_%s_%03d', concentration,video(k)));

end


