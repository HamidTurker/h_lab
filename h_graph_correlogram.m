function  [pop, correlogram] = h_graph_correlogram(events, spikes, bin_bounds, do_Hz, do_average, do_plot, do_smooth, kernelbw)

%  Make a correlogram
%  
%    Args:
%      events (vector)          : Vector of event times (numeric) for events of interest.
%      spikes (vector)          : Vector of time points where spikes occurred (numeric, in seconds).
%      bin_bounds (vector)      : Vector of time bin time points within which the population code will be computed and between which the
%                                 correlations will be computed (relative to events, in seconds).
%     
%    Return:
%      A dataframe with the moment-to-moment population activity and the
%      corresponding correlogram (i.e., correlations between the
%      moment-to-moment population activity).

%% Initialize
n_events = length(events);
n_cells = length(spikes);
n_bins = length(bin_bounds)-1;
pop = zeros(n_cells, n_bins, n_events); % Final frame for population acitivity
eventframe = events+bin_bounds; % All the bin boundaries, surrounding each of the events-of-interest

%% Mark spikes in the raster surrounding each event
for i = 1:n_events
    for c = 1:n_cells

        % Find all spikes between first and last bin bound
        thistrial = spikes{c}(spikes{c} >= eventframe(i,1) & spikes{c} <= eventframe(i,end));

        % If there are spikes on this trial, sort the spikes.
        % Otherwise, just skip to the next trial.
        if ~isempty(thistrial)

            % Sort the spikes into the bins
            for b = 1:n_bins
                pop(c,b,i) = length(spikes{c}(spikes{c} >= eventframe(i,b) & spikes{c} <= eventframe(i,b+1)));
            end
        end
    end
end

% Smooth the population code?
if (do_smooth)
    for i = 1:n_events
        for c = 1:n_cells
            pop(c,:,i) = filter(gausswin(kernelbw),1,pop(c,:,i));
        end
    end
end

% Convert spike counts to Hz?
if (do_Hz)
    for i = 1:n_events
        for b = 1:n_bins
            timediff = eventframe(i,b+1)-eventframe(i,b);
            pop(:,b,i) = pop(:,b,i) / timediff;
        end
    end
end

% Correlogram
correlogram = zeros(n_bins, n_bins, n_events);
for i = 1:n_events
    for b1 = 1:n_bins
        for b2 = 1:n_bins
            correlogram(b1,b2,i) = corr(pop(:,b1,i),pop(:,b2,i));
        end
    end
end

% Average?
if (do_average)
    correlogram = tanh(mean(atanh(correlogram),3,'omitnan'));
end

% Plot correlogram?
if (do_plot)
    figure
    h = heatmap(correlogram,'Colormap',summer);
    h.NodeChildren(3).YDir='normal'; 
end

end