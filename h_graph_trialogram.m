function  [pop, trialogram] = h_graph_trialogram(events, spikes, bin_bounds, do_Hz, do_plot, do_smooth, kernelbw)

%  Make a correlogram
%  
%    Args:
%      events (vector)          : Vector of event times (numeric) for events of interest.
%      spikes (vector)          : Vector of time points where spikes occurred (numeric, in seconds).
%      bin_bounds (vector)      : Vector of time bin time points within which the population code will be computed and between which the
%                                 correlations will be computed (relative to events, in seconds).
%      do_Hz (bool)             : Convert data to Hz? Recommended
%      do_plot (bool)           : Create a plot of the trialogram?
%      do_smooth (bool)         : Apply a Gaussian filter to the binned population vectors?
%      kernelbw (scalar)        : Bandwidth of the Gaussian kernel
%     
%    Return:
%      A dataframe with the moment-to-moment population activity and the
%      corresponding trialogram (i.e., correlations between the
%      moment-to-moment population activity /across/ trials, as in Bulkin et al., 
%      as opposed to auto-correlation across time).

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

% Trialogram
trialhold = zeros(n_bins, n_bins, n_events-1);
trialhold2 = zeros(n_bins, n_bins, n_events);
triallist = 1:n_events;
for t1 = 1:n_events
    allothertrials = pop(:,:,setdiff(triallist,t1));
    for b1 = 1:n_bins
        PVb = pop(:,b1,t1);
        for t2 = 1:(n_events-1)
            for b2 = 1:n_bins
                trialhold(b1,b2,t2) = corr(PVb, allothertrials(:,b2,t2));
            end
        end
    end
    trialhold2(:,:,t1) = tanh(mean(atanh(trialhold),3,'omitnan'));
end
trialogram = tanh(mean(atanh(trialhold2),3,'omitnan'));

% Plot trialogram?
if (do_plot)
    figure
    h = heatmap(trialogram,'Colormap',summer);
    h.NodeChildren(3).YDir='normal'; 
end

end


% 
% 
% events=df.spikes{s}.vid_arriv(trials.Mcorr);
% spikes=df.spikes{s}.spike;
% do_Hz=1;
% do_average=1;
% do_plot=1;
% do_smooth=1;
% kernelbw=3;
% 
% 
% % Trialogram
% trialogram = zeros(n_bins, n_bins);
% trialhold = zeros(n_bins, n_bins, n_events-1);
% trialhold2 = zeros(n_bins, n_bins, n_events);
% trialcorr = zeros(n_bins,n_bins);
% triallist = 1:n_events;
% 
% for t1 = 1:n_events
%     allothertrials = pop(:,:,setdiff(triallist,t1));
%     for b1 = 1:n_bins
%         PVb = pop(:,b1,t1);
%         for t2 = 1:(n_events-1)
%             for b2 = 1:n_bins
%                 trialhold(b1,b2,t2) = corr(PVb, allothertrials(:,b2,t2));
%             end
%         end
%     end
%     trialhold2(:,:,t1) = tanh(mean(atanh(trialhold),3,'omitnan'));
% end
% 
% 
% figure
% a=tanh(mean(atanh(trialhold2),3,'omitnan'));
% h = heatmap(a,'Colormap',summer);
% h.NodeChildren(3).YDir='normal';
% 
% 
% 
% for t1 = 1:n_events
%     for b1 = 1:n_bins
%         PVb = pop(:,b1,t1);
%         curr_triallist = setdiff(triallist,t1);
%         allothertrials = pop(:,:,curr_triallist);
%         for t2 = 1:(n_events-1)
%             for b2 = 1:n_bins
%                 trialhold(b1,b2,t2) = corr(PVb, allothertrials(:,b2,t2));
%                 trialhold(b2,b1,t2) = trialhold(b1,b2,t2);
%             end
%         end
%         %trialhold2(t1,:) = tanh(mean(atanh(trialhold),'omitnan'));
%     end
% end
