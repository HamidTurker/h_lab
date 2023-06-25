function  raster_frame = h_graph_raster(events, spikes, per, group, pre, post, bin, do_count, do_Hz, do_plot)

%  Make a raster around events (events, in seconds into global time),
%  with a pre-event and post-event window, with a given bin size.
%  
%    Args:
%      events (vector)         : Vector of event times (numeric) for events of interest.
%      spikes (vector)         : Vector of time points where spikes occurred (numeric, in seconds).
%      per (vector)            : Plot the raster 'per' what (on the y-axis)? Typically,
%                                there is one row 'per' trial, in which case you should provide
%                                a vector of trial ID numbers of each event time.
%      group (vector)          : Vector where each entry notes the group label of the corresponding entry in the event vector.
%      pre (num)               : Pre-event window (in seconds).
%      post (num)              : Post-event window (in seconds).
%      plot (bool)             : Make a plot (TRUE) or just return the final data frame (FALSE)?
%      
%     
%    Return:
%      A raster plot or data frame.


%% Initialize
pre = -(abs(pre));
n_events = length(events);
bin_seq = pre:bin:post;
raster_frame = zeros(length(per), length(bin_seq));

%% Mark spikes in the raster surrounding each event
for i = 1:n_events
    
    % Which per entry?
    per_entry = per(i);

    % Find the trial section
    starttime = events(i)+pre; endtime = events(i)+post+bin;
    thistrial = spikes(spikes >= starttime & spikes <= endtime);

    % If there are spikes on this trial
    if ~isempty(thistrial)

        % Spike times in relation to the start of this trial
        time_corrected_spikes = thistrial - starttime;

        % Modulate out bin size and correct time for pre-window size
        % This gives us the bin entries (i.e. bins that need to be marked)
        bin_times = time_corrected_spikes-mod(time_corrected_spikes,bin);
        bin_entries = round((1/bin) * bin_times + 1);

        % Mark the bin_entries in the raster_frame
        if (do_count)
            for j = 1:length(bin_entries)
                raster_frame(per_entry,bin_entries(j)) = raster_frame(per_entry,bin_entries(j)) + 1;
            end
        else
            raster_frame(per_entry,bin_entries) = 1;
        end

    end
end

% Convert counts to Hz? Only possible if we did counts in the previous
% step. Otherwise, the results won't make sense.
if (do_count && do_Hz)
    raster_frame = raster_frame / bin;
end
