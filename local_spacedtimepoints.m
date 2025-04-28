function [df_timepoints, tvec_leftedge2firstevent] = local_spacedtimepoints(df_spikes, perc, edges, df_task, condition)
    
    % Initialize outs
    df_timepoints = {};
    tvec_leftedge2firstevent = linspace(edges(1), 0, perc);
    
    % For each trial, get the EoI points and linearly space out time points
    hold = {};
    for i = 1:length(df_task)
        for t = 1:length(df_spikes{i}.vid_arriv)
    
            % Do we have digging on this trial?
            if ~isnan(df_spikes{i}.vid_digend(t)) % We have a dig
    
                digtime = df_spikes{i}.vid_digend(t);
                leavetime = digtime+df_spikes{i}.vid_leave(t);

                bin_bounds = horzcat([linspace(edges(1),  0,                    perc)], ...   % e.g. 5 linear seconds prior to arrival
                                     [linspace(0,         digtime,              perc+1)], ... % Variable period to dig_end
                                     [linspace(digtime,   leavetime,            perc+1)], ... % Variable period to leaving
                                     [linspace(leavetime, leavetime+edges(2),   perc+1)]);    % e.g. 5 linear seconds after leaving
                bin_bounds = unique(bin_bounds); % The previous line creates some duplicates, which we remove here
    
            else                         % No digging
    
                leavetime = df_spikes{i}.vid_leave(t);
                
                bin_bounds = horzcat([linspace(edges(1),  0,                    perc)], ...   % 5 seconds prior to arrival
                                     [linspace(0,         leavetime,            perc+1)], ... % Variable period to leaving
                                     [linspace(leavetime, leavetime+edges(2),   perc+1)]);    % 5 seconds after leaving
                bin_bounds = unique(bin_bounds); % The previous line creates some duplicates, which we remove here

            end

            % Save spacings
            hold{i,1}{t,1} = bin_bounds;

        end
    end

    % Subset only to those trials part of the requested condition
    tcount = 0;
    for i = 1:length(df_task)
        toi = logical(getfield(df_task{i},condition));
        for t = 1:length(toi)
            if toi(t)
                tcount=tcount+1;
                df_timepoints{tcount,1} = hold{i}{t};
            end
        end
    end

    return
