function [epochs] = local_trialify(lfps, tvecs, events, subIDs, window, buffer, interpolation, Fs)
    
    % Initialize
    epochs = {};

    % Get the trial epochs for all conditions in each session
    for s = 1:length(subIDs) % ..for each subject's session

        % Trial-specific tvecs (from the start of the trial window to the
        % end of the trial window, including buffer time)
        for t = 1:length(events{s}.arriv) % ..for each trial
            trialtvecs(t,:) = events{s}.arriv(t)+(window(1)-buffer):1/Fs:events{s}.arriv(t)+(window(2)+buffer);
        end

        % Epoched, interpolated, buffered trial-level LFPs
        trialdata = interp1(tvecs{s},lfps{s}.data,trialtvecs,interpolation);

        % Create a data structure as expected by Fieldtrip
        for t = 1:length(events{s}.arriv) % ..for each trial
            epochs{s}.trial{t} = trialdata(t,:);
            epochs{s}.time{t} = (window(1)-buffer):1/Fs:(window(2)+buffer);
        end
    end


end