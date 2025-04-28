function [ft_df] = local_combine_sess(trialified_lfp,roi_label,Fs)
    
    % Initialize
    ft_df = {};

    % Get the trial epochs for all conditions in each session
    count=0;
    for s = 1:length(trialified_lfp) % ..for each session
        for t = 1:length(trialified_lfp{s}.trial) % ..for each trial in that session
            
            % Update counter
            count=count+1;

            % Transfer to the final data frame
            ft_df.trial{1,count} = trialified_lfp{s}.trial{t};
            ft_df.time{1,count} = trialified_lfp{s}.time{t};
        end
    end

    % Add additional meta info that FT wants
    ft_df.label = {roi_label};
    ft_df.fsample = Fs;

end