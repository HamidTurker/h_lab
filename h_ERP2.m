function [ERPs, ERPs_ft, tvec] = h_ERP2(epoched_lfps, df_task, subIDs, conditions, normalization, baseline)

    % Initialize
    ERPs = {}; % Basic structure
    ERPs_ft = {}; % Structure in Fieldtrip format

    % Global tvec
    tvec = epoched_lfps{1}.time{1};

    % Baseline bins
    [~, idx_start] = min(abs(tvec-baseline(1))); % Start
    [~, idx_end] = min(abs(tvec-baseline(2))); % End

    % Create the ERPs for each requested condition
    for c = conditions

        % Set the fields
        ERPs = setfield(ERPs,c,{});
        ERPs_ft = setfield(ERPs_ft,c,{});
        ERPs.(c).tvec = tvec;

        % Subset the trialified data to include only those trials we're
        % currently interested in (trials of interest)
        for s = 1:length(subIDs) % ..for each subject's session

            % Trials of Interest
            toi = getfield(df_task{s},c);

            chosen_trials = epoched_lfps{s}.trial(:,toi);
            chosen_times = epoched_lfps{s}.time(:,toi);

            for t = 1:sum(toi)
                % Baseline
                base_avg = mean(chosen_trials{t}(idx_start:idx_end));

                % Basic formatting
                ERPs.(c).sess_trials{s}(t,:) = chosen_trials{t} - base_avg;

                % Fieldtrip formatting
                ERPs_ft.(c){s}.trial{t} = chosen_trials{t} - base_avg;
                ERPs_ft.(c){s}.time{t} = chosen_times{t};
            end

            % Session average
            if (normalization == "z-score")
                ERPs.(c).sess_avg(s,:) = zscore(mean(ERPs.(c).sess_trials{s}));
            elseif (normalization == "range")
                ERPs.(c).sess_avg(s,:) = normalize(mean(ERPs.(c).sess_trials{s}), "range");
            else
                ERPs.(c).sess_avg(s,:) = mean(ERPs.(c).sess_trials{s});
            end

        end

        % Group-level average (unweighted)
        ERPs.(c).unweighted_grp_avg = mean(ERPs.(c).sess_avg);

        % Group-level average (weighted by subject)
        unique_subs = unique(subIDs);
        for u = 1:length(unique_subs)
            sub_sess = subIDs == unique_subs(u);
            ERPs.(c).sub_avg(u,:) = mean(ERPs.(c).sess_avg(sub_sess,:));
        end
        ERPs.(c).subweighted_grp_avg = mean(ERPs.(c).sub_avg);

    end

end





 

