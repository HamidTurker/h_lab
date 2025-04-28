function [ITPC, ratioITPC, stats_out, test_out] = h_ITPC(ft_trialified, tvec_in, tvec_out, frequencies, samplerate, wavelet_cycles, baseline, subsample, save_stats, cfg_test)

    % Compute Inter-Trial Phase Clustering
    fprintf('\n\tStarting Inter-Trial Phase Clustering analysis.. \n')

    % Subsampling the number of trials to match another condition?
    if ~isempty(subsample)
        % Subsample of trial IDs
        n_trials = length(ft_trialified.trial);
        chosen = randsample(n_trials,subsample)';
        
        % Talk
        fprintf('\n\tSubsampling %d out of %d total number of trials available!',subsample,n_trials)

        % Select the sampled trials
        hold = {};
        for t = 1:length(chosen)
            hold.trial{t} = ft_trialified.trial{chosen(t)};
            hold.time{t} = ft_trialified.time{chosen(t)};
        end

        % Replace in the main data frame
        ft_trialified.trial = hold.trial;
        ft_trialified.time = hold.time;

        % Save_stats aren't provided if we're subsampling
        if ~isempty(save_stats)
            save_stats = [];
        end
    end

    % Params
    n_trials = length(ft_trialified.trial);
    n_timepoints_trial = length(tvec_in);
    n_wavelet = n_timepoints_trial;
    n_timepoints_alldata = n_timepoints_trial * n_trials;
    n_frequencies = length(frequencies);
    n_convolution = n_wavelet+n_timepoints_alldata-1;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    wavelet_time = -n_timepoints_trial/samplerate/2 : 1/samplerate : n_timepoints_trial/samplerate/2-1/samplerate;
    gausswindow = wavelet_cycles ./ (2*pi*frequencies);
    timepoints_of_interest = dsearchn(tvec_in',tvec_out');

    % Initialize
    ITPC = zeros(n_frequencies,length(timepoints_of_interest));

    % Vectorize data (put all trials on a single row) and take the fft()
    fprintf('\n\n\tVectorizing the data and performing the Fourier transform..\n')
    fft_data = [];
    for t = 1:n_trials
        fft_data = horzcat(fft_data,ft_trialified.trial{t});
    end
    fft_data = fft(fft_data,n_conv_pow2);

    % Are we saving any stats? Initialize stats_out
    if ~isempty(save_stats)

        % Saving at the session level
        n_stat_windows = length(save_stats.window(:,1));
        n_stat_sessions = length(save_stats.size_sess(:,1));
        stats_out = {};
        stats_out.raw = NaN(n_stat_sessions, n_frequencies, n_stat_windows);

    else
        stats_out = [];
    end

    % Compute the ITPC
    for fi = 1:n_frequencies
        fprintf('\t- Computing ITPC for frequency %d of %d\n', fi, n_frequencies)

        % Wavelet
        wavelet = exp(2*1i*pi*frequencies(fi).*wavelet_time) .* exp(-wavelet_time.^2./(2*(gausswindow(fi)^2)))/frequencies(fi);

        % Convolution
        conv_data = ifft(fft(wavelet,n_conv_pow2) .* fft_data);
        conv_data = conv_data(1:n_convolution);
        conv_data = reshape(conv_data(floor((n_timepoints_trial-1)/2):end-1-ceil((n_timepoints_trial-1)/2)),n_timepoints_trial,n_trials);
        % Reshape into a complex double (row = time point, column = trial, data = wavelet-convolved data at the current frequency)

        % Extract ITPC
        hold = abs(mean(exp(1i*angle(conv_data)),2));
        ITPC(fi,:) = hold(timepoints_of_interest);

        % Stats_out, if requested
        if ~isempty(save_stats)

            % Indices of the which trials to select for a given session
            t_idx = NaN(n_stat_sessions,2);
            t_idx(1,1) = 1;
            count = 0;
            for s = 1:n_stat_sessions
                count = count+save_stats.size_sess(s);
                t_idx(s,2) = count;
            end
            for s = 2:n_stat_sessions
                t_idx(s,1) = t_idx(s-1,2)+1;
            end

            % For each requested stats window
            for w = 1:n_stat_windows
                win_points = dsearchn(tvec_in',save_stats.window(w,:)'); % Indices of the stats window start and end

                for s=1:n_stat_sessions
                    stats_out.raw(s,fi,w) = mean(abs(mean(exp(1i*angle(conv_data(win_points(1):win_points(2),t_idx(s,1):t_idx(s,2)))),2)));
                end
            end
        end

    end

    % Baseline correction?
    if ~isempty(baseline)
        fprintf('\n\n\tPerforming baseline correction..\n')
        ratioITPC = zeros(n_frequencies, length(timepoints_of_interest));
        bl_points = dsearchn(tvec_out',baseline'); % Indices of the baseline start and end

        for f=1:n_frequencies
            bl_avg  = mean(ITPC(f,bl_points(1):bl_points(2)));
            ratioITPC(f,:) = ITPC(f,:) / bl_avg;
        end
    end

    % Non-parametric permutation testing?
    if ~isempty(cfg_test)
        fprintf('\n\n\tCreating a shuffled distribution..\n')

        % Permutation testing by computing the ratio with shuffled baseline periods
        if any(cfg_test.method == 'ratio')
            
            shuff = zeros(n_frequencies, length(timepoints_of_interest), cfg_test.iters); % Shuffling data frame
            bl_dur = abs(diff(baseline)); % What's the total duration of the baseline? We'll match its length
            bl_starts = tvec_out(tvec_out <= tvec_out(end)-bl_dur); % All possible random baseline starting points
            
            for i = 1:cfg_test.iters
                start = bl_starts(randi(length(bl_starts))); bl_shuff = [start start+bl_dur];
                bl_points = dsearchn(tvec_out',bl_shuff'); % Random baseline within trial window, same duration as original baseline

                for f=1:n_frequencies
                    bl_avg  = mean(ITPC(f,bl_points(1):bl_points(2)));
                    shuff(f,:,i) = ITPC(f,:) / bl_avg;
                end

            end

            % Permutation testing out
            test_out = mean(shuff,3);

        end
    end

    % Finished!
    fprintf('\n\n\tFinished!\n')
end

