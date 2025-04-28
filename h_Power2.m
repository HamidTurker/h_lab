function [Power, ratioPower, stats_out, test_out] = h_Power2(ft_trialified, tvec_in, tvec_leftedge2firstevent, frequencies, samplerate, wavelet_cycles, baseline, df_timepoints, mintimes, subsample, save_stats, cfg_test)

    % Compute power at linearly spaced time points that can be unique for each trial
    fprintf('\n\tStarting Power analysis!\n')


    % Subsampling the number trials available?
    if ~isempty(subsample)

        % Subsample of trial IDs
        n_trials = length(ft_trialified.trial);
        vec_trials = randperm(length(1:n_trials));
        chosen = vec_trials(1:subsample)';
        
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

        % Also subsample the df_timepoints
        hold={};
        for t = 1:length(chosen)
            hold{t} = df_timepoints{chosen(t)};
        end
        % Replace in the main data frame
        df_timepoints = hold';
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

    % Initialize outs
    Power = zeros(n_frequencies,mintimes);
    ratioPower = zeros(n_frequencies,mintimes);
    stats_out = [];
    test_out = [];

    % Vectorize data (put all trials on a single row) and take the fft()
    fprintf('\n\n\tVectorizing the data and performing the Fourier transform..\n')
    fft_data = [];
    for t = 1:n_trials
        fft_data = horzcat(fft_data,ft_trialified.trial{t});
    end
    fft_data = fft(fft_data,n_conv_pow2);

    % % Compute the Power
    for fi = 1:n_frequencies
        fprintf('\t- Computing Power for frequency %d of %d\n', fi, n_frequencies)

        % Wavelet
        wavelet = exp(2*1i*pi*frequencies(fi).*wavelet_time) .* exp(-wavelet_time.^2./(2*(gausswindow(fi)^2)))/frequencies(fi);

        % Convolution
        conv_data = ifft(fft(wavelet,n_conv_pow2) .* fft_data);
        conv_data = conv_data(1:n_convolution);
        conv_data = reshape(conv_data(floor((n_timepoints_trial-1)/2):end-1-ceil((n_timepoints_trial-1)/2)),n_timepoints_trial,n_trials);
        % Reshape into a complex double (row = time point, column = trial, data = wavelet-convolved data at the current frequency)

        % Extract power for subsampled time points
        spaced_data = []; tcount = 0; error_count = 0;
        for ti = 1:n_trials
            if length(df_timepoints{ti}) == mintimes
                tcount=tcount+1;
                timepoints_of_interest = dsearchn(tvec_in',df_timepoints{ti}');
                spaced_data(tcount,:) = conv_data(timepoints_of_interest,ti);
            else
                error_count = error_count+1;
                fprintf('Warning: The current df.timepoints entry ( trial %d ) did not have the same number of time points as mintimes! This has now happened %d time(s)..\n', ti, error_count)
            end
        end
        Power(fi,:) = mean(abs(spaced_data).^2);
    end

    % Baseline correction?
    if ~isempty(baseline)
        fprintf('\n\n\tPerforming baseline correction..\n')
        ratioPower = zeros(n_frequencies, mintimes);
        bl_points = dsearchn(tvec_leftedge2firstevent',baseline'); % Indices of the baseline start and end

        for fi=1:n_frequencies
            bl_avg  = mean(Power(fi,bl_points(1):bl_points(2)));
            ratioPower(fi,:) = Power(fi,:) / bl_avg;
        end
    end

    % Non-parametric permutation testing?
    if ~isempty(cfg_test)
        fprintf('\n\n\tCreating a shuffled distribution..\n')

        % Permutation testing by computing the ratio with shuffled baseline periods
        if any(cfg_test.method == 'ratio')
            
            shuff = zeros(n_frequencies, length(timepoints_of_interest), cfg_test.iters); % Shuffling data frame
            bl_dur = abs(diff(baseline)); % What's the total duration of the baseline? We'll match its length
            bl_starts = tvec_leftedge2firstevent(tvec_leftedge2firstevent <= tvec_leftedge2firstevent(end)-bl_dur); % All possible random baseline starting points
            
            for i = 1:cfg_test.iters
                start = bl_starts(randi(length(bl_starts))); bl_shuff = [start start+bl_dur];
                bl_points = dsearchn(tvec_leftedge2firstevent',bl_shuff'); % Random baseline within trial window, same duration as original baseline

                for f=1:n_frequencies
                    bl_avg  = mean(Power(f,bl_points(1):bl_points(2)));
                    shuff(f,:,i) = Power(f,:) / bl_avg;
                end

            end

            % Permutation testing out
            test_out = mean(shuff,3);

        end
    end

    % Finished!
    fprintf('\n\n\tFinished!\n')
end

