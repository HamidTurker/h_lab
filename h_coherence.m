function [coherence, ratioCoh, test_out] = h_coherence(ft_trialifiedX, ft_trialifiedY, tvec_in, tvec_out, frequencies, samplerate, wavelet_cycles, baseline, subsample, cfg_test)

%ft_trialifiedX=coh.PFC.match_corr;
%ft_trialifiedY=coh.HPC.match_corr;
%tvec_in=tvec_unbuffed;
%frequencies=param.FOIs;
%samplerate=param.Fs;
%wavelet_cycles=param.wavelet_cycles;
%baseline=param.baseline;
%subsample = 100;

    % Compute spectral coherence
    fprintf('\n\tStarting spectral coherence analysis.. \n')

    % Subsampling the number trials available?
    if ~isempty(subsample)

        % Subsample of trial IDs
        n_trials = length(ft_trialifiedX.trial);
        vec_trials = randperm(length(1:n_trials));
        chosen = vec_trials(1:subsample)';

        % Talk
        fprintf('\n\tSubsampling %d out of %d total number of trials available!',subsample,n_trials)
        
        % Select the sampled trials
        hold = {}; % X data
        for t = 1:length(chosen)
            hold.trial{t} = ft_trialifiedX.trial{chosen(t)};
            hold.time{t} = ft_trialifiedX.time{chosen(t)};
        end
        % Replace in the main data frame
        ft_trialifiedX.trial = hold.trial;
        ft_trialifiedX.time = hold.time;

        hold = {}; % Y data
        for t = 1:length(chosen)
            hold.trial{t} = ft_trialifiedY.trial{chosen(t)};
            hold.time{t} = ft_trialifiedY.time{chosen(t)};
        end
        % Replace in the main data frame
        ft_trialifiedY.trial = hold.trial;
        ft_trialifiedY.time = hold.time;
        
    else
        fprintf('\tNo subsampling..\n')
    end

    % Params
    n_trials = length(ft_trialifiedX.trial);
    n_timepoints_trial = length(tvec_in);
    n_wavelet = n_timepoints_trial;
    n_timepoints_alldata = n_timepoints_trial * n_trials;
    n_frequencies = length(frequencies);
    n_convolution = n_wavelet+n_timepoints_alldata-1;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    wavelet_time = -n_timepoints_trial/samplerate/2 : 1/samplerate : n_timepoints_trial/samplerate/2-1/samplerate;
    gausswindow = wavelet_cycles ./ (2*pi*frequencies);
    timepoints_of_interest = dsearchn(tvec_in',tvec_out');
    cycle_width = linspace(1.5,3,n_frequencies); % number of cycles on either end of the center point (1.5 means a total of 3 cycles))

    % Initialize
    coherence = zeros(n_frequencies,length(timepoints_of_interest));
    ratioCoh = [];
    test_out = [];
    
    % Vectorize data (put all trials on a single row) and take the fft()
    fprintf('\n\n\tVectorizing the data and performing the Fourier transform..\n')
    fft_sigX = []; fft_sigY = [];
    for t = 1:n_trials
        fft_sigX = horzcat(fft_sigX,ft_trialifiedX.trial{t});
        fft_sigY = horzcat(fft_sigY,ft_trialifiedY.trial{t});
    end
    fft_sigX = fft(fft_sigX,n_conv_pow2);
    fft_sigY = fft(fft_sigY,n_conv_pow2);

    % Compute spectral coherence
    for fi = 1:n_frequencies
        fprintf('\t- Computing spectral coherence for frequency %d of %d\n', fi, n_frequencies)

        % Wavelet
        wavelet = exp(2*1i*pi*frequencies(fi).*wavelet_time) .* exp(-wavelet_time.^2./(2*(gausswindow(fi)^2)))/frequencies(fi);

        % Convolution
        conv_sigX = ifft(fft(wavelet,n_conv_pow2) .* fft_sigX);
        conv_sigX = conv_sigX(1:n_convolution);
        conv_sigX = reshape(conv_sigX(floor((n_timepoints_trial-1)/2):end-1-ceil((n_timepoints_trial-1)/2)),n_timepoints_trial,n_trials);

        conv_sigY = ifft(fft(wavelet,n_conv_pow2) .* fft_sigY);
        conv_sigY = conv_sigY(1:n_convolution);
        conv_sigY = reshape(conv_sigY(floor((n_timepoints_trial-1)/2):end-1-ceil((n_timepoints_trial-1)/2)),n_timepoints_trial,n_trials);

        % Computer power & cross-spectral power
        specXX = mean(conv_sigX.*conj(conv_sigX),2);
        specYY = mean(conv_sigY.*conj(conv_sigY),2);
        specXY = abs(mean(conv_sigX.*conj(conv_sigY),2)).^2;

        % Spectral coherence
        coherence(fi,:) = specXY(timepoints_of_interest) ./ (specXX(timepoints_of_interest).*specYY(timepoints_of_interest));
    end

    % Baseline correction?
    if ~isempty(baseline)
        fprintf('\n\n\tPerforming baseline correction..\n')
        ratioCoh = zeros(n_frequencies, length(timepoints_of_interest));
        baseline_points = dsearchn(tvec_out',baseline'); % Indices of the baseline start and end

        for f=1:n_frequencies
            bl_avg  = mean(coherence(f,baseline_points(1):baseline_points(2)));
            ratioCoh(f,:) = coherence(f,:) / bl_avg;
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
                    bl_avg  = mean(coherence(f,bl_points(1):bl_points(2)));
                    shuff(f,:,i) = coherence(f,:) / bl_avg;
                end

            end

            % Permutation testing out
            test_out = mean(shuff,3);

        end
    end

    % Finished!
    fprintf('\n\n\tFinished!\n')

end