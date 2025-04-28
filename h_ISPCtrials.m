function [ISPC, ratioISPC, ps] = h_ISPCtrials(ft_trialifiedX, ft_trialifiedY, tvec_in, tvec_out, frequencies, samplerate, wavelet_cycles, baseline)

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
    ISPC = zeros(n_frequencies,length(timepoints_of_interest));
    PS = zeros(n_frequencies,length(timepoints_of_interest));
    ps = zeros(n_frequencies,length(timepoints_of_interest));

    % Vectorize data (put all trials on a single row) and take the fft()
    fprintf('\n\n\tVectorizing the data and performing the Fourier transform..\n')
    fft_sigX = []; fft_sigY = [];
    for t = 1:n_trials
        fft_sigX = horzcat(fft_sigX,ft_trialifiedX.trial{t});
        fft_sigY = horzcat(fft_sigY,ft_trialifiedY.trial{t});
    end
    fft_sigX = fft(fft_sigX,n_conv_pow2);
    fft_sigY = fft(fft_sigY,n_conv_pow2);

    % Compute the ISPC
    for fi = 1:n_frequencies
        fprintf('\t- Computing ISPC-trial for frequency %d of %d\n', fi, n_frequencies)

        % Wavelet
        wavelet = exp(2*1i*pi*frequencies(fi).*wavelet_time) .* exp(-wavelet_time.^2./(2*(gausswindow(fi)^2)))/frequencies(fi);

        % Convolution
        conv_sigX = ifft(fft(wavelet,n_conv_pow2) .* fft_sigX);
        conv_sigX = conv_sigX(1:n_convolution);
        conv_sigX = reshape(conv_sigX(floor((n_timepoints_trial-1)/2):end-1-ceil((n_timepoints_trial-1)/2)),n_timepoints_trial,n_trials);
        phase_sigX = angle(conv_sigX);

        conv_sigY = ifft(fft(wavelet,n_conv_pow2) .* fft_sigY);
        conv_sigY = conv_sigY(1:n_convolution);
        conv_sigY = reshape(conv_sigY(floor((n_timepoints_trial-1)/2):end-1-ceil((n_timepoints_trial-1)/2)),n_timepoints_trial,n_trials);
        phase_sigY = angle(conv_sigY);

        % Phase angle differences
        phase_diffs = phase_sigX-phase_sigY;

        % Compute phase synchrony
        hold = abs(mean(exp(1i*phase_diffs),2));
        ps(fi,:) = hold(timepoints_of_interest);

        % Cycle window in indices for this frequency
        cyclewindow_idx = round((1000/frequencies(fi))*cycle_width(fi)/(1000/samplerate));

        % Phase synchronization
        for ti=1:length(timepoints_of_interest)
            phasesynch = abs(mean(exp(1i*phase_diffs(timepoints_of_interest(ti)-cyclewindow_idx:timepoints_of_interest(ti)+cyclewindow_idx,:)),1));
            ISPC(fi,ti) = mean(phasesynch);
        end

    end

    % Baseline correction?
    if ~isempty(baseline)
        fprintf('\n\n\tPerforming baseline correction..\n')
        ratioISPC = zeros(n_frequencies, length(timepoints_of_interest));
        baseline_points = dsearchn(tvec_out',baseline'); % Indices of the baseline start and end

        for f=1:n_frequencies
            baseline  = mean(ISPC(f,baseline_points(1):baseline_points(2)));
            ratioISPC(f,:) = ISPC(f,:) / baseline;
        end
    end

    % Finished!
    fprintf('\n\n\tFinished!\n')
end


%ft_trialifiedX=ft_itpc.PFC.match_corr;
%ft_trialifiedY=ft_itpc.HPC.match_corr;
%tvec_in=tvec_unbuffed;
%frequencies=1:1:140;
%samplerate=1000;
%wavelet_cycles=2;
%tvec_out=-3:.05:8;
%baseline=[-3 -2];