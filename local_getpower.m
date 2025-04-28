function [power] = local_getpower(ft_spect, win, freqs, baseline, level)

    %level = "session";
    ft_spect = ft.PFC.unbuffed_epochs;
    %freqs = 20:1:30;
    %baseline=param.baseline;
    %win=[0 2];
    %cycles=param.wavelet_cycles;

    % Global tvec
    tvec = ft_spect{1}.time{1};

    % Feature-window and baseline bins
    [~, w_start] = min(abs(tvec-win(1))); % Start
    [~, w_end] = min(abs(tvec-win(2))); % End
    [~, b_start] = min(abs(tvec-baseline(1))); % Start
    [~, b_end] = min(abs(tvec-baseline(2))); % End

    n_sess = length(ft_spect);

    for s = 1:n_sess

            % Run FT on each session
            bline=[-3 -2];

            ft_spect{s}.label = {'PFC'};
            ft_spect{s}.fsample = 1000;

            spect = ft_spect{s};
            cfg              = [];
            cfg.output       = 'pow';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.foi          = freqs;
            cfg.toi          = tvec_unbuffed;
            %cfg.t_ftimwin    = cycles./cfg.foi;  % Cycles per time window
            cfg.t_ftimwin(1:length(cfg.foi)) = 0.2;
            cfg.tapsmofrq(1:length(cfg.foi)) = 10;
            %cfg.keeptrials   = 'yes'; % Need this for stats later
            powTFR = ft_freqanalysis(cfg, spect);

            powTFR.powspctrm = squeeze(powTFR.powspctrm);
            power(s) = mean(mean(powTFR.powspctrm(:,w_start:w_end)));

    end

end

power_nonmatchincr=power;
scatter(power_nonmatchincr,perf);

scatter(power_sess_nonmatchcorr,perf);

