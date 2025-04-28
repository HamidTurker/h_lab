function [XtoY, YtoX, ratioXtoY, ratioYtoX, percXtoY, percYtoX, zXtoY, zYtoX, rangeXtoY, rangeYtoX, stats_out, test_outXtoY, test_outYtoX] = h_Granger2(ft_trialifiedX, ft_trialifiedY, order, tvec_in, tvec_in_samplerate, ...
                                                                                                                                                        tvec_out, timewindow, frequencies, baseline, subtractERP, ...
                                                                                                                                                        tvec_leftedge2firstevent, df_timepoints, mintimes, subsample, save_stats, cfg_test)
    % Compute Granger at linearly spaced time points that can be unique for each trial
    fprintf('\n\tStarting Granger causality analysis.. \n')


%ft_trialifiedX=granger.PFC.match_corr;
%ft_trialifiedY=granger.HPC.match_corr;
%order=granger.order;
%tvec_in=tvec_buffed;
%tvec_in_samplerate=param.Fs; 
%timewindow=tsegment;
%frequencies=param.FOIs;
%baseline=param.baseline;
%subtractERP=true;
%tvec_leftedge2firstevent;
%df_timepoints=df.timepoints;
%subsample=[];
%save_stats=[];
%toi=287;
%mintimes=400;

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

        % Also subsample the df_timepoints
        hold={};
        for t = 1:length(chosen)
            hold{t} = df_timepoints{chosen(t)};
        end
        % Replace in the main data frame
        df_timepoints = hold';

    else
        fprintf('\tNo subsampling..\n')
    end

    % Params
    n_trials = length(ft_trialifiedX.trial);
    n_timepoints_data = length(tvec_in);
    %n_timepoints_out = length(tvec_out);
    n_timepoints_out = mintimes;
    n_frequencies = length(frequencies);
    n_timewin_points = round(timewindow*tvec_in_samplerate);
    %n_order_points = round(order*tvec_in_samplerate);

    % Initialize
    XtoY = zeros(n_frequencies, n_timepoints_out); % Granger prediction X -> Y
    YtoX = zeros(n_frequencies, n_timepoints_out); % Granger prediction Y -> X

    zXtoY=[]; zYtoX=[];
    ratioXtoY=[]; ratioYtoX=[];
    percXtoY=[]; percYtoX=[];
    rangeXtoY=[]; rangeYtoX=[];

    stats_out = [];
    test_outXtoY = [];
    test_outYtoX = [];

    % Are we subtracting the ERP?
    if (subtractERP)
        
        % Talk
        fprintf('\tSubtracting ERP..\n')

        % Format for ERP computation
        df_ERP_X = zeros(n_trials,n_timepoints_data); df_ERP_Y = zeros(n_trials,n_timepoints_data);
        for t = 1:n_trials
            df_ERP_X(t,:) = ft_trialifiedX.trial{t};
            df_ERP_Y(t,:) = ft_trialifiedY.trial{t};
        end
        
        % ERP
        ERP_X = mean(df_ERP_X);
        ERP_Y = mean(df_ERP_Y);

        % Subtract ERP
        for t = 1:n_trials
            ft_trialifiedX.trial{t} = ft_trialifiedX.trial{t} - ERP_X;
            ft_trialifiedY.trial{t} = ft_trialifiedY.trial{t} - ERP_Y;
        end
    else
        % Talk
        fprintf('\tNot subtracting ERP..\n')
    end

    % Compute the Granger prediction score for each timepoint of interest
    for toi = 1:n_timepoints_out

        % Talk
        fprintf('\n\tSetting up for time point %d/%d.. \n',toi,n_timepoints_out)

        % We're gonna grab the time window from each trial and line them all up, 
        % as though they are all one continuous signal (row = X,Y ; column = time point).
        % Note that we are also detrending and z-scoring the time window.
        snipX = []; snipY = [];
        for t = 1:n_trials

            % This particular trial's time points of interest
            timepoints_of_interest = dsearchn(tvec_in',df_timepoints{t}');

            % Current TOI and time window indices
            idx_TOI = timepoints_of_interest(toi);                                   % The index of the time point for the current TOI
            idx_start = idx_TOI - floor(n_timewin_points/2);                         % The index of the start of the sliding window surrounding the TOI
            idx_end = idx_TOI + floor(n_timewin_points/2)-mod(n_timewin_points+1,2); % The index of the end of the sliding window surrounding the TOI

            % Snippet of data (the size of timewindow) surrounding the current
            % TOI, across all trials and for both X and Y time series. We also
            % need to do some reformatting for the bivariate autoregression.
            snipX = horzcat(snipX, zscore(detrend(ft_trialifiedX.trial{t}(idx_start:idx_end))));
            snipY = horzcat(snipY, zscore(detrend(ft_trialifiedY.trial{t}(idx_start:idx_end))));
        end
        % Top row is all the snippets from the X signal trials, lined up horizontally.
        % The bottom row is all the snippets from the Y signal, also concattenated horizontally.
        snipXY = [snipX; snipY];

        % Fit AR models
        %[Ax,Ex] = armorf(snipXY(1,:), n_trials, n_timewin_points, order);
        %[Ay,Ey] = armorf(snipXY(2,:), n_trials, n_timewin_points, order);
        [Axy,E] = armorf(snipXY,      n_trials, n_timewin_points, order);

        % Code below is from the adapted version (Cohen, 2014) of the BSMART toolbox function pwcausal.m (Cui et al., 2008)
        % corrected covariance
        eyx = E(2,2) - E(1,2)^2/E(1,1);
        exy = E(1,1) - E(2,1)^2/E(2,2);
        N = size(E,1);

        for fi=1:n_frequencies

            % Talk
            fprintf('\n\tComputing Granger score at time point %d/%d for frequency %d/%d..',toi,n_timepoints_out,fi,n_frequencies)

            % Transfer matrix (note the similarity to Fourier transform)
            H = eye(N);
            for m = 1:order
                H = H + Axy(:, (m-1)*N+1:m*N) * exp(-1i*m*2*pi*frequencies(fi) / tvec_in_samplerate);
            end

            Hi = inv(H);
            S  = H\E*Hi'/tvec_in_samplerate;

            % Granger prediction per frequency
            XtoY(fi,toi) = log( abs(S(2,2)) / abs(S(2,2)-(Hi(2,1) * exy * conj(Hi(2,1))) / tvec_in_samplerate) );
            YtoX(fi,toi) = log( abs(S(1,1)) / abs(S(1,1)-(Hi(1,2) * eyx * conj(Hi(1,2))) / tvec_in_samplerate) );
        end

    end

    % Normalizations
    zXtoY = zeros(n_frequencies, n_timepoints_out); rangeXtoY = zeros(n_frequencies, n_timepoints_out);
    zYtoX = zeros(n_frequencies, n_timepoints_out); rangeYtoX = zeros(n_frequencies, n_timepoints_out);
    for f=1:n_frequencies
        zXtoY(f,:)=zscore(XtoY(f,:));
        zYtoX(f,:)=zscore(YtoX(f,:));
        rangeXtoY(f,:)=normalize(XtoY(f,:),'range');
        rangeYtoX(f,:)=normalize(YtoX(f,:),'range');
    end

    % Baseline correction?
    if ~isempty(baseline)
        
        ratioXtoY = zeros(n_frequencies, n_timepoints_out); percXtoY = zeros(n_frequencies, n_timepoints_out);
        ratioYtoX = zeros(n_frequencies, n_timepoints_out); percYtoX = zeros(n_frequencies, n_timepoints_out);
        baseline_points = dsearchn(tvec_leftedge2firstevent',baseline'); % Indices of the baseline start and end

        for f=1:n_frequencies
            baseline_XtoY  = mean(XtoY(f,baseline_points(1):baseline_points(2)));
            ratioXtoY(f,:) = XtoY(f,:) / baseline_XtoY;
            percXtoY(f,:) = (XtoY(f,:) / baseline_XtoY * 100) - 100;
                 
            baseline_YtoX = mean(YtoX(f,baseline_points(1):baseline_points(2)));
            ratioYtoX(f,:) = YtoX(f,:) / baseline_YtoX;
            percYtoX(f,:) = (YtoX(f,:) / baseline_YtoX * 100) - 100;
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
                    bl_avg  = mean(XtoY(f,bl_points(1):bl_points(2)));
                    test_outXtoY(f,:) = XtoY(f,:) / bl_avg;

                    bl_avg  = mean(YtoX(f,bl_points(1):bl_points(2)));
                    test_outYtoX(f,:) = YtoX(f,:) / bl_avg;
                end

            end
        end
    end

end


function varargout = armorf(x,Nr,Nl,p)
%ARMORF   AR parameter estimation via LWR method by Morf modified.
%   x is a matrix whose every row is one variable's time series
%   Nr is the number of realizations, Nl is the length of every realization
%   If the time series are stationary long, just let Nr=1, Nl=length(x)
%   p is the order of AR model
%
%   A = ARMORF(X,NR,NL,P) returns the polynomial coefficients A corresponding to
%     the AR model estimate of matrix X using Morf's method.
%
%   [A,E] = ARMORF(...) returns the final prediction error E (the
%   covariance matrix of the white noise of the AR model).
%
%   [A,E,K] = ARMORF(...) returns the vector K of reflection
%     coefficients (parcor coefficients).
%
%   Ref: M. Morf, etal, Recursive Multichannel Maximum Entropy Spectral Estimation,
%              IEEE trans. GeoSci. Elec., 1978, Vol.GE-16, No.2, pp85-94.
%        S. Haykin, Nonlinear Methods of Spectral Analysis, 2nd Ed.
%              Springer-Verlag, 1983, Chapter 2
%
%   finished on Aug.9, 2002 by Yonghong Chen

% Initialization
[L,N]=size(x);
R0=zeros(L,L);
R0f=R0;
R0b=R0;
pf=R0;
pb=R0;
pfb=R0;
ap(:,:,1)=R0;
bp(:,:,1)=R0;
En=R0;
for i=1:Nr
    En=En+x(:,(i-1)*Nl+1:i*Nl)*x(:,(i-1)*Nl+1:i*Nl)';
    ap(:,:,1)=ap(:,:,1)+x(:,(i-1)*Nl+2:i*Nl)*x(:,(i-1)*Nl+2:i*Nl)';
    bp(:,:,1)=bp(:,:,1)+x(:,(i-1)*Nl+1:i*Nl-1)*x(:,(i-1)*Nl+1:i*Nl-1)';
end
ap(:,:,1) = inv((chol(ap(:,:,1)/Nr*(Nl-1)))');
bp(:,:,1) = inv((chol(bp(:,:,1)/Nr*(Nl-1)))');
for i=1:Nr
    efp = ap(:,:,1)*x(:,(i-1)*Nl+2:i*Nl);
    ebp = bp(:,:,1)*x(:,(i-1)*Nl+1:i*Nl-1);
    pf = pf + efp*efp';
    pb = pb + ebp*ebp';
    pfb = pfb + efp*ebp';
end
En = chol(En/N)'; % Covariance of the noise

% Initial output variables
coeff = [];%  Coefficient matrices of the AR model
kr=[];  % reflection coefficients

for m=1:p
    % Calculate the next order reflection (parcor) coefficient
    ck = inv((chol(pf))')*pfb*inv(chol(pb));
    kr=[kr,ck];
    % Update the forward and backward prediction errors
    ef = eye(L)- ck*ck';
    eb = eye(L)- ck'*ck;

    % Update the prediction error
    En = En*chol(ef)';
    E = (ef+eb)./2;

    % Update the coefficients of the forward and backward prediction errors
    ap(:,:,m+1) = zeros(L);
    bp(:,:,m+1) = zeros(L);
    pf = zeros(L);
    pb = zeros(L);
    pfb = zeros(L);

    for i=1:m+1
        a(:,:,i) = inv((chol(ef))')*(ap(:,:,i)-ck*bp(:,:,m+2-i));
        b(:,:,i) = inv((chol(eb))')*(bp(:,:,i)-ck'*ap(:,:,m+2-i));
    end
    for k=1:Nr
        efp = zeros(L,Nl-m-1);
        ebp = zeros(L,Nl-m-1);
        for i=1:m+1
            k1=m+2-i+(k-1)*Nl+1;
            k2=Nl-i+1+(k-1)*Nl;
            efp = efp+a(:,:,i)*x(:,k1:k2);
            ebp = ebp+b(:,:,m+2-i)*x(:,k1-1:k2-1);
        end
        pf = pf + efp*efp';
        pb = pb + ebp*ebp';
        pfb = pfb + efp*ebp';
    end
    ap = a;
    bp = b;
end
for j=1:p
    coeff = [coeff,inv(a(:,:,1))*a(:,:,j+1)];
end

varargout{1} = -coeff;
if nargout >= 2
    varargout{2} = En*En';
end
if nargout >= 3
    varargout{3} = kr;
end
end
    