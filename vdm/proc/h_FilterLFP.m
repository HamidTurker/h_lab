function lfp_tsd = h_FilterLFP(cfg_in,lfp_tsd)
% Filter the LFP data
%
% INPUT:
% tsd with LFP data
%
% OUTPUT:
% tsd with filtered LFP data
%
% CFG OPTIONS with defaults:
%
% cfg.Type = 'butter'; % {'cheby1','butter','fdesign'} -- type of filter to be used
% cfg.Order = 4; % filter order
% cfg.Display_filter = 0; % show output of fvtool on filter
% cfg.Bandtype = 'bandpass'; % 'highpass', 'lowpass'
% cfg.Ripple = 0.5; % passband ripple (in dB) for Chebyshev filters only
% cfg.Filter = [6 10]; filter range to use (in Hz)
% cfg.Verbose = 1; If 1 display helpful text in command window, if 0 don't
%
% HBT 2023 Apr 23

%% Initialize default parameters
cfg_def.Type = 'butter';
cfg_def.Order = 4;
cfg_def.Display_filter = 0;
cfg_def.Bandtype = 'bandpass'; % only implemented for butterworth
cfg_def.Ripple = 0.5; % passband ripple (in dB) for Chebyshev filters only
cfg_def.Filter = [1 500];
cfg_def.Verbose = 1;

mfun = mfilename;
cfg = h_ProcessConfig(cfg_def,cfg_in); % this takes fields from cfg_in and puts them into cfg


%% Run some checks on the Fs
if ~CheckTSD(lfp_tsd)
   error('h_FilterLFP.m: CheckTSD failed.'); 
end

% Check reported Fs in headers
nSignals = length(lfp_tsd.cfg.Header);
for iS = nSignals:-1:1    
   if isfield(lfp_tsd.cfg.Header{iS},'SamplingFrequency')
       reported_Fs(iS) = lfp_tsd.cfg.Header{iS}.SamplingFrequency;
   elseif isfield(lfp_tsd.cfg.Header{iS},'Fs')
       reported_Fs(iS) = lfp_tsd.cfg.Header{iS}.Fs;
   else
      error('Unknown Fs.'); 
   end  
end

reported_Fs = unique(reported_Fs);
if length(reported_Fs) > 1
   error('h_FilterLFP.m: multiple sampling frequencies in header.'); 
end

tvec_diffs = diff(lfp_tsd.tvec);
median_Fs = 1./median(tvec_diffs);

if cfg.Verbose; fprintf('h_FilterLFP.m: reported Fs %.2f, median tvec Fs %.2f.\n',reported_Fs,median_Fs); end

Fs = reported_Fs; % Reported Fs checks out. Proceed..


%% Construct the filter
Wn = cfg.Filter ./ (Fs/2); % Convert to units of half sample rate
switch cfg.Type
   
    case 'butter'
        
        %[z,p,k] = butter(cfg.order,Wn);
        switch cfg.Bandtype
            case 'bandpass'
                [b,a] = butter(cfg.Order,Wn); 
            case 'highpass'
                [b,a] = butter(cfg.Order,Wn,'high'); 
            case 'lowpass'
                [b,a] = butter(cfg.Order,Wn,'low'); 
        end
        
    case 'cheby1'
        
        %[z,p,k] = cheby1(cfg.Order,cfg.R,Wn);
        [b,a] = cheby1(cfg.Order,cfg.R,Wn);
    
    case 'fdesign'
        
        d = fdesign.bandpass('N,F3dB1,F3dB2',cfg.Order,cfg.Filter(1),cfg.Filter(2),Fs);
        Hd = design(d,'butter');
        b = Hd.sosMatrix; a = Hd.scaleValues;
end

% Convert to SOS format
%[sos,g] = zp2sos(z,p,k);

% Display if requested
if cfg.Display_filter
    %h = dfilt.df2sos(sos,g);
    %fvtool(h);
    if exist('Hd','var')
        fvtool(Hd);
    else
        fvtool(b,a);
    end
    fprintf('h_FilterLFP.m: Paused to display filter, press a key to continue...\n');
    pause;
end


%% Process signals
for iS = 1:nSignals
    
    if cfg.Verbose; fprintf('h_FilterLFP.m: Filtering signal %d/%d...\n',iS,nSignals); end
    
    % Check for NaNs in the data. If there are, issue a warning and replace with zeros
    temp_sig = lfp_tsd.data(iS,:);
    
    nan_idx = find(isnan(temp_sig));
    
    if ~isempty(nan_idx)
        fprintf('WARNING: h_FilterLFP.m: Signal %d contains NaNs (%d).\n',iS,length(nan_idx));
        temp_sig(nan_idx) = 0;
    end
    
    % Filter (finally)
    %temp_sig = filtfilt(sos,g,temp_sig);
    temp_sig = filtfilt(b,a,temp_sig);
    
    % Reinstate NaNs and put signal back into tsd
    temp_sig(nan_idx) = NaN;
    lfp_tsd.data(iS,:) = temp_sig;
end

%% Housekeeping
lfp_tsd.cfg.history.mfun = cat(1,lfp_tsd.cfg.history.mfun,mfun);
lfp_tsd.cfg.history.cfg = cat(1,lfp_tsd.cfg.history.cfg,{cfg});