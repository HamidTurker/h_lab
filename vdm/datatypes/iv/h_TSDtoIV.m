function [iv_out,thr_out] = h_TSDtoIV(cfg_in,tsd_in)
% Create interval data from tsd by tresholding
%
% INPUTS:
%
% tsd_in: input tsd
%
% CFG OPTIONS:
% cfg.Method = 'zscore';
% cfg.Threshold = 5;
% cfg.Operation =  '>'; % '<', '>'
% cfg.Merge_thr = 0.05; % merge events closer than this
% cfg.Target = []; % which data (label) to use
% cfg.Minlen = 0.05; % minimum interval length
% cfg.Verbose = 1; 1 display command window text, 0 don't
%
% OUTPUTS:
%
% iv_out: output interval data
% thr_out: threshold used (in same units as input; helpful if z-scored
%   threshold is to be applied later to other data)
%
% HBT 2023 Apr 23

%% Define default parameters
cfg_def.Method = 'zscore';
cfg_def.Threshold = 5;
cfg_def.Operation = '>'; % return intervals where threshold is exceeded
cfg_def.Merge_thr = 0.05; % merge events closer than this
cfg_def.Target = [];
cfg_def.Minlen = 0.05; % minimum interval length
cfg_def.Verbose = 1;

mfun = mfilename;
cfg = h_ProcessConfig(cfg_def,cfg_in); % should take whatever is in cfg_in and put it into cfg!

iv_out = h_iv; % initialize new iv struct

%% Check if conditions are in place
nData = size(tsd_in.data,1);
if nData > 1
    if ~isempty(cfg.Target)
        temp_data = getd(tsd_in,cfg.Target);
    else
       error('Multiple data dimensions exist but no label is specified.'); 
    end
else
    temp_data = tsd_in.data;
end

%% Apply transform to data if requested
thr_out = cfg.Threshold;
switch cfg.Method
    case 'zscore'
        [temp_data,mu,sigma] = zscore(temp_data);
        switch cfg.Operation
            case '>' 
                thr_out = mu + cfg.Threshold.*sigma; % store threshold for output argument
            case '<'
                thr_out = mu - cfg.Threshold.*sigma;
        end
    case 'percentile'
        temp_data_sorted = sort(temp_data,'ascend');
        thr_out_idx = round(cfg.Threshold.*length(temp_data));
        thr_out = temp_data_sorted(thr_out_idx);
        
        temp_data = tiedrank(temp_data)./length(temp_data); % assign percentile to each data point       
end

%% Detect crossings
switch cfg.Operation
    case '>'
        detec = temp_data > cfg.Threshold;
    case '>='
        detec = temp_data >= cfg.Threshold;
    case '<'
        detec = temp_data < cfg.Threshold;
    case '<='
        detec = temp_data <= cfg.Threshold;
    case 'range'
        detec = temp_data > cfg.Threshold(1) & temp_data < cfg.Threshold(2);
end

% Pad the detection so we can deal with cases where first or last samples are detects
detec = cat(2,0,detec,0);

dfs = diff(detec);
up_idx = find(dfs == 1);
down_idx = find(dfs == -1) - 1;

up_t = tsd_in.tvec(up_idx);
down_t = tsd_in.tvec(down_idx);

% remove doubles
d = up_t(2:end)-down_t(1:end-1);

merge_idx = find(d < cfg.Merge_thr);
down_t(merge_idx) = [];
up_t(merge_idx+1) = [];

% remove too short
iv_len = down_t-up_t;
keep_idx = iv_len > cfg.Minlen;

iv_out.tstart = up_t(keep_idx);
iv_out.tend = down_t(keep_idx);

% ensure column vectors
if ~iscolumn(iv_out.tstart)
    iv_out.tstart = iv_out.tstart';
end

if ~iscolumn(iv_out.tend)
    iv_out.tend = iv_out.tend';
end

if cfg.Verbose
    disp([mfun,': ',num2str(length(iv_out.tstart)),' intervals found.'])
end

%% Housekeeping
iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});