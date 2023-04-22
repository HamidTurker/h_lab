function csc_tsd = h_LoadCSC(cfg_in)
% Load Neuralynx .ncs files, handling missing samples/discontinuities automatically.
%
% INPUTS:   A single config structure. The structure can have the
%           following fields. If left empty, these default values specified
%           here will be used.
%
%   cfg.Flist (cell array):             Cell array containing filenames to load
%                                       If no file list field is specified, load all *.Ncs files in current dir
%   cfg.TimeConvFactor = 10^-6;         From nlx units to seconds
%   cfg.VoltageConvFactor = 1;          If 0, return raw AD bits.
%                                       If 1, output in Volts, 1000 in mV, 10^6 in uV
%   cfg.DecimateByFactor = [];          If non-empty, decimate data by specified factor
%   cfg.Verbose = 1;                    Allow or suppress displaying of command window text
%
%
% OUTPUTS:
%
%   csc_tsd: CSC data tsd struct
%
% HBT 2023 Apr 20

%% Initialize some default parameters, if not specified by user
cfg_def.Flist = {};
cfg_def.TimeConvFactor = 10^-6; % 10^-6 means convert nlx units to seconds
cfg_def.VoltageConvFactor = 1;
cfg_def.DecimateByFactor = [];
cfg_def.Verbose = 1;

mfun = mfilename;

%% Overwrite default cfg parameters with user-defined parameters
cfg = h_ProcessConfig(cfg_def,cfg_in);

%% Check if we can load the requested file(s) and set some more fields
% No file list provided, load everything
if isempty(cfg.Flist)
    cfg.Flist = FindFiles('*.Ncs');    
else    
    if ~isa(cfg.Flist,'cell')
        error('cfg.Flist should be a cell array.');
    end   
end

% File doesn't exist
if isempty(FindFiles(cfg.Flist{1}))
    error('File does not exist in directory. Check spelling.');
end

% Initialize final output data as Time-Stamped Datatype (TSD)
csc_tsd = h_tsd;

% Write in the unit
if cfg.VoltageConvFactor == 1
    csc_tsd.units = 'V';
elseif cfg.VoltageConvFactor == 1000
    csc_tsd.units = 'mV';
elseif cfg.VoltageConvFactor == 10^6
    csc_tsd.units = 'uV';
elseif cfg.VoltageConvFactor == 0
    csc_tsd.units = 'ADbits';
else
    error('Input voltage conversion factor is not a valid value.')
end

%% Proceed to load
Flist = sort(cfg.Flist);
nFiles = length(Flist);

% Track sample counts for different files, will throw error if not equal..
sample_count_tvec = nan;
sample_count_data = nan;

% Talk to me
if cfg.Verbose; fprintf('%s: Loading %d file(s)...\n',mfun,nFiles); end

% Loop through the list of all requested files and do the thing..
for iF = 1:nFiles
    
    fname = Flist{iF};
    
    % Load raw data
    %[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);
    %hold_ft = ft_read_neuralynx_interp({fname});
    hold = read_neuralynx_ncs(fname);
    Timestamps = double(hold.TimeStamp);
    SampleFrequencies = hold.SampFreq;
    NumberOfValidSamples = hold.NumValidSamp;
    Samples = hold.dat;
    Header = hold.hdr;

    % Talk
    if cfg.Verbose; fprintf('Any discontinuities in the recording will be corrected..\n'); end

    % Disabled channels cannot be loaded
    if Timestamps == 0 
        error(['No csc data (disabled tetrode channel). Consider deleting ',fname,'.']);
    end
    
    % Check for constant sampling frequency
    Fs = unique(SampleFrequencies);
    if length(Fs) ~= 1
        error('More than one sampling frequency found!');
    end
       
    % data format is n columns (blocks) of 512 samples, with a timestamp for each first
    % sample per block
    %
    % b1_1    b2_1    b3_1    ...   bn_1
    % b1_2    b2_2    b3_2    ...   bn_2
    % .       .       .       ...   .
    % .       .       .       ...   .
    % b1_512  b2_512  b3_512  ...   bn_512
    %
    % But not all blocks have 512 good samples, some have fewer than that.
    % This number is given in NumberOfValidSamples.
    
    % Convert units
    Timestamps = Timestamps .* cfg.TimeConvFactor;
    if cfg.VoltageConvFactor ~= 0
        Samples = Samples .* cfg.VoltageConvFactor .* Header.ADBitVolts;
    end
    
    % Construct within-block time vector
    nSamplesPerBlock = size(Samples,1);
    block_tvec = 0:1./Fs:(nSamplesPerBlock-1)./Fs;
    
    % Slow but tractable: loop over blocks remembering to initialize variables first   
    data = nan(numel(Samples),1); % allocate memory and then fill; can trim later
    tvec = nan(numel(Samples),1);
    
    nBlocks = length(Timestamps);
    badBlocks = 0; % Bad blocks counter
    idx = 1; % Index for valid samples (to be updated with NumberOfValidSamples
    for iB = 1:nBlocks
        
        nvs = NumberOfValidSamples(iB);
        if nvs ~= 512, badBlocks = badBlocks + 1; end
        
        currentData = Samples(1:nvs,iB);
        currentTime = block_tvec(1:nvs)+Timestamps(iB);
        
        data(idx:idx+nvs-1) = currentData;
        tvec(idx:idx+nvs-1) = currentTime;
        
        idx = idx + nvs;
        
    end % of block loop
    
    % How many bad blocks?
    cfg.badBlocks = badBlocks;
    if cfg.Verbose; fprintf('%s: %s %d/%d bad blocks found (%.2f%%).\n',mfun,fname,badBlocks,nBlocks,(badBlocks./nBlocks).*100); end
    
    % Remove nans
    data = data(~isnan(data));
    tvec = tvec(~isnan(tvec));
    
    % Note the vector sizes
    sample_count_tvec(iF) = length(tvec);
    sample_count_data(iF) = length(data);
    
    % Check if the data is the same length for each channel.  
    if iF >1 && length(data) ~= length(csc_tsd.data(iF-1,:))
        error('Data lengths differ across channels.');
    end
    
    % Decimate data, if specified..
    if ~isempty(cfg.DecimateByFactor)
        
        fprintf('%s: Decimating by factor %d...\n',mfun,cfg.DecimateByFactor)
        data = decimate(data,cfg.DecimateByFactor);
        tvec = tvec(1:cfg.DecimateByFactor:end);
        Header.SamplingFrequency = Header.SamplingFrequency./cfg.DecimateByFactor;
        
    end

    % Done, add to tsd
    csc_tsd.tvec = tvec;
    csc_tsd.data(iF,:) = data;
    
    [~,fn,fe] = fileparts(fname);
    csc_tsd.label{iF} = cat(2,fn,fe);

    csc_tsd.cfg.Header{iF} = Header;
    
end

% Check if anything unequal --> error
if numel(unique(sample_count_data)) > 1
    error('Data sizes unequal.');
end

if numel(unique(sample_count_tvec)) > 1
   error('tvec sizes unequal.');
end

% Add sessionID
csc_tsd.cfg.SessionID = cfg.Flist{1};

%% Housekeeping
csc_tsd.cfg.history.mfun = cat(1,csc_tsd.cfg.history.mfun,mfun);
csc_tsd.cfg.history.cfg = cat(1,csc_tsd.cfg.history.cfg,{cfg});
