%% Vandermeer: LFP analysis with Neuralynx data, from start to finish
% https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:nsb2015:week0
% Note: the vandermeerlab (VDM) codebase functions have been adapted here into
%       h_lab functions
% Note: we'll be using Fieldtrip functions along the way

clear all; clc

%% Vandermeer codebase setup
%addpath '/Users/hamid/Neuralynx/vdm/'
%addpath(genpath('/Users/hamid/matlab/vandermeerlab/code-matlab/shared/'))
addpath(genpath('/Users/hamid/matlab/matlab/hlab_matlab/vdm/'))

%% Global directories and parameters
path.tutorial = '/Users/hamid/Desktop/Research/+Github/vandermeerlab/';
path.ratfolder = '/Users/hamid/Desktop/Research/+Github/vandermeerlab/data/MotivationalT/R042/R042-2013-08-18/';
path.cscdata = '/Users/hamid/Desktop/Research/+Github/vandermeerlab/data/MotivationalT/R042/R042-2013-08-18/R042-2013-08-18-CSC03a.ncs';
path.eventdata = '/Users/hamid/Desktop/Research/+Github/vandermeerlab/data/MotivationalT/R042/R042-2013-08-18/R042-2013-08-18-Events.nev';
path.io = '/Users/hamid/matlab/fieldtrip-20230418/fileio/private';


%% Load in CSC data
% Since we're loading in a .ncs file, we need Fieldtrip's
% read_neuralynx_ncs. We'll cd to the directory that has the input/output
% scripts needed to read in the data. But, we're running the adapted
% LoadCSC function from the VDM lab. This is done by setting up a so-called
% configuration (cfg) structure and giving that to the function.
cd(path.io)

% Set up the cfg structure. Provide the path to the filename of interest
% and ensure that it's in a cell array.
cfg = [];

% Load CSC data from its location
cfg.Flist = {path.cscdata};
LFP = h_LoadCSC(cfg);

% Load Event data from its location
cfg.Elist = {path.eventdata};
events = h_LoadEvents(cfg);

%% Filter the data in the ripple band (150-220Hz)
cd(path.ratfolder)
cfg = []; cfg.Filter = [150 220]; % specify passband
LFPfilt = h_FilterLFP(cfg,LFP);

% extract the signal envelope by Hilbert transform

%% Detect intervals that pass a threshold
cfg = []; cfg.Method = 'zscore';
cfg.Threshold = 5; cfg.Operation =  '>'; % return intervals where threshold is exceeded
SWR = h_TSDtoIV(cfg,LFPfilt); % make intervals (corresponding to SWR events)





%%%%%%%%%%%%



%% Global directories and parameters
path.base = '/Users/hamid/Desktop/Research/+Projects/Smith Lab/N-back_interference/+analyses/';
path.ratfolder = '/Users/hamid/Desktop/Research/+Projects/Smith Lab/N-back_interference/CMTS_recording/R1813/R1813_2014-04-13/';
path.cscdata = '/Users/hamid/Desktop/Research/+Projects/Smith Lab/N-back_interference/CMTS_recording/R1813/R1813_2014-04-13/CSC7.ncs';
path.eventdata = '/Users/hamid/Desktop/Research/+Projects/Smith Lab/N-back_interference/CMTS_recording/R1813/R1813_2014-04-13/Events.nev';
path.io = '/Users/hamid/matlab/fieldtrip-20230418/fileio/private';


%% Load in CSC data
% Since we're loading in a .ncs file, we need Fieldtrip's
% read_neuralynx_ncs. We'll cd to the directory that has the input/output
% scripts needed to read in the data. But, we're running the adapted
% LoadCSC function from the VDM lab. This is done by setting up a so-called
% configuration (cfg) structure and giving that to the function.
cd(path.io)

% Set up the cfg structure. Provide the path to the filename of interest
% and ensure that it's in a cell array. Our original sampling rate is 32k
% Hz. So, taking every 16th sample means we end up with 2k Hz.
cfg = []; cfg.DecimateByFactor = 16;

% Load CSC data from its location
cfg.Flist = {path.cscdata};
LFP = h_LoadCSC(cfg);

% Load Event data from its location
cfg.Elist = {path.eventdata};
events = h_LoadEvents(cfg);

%% Filter the data in the ripple band (150-220Hz)
cd(path.ratfolder)
cfg = []; cfg.Filter = [150 220]; % specify passband
LFPfilt = h_FilterLFP(cfg,LFP);

% extract the signal envelope by Hilbert transform

%% Detect intervals that pass a threshold
cfg = []; cfg.Method = 'zscore';
cfg.Threshold = 5; cfg.Operation =  '>'; % return intervals where threshold is exceeded
SWR = h_TSDtoIV(cfg,LFPfilt); % make intervals (corresponding to SWR events)











