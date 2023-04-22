function tsd_out = h_tsd(varargin)
% TSD Time-stamped data datatype constructor
%   Time-stamped data, or tsd, is one of the main data types. Each tsd
%   struct contains a time vector and corresponding data (such as voltages 
%   or position). 
%
%   function tsd_out = h_tsd(tvec,data)
%
%   function tsd_out = h_tsd(tvec,data,'label')
%
%   INPUTS:
%      tvec: time vector (time stamps in ascending order)
%      data: measurements taken at each time stamp
%     label: unique identifier
%     units: measurement units
%
%   OUTPUTS
%      iv_out: iv struct with fields:
%           .type   - 'tsd'; datatype identification
%           .units  - 'au'; unit marker, e.g., 'cm','px','mV'
%           .tvec   - [nx1] double, where n is the number of time samples taken
%           .data   - [1xn] double containing measures taken at each sample
%                     point
%           .label  - Unique identifier for the data contained with the tsd,
%                     such as a record of the data source (in the case of a 
%                     CSC, the label corresponds to the tetrode that recorded 
%                     the signal)
%           .cfg    - record of the data's history including config 
%                     parameters and functions visited
%
% HBT 2023 Apr 21

%% Initialize some default parameters, if not specified by user
tsd_out.type = 'tsd';
tsd_out.units = 'au';
tsd_out.tvec = [];
tsd_out.data = [];
tsd_out.label = {};

%% Check user-defined args
if nargin == 2
    tsd_out.tvec = varargin{1};
    tsd_out.data = varargin{2};
elseif nargin == 3
    tsd_out.tvec = varargin{1};
    tsd_out.data = varargin{2};
    tsd_out.label = varargin{3};
elseif nargin == 4
    tsd_out.tvec = varargin{1};
    tsd_out.data = varargin{2};
    tsd_out.label = varargin{3};
    tsd_out.units = varargin{4};
end

%% Housekeeping
tsd_out.cfg.history.mfun{1} = mfilename;
tsd_out.cfg.history.cfg{1} = [];