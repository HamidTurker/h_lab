function events_ts = h_LoadEvents(cfg_in)
% Load Neuralynx .nev files
%
% INPUTS:
%
%   cfg.eventList:    cell array with event strings to process
%                     if not specified, load all
%   cfg.eventLabel:   cell array with strings to use as labels
%                     (should match eventList in size)
%
% OUTPUTS:
%
%   events_ts:        ts object with event timestamps
%
% HBT 2023 Apr 21

%% Initialize some parameters
cfg_def = [];
cfg = h_ProcessConfig(cfg_def,cfg_in); % this takes fields from cfg_in and puts them into cfg

mfun = mfilename;

%% Load events file
fname = cfg_in.Elist{1};
if isempty(fname)
   error('No events file provided in cfg.Elist{1}.'); 
end

% Initialize frames and read in event file
EVTimeStamps = []; EventIDs = [];
TTLs = []; EventStrings = {}; EventNumber = [];
hold = read_neuralynx_nev(fname);
for ii = 1:length(hold)
    EVTimeStamps(ii,1)    = double(hold(ii).TimeStamp);
    EventIDs(ii,1)        = double(hold(ii).EventId);
    TTLs(ii,1)            = double(hold(ii).TTLValue);
    EventStrings(ii,1)    = {hold(ii).EventString};
    EventNumber(ii,1)     = double(hold(ii).EventNumber);
end

% Obtain event LIST to process
if ~isfield(cfg,'eventList')      
        cfg.eventList = unique(EventStrings);        
else     
        if ~isa(cfg.eventList,'cell')
            error('LoadEvents: cfg.eventList should be a cell array.');
        end     
end

% Obtain event LABEL to process 
if ~isfield(cfg,'eventLabel')     
        cfg.eventLabel = cfg.eventList; 
end

%% Get times for each event label
events_ts = h_ts;

for iEvent = 1:length(cfg.eventList)
   
    % Which event?
    ev_string = cfg.eventList{iEvent};

    % Find where those events are and correct their time stamp
    ev_id = strmatch(ev_string, EventStrings, 'exact');
    ev_t = EVTimeStamps(ev_id)*10^-6;
    
    % Check if this eventLabel already exists, if so append (not create new)
    label_idx = strmatch(cfg.eventLabel{iEvent}, events_ts.label, 'exact');
    
    if isempty(label_idx) % label doesn't exist, create new        
        events_ts.t{iEvent} = ev_t;
        events_ts.label{iEvent} = cfg.eventLabel{iEvent};
        
        if ~iscolumn(events_ts.t{iEvent})
            events_ts.t{iEvent} = events_ts.t{iEvent}';
        end
        
    else % label exists, append and sort (sort is helpful for later use with nearest_idx3)        
        events_ts.t{label_idx} = sort(cat(1, events_ts.t{label_idx}, ev_t));
                
    end    
end

% Add session info
events_ts.cfg.SessionID = cfg.Flist{1};
events_ts.cfg.filename = fname;

%% Housekeeping
events_ts.cfg.history.mfun = cat(1,events_ts.cfg.history.mfun,mfun);
events_ts.cfg.history.cfg = cat(1,events_ts.cfg.history.cfg,{cfg});
