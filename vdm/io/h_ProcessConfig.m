function [cfg_out] = h_ProcessConfig(cfg_default,cfg_in)
% Processes the cfg file and the cfg_in
%
%   INPUTS:
%       cfg_default: default parameters stored as a struct within the containing function
%       cfg_in: input parameters from the input cfg file
%
%   OUTPUT:
%       Returns a new cfg file containing the default parameters and any updated cfg
%       parameters given the input file.
% 
% HBT 2023 Apr 21


%% Set output cfg to cfg structure with the default values
cfg_out = cfg_default;

%% Process input cfg parameters into workspace
% NOTE - The two processes (replace and add) are separated only for clarity 
% and ease of debugging. It can easily be compressed to half the amount of code (and
% therefore time).

if ~isempty(cfg_in)
    
    % If there are default parameters,replace them with cfg_in input
    if ~isempty(cfg_default)
        default_F = fieldnames(cfg_default);
        nFields = length(default_F);

        for i = 1:nFields
            idF = default_F{i};
            if isfield(cfg_in,idF)
                cfg_out.(idF) = cfg_in.(idF);
            end % set cfg_out fields    
        end % iterate default cfg fields
    end

    % Add input parameters
    new_F = fieldnames(cfg_in);
    nFields = length(new_F);
    
    for i = 1:nFields
        inF = new_F{i};
        if ~isfield(cfg_default,inF)
            cfg_out.(inF) = cfg_in.(inF);
        end % set cfg_out fields    
    end % iterate extra cfg fields
    
end
