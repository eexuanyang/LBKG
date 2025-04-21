% hNewHARQProcess Initialize HARQ processes
%
%   HARQPROCESSES = hNewHARQProcess(ENBUE) is a structure containing
%   the initialized parameters of hybrid automatic repeat request processes
%   given a downlink, uplink or sidelink configuration ENBUE.
%
%   HARQPROCESSES = hNewHARQProcess(ENB,PDSCH) initializes HARQ processes
%   for a given ENB structure and PDSCH configuration structure. 

%   Copyright 2015-2018 The MathWorks, Inc.

function harqProcesses = hNewHARQProcess(enbue,varargin)
    
    if(nargin>1)
        enbue.PDSCH = varargin{1};
    end
    % Identify the structure containing HARQ entity parameters (PDSCH in
    % downlink, PUSCH in uplink, UE in sidelink)
    if isfield(enbue,'PDSCH')
        harqstruct = enbue.PDSCH;
    elseif isfield(enbue,'PUSCH')
        harqstruct = enbue.PUSCH;
    else
        harqstruct = enbue; % Input structure is 'UE'
    end
    
    % If Modulation is a cell array multiple codewords are used 
    if (iscell(harqstruct.Modulation)) || isstring(harqstruct.Modulation)
        ncw = length(harqstruct.Modulation);
    else
        ncw = 1;
    end
    harqProcess.data = {};                     % Initialize data
    harqProcess.iniConfig = harqstruct;        % Initialization config
    harqProcess.iniConfig.RVIdx = ones(1,ncw); % Add RVIdx to config
    harqProcess.blkerr = zeros(1,ncw);         % Initialize block errors
    harqProcess.ncw = ncw;                     % Set number of codewords
    
    % Initialize transmit configuration structure
    harqProcess.txConfig = harqProcess.iniConfig;
    
    % Create an empty structure for each codeword
    harqProcess.decState(1:ncw) = deal(struct());
    
    % Create HARQ processes as indicated by NHARQProcesses
    harqProcesses(1:harqstruct.NHARQProcesses) = harqProcess;
    
end