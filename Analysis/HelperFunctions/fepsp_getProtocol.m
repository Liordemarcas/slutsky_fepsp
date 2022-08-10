function protocol_info = fepsp_getProtocol(varargin)

% Load / create protocol info. If sampling frequency & dead time are
% provided, convert time into sample numbers in important windows.
%
% INPUT (required):
%   protocol_id - string or char. Can be id (case sensitive) of any
%                 protocol saved before (default options are 'io' & 'stp').
%                 Can also be "custom" (lower case) to save a new protocol,
%                 or "display" (lower case) in order to get all saved protocols.

% INPUT (optional):
%   fs         - numeric scalar. Sampling frequency of the 
%                recorded data [Hz]. 
%   dt         - non-negative scalar. Dead time between
%                stimulus onset & earliest possible response [ms].
%                Used to omit stimulus artifact from analysis, by creating
%                windows distanced from each stimuli start & end.
%                Default: 2.
%
% OUTPUT:
%   struct "protocol_info" with the following fields:
%
%   protocol_id - char vector. protocol identifier.
%   nStim       - numeric scalar. Number of stimuli given in protocol.
%   stim_times  - numeric column vector. Start time of each stimuli [ms].
%   stim_length - numeric scalar. How long stimuli are [us].
%   rec_length  - numeric scalar. How long each trace recording is [ms].
%   pstim_win   - numeric n x 2 mat of stimulus number (row) x [start, end]
% %               of window (column). Area after each stimulation, from
%                 current stimulus end to the start of the next one, or 30
%                 millisecond after the last stimulation.
%
%   the following fields will be empty (in subfields if exist) if fs is
%   empty:
%
%   response    - struct. Area of the response, with fields:
%       base    - numeric mat, stimulus number (row) x start / end of window
%                 (column). Basis of response window, from dt after
%                 pstim_win starts to dt before it ends, in sample numbers.
%       win     - numeric row vector. All the samples inside the limits
%                 defined by base of response window (above).
%   traces_xlim - numeric 2 elements. Time span relative to stimulus onset
%                 for displaying the trace [ms]. From 1 dt before first
%                 stimulus, to end of pstim_win (30 ms after last stimulus).
%   baseline    - numeric row vector. All the samples in baseline window -
%                 until 1 dt before first stimulus.
%   Tstamps     - numeric column vector. Time stamps matching the protocol,
%                 centred on the first stimulus - its time stamp is 0 [ms].
%
%   if protocol_id is "display", output will be a cell vector with all
%   existing protocols ID.
%
% CALL:
%   protocol_info = get_protocol("protocol_id",<string
%   scalar>,"fs",<numeric scalar>,"dt",<numeric non-negative scalar>)
%
%   OR (for saved protocols display)
%
%   protocol_list = get_protocol("protocol_id","display")

% Note1: add validation that protocol has needed fields & fill user-defined
%        fields if needed when loading, in order to let user make his own
%        protocols in a diffrent function and load them using get_protocol?
%        seem unnessary.
% Note2: convert so "custom" is not nessary - simply if user gave non saved
%        protocol id, start creating it? easy way for duplicated protocols
%        due to spelling errors...
% Note3: add Tstamps when calling using fs !!!!Done!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('protocol_id',   [], @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('fs',            [], @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('dt',            2,  @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))

p.parse(varargin{:});

protocol_id     = p.Results.protocol_id;
fs              = p.Results.fs;
dt              = p.Results.dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load or create protocol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get existing protocols form pre-defined folder inside package. We get the
% path to the package from the path to "get_protocol" file.
function_path = mfilename('fullpath');
[function_folder, ~] = fileparts(function_path);
[package_folder, ~] = fileparts(function_folder);
saved_protocols = dir(fullfile(package_folder,'stimProtocols','*.mat'));
[~,IDs] = cellfun(@fileparts,{saved_protocols.name},'UniformOutput',false);

% force IDs to be cell (if there is only 1 protocol, it will be char instead)
if ~iscell(IDs)
    IDs = {IDs};
end

switch protocol_id
    case "custom"
        % create new protocol from user inputs
        
        % collect protocol info from user
        AllPrompets = {sprintf('Protocol ID (leave empty to avoid saving protocol)\n(can''t be reserved words "custom" or "display"):'),...
            sprintf('Stimulation times [ms] (relative to recording start)\n(use spaces to separate the times):'),...
            'Stimulus length [us]:',...
            'Recording length [ms]'};
        user_inputs = inputdlg(AllPrompets,'Make Custom Protocol');
        protocol_info = struct();
        protocol_info.protocol_id = user_inputs{1};
        protocol_info.stim_times = sort(str2double(split(user_inputs{2})));%[ms] %always convert to col vec
        protocol_info.stim_length = str2double(user_inputs{3}); %[us]
        protocol_info.rec_length = str2double(user_inputs{4}); %[ms]
        
        % validate numeric data - if NaN they were not parse right, they
        % cannot be negative
        numeric_fields = {'stim_times','stim_length','rec_length'};
        for iField = 1:3
            if any(isnan(protocol_info.(numeric_fields{iField}))) || any(protocol_info.(numeric_fields{iField}) < 0)
                error('Error parsing the value in field %s',numeric_fields{iField})
            end
        end
        
        % create general helper fields from user inputs:
        % number of stimuli given
        protocol_info.nStim = numel(protocol_info.stim_times);
        % post stimulation windows
        stim_ends = protocol_info.stim_times+protocol_info.stim_length/1000;
        if protocol_info.nStim > 1
            % from each stimulus end to the start of the next one. After
            % the last one, give a window of 30 ms (arbitrary but usally enough)
            next_stim_starts = protocol_info.stim_times(2:end);
            next_stim_starts(end+1) = stim_ends(end) + min(30,protocol_info.rec_length);
        else
            % give an window of 30 ms after stimulus ends (arbitrary but usally enough)
            next_stim_starts = stim_ends + min(30,protocol_info.rec_length);
        end
        protocol_info.pstim_win = [stim_ends next_stim_starts];     
        
        % save protocol in the expected folder. Make it if needed
        if ~isempty(protocol_info.protocol_id)
            % prevent naming protocols "custom" or "display" - they won't load
            if ismember(protocol_info.protocol_id,{'custom','display'})
                error('Protocol ID cannot be reserved words "custom" or "display"')
            end
            
            % warn about overwriting existing protocols
            if ismember(protocol_info.protocol_id,IDs)
                answer = questdlg({sprintf('Protocol id %s already exist.',protocol_info.protocol_id)...
                    'Do you wish to overwrite it, or cancel?',...
                    'This cannot be undone'},...
                    'Overwrite?','Overwrite','Cancel','Cancel');
                if isempty(answer) || strcmp(answer,'Cancel Save')
                    error(['User canceled due to existing protocol with the same name {%s}.',...
                        ' Use get_protocol("protocol_id","display") to see all existing protocols ID'],...
                        protocol_info.protocol_id)
                end
            end
            
            % make dir if needed and save protocol
            if ~exist(fullfile(package_folder,'stimProtocols'),'dir')
                mkdir(fullfile(package_folder,'stimProtocols'))
            end
            save(fullfile(package_folder,'stimProtocols',[protocol_info.protocol_id '.mat']),'protocol_info')
        end
        
        
    case "display"
        % return existing protocols ID as output
        protocol_info = IDs;
        return
        
        
    otherwise
        % load existing protocol. No need to varify it - it was checked on creation
        if ~ismember(protocol_id,IDs)
            error(['Cannot find protocol ID %s. Known protocols ID are: {%s}. ',...
                'Use protocol_id "custom" to add a new protocol'],protocol_id,strjoin(IDs,','))
        end
        filename = fullfile(package_folder,'stimProtocols',join([protocol_id '.mat'],''));
        load(filename,'protocol_info')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert from times to sample numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initaite empty fields
protocol_info.response      = struct('base',[],'win',[]);
protocol_info.traces_xlim   = [];
protocol_info.baseline      = [];
protocol_info.Tstamps       = [];

if ~isempty(fs) 
    % if user inserted fs (& dt), use them to convert from times to sample
    % numbers (that can be used as index) in important windows. Else leave
    % fields empty
    
    % define a window of expected response - 
    % sample numbers that are far by dt from each stimulus end / start
    protocol_info.response.base = round([protocol_info.pstim_win(:,1)+dt,protocol_info.pstim_win(:,2)-dt] * fs / 1000);
    
    % create a vector of sample numbers capturing all response windows
    for iStim = 1:protocol_info.nStim
        protocol_info.response.win = [protocol_info.response.win, protocol_info.response.base(iStim,1):protocol_info.response.base(iStim,2)];
    end
    
    % create window to display (xlim) - from 1 dt before first stimulus
    % starts until dt after last stimulus ends
%     protocol_info.traces_xlim = [-dt, protocol_info.pstim_win(end,2)-protocol_info.stim_times(1)];
    protocol_info.traces_xlim = [-dt, protocol_info.pstim_win(end,2)];
    
    % get all the sample numbers for the baseline - 1 to dt before first stimulation
    protocol_info.baseline = 1:round((protocol_info.stim_times(1)-dt)* fs / 1000);
    
    % create time stamps
    protocol_info.Tstamps = -protocol_info.stim_times(1) :...
        (1 / fs) * 1000 : (protocol_info.rec_length - protocol_info.stim_times(1));
    protocol_info.Tstamps = protocol_info.Tstamps(1 : end - 1)';
end

end

% EOF