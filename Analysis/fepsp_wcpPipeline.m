function lfp = fepsp_wcpPipeline(varargin)
% Read wcp files and pass them through the whole fepsp pipline.
% 
% INPUT
%   basepath        char. dir with .wcp files.
%                   Defualt: pwd.
%   wcpfiles        wcp filenames to analyze. 
%                   string vector or a cell array of chars.
%                   Defualt: pick by uigetfile.
%   fepsp_protocol  char. see 'fepsp_getProtocol.m'.
%                   Defualt: 'io'.
%   intens          numeric vector. intensity values of stimulations [uA]
%                   Defualt: 1 to number of files.
%   out_name        Text scalar, suffix to add to saved lfp files & figures saved.
%                   Defualt: ''.
%   fsOut           numeric. requested sampling frequency. if empty will not downsample.
%                   Defualt: 5000.
%   cf              numeric, cut off frequency for sync filter (will filter
%                   each trace sepertly, may be bad filtering due to very short data...).
%                   If empty won't filter.
%                   Defualt: [].
%   max_jitter      numeric, see 'fepsp_markings.m'.
%                   Defualt: 0.5.
%   plot_summary    logical, to plot & save summary plot or not. Does not
%                   support stability. See 'fepsp_summaryPlot.m'.
%                   Defualt: true.
%   
% 
% DEPENDENCIES
%   getLFP
%   specBand
%   slutsky_fepsp (package)
%   import_wcp (external)
%   iosr.dsp
% 
% 14 Jun 22 LdM (mostly cut & rearrange of legacy function)
% 17 mar 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'wcpfiles', []);
addOptional(p, 'fepsp_protocol', 'io', @ischar);
addOptional(p, 'intens', [], @isnumeric);
addOptional(p, 'out_name', '', @mustBeTextScalar);
addOptional(p, 'fsOut', 5000, @isnumeric);
addOptional(p, 'cf', [], @isnumeric);
addOptional(p, 'max_jitter', 0.5, @isnumeric);
addOptional(p, 'plot_summary', true, @(x) validateattributes(x, {'logical','numeric'}, {'binary','scalar'}))

parse(p, varargin{:})
basepath        = p.Results.basepath;
wcpfiles        = p.Results.wcpfiles;
fepsp_protocol  = p.Results.fepsp_protocol;
intens          = p.Results.intens;
out_name        = p.Results.out_name;
fsOut           = p.Results.fsOut;
cf              = p.Results.cf;
max_jitter      = p.Results.max_jitter;
plot_summary    = p.Results.plot_summary;

if isempty(wcpfiles)
    [files_names,files_paths] = uigetfile(join([basepath,filesep,'*.wcp'],''));
    wcpfiles = fullfile(files_paths,files_names);
end

if isempty(intens)
    intens = 1 : length(wcpfiles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize protocol info
protocol_info = fepsp_getProtocol('protocol_id', fepsp_protocol);
stim_times = protocol_info.stim_times / 1000;

% initialize
nfiles = length(wcpfiles);
cntdata = [];
stim_locs = cell(1, nfiles);
stim_start = 1;
filelength = nan(1, nfiles);

% load and organize data
for ifile = 1 : nfiles
    
    % remove extention if any
    if contains(wcpfiles{ifile},'.wcp')
        basename = extractBefore(wcpfiles{ifile},'.wcp');
    else
        basename = wcpfiles{ifile};
    end

    % load lfp
    lfp = getLFP('basepath', basepath, 'basename', basename,...
        'ch', 1, 'chavg', {},...
        'fs', fsOut, 'interval', [0 inf], 'extension', 'wcp',...
        'savevar', false, 'forceL', true, 'cf', cf);
    fs = lfp.fs;
    [nsamps, ntraces] = size(lfp.data);

    % organize according to protocol
    switch fepsp_protocol
        
        case 'freerun'
            % remove incomplete data from last tract
            tmpdata = lfp.data(:);
            rmidx = find(movmax(diff(tmpdata), [0, 100]) == 0);
            if max(diff(rmidx)) > 1 
                warning('check')
            end
            if ~isempty(rmidx)
                tmpdata(rmidx(1) : end) = [];
            end
            cntdata = [cntdata; tmpdata];
            filelength(ifile) = length(tmpdata);            

        otherwise
            % cat data and create stim indices
            cntdata = [cntdata; lfp.data(:)];
            stim_locs{ifile} = stim_start + [stim_times(1) * fs :...
                size(lfp.data, 1) : length(lfp.data(:))];
            stim_start = length(cntdata);
    end
end

% find duplicated intens, and merge them
[unique_intens, ~, ic] = unique(intens,'stable');
idx_2_del = [];
for iDup = unique(ic)'
    dup_idx = find(ic == iDup);
    if numel(dup_idx) > 1
        fprintf('Merging [%s] due to same intensity\n',join(wcpfiles(dup_idx),', '))
        stim_locs{dup_idx(1)} = [stim_locs{dup_idx}];
        wcpfiles(dup_idx(1)) = join(wcpfiles(dup_idx)," + ");
        idx_2_del = [idx_2_del dup_idx(2:end)]; %#ok<AGROW> Very Small growth
    end
end
stim_locs(idx_2_del) = [];
wcpfiles(idx_2_del) = [];
intens = unique_intens;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize lfp struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add info to lfp struct
lfp.files       = wcpfiles;
lfp.protocol_id = fepsp_protocol;
lfp.intens      = intens;
lfp.stim_locs   = stim_locs;
lfp.cf          = cf;
lfp.max_jitter  = max_jitter;
lfp.slope_area  = [0.2 0.8];
lfp.data_in     = cntdata;

% save
recdir = fullfile(basepath);
warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir(recdir);
warning('on','MATLAB:MKDIR:DirectoryExists')
save(fullfile(recdir, [out_name, '.lfp.mat']), 'lfp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% continue processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orginaze traces
lfp.traces  = fepsp_org2traces(lfp,"graphics",false,"save_var",false);

% mark points
marking_win = fepsp_markings(lfp);
% load markings and updated traces
waitfor(marking_win)
[~, basename] = fileparts(basepath);
marking_files = fullfile(basepath,[basename '_fepsp_markings.mat']);
load(marking_files, "markings");
lfp.markings = markings;
load(marking_files, "traces");
lfp.traces_pre_del  = lfp.traces;
lfp.traces  = traces;
delete(marking_files)

% calculate results from points
lfp.results = fepsp_analyse(lfp,"save_var",false);

% save
recdir = fullfile(basepath);
warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir(recdir);
warning('on','MATLAB:MKDIR:DirectoryExists')
save(fullfile(recdir, [out_name, '.lfp.mat']), 'lfp')

% create summary plots
if plot_summary
    warning('off','MATLAB:MKDIR:DirectoryExists')
    fepsp_summaryPlot(lfp);
    warning('on','MATLAB:MKDIR:DirectoryExists')
    saved_name = fullfile(basepath, 'graphics', 'fepsp',sprintf('%s_fepsp_results.tif', basename));
    new_name   = fullfile(basepath, 'graphics', 'fepsp',sprintf('%s_%s_fepsp_results.tif', basename, out_name));
    movefile(saved_name,new_name);
end

% EOF