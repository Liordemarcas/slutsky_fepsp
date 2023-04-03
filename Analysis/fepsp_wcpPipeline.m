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
%   inverse_all     Logical flag. If true, inverse all input traces.
%                   Defualt: false.
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
addOptional(p, 'plot_summary', true, @(x) validateattributes(x, {'logical','numeric'}, {'binary','scalar'}))
addOptional(p, 'add_fields',cell.empty([0,2]),@(x) validateattributes(x,{'cell'},{'ncols',2}))
addOptional(p, 'inverse_all', false, @(x) validateattributes(x, {'logical','numeric'}, {'binary','scalar'}))
addOptional(p, 'overwrite_existing', 0, @(x) validateattributes(x, {'numeric'}, {'scalar'}))

parse(p, varargin{:})
basepath        = p.Results.basepath;
wcpfiles        = p.Results.wcpfiles;
fepsp_protocol  = p.Results.fepsp_protocol;
intens          = p.Results.intens;
out_name        = p.Results.out_name;
fsOut           = p.Results.fsOut;
cf              = p.Results.cf;
plot_summary    = p.Results.plot_summary;
add_fields      = p.Results.add_fields;
inverse_all     = p.Results.inverse_all;
overwrite_existing = p.Results.overwrite_existing;

if isempty(wcpfiles)
    [files_names,files_paths] = uigetfile(join([basepath,filesep,'*.wcp'],''),'MultiSelect','on');
    % return if no file selected
    if isnumeric(files_names)
        err = MException('fepsp_wcpPipeline:NoFileSelected','No file selected');
        throw(err)
    end

    wcpfiles = fullfile(files_paths,files_names);
    % fix for uigetfile return 1 file as char
    if ischar(wcpfiles)
        wcpfiles = {wcpfiles};
    end
end

if isempty(intens)
    intens = 1 : length(wcpfiles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if file exist, and if to overwrite 
% (overwrite_existing == 1 -> always, overwrite_existing == 0 -> ask, overwrite_existing == -1 -> never)
recdir = fullfile(basepath);
if exist(fullfile(recdir, join([out_name, '.lfp.mat'],'')),"file") && ~overwrite_existing
    answer = questdlg(sprintf('File named %s already exist. Overwrite it?',[out_name, '.lfp.mat']),'Overwrite','Yes','No','Yes');
    if isempty(answer) || strcmp(answer,'No')
        err = MException('fepsp_wcpPipeline:UserNotOverWrite','Aborting, file exist, user choose not to overwrite');
        throw(err)
    end
elseif exist(fullfile(recdir, join([out_name, '.lfp.mat'],'')),"file") && overwrite_existing == -1
    err = MException('fepsp_wcpPipeline:UserNotOverWrite','Aborting, file exist, user choose not to overwrite');
    throw(err)
end

% organize protocol info
protocol_info = fepsp_getProtocol('protocol_id', fepsp_protocol);
stim_times = protocol_info.stim_times / 1000;

% initialize
nfiles = length(wcpfiles);
cntdata = [];
traces_T = {};
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
    lfp_tmp = getLFP('basepath', basepath, 'basename', basename,...
        'ch', 1, 'chavg', {},...
        'fs', fsOut, 'interval', [0 inf], 'extension', 'wcp',...
        'savevar', false, 'forceL', true, 'cf', cf);
    fs = lfp_tmp.fs;
    [nsamps, ntraces] = size(lfp_tmp.data);

    % organize according to protocol
    switch fepsp_protocol
        
        case 'freerun'
            % remove incomplete data from last tract
            tmpdata = lfp_tmp.data(:);
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
            cntdata = [cntdata; lfp_tmp.data(:)];
            stim_locs{ifile} = stim_start + [stim_times(1) * fs :...
                size(lfp_tmp.data, 1) : length(lfp_tmp.data(:))];
            stim_start = length(cntdata);
            traces_T = [traces_T,lfp_tmp.trace_T];
    end
end

% find duplicated intens, and merge them
[unique_intens, ~, ic] = unique(intens,'stable');
idx_2_del = [];
for iDup = unique(ic)'
    dup_idx = find(ic == iDup);
    if numel(dup_idx) > 1
        fprintf('Merging [%s] due to same intensity\n',string(join(wcpfiles(dup_idx),', ')))
        stim_locs{dup_idx(1)} = [stim_locs{dup_idx}];
        wcpfiles(dup_idx(1)) = join(wcpfiles(dup_idx)," + ");
        idx_2_del = [idx_2_del dup_idx(2:end)]; %#ok<AGROW> Very Small growth
    end
end
stim_locs(idx_2_del) = [];
wcpfiles(idx_2_del) = [];
intens = unique_intens;

% inverse all traces if requested
if inverse_all
    cntdata = -cntdata;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize lfp struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add info to lfp struct
lfp.files       = wcpfiles;
lfp.fs          = lfp_tmp.fs;
lfp.fs_orig     = lfp_tmp.fs_orig;
lfp.protocol_id = fepsp_protocol;
lfp.intens      = intens;
lfp.stim_locs   = stim_locs;
lfp.cf          = cf;
lfp.data_in     = cntdata;
lfp.title       = out_name;
lfp.traces_T    = traces_T; % time each trace was recorded
for iField = 1:size(add_fields,1)
    lfp.(add_fields{iField,1}) = add_fields{iField,2};
end


% save
warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir(recdir);
warning('on','MATLAB:MKDIR:DirectoryExists')
save(fullfile(recdir, join([out_name, '.lfp.mat'],'')), 'lfp','-v7.3')


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

% fit models to responce in case of IO
for iCh = 1 : size(lfp.traces,1)
    for iStim = 1:protocol_info.nStim
        lfp.Amp_mdl{iCh}{iStim}   = fepsp_fitIO(lfp,"resp_type","Amp","stim2work",iStim,"Chan",iCh);
        lfp.Slope_mdl{iCh}{iStim} = fepsp_fitIO(lfp,"resp_type","Slope","stim2work",iStim,"Chan",iCh);
    end
end

% save
recdir = fullfile(basepath);
warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir(recdir);
warning('on','MATLAB:MKDIR:DirectoryExists')
save(fullfile(recdir, join([out_name, '.lfp.mat'],'')), 'lfp',"-v7.3")

% create summary plots
if plot_summary
    warning('off','MATLAB:MKDIR:DirectoryExists')
    fepsp_summaryPlot(lfp);
    warning('on','MATLAB:MKDIR:DirectoryExists')
    if ~isfield(lfp,'saveFig') || lfp.saveFig
        saved_name = fullfile(basepath, 'graphics', 'fepsp',sprintf('%s_fepsp_results.tif', basename));
        new_name   = fullfile(basepath, 'graphics', 'fepsp',sprintf('%s_%s_fepsp_results.tif', basename, out_name));
        movefile(saved_name,new_name);
    end
end

% EOF