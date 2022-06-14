function Legacy_fepsp_wcpPipeline(varargin)

% organizes recordings from wcp. assumes all .wcp files are in basepath and
% are named with unique numbers as suffixes (e.g. "xxx_034"). user selects
% specific files by specifying the last 2-3 digits of the filename. the
% function load the data, organizes it in a separate folder and calls the
% slutsky_fepsp functions.
% 
% INPUT
%   basepath        char. dir with .wcp files
%   wcpfiles        wcp filenames to analyze. can be numeric (last digits)
%                   or a cell array of chars with the filenames not including
%                   the extension .wcp. for example {'xxx.034'}
%   fepsp_protocol  char. can be 'freerun', 'io', or 'stp'
%   recname         char. name of output folder
%   fsOut           numeric. requested sampling frequency. if empty will
%                   not downsample
%   intens          numeric. intensity values of stimulations [uA]
% 
% DEPENDENCIES
%   getLFP
%   specBand
%   slutsky_fepsp (package)
%   import_wcp (external)
%   iosr.dsp
% 
% 17 mar 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'wcpfiles', []);
addOptional(p, 'fepsp_protocol', 'io', @ischar);
addOptional(p, 'recname', '', @ischar);
addOptional(p, 'fsOut', 1250, @isnumeric);
addOptional(p, 'intens', [], @isnumeric);

parse(p, varargin{:})
basepath        = p.Results.basepath;
wcpfiles        = p.Results.wcpfiles;
fepsp_protocol  = p.Results.fepsp_protocol;
recname         = p.Results.recname;
fsOut           = p.Results.fsOut;
intens          = p.Results.intens;

if isempty(intens)
    intens = 1 : length(wcpfiles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find requested .wcp files
cd(basepath)
rawfiles = dir('*.wcp');
rawnames = {rawfiles.name};
rawsuffix = cellfun(@(x) x(end - 6 : end - 4) + ".wcp", rawnames, 'uni', false);

% get all .wcp files in folder if non selected
if isempty(wcpfiles)
    wcpfiles = rawnames;
end
nfiles = length(wcpfiles);

% convert selected numbers to wcp filenames
if isnumeric(wcpfiles)
    for ifile = 1 : nfiles
        fileidx = contains(rawnames, rawsuffix{wcpfiles(ifile)});
        [~, basename{ifile}] = fileparts(rawfiles(fileidx).name);
    end
else
    basename = wcpfiles;
end

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

    % load lfp
    lfp = getLFP('basepath', basepath, 'basename', basename{ifile},...
        'ch', 1, 'chavg', {},...
        'fs', [], 'interval', [0 inf], 'extension', 'wcp',...
        'savevar', false, 'forceL', true, 'cf', []);
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

        case {'io', 'stp', 'stp5'}
            % cat data and create stim indices
            cntdata = [cntdata; lfp.data(:)];
            stim_locs{ifile} = stim_start + [stim_times(1) * fs :...
                size(lfp.data, 1) : length(lfp.data(:))];
            stim_start = length(cntdata);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter and downsample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter params
import iosr.dsp.*       % import sync filter toolbox
cf = 450;               % low pass cutoff
fsIn = fs; 
if isempty(fsOut)
    fsOut = fsIn;
end
fsRatio = (fsIn / fsOut);
if cf > fsOut / 2
    warning('low pass cutoff beyond nyquist')
end
filtRatio = cf / (fsIn / 2);
ntbuff = 525;           % default filter size in iosr toolbox
if mod(ntbuff, fsRatio) ~= 0
    ntbuff = round(ntbuff + fsRatio - mod(ntbuff, fsRatio));
end

% do the filtering
lfp.data = [iosr.dsp.sincFilter(cntdata, filtRatio)]';

% downsample
if fsRatio ~= 1
    lfp.data = real(lfp.data(:, fsRatio : fsRatio :...
        length(lfp.data) - ntbuff));
end
lfp.fs = fsOut;
lfp.filelength = filelength / fsIn;       % [s]

% downsample stimulus indices
stim_locs = cellfun(@(x) round(x / fsRatio), stim_locs, 'uni', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize lfp struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add info to lfp struct
lfp.files = wcpfiles;
lfp.fepsp_protocol = fepsp_protocol;
lfp.intens = intens;

% save
recdir = fullfile(basepath, recname);
mkdir(recdir);
save(fullfile(recdir, [recname, '.lfp.mat']), 'lfp')

% debugging for cntdata
dbflag = true;
if dbflag
    fh = figure;
    plot([1 : length(cntdata)] / fsRatio, cntdata)
    hold on
    plot([[stim_locs{:}]; [stim_locs{:}]], ylim, '--k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% continue processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = recdir;
cd(basepath)
[~, basename] = fileparts(basepath);
fs = lfp.fs;

% fepsp --------------------------------------------------------------------

% organize in fepsp cell format
[traces] = fepsp_org2traces('data_in', lfp.data',...
    'basepath', basepath, 'fs', fs, 'protocol_id', fepsp_protocol,...
    'stim_locs', stim_locs);

% mark traces via gui
marking_win = fepsp_markings("traces", traces, "fs", fs,...
    "protocol_id", fepsp_protocol, "base_path", basepath,...
    "intens", intens, "traces_Xlimit", [], "traces_Ylimit", [],...
    "dt", 2, "max_jitter", 0.5, "fast_mark", false);

% load markings and updated traces
waitfor(marking_win)
load([basename, '_fepsp_markings.mat'], "markings")
load([basename, '_fepsp_traces.mat'], "traces")

% analyze traces according to the manual markings
results = fepsp_analyse("traces", traces, "fs", fs,...
    "protocol_id", fepsp_protocol, "markings", markings,...
    "base_path", basepath, "save_var", true, "slope_area", [0.2 0.9]);

% step 4        dispalys the results
analysed_fepsp = fepsp_summaryPlot("traces", traces, "fs", fs,...
    "protocol_id", fepsp_protocol, "markings", markings, "results", results,...
    "base_path", basepath, "intens", intens);


% freerun -----------------------------------------------------------------

if strcmp(fepsp_protocol, 'freerun')
    % spectrogram
    [s, tstamps, freq] = specBand('basepath', basepath, 'sig', lfp.data,...
        'fs', lfp.fs, 'graphics', true, 'logfreq', true);
    yLimit = ylim;
    hold on
    tidx = [cumsum(lfp.filelength) / 60 / 60]';
    tidx = [1; tidx(1 : end - 1)];
    plot([tidx, tidx], ylim, '--k', 'LineWidth', 2)
    hold on
    text(tidx, repmat(yLimit(2) + 0.01 * yLimit(2), length(tidx), 1), string(wcpfiles))
end

end

% EOF
