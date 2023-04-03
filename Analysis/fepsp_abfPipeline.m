function lfp = fepsp_abfPipeline(varargin)
% load and analyse abf files in the fepsp pipeline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'abffiles', []);
addOptional(p, 'fepsp_protocol', 'io', @ischar);
addOptional(p, 'intens', [], @isnumeric);
addOptional(p, 'out_name', '', @mustBeTextScalar);
%addOptional(p, 'fsOut', 5000, @isnumeric);
addOptional(p, 'Ch', 1, @isnumeric)
addOptional(p, 'overwrite_existing', 0, @(x) validateattributes(x, {'numeric'}, {'scalar'}))
addOptional(p, 'plot_summary', true, @(x) validateattributes(x, {'logical','numeric'}, {'binary','scalar'}))
addOptional(p, 'add_fields',cell.empty([0,2]),@(x) validateattributes(x,{'cell'},{'ncols',2}))
p.parse(varargin{:})

basepath = p.Results.basepath;
abffiles = string(p.Results.abffiles);
fepsp_protocol = p.Results.fepsp_protocol;
intens = p.Results.intens;
out_name = p.Results.out_name;
% fsOut = p.Results.fsOut;
Ch = p.Results.Ch;
overwrite_existing = p.Results.overwrite_existing;
plot_summary = p.Results.plot_summary;
add_fields = p.Results.add_fields;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

for iFile = numel(abffiles):-1:1
    [tmp_d,si(iFile)] = abfload2(abffiles(iFile));
    for iCh = size(tmp_d,2):-1:1
        traces{iCh,iFile} = squeeze(tmp_d(:,iCh,:));
    end
end
traces = traces(Ch,:);

if ~all(si(1) == si)
    error('Files have diffrent si from each other')
end

% find duplicated intens, and merge them
[unique_intens, ~, ic] = unique(intens,'stable');
idx_2_del = [];
for iDup = unique(ic)'
    dup_idx = find(ic == iDup);
    if numel(dup_idx) > 1
        fprintf('Merging [%s] due to same intensity\n',string(join(abffiles(dup_idx),', ')))
        abffiles(dup_idx(1)) = join(abffiles(dup_idx)," + ");
        for iCh = 1 : size(traces,1)
            traces{iCh,dup_idx(1)} = [traces{iCh,dup_idx}];
        end
        idx_2_del = [idx_2_del dup_idx(2:end)]; %#ok<AGROW> Very Small growth
    end
end
abffiles(idx_2_del) = [];
traces(:,idx_2_del) = [];
intens = unique_intens;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize lfp struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add info to lfp struct
lfp.files       = abffiles;
lfp.protocol_id = fepsp_protocol;
lfp.intens      = intens;
% lfp.stim_locs   = stim_locs;
% lfp.cf          = cf;
lfp.traces      = traces;
lfp.title       = out_name;
lfp.fs          = 1 / (si(1)/1e6); % us to ms
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
% lfp.traces  = fepsp_org2traces(lfp,"graphics",false,"save_var",false);

% detrend
% build regressor to find the trend of a trace using the first and last 5
% ms. This is to exclude the evoked response and stimulus artifact from the
% linear fit. Based on matrix form of polynomial regression.
lip = ceil(5 * lfp.fs / 1000);
reg_len = 2 * lip;
trace_len = length(lfp.traces{1}(:,1));

use_idx = [1 : lip, trace_len : -1 : trace_len - lip + 1];
w = [0 : reg_len - 1] / (reg_len - 1);
w = [w', ones(length(w), 1)];
[q, r] = qr(w, 0);
w_full = [0 : trace_len - 1] / (trace_len - 1);
w_full = [reshape(w_full, trace_len, []), ones(trace_len, 1)];

for iIntens = 1 : size(lfp.traces,2)
    for iCh = 1 : size(lfp.traces,1)
        trend_coeff = (r \ q' * lfp.traces{iCh, iIntens}(use_idx, :));
        trend{iCh, iIntens} = w_full * trend_coeff;
        lfp.traces{iCh, iIntens} = lfp.traces{iCh, iIntens} - trend{iCh, iIntens};
    end
end


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
protocol_info = fepsp_getProtocol("protocol_id",lfp.protocol_id);

% fit models to responce in case of IO
% if protocol_info.nStim == 1 && numel(lfp.intens) > 1
for iCh = 1 : size(lfp.traces,1)
    for iStim = 1:protocol_info.nStim
        lfp.Amp_mdl{iCh}{iStim}   = fepsp_fitIO(lfp,"resp_type","Amp","stim2work",iStim);
        lfp.Slope_mdl{iCh}{iStim} = fepsp_fitIO(lfp,"resp_type","Slope","stim2work",iStim);
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