function [results]= fepsp_analyse(varargin)

% Calculate amplitude and slope for each trace (repetition) according to
% the start and peak markings of each response. the area of the slop can be
% determined by the user (slope_area). 
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   markings    - struct. see "epsp_markings.m.
%   fs          - numeric scalar. sampling frequency [Hz].
%   protocol_id - string or char. ID of stimulation protocol.
%                 e.g. "io", "stp" or "custom". See "fepsp_getProtocol.m".

% INPUT (optional):
%   slope_area  - optional. Numeric 2 elements between 0-1. Part of
%                 amplitude to measure the slope between. For example, if
%                 slope_area is [0.2 0.9] then slope will be measured in
%                 the response between the points that match 20% & 90% of
%                 the amplitude. Order does not matter - slope between 20%
%                 & 90% of amplitude is the same as slope between 90% & 20%.
%                 Default: [0.2 0.8]
%   base_path -   string or char. Full path to where the output should be
%                 saved. The name of the last folder in base_path will be
%                 the prefix for all saved data. e.g. base_path =
%                 'lh85_211012_132000' than the output of fepsp_markings
%                 will be 'lh85_211012_132000_fepsp_markings'.
%                 Default: pwd.
%   save_var    - logical flag. Save output in base_path as
%                 mat file named "<base_name>_fepsp_results".
%                 Default: true
%
% OUTPUT:
%   results     - struct scalar. Holds all of the results & analysis info.
%                 Has fields:
%       slope_area - the same as input slope_area, see above. Between what
%                    parts of response amplitude slope was measured.
%       all_traces - struct. Analysis results for all the traces in
%                    every trace group. Has fields:
%           Amp         - standard cell (see README). Each value is the
%                         amplitude of matching response.
%           Slope       - standard cell (see README). Each value is the
%                         slope of matching response.
%           slope_err   - 2d cell array of channel (row) x intensity
%                         (column), with each cell containing a 2d struct
%                         array of stimulus number (row) x repetition (column).
%                         Error estimation structure generated by polyfit,
%                         for each calculated slope. 
%                         See "polyfit.m" for more info.
%           slope_win -   2d cell array of channel (row) x intensity
%                         (column), with each cell containing a 3d numeric
%                         mat, slope_area element (dim 1, always size 2) x
%                         stimulus number (dim 2) x repetition (dim 3).
%                         The sample number that define the window,
%                         corresponding element in slope_area. For example
%                         if slope are is [0.5 0.1], the first element in
%                         slope_window is the sample that match 50% of
%                         amplitude inside of response, and the second is
%                         the sample number that match 10%.
%       avg_traces - struct. Analysis results for the average traces
%                    matching each traces group. Has the same fields as
%                    all_traces, but inside each cell array the size of the
%                    dimension of "repetition" is 1 (as there is only 1
%                    average trace). 
%                    Each field have the same measurements, but for the
%                    average trace of each traces group.
%
% CALL:
%   results = fepsp_analyse("traces", <cell array>, "fs", <numeric scalar>,
%   "protocol_id", <string scalar>, "markings", <struct scalar>,
%   "base_folder", <folder path>, "save_var", <logical flag>, "slope_area",
%   <numeric 2 elements 0<= & >=1>)


% Note1: Change polyfit to fit. Problem - too few points case.
% Note2: Path conflict with the fepsp_analysis in slutskycode.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('traces',        [],             @(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',            [],             @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',   [],             @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('markings',      [],             @(x) validateattributes(x,{'struct'},{'scalar'}))
p.addParameter('base_path',     pwd,            @isfolder)
p.addParameter('save_var',      true,           @(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))
p.addParameter('slope_area',    [0.20 0.80],    @(x) validateattributes(x,{'numeric'},{'numel',2,'>=',0,'<=',1}))

p.parse(varargin{:})

traces          = p.Results.traces;
fs              = p.Results.fs;
protocol_id     = p.Results.protocol_id;
markings        = p.Results.markings;
base_path       = p.Results.base_path;
save_var        = p.Results.save_var;
slope_area      = p.Results.slope_area(:)'; % make sure it is row vec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get protocol info
protocol_info = fepsp_getProtocol("protocol_id",protocol_id,"fs",fs);
% we only need the Tstamps from the parameters that require fs, so dt does
% not have any effect - we can use the default

% all traces outputs:
[Amp,Slope]         = deal(cellfun(@(x) nan(protocol_info.nStim,size(x,2)),traces,'UniformOutput',0));
fit_err_struct      = struct('R',[],'df',[],'normr',[]);
slope_err           = cellfun(@(x) repmat(fit_err_struct,protocol_info.nStim,size(x,2)),traces,'UniformOutput',0);
slope_win           = cellfun(@(x) nan(2,protocol_info.nStim,size(x,2)),traces,'UniformOutput',0);

% average traces outputs:
[Amp_avg,slope_avg] = deal(repmat({nan(protocol_info.nStim,1)},size(traces)));
slope_err_avg       = repmat({repmat(fit_err_struct,protocol_info.nStim,1)},size(traces));
slope_win_avg       = repmat({nan(2,protocol_info.nStim)},size(traces));

% warning msg base, instead of "Polynomial not unique"
warn_base = ['As a results, relatively big fit error may occur. '...
    'This may be due to start & peak marker too close, very small response amplitude, or a biphasic event between markers. ' ...
    'If slope is expected to be close to 0 in this trace, you can ignore. '...
    'Check matching slope_err at the for more info about the fit error'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% polynomial not unique might occur when using polyfit. As we are using
% polyfit multiple times, it will be uninformative and clutter the command
% window. We will silence it, and present a clear warning instead when it
% rises (see later)
warning('off','MATLAB:polyfit:PolyNotUnique')

for iTraces_group = numel(traces):-1:1 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract data of all the traces of this trace group
    loop_starts = markings.starts{iTraces_group};
    loop_peaks = markings.peaks{iTraces_group};
    loop_traces_group = traces{iTraces_group};
    % fix deleted traces to NaN. Identify by columns that are NaN in markings
    deleted_traces = all(isnan(loop_starts),1) | all(isnan(loop_peaks),1);
    loop_traces_group(:,deleted_traces) = nan;
    % convert the marks to 1 so they can be used as index. Because
    % loop_traces_group for traces that were deleted are now NaN, amplitude
    % will also be NaN.
    loop_starts(:,all(isnan(loop_traces_group),1)) = 1;
    loop_peaks(:,all(isnan(loop_traces_group),1)) = 1;
    
    % extract data of the average trace of the same traces group
    [iChan,iInten] = ind2sub(size(traces),iTraces_group);
    loop_avg_trace = mean(loop_traces_group,2,'omitnan');
    loop_starts_avg = markings.starts_avg{iTraces_group};
    loop_peaks_avg = markings.peaks_avg{iTraces_group};
    
    % calculate amplitude for all the traces:
    % convert indices from subscripts to linear, in order to maintain 1
    % index refer 1 sample.
    linear_starts = sub2ind(size(loop_traces_group),loop_starts,repmat(1:size(loop_traces_group,2),protocol_info.nStim,1));
    linear_ends = sub2ind(size(loop_traces_group),loop_peaks,repmat(1:size(loop_traces_group,2),protocol_info.nStim,1));
    % amplitude is simply absolute of: start of response minus end of response
    Amp{iTraces_group} = abs(loop_traces_group(linear_starts) - loop_traces_group(linear_ends));
    
    % redo for average trace of the same traces group
    Amp_avg{iTraces_group} = abs(loop_avg_trace(loop_starts_avg) - loop_avg_trace(loop_peaks_avg));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate slope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for iStim = protocol_info.nStim:-1:1
        for iTrace = size(loop_traces_group,2):-1:1
            % skip deleted traces - their amplitude will be NaN
            if isnan(Amp{iTraces_group}(iStim,iTrace))
                continue
            end
            
            % create a window between response start & peak
            s2p_win = loop_starts(iStim,iTrace):loop_peaks(iStim,iTrace);
            
            % calculate Slope:
            [slope_win{iTraces_group}(:,iStim,iTrace),Slope{iTraces_group}(iStim,iTrace),slope_err{iTraces_group}(iStim,iTrace)] = ...
                calculate_slope(protocol_info.Tstamps,loop_traces_group(:,iTrace),s2p_win,Amp{iTraces_group}(iStim,iTrace),slope_area);
            
            % there are a few cases when there are not enough points
            % between percentiles - usually when response amplitude is very
            % small. We silenced those warnings earlier, and we will
            % produce informative warning instead
            if strcmp(lastwarn,'Polynomial is not unique; degree >= number of data points.')
                lastwarn('') % clear lastwarn to prevent warning reappear next loop
                warning(['Its seems there are too few points when searching slope, in channel %d, intensity %d, stimulus number %d in trace %d. '...
                    warn_base],iChan,iInten,iStim,iTrace)
            end
        end
        
        % redo for average trace of the same traces group:
        
        % create a window
        s2p_win = loop_starts_avg(iStim):loop_peaks_avg(iStim);
        
        % Calculate Slope:
        [slope_win_avg{iTraces_group}(:,iStim),slope_avg{iTraces_group}(iStim),slope_err_avg{iTraces_group}(iStim)] = ...
            calculate_slope(protocol_info.Tstamps,loop_avg_trace,s2p_win,Amp_avg{iTraces_group}(iStim),slope_area);
        
        % Give warning if needed
        if strcmp(lastwarn,'Polynomial is not unique; degree >= number of data points.')
            lastwarn('')
            warning(['Its seems there are too few points when searching slope, in channel %d, intensity %d, stimulus number %d in the average trace. '...
                warn_base],iChan,iInten,iStim)
        end
    end
end

% restore the polyfit warning
warning('on','MATLAB:polyfit:PolyNotUnique')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize outputs & save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pack all outputs in struct
results.all_traces = struct('Amp',{Amp},'Slope', {Slope},'slope_err',{slope_err},'slope_win',{slope_win});
results.avg_traces = struct('Amp',{Amp_avg},'Slope',{slope_avg},'slope_err',{slope_err_avg},'slope_win',{slope_win_avg});
results.slope_area = slope_area;

% info
results.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
results.info.protocol_id = protocol_id;

% save
if save_var
    % create file name from base_folder
    [~,base_name] = fileparts(base_path);
    results_file = [base_path filesep base_name '_fepsp_results.mat'];
    save(results_file,'results')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [slope_win,Slope,slope_err] = calculate_slope(Tstamps,trace,s2p_win,Amp,slope_area)
% Calculate slope in a single trace between 2 point inside of response_win
% that match slope_area part of amplitude.
% Algorithm: Calculate all of the amplitudes of the points in
% s2p_win (absolute difference from first point in s2p_win). Find
% the points whose amplitude value is the closest to the part of Amp
% specified by slope_area. Fit a line to all the points between the two
% points found - take this line slope.
%
% INPUT
%   Tstamps     - numeric vector. Time [ms] that match the trace.
%   trace       - numeric vector. Single data sample repetition in [mV].
%   s2p_win     - positive ascending integer vector. Samples number between start & peak markers
%   Amp         - numeric scalar, calculated response amplitude for this trace [mV].
%   slope_area  - Numeric 2 elements row vectoe between 0 & 1. Part of
%                 amplitude to measure the slope between.  
%
% OUTPUT
%   slope_win   - Numeric 2 elements positive integer.
%                 The sample number that define the window. Each element
%                 matches the corresponding element in slope_area. For
%                 example if slope are is [0.5 0.1], the first element in
%                 slope_window is the sample that match 50% of amplitude
%                 inside of response, and the second is the sample number
%                 of 10%.
%   Slope       - Numeric scalar, calculated response slope for this trace [mV/ms].
%   slope_err   - Struct scalar. Error estimation structure generated by
%                 polyfit. See polyfit for more info. 

% find all the responses points amplitude - absolute difference from
% response start data sample
points_amp = abs(trace(s2p_win(1)) - trace(s2p_win));

% create a vector of "how much are slope_area's percentiles of calculated
% amplitude?"
p_vec = Amp*slope_area;

% using implicit expansion, find data samples with the closest
% amplitude to wanted percentile of general amplitude
[~,slope_win] = min(abs(points_amp-p_vec),[],1);

% convert slope_window to sample numbers
slope_win = slope_win + s2p_win(1) - 1;

% if there is a significant biphasic event inside the
% response, the lower percentage might be found after the
% higher one. Will break colon operator later. So sort them.
slope_win_sorted = sort(slope_win,'ascend');
% notice that we keep output unsorted, so the first the 1st element in
% slope_window matching the first percentile in slope_area

% calculate Slope by fitting a line. Save results & fit error to output:
[fit_params,slope_err] = ...
    polyfit(Tstamps(slope_win_sorted(1):slope_win_sorted(2)), trace(slope_win_sorted(1):slope_win_sorted(2)), 1);

% polyfit first output is in descending power order. We don't
% need the intercept (power == 0) value.
Slope = fit_params(1);
end

% EOF