function [sumPlot] = fepsp_plot_avgTraces(varargin)
% Create a figure for each channel, with simply the average traces of each
% intensity
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - Numeric scalar. sampling frequency [Hz].
%   protocol_id - String or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "fepsp_getProtocol.m" for
%                 more info.
% INPUT (optional):
%   intens      - numeric vector of stimulus intensities used for each
%                 column of traces. Typically provided in units of uA. if
%                 empty will be set as a unit increasing vector.
%                 Default: [].
%   traces_xlim - numeric vector of 2 elements. Time span relative
%                 to stimulus onset for displaying the trace [ms]. For
%                 example, xLimit = [-1 30] will show the trace from 1 ms
%                 before the stimulus until 30 ms after the stimulus. If
%                 empty will use protocol default values. See
%                 "fepsp_getProtocol.m" for more info. 
%                 Default: [].
%   traces_ylim - numeric vector with 2 elements. voltage limits for
%                 displaying the trace [mV]. These y-limits remain constant
%                 throughout all intensities. If empty will be set
%                 according to the max range in each channel (excluding
%                 stimulus artifact). 
%                 Default: [].
%   dt          - non-negative scalar. Dead time between
%                 stimulus onset and earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See "fepsp_getProtocol.m" for more info.
%                 Default: 2.
%   saveFig     - logical. save figure in basepath {true}
%
% OUTPUT:
%   AvgPlot     - figure handles array. Handles to the averge traces figures, one
%                 for each channel.

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;
p.addParameter('traces',        [], @(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',            [], @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',   [], @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('intens',        [], @(x) (isnumeric(x) && isvector(x)) || isempty(x))
p.addParameter('traces_xlim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('traces_ylim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('dt',            2,  @(x) validateattributes(x,{'numeric'},{'vector','nonnegative'}))
p.addParameter('saveFig',       true, @islogical)

parse(p, varargin{:})

traces          = p.Results.traces;
fs              = p.Results.fs;
protocol_id     = p.Results.protocol_id;
intens          = p.Results.intens;
traces_xlim     = sort(p.Results.traces_xlim);
traces_ylim     = sort(p.Results.traces_ylim);
dt              = p.Results.dt;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate intens
if isempty(intens)
    intens = -(size(traces, 2): -1 : 1);
elseif numel(intens) ~= size(traces, 2)
    error('discripancy between the length of intens and the size of traces')
end

% protocol info
protocol_info = fepsp_getProtocol("protocol_id",protocol_id,"fs",fs,"dt",dt);

% xlim for plot
if ~isempty(traces_xlim)
    protocol_info.traces_xlim = traces_xlim;
end

% stim params 
nIntens = length(intens);
[intens_sorted, intens_order] = sort(intens, 'ascend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iChan = size(traces, 1) : -1 : 1
 
    % open figure 
    sumPlot(iChan) = figure();
    sumPlot(iChan) = fepsp_graphics(sumPlot(iChan));          % set graphics
    title(sprintf('Channel %d - %s', iChan,upper(protocol_info.protocol_id)),...
        'FontSize', 28, 'FontWeight', 'bold', 'FontName', 'FixedWidth','Interpreter','none')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot average traces for all intensities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % avg trace for each intensity
    traces_avg = cell2mat(cellfun(@(x) mean(x, 2, 'omitnan'),...
        traces(iChan, :), 'UniformOutput', false));
    
    % plot and change data tips
    if nIntens > 1
%         colororder(green_magenta_Cmap(nIntens))
%         colororder(distinguishable_colors(nIntens))
        base_colors = parula(nIntens);
        base_colors = rgb2hsv(base_colors);
        no_more_down = base_colors(:,3) <= 0.3;
        base_colors(~no_more_down,3) = base_colors(~no_more_down,3) - 0.2;
        base_colors = hsv2rgb(base_colors);
        colororder(base_colors)
    end
    traces_avg_h = plot(protocol_info.Tstamps, traces_avg, 'LineWidth', 2);
    for iIntens = 1 : nIntens
        traces_avg_h(iIntens).Tag = num2str(intens(iIntens));
        r = dataTipTextRow('intensity', ones(size(traces_avg_h(iIntens).XData)) * intens(iIntens));
        traces_avg_h(iIntens).DataTipTemplate.DataTipRows(end+1) = r;
    end

    % axes limits
    xlim(protocol_info.traces_xlim)
    if isempty(traces_ylim)
        ylim([min(traces_avg(protocol_info.response.win, :), [], 'all')...
            max(traces_avg(protocol_info.response.win,:), [], 'all')] .* 1.1)
    else
        ylim(traces_ylim)
    end

    % legend
    legend(traces_avg_h(intens_order),string(intens_sorted),'Location','southeast','NumColumns', 2)

    % labels
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    title('Average Traces')
    fepsp_graphics(gca);          % set graphics
end

% save
if saveFig
    fprintf(['This is the time to make edits (such as moving the legends) before figure export.'...
        '\ncontinure running the code (by F5) ones you are ready to export!\n'])
    keyboard
    basepath = pwd;
    [~, basename] = fileparts(basepath);
    figpath = fullfile(basepath, 'graphics', 'fepsp');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_traces_only_fepsp_results', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF