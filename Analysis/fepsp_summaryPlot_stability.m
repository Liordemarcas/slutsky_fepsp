function [sumPlot] = fepsp_summaryPlot_stability(varargin)
% Create a figure for each channel, summarising the stability of i/o over
% time. Right now do not support more then 1 intensity in the stability
% analysis, of protocols with more then 1 simulation.
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - Numeric scalar. sampling frequency [Hz].
%   protocol_id - String or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "fepsp_getProtocol.m" for
%                 more info.
%   results     - Struct. See "fepsp_analyse.m" for more info.

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
%   slidingwin_size
%               - positive interger, window size for sliding window at
%                 trend line.
%                 Default: 10.
%   dt          - non-negative scalar. Dead time between
%                 stimulus onset and earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See "fepsp_getProtocol.m" for more info.
%                 Default: 2.
%   Trace2Time_factor
%               - positive numeric, how much time each trace takes,
%                 including waiting time to next trace, in secounds. for
%                 example if recourding length is 148.47ms & there is 15
%                 secounds between traces (io protocol @ 0.06Hz), then each
%                 trace takes 15.14748 (it may not be exsact to get general
%                 time scale).
%                 Default: 15.14748.
%   saveFig     - logical. save figure in basepath {true}
%
% OUTPUT:
%   sumPlot     - figure handles array. Handles to the summary figures, one
%                 for each channel.
%
% CALL:
%   sumPlot = fepsp_show("traces", <cell array>, "fs", <numeric
%   scalar>, "protocol_id", <string scalar>, "markings", <struct
%   scalar>, "results", <struct scalar or empty>, "intens", <numeric
%   vector>, "traces_Xlimit", <numeric 2 elements or empty>,
%   "traces_Ylimit", <numeric 2 elements or empty>, "dt", <numeric
%   non-negative scalar>)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('traces',           [],       @(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',               [],       @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',      [],       @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('results',          [],       @(x) validateattributes(x,{'struct'},{'scalar'}))
p.addParameter('intens',           [],       @(x) (isnumeric(x) && isvector(x)) || isempty(x))
p.addParameter('traces_xlim',      [],       @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('traces_ylim',      [],       @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('slidingwin_size',  20,       @(x) validateattributes(x,{'numeric'},{'scalar','positive','integer'}))
p.addParameter('dt',               2,        @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}))
p.addParameter('Trace2Time_factor',15.14748, @(x) validateattributes(x,{'numeric'},{'scalar','positive'}))
p.addParameter('saveFig',          true,     @islogical)

parse(p, varargin{:})

traces            = p.Results.traces;
fs                = p.Results.fs;
protocol_id       = p.Results.protocol_id;
results           = p.Results.results;
intens            = p.Results.intens;
traces_xlim       = sort(p.Results.traces_xlim);
traces_ylim       = sort(p.Results.traces_ylim);
slidingwin_size   = p.Results.slidingwin_size;
dt                = p.Results.dt;
Trace2Time_factor = p.Results.Trace2Time_factor;
saveFig           = p.Results.saveFig;

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
if nIntens > 1
    error('Right now more then 1 intens is not supported for stability')
end
nStim = protocol_info.nStim;
if nStim > 1
    error('Right now more then 1 stim is not supported for stability')
end

% results data
Amp = results.all_traces.Amp;
Slope =  results.all_traces.Slope;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iChan = size(traces, 1) : -1 : 1
 
    % open figure 
    sumPlot(iChan) = figure();
    sumPlot(iChan) = fepsp_graphics(sumPlot(iChan));          % set graphics
    sgtitle(sprintf('Channel %d - Stability', iChan),...
        'FontSize', 28, 'FontWeight', 'bold', 'FontName', 'FixedWidth')
    
    
    sb1 = subplot(2, 2, [1, 2]);

    nTraces = size(traces{iChan,1},2);
    cMap = green_magenta_Cmap(nTraces);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot traces through time, on top of each other
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    colororder(cMap)
    traces_plot_h = plot(protocol_info.Tstamps,traces{iChan,1},'LineWidth',1);
    for iTrace = 1 : nTraces
        traces_plot_h(iTrace).Color(4) = 0.5;
        traces_plot_h(iTrace).Tag = num2str(iTrace);
        r = dataTipTextRow('Trace Number', ones(size(traces_plot_h(iTrace).XData)) * iTrace);
        traces_plot_h(iTrace).DataTipTemplate.DataTipRows(end+1) = r;
    end
    xlim(protocol_info.traces_xlim)
    if isempty(traces_ylim)
        ylim([min(traces{iChan,1}(protocol_info.response.win, :), [], 'all')...
            max(traces{iChan,1}(protocol_info.response.win,:), [], 'all')] .* 1.1)
    else
        ylim(traces_ylim)
    end
    sb1.Colormap = cMap;
    cBar = colorbar();
    caxis([1 nTraces])
    clim([1 nTraces])
    ylabel(cBar,'Trace Number [1 is oldest]')
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    fepsp_graphics(sb1);          % set graphics
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot analysis results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    TimeFrame = (1:nTraces).*Trace2Time_factor;
    
    % scatter amplitutde
    sb2 = subplot(2, 2, 3);
    scatter(TimeFrame,Amp{iChan,1},30,'k','o','filled');
    hold on
    % regression line
    p = polyfit(TimeFrame(~isnan(Amp{iChan,1})),Amp{iChan,1}(~isnan(Amp{iChan,1})),1);
    p1 = plot(TimeFrame,polyval(p,TimeFrame),'g-','LineWidth',1.5,"DisplayName",sprintf('regression: y=(%.2g)x%+.2f',p(1),p(2)));
    % trend line
    avg_data = movmean(Amp{iChan,1},slidingwin_size,'omitnan');
    std_data = movstd(Amp{iChan,1},slidingwin_size,'omitnan');
    p2 = plot(TimeFrame,avg_data,'Color',[1,0,1],'LineWidth',2,"DisplayName",sprintf('Trend Line (win size %d)',slidingwin_size));
    UL = avg_data + std_data;
    LL = avg_data - std_data;
    p3 = patch([TimeFrame,fliplr(TimeFrame)],[LL,fliplr(UL)],[1,0,1],'FaceAlpha',0.2,'EdgeColor','r','DisplayName','1 std');
    uistack(p3,"bottom");
    xlabel('Time [s]')
    ylabel('Amplidute [mV]')
    legend([p1,p2,p3],"Location","best")
    fepsp_graphics(sb2);
    
    % scatter slope
    sb3 = subplot(2, 2, 4);
    scatter(TimeFrame,Slope{iChan,1},30,'k','o','filled');
    hold on
    % regression line
    p = polyfit(TimeFrame(~isnan(Slope{iChan,1})),Slope{iChan,1}(~isnan(Slope{iChan,1})),1);
    p1 = plot(TimeFrame,polyval(p,TimeFrame),'g-','LineWidth',1.5,"DisplayName",sprintf('regression: y=(%.2g)x%+.2f',p(1),p(2)));
    % trend line
    avg_data = movmean(Slope{iChan,1},slidingwin_size,'omitnan');
    std_data = movstd(Slope{iChan,1},slidingwin_size,'omitnan');
    p2 = plot(TimeFrame,avg_data,'Color',[1,0,1],'LineWidth',2,"DisplayName",sprintf('Trend Line (win size %d)',slidingwin_size));
    UL = avg_data + std_data;
    LL = avg_data - std_data;
    p3 = patch([TimeFrame,fliplr(TimeFrame)],[LL,fliplr(UL)],[1,0,1],'FaceAlpha',0.2,'EdgeColor','r','DisplayName','1 std');
    uistack(p3,"bottom");
    xlabel('Time [s]')
    ylabel('Slope [mV/ms]')
    legend([p1,p2,p3],"Location","best")
    fepsp_graphics(sb3);
    
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
    figname = fullfile(figpath, sprintf('%s_fepsp_results', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF