function [sumPlot] = fepsp_summaryPlot(varargin)
% Create a figure for each channel, summarising the fepsp experiment.
% Will create input-output graph if protocol has only 1 stimulus, and
% facilitation graph if it has more.
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - Numeric scalar. sampling frequency [Hz].
%   protocol_id - String or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "fepsp_getProtocol.m" for
%                 more info.
%   markings    - Struct. See "fepsp_markings.m" for more info.
%   results     - Struct. See "fepsp_analyse.m" for more info.
%
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
%   dt          - non-negative vector. For each stimuli, Dead time between
%                 stimulus onset and earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See "fepsp_getProtocol.m" for more info.
%                 Default: 2.
%   stim2plot   - numeric 2d of intergers, X on 2. 
%                 Fist col is the dividend stimuli, secound col is the divisor stimuli.
%                 What stimuli to plot (in 1st col), and for each stimuli
%                 if to divide other other stimuli by it. Placing NaN in
%                 the secound col will tell the function not to divide, and
%                 to plot IO graph instead.
%                 For example,[1 nan ; 2 1 ; 3 1 ; 2 3] will plot IO for
%                 1st stimuli, fasilitation between first & secound, 
%                 fasilitation between first & third, and fasilitation between secound & third.
%                 Default: [] in this case:
%                   if only 1 stimulus will plot IO curve - simullar to [1 nan]
%                   if more than 1, will plot fasilitation to 1st stimulimi
%                   - simular to [2 1 ; 3 1; 4 1 ; ...]
%   saveFig     - logical. save figure in basepath {true}
%   Slope_mdl   - cell vector (1 elem per channel), contain cell vectors (1 elem per intens) containing
%                 LinearModel or NonLinearModel. 
%                 Model fitted to IO curve, based on slope responce.
%   Amp_mdl     - cell vector (1 elem per channel), contain cell vectors (1 elem per intens) containing
%                 LinearModel or NonLinearModel. 
%                 Model fitted to IO curve, based on amp responce.
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

p.addParameter('traces',        [], @(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',            [], @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',   [], @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('markings',      [], @(x) validateattributes(x,{'struct'},{'scalar'}))
p.addParameter('results',       [], @(x) validateattributes(x,{'struct'},{'scalar'}))
p.addParameter('intens',        [], @(x) (isnumeric(x) && isvector(x)) || isempty(x))
p.addParameter('traces_xlim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('traces_ylim',   [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('dt',            2,  @(x) validateattributes(x,{'numeric'},{'vector','nonnegative'}))
p.addParameter('stim2plot',     [], @(x) validateattributes(x,{'numeric'},{'2d','positive'}))
p.addParameter('saveFig',       true, @(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))
p.addParameter('Slope_mdl',     [], @(x) mustBeA(x,{'cell'}))
p.addParameter('Amp_mdl',       [], @(x) mustBeA(x,{'cell'}))

parse(p, varargin{:})

traces          = p.Results.traces;
fs              = p.Results.fs;
protocol_id     = p.Results.protocol_id;
markings        = p.Results.markings;
results         = p.Results.results;
intens          = p.Results.intens;
traces_xlim     = sort(p.Results.traces_xlim);
traces_ylim     = sort(p.Results.traces_ylim);
dt              = p.Results.dt;
stim2plot       = p.Results.stim2plot;
saveFig         = p.Results.saveFig;
Slope_mdl       = p.Results.Slope_mdl;
Amp_mdl         = p.Results.Amp_mdl;

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
nStim = protocol_info.nStim;
[intens_sorted, intens_order] = sort(intens, 'ascend');

% slop params
slope_window_avg = results.avg_traces.slope_win;
slope_area_label = sprintf('%.3g%% to %.3g%%', results.slope_area * 100);

% results data
Amp = results.all_traces.Amp;
Slope =  results.all_traces.Slope;

% create default for stim2plot, understand how many plots should be
if isempty(stim2plot)
    if nStim == 1
        stim2plot = [1 nan];
    else
        stim2plot = [(2:nStim)',ones(nStim-1,1)];
    end
end
nAn_rows = size(stim2plot,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iChan = size(traces, 1) : -1 : 1
 
    % open figure 
    sumPlot(iChan) = figure('WindowState','maximized');
    sumPlot(iChan) = fepsp_graphics(sumPlot(iChan));          % set graphics
    sgtitle(sprintf('Channel %d - %s', iChan,upper(protocol_info.protocol_id)),...
        'FontSize', 28, 'FontWeight', 'bold', 'FontName', 'FixedWidth','Interpreter','none')
    
    sb1 = subplot(nAn_rows + 1, 2, [1 2]);
    %{
    if nStim > 1
        if ~stp_legacy
            sb1 = subplot(nStim, 2, [1 2]); % New option: 2 subplots for each stim
        else
            sb1 = subplot(length(intens) + 1, 1, 1); % Old option: 1 subplot for each intens
        end
        
    else
        sb1 = subplot(2, 2, [1, 2]);
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot average traces for all intensities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % avg trace for each intensity
    traces_avg = cell2mat(cellfun(@(x) mean(x, 2, 'omitnan'),...
        traces(iChan, :), 'UniformOutput', false));
    
    % plot and change data tips
    if nIntens > 1
        colororder(green_magenta_Cmap(nIntens))
%         colororder(distinguishable_colors(nIntens))
    end
    traces_avg_h = plot(protocol_info.Tstamps, traces_avg, 'LineWidth', 2);
    for iIntens = 1 : nIntens
        traces_avg_h(iIntens).Tag = num2str(intens(iIntens));
        r = dataTipTextRow('intensity', ones(size(traces_avg_h(iIntens).XData)) * intens(iIntens));
        traces_avg_h(iIntens).DataTipTemplate.DataTipRows(end+1) = r;
    end
    hold on
    
    % plot start markings
    starts_avg = [markings.starts_avg{iChan,:}];
    starts_lin = sub2ind(size(traces_avg), starts_avg,...
        repmat(1 : nIntens, nStim, 1));
    starts_h = plot(protocol_info.Tstamps(starts_avg), traces_avg(starts_lin), 'k*');
    
    % plot peak markings
    peaks_avg = [markings.peaks_avg{iChan,:}];
    peaks_lin = sub2ind(size(traces_avg), peaks_avg,...
        repmat(1:nIntens, nStim, 1));
    peaks_h = plot(protocol_info.Tstamps(peaks_avg), traces_avg(peaks_lin), 'kO');
    
    % color slope area in each trace
    for iStim = nStim : -1 : 1
        for iIntens = nIntens : -1 : 1
            % get slope window & plot on top of samples matching it
            loop_slope_window = sort(slope_window_avg{iChan, iIntens}(:, iStim)); 
            slope_area_mark(iIntens) =...
                plot(protocol_info.Tstamps(loop_slope_window(1) : loop_slope_window(2)),...
                traces_avg(loop_slope_window(1) : loop_slope_window(2), iIntens), 'c', 'LineWidth', 3);
            % make transparent
            slope_area_mark(iIntens).Color(4) = 0.5;
            % remove from legend
            slope_area_mark(iIntens).Annotation.LegendInformation.IconDisplayStyle = 'off';
            % remove from data tips selection
            slope_area_mark(iIntens).PickableParts = 'none';
        end
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
    legend([traces_avg_h(intens_order); starts_h(1); peaks_h(1)],...
        [string(intens_sorted), {'response Start', 'response Peak'}],...
        'Location','southeast','NumColumns', 2)
    
    % show one reprasentative slope
    slope_area_mark(1).Annotation.LegendInformation.IconDisplayStyle = 'on';
    slope_area_mark(1).DisplayName = ['Slope: ' slope_area_label];
    
    % labels
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    title('Average Traces')
    sb1 = fepsp_graphics(sb1);          % set graphics

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot analysis results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_all_intens = [Amp{iChan,:}];
    slope_all_intens = [Slope{iChan,:}];
    % get how many traces are in each intens
    nTraces = cellfun(@(x) size(x,2),traces(iChan,:));
    intens_group = repelem(intens,nTraces);
    for iDivision = 1:nAn_rows
        dividend_stim = stim2plot(iDivision,1);
        divisor_stim  = stim2plot(iDivision,2);
        if ~isnan(divisor_stim)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2 plots per stimulation, simillar to input-output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % collect data & normalize it
            amp_fac_all    = amp_all_intens./amp_all_intens(divisor_stim,:);
            slope_fac_all    = slope_all_intens./slope_all_intens(divisor_stim,:);

            % create axes
            amp_subplot = subplot(nAn_rows + 1, 2, 2*iDivision + 1);
            slope_subplot = subplot(nAn_rows + 1, 2, 2*iDivision + 2);
            % If we number each element in an increasing rowwise manner, in a row-by-2 matrix:
            % rows end in nRow*2, therefore they start with nRow*2-1
            %   1 | 1, 2    ->  1*2-1, 1*2
            %   2 | 3, 4    ->  2*2-1, 2*2
            %   3 | 5, 6    ->  2*3-1, 3*2
            % here I add +2 to make it start at secound row

            % plot "input - output" for amplitude
            %                 boxes = boxchart(loop_subplot,intens_group,slope_fac_all(iStim,:),"GroupByColor",intens_group);
            boxplot(amp_subplot,amp_fac_all(dividend_stim,:),intens_group,'GroupOrder',string(intens_sorted),'Positions',intens_sorted)

            % plot "input - output" for slope
            %                 boxes = boxchart(loop_subplot,intens_group,slope_fac_all(iStim,:),"GroupByColor",intens_group);
            boxplot(slope_subplot,slope_fac_all(dividend_stim,:),intens_group,'GroupOrder',string(intens_sorted),'Positions',intens_sorted)

            % color boxplots 2 match avg_traces_plot colors
            amp_box_edges = findobj(amp_subplot,'Tag','Box');
            slope_box_edges = findobj(slope_subplot,'Tag','Box');
            for iBox = nIntens:-1:1
                iIntens = intens_order(end-iBox+1); % intensity is inversed & sorted between box & input intens order
                patch(amp_subplot,get(amp_box_edges(iBox),'XData'),get(amp_box_edges(iBox),'YData'),traces_avg_h(iIntens).Color,'FaceAlpha',.5,'PickableParts','none');
                patch(slope_subplot,get(slope_box_edges(iBox),'XData'),get(slope_box_edges(iBox),'YData'),traces_avg_h(iIntens).Color,'FaceAlpha',.5,'PickableParts','none');
            end

            % arrange graphics
            amp_subplot.PositionConstraint = 'innerposition';      % Let the graphic part take more space
            slope_subplot.PositionConstraint = 'innerposition';    % Let the graphic part take more space
            xticklabels(amp_subplot,string(intens_sorted))   % for clearer view
            xticklabels(slope_subplot,string(intens_sorted)) % for clearer view
            ylabel(amp_subplot,{'Facilitation' sprintf('$\\frac{stim %d}{stim %d}$',dividend_stim,divisor_stim)},'Interpreter','latex') % 1 label is enough
            fepsp_graphics(amp_subplot);
            fepsp_graphics(slope_subplot);
            
            if iDivision == nAn_rows
                xlabel(amp_subplot,sprintf('\\fontsize{11}Intensity\n\\fontsize{16}Amplitude'))
                xlabel(slope_subplot,sprintf('\\fontsize{11}Intensity\n\\fontsize{16}Slope'))
            end
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % old option: create 1 subplot for each intensity, with both stimulations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % multi stimuli protocol - show facilitation.
            % for each intensity create 1 box plot for amplitude & 1 for slope
            for iIntens = 1 : nIntens
                loop_subplot = subplot(nIntens + 1, 1, iIntens + 1);

                % normalised amplitude & slope
                loop_amp_norm = Amp{iChan, iIntens}./Amp{iChan, iIntens}(1,:);
                loop_slope_norm = Slope{iChan, iIntens}./Slope{iChan, iIntens}(1,:);

                % order data so amplitude & slope are next to each other in
                % every stimulus number
                cat_data = nan(size(traces{iChan, iIntens},2),2*nStim);
                col2fill = mat2cell(1:2*nStim,1,ones(1,nStim)*2);
                for iStim = nStim:-1:1
                    cat_data(:,col2fill{iStim}) = [loop_amp_norm(iStim,:)' loop_slope_norm(iStim,:)'];
                end

                % lazy fix for 1 trace issue - cat_data will be a vector,
                % boxplot will create only 1 box. In general, using only 1
                % trace is not recommended
                if isvector(cat_data) && size(traces{iChan, iIntens},2) == 1
                    % place nans on top, just to plot a vertical bar at each value
                    cat_data = [nan(size(cat_data));cat_data]; %#ok<AGROW> % not really, is predefined earlier
                end

                % plot boxplots
                boxplot(cat_data(:,3:end));

                % colour boxplots by type
                colors = repmat({'m','g'},1,nStim);
                box_edges = findobj(loop_subplot,'Tag','Box');
                for iBox = (nStim*2-2):-1:1
                    color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),colors{iBox},'FaceAlpha',.5,'PickableParts','none');
                end

                % add separating lines between stimulations
                xline(loop_subplot, (2:(nStim-1)) +0.5,'--k','Alpha',0.5)

                % finishing touch:
                % fix ticks to be centred between boxes - make it seem like
                % the number refer both of the boxes above it
                xticks(loop_subplot,1.5:2:(nStim*2-2));
                text(max(xlim()).*0.95,max(ylim()).*0.90,sprintf('\\bf{%guA}',intens(iIntens)))
                if iIntens == ceil(nIntens/2)
                    % This is the middle graph, place there the ylabel & legend
                    legend(color_box(1:2),['Slope: ' slope_area_label],'Amplitude','Location','best')
                    ylabel({'Mean Part of' 'First Stimulation' '[stim/(stim num 1)]'})
                end
                if iIntens == nIntens
                    % This is the bottom graph. Place there the xlabel & xticklabels
                    xlabel('Number of Stimulation [#]')
                    xticklabels(string(2:nStim))
                else
                    xticklabels({})
                end
                loop_subplot.PositionConstraint = 'innerposition';    % Let the graphic part take more space
                loop_subplot = fepsp_graphics(loop_subplot);          % set graphics
            end
            %}


        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % only 1 stimulation in protocol - show input-output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % plot input-output for amplitude
            loop_subplot = subplot(nAn_rows + 1, 2, 2*iDivision + 1);
            boxplot(loop_subplot,amp_all_intens(dividend_stim,:),string(intens_group),'GroupOrder',string(intens_sorted),'Positions',intens_sorted)
            xticklabels(string(intens_sorted))
            xlabel('Intensity [uA]')
            ylabel({sprintf('Stim %d',dividend_stim),'Amplidute [mV]'})
%             title('Input/Output (Amplidute)')

            % color boxplots 2 match avg_traces_plot colors
            box_edges = findobj(loop_subplot,'Tag','Box');
            for iBox = nIntens:-1:1
                iIntens = intens_order(end-iBox+1); % intensity is inversed & sorted between box & input intens order
                color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),traces_avg_h(iIntens).Color,'FaceAlpha',.5,'PickableParts','none');
            end

            % plot the model below everything, if it exist
            if ~isempty(Amp_mdl{iChan}{dividend_stim})
                hold(loop_subplot,"on")
                [p_res,p_res_CI] = predict(Amp_mdl{iChan}{dividend_stim},intens_sorted');
                mdl_lines = plot(loop_subplot,intens_sorted,p_res,'b-',intens_sorted',p_res_CI,'b--');
                uistack(mdl_lines,"bottom")
            end
            loop_subplot = fepsp_graphics(loop_subplot);          % set graphics


            % plot input-output for slope
            loop_subplot = subplot(nAn_rows + 1, 2, 2*iDivision + 2);
            boxplot(loop_subplot,slope_all_intens(dividend_stim,:),string(intens_group),'GroupOrder',string(intens_sorted),'Positions',intens_sorted)
            xlabel('Intensity [uA]')
            ylabel('Slope [mV/ms]')
%             title(['Input/Output (Slope: ' slope_area_label ')'])

            % color boxplots 2 match avg_traces_plot colors
            box_edges = findobj(loop_subplot,'Tag','Box');
            for iBox = nIntens:-1:1
                iIntens = intens_order(end-iBox+1); % intensity is inversed & sorted between box & input intens order
                color_box(iBox) = patch(get(box_edges(iBox),'XData'),get(box_edges(iBox),'YData'),traces_avg_h(iIntens).Color,'FaceAlpha',.5,'PickableParts','none');
            end

            % plot the model below everything, if it exist
            if ~isempty(Slope_mdl{iChan}{dividend_stim})
                hold(loop_subplot,"on")
                [p_res,p_res_CI] = predict(Slope_mdl{iChan}{dividend_stim},intens_sorted');
                mdl_lines = plot(loop_subplot,intens_sorted,p_res,'b-',intens_sorted,p_res_CI,'b--');
                uistack(mdl_lines,"bottom")
            end

            loop_subplot = fepsp_graphics(loop_subplot);          % set graphics
        end
    end
    
    % let user open graphs in a new fig, enlarged
    all_ax = findobj(sumPlot(iChan),'Type','axes');
    cm = uicontextmenu(sumPlot(iChan));
    uimenu(cm,"Text","Open in new fig","MenuSelectedFcn",@full_fig);
    
    [all_ax.ContextMenu] = deal(cm);

    % save
    if saveFig
        fprintf(['This is the time to make edits (such as moving the legends) before figure export.'...
            '\ncontinue running the code (by F5) ones you are ready to export!\n'])
        keyboard
        basepath = pwd;
        [~, basename] = fileparts(basepath);
        figpath = fullfile(basepath, 'graphics', 'fepsp');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_fepsp_results', basename));
        export_fig(figname, '-tif', '-transparent', '-r300')
        fprintf('Saved!\n')
    end

end

end

function full_fig(~,evt)
    ax = evt.Source.Parent.Parent.CurrentAxes;
    f = figure();
    new_ax = copyobj(ax,f);
    drawnow;
    subplot(1,1,1,new_ax)
    fepsp_graphics(f);
    fepsp_graphics(new_ax);
end
% EOF