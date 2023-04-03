function marking_fig = fepsp_markings(varargin)

% creates gui for user to mark start (green) and end (red) points of a
% response. user can choose to go through all intensities or just one
% intensity per channel (see fast_mark). the number of marker pairs in a
% trace is determined by the number of stimulations derived from the
% stimulus protocol. provides easy navigation between all channels and
% intensities. the output consists of markings for each trace separately
% and markings for the average trace. allows removal or inversion of
% individual traces (right click on plot). compensates for small
% inaccuracies in the manual markings by jittering in time the markers
% closer to the peak (max_jitter).
%
% INPUT (required):
%   traces      - 2d cell array. see fepsp_org2traces.m. 
%   fs          - numeric scalar. sampling frequency [Hz].
%   protocol_id - string or char. ID of stimulation protocol.
%                 e.g. "io","stp" or "custom". See "fepsp_getProtocol.m" for
%                 more info.

% INPUT (optional):
%   intens      - numeric vector of stimulus intensities used for each
%                 column of traces. Typically provided in units of uA. if
%                 empty will be set as a unit increasing vector.
%   base_path -   string or char. Full path to where the output should be
%                 saved. The name of the last folder in base_path will be
%                 the prefix for all saved data. e.g. base_path =
%                 'lh85_211012_132000' than the output of fepsp_markings
%                 will be 'lh85_211012_132000_fepsp_markings'.
%                 Default: pwd.
%   traces_xlim - numeric vector of 2 elements. Time span relative
%                 to stimulus onset for displaying the trace [ms]. For
%                 example, xLimit = [-1 30] will show the trace from 1 ms
%                 before the stimulus until 30 ms after the stimulus. If
%                 empty will use protocol default values. See
%                 "fepsp_getProtocol.m" for more info. 
%                 Default: [].
%   traces_ylim - numeric vector with 2 elements. voltage limits for
%                 displaying the trace [mV]. these y-limits remain constant
%                 throughout all intensities. If empty will be set
%                 according to the range of traces in highest intensity
%                 (excluding stimulus artifact).
%                 Default: [].
%   dt          - non-negative scalar. Dead time between
%                 stimulus onset & earliest possible response [ms].
%                 Used to omit stimulus artifact from analysis. 
%                 See "fepsp_getProtocol.m" for more info.
%                 Default: 2.
%   max_jitter  - non-negative scalar. Time tolerance for selecting
%                 the response boundaries [ms]. I.e., the precise start
%                 / end points of a response will be the max or min value
%                 within +/- max_jitter of the points selected by the user.
%                 This jitter compensates for small inaccuracies when
%                 marking the peak / trough.
%                 Default: 0.5.
%   fast_mark   - logical flag. If true, user will be asked to
%                 mark response boundaries for only the largest intensity.
%                 Those boundaries will be applied to all intensities
%                 without manual inspection. Note that for each intensity,
%                 the precise points of end / start will still be within
%                 +/- max_jitter of the manual selection. 
%                 Default: false.
%   hold_dt     - logical flag. If true, mark bounderies will have a
%                 fixed time distance between them (before jitter!). 
%                 Time can be changed in the GUI, defualts at 5ms. 
%                 It is recommended to use this option with max_jitter = 0.
%                 Default: false.
%
% OUTPUTS:
%   while gui is open the function will return marking_fig which is a
%   figure handle to the gui.
%   
%   after the gui window is closed and/or the user clicks "export", the
%   output will be saved (see base_path for location & name). Output will
%   consist of the variable "traces", same as input but with all the
%   deletions (as NaN) & inversions happened during marking, & struct
%   "markings" with the following fields:
% 
%   peaks       - standard cell. 2d cell array of channel (row) x intensity
%                 (column), with each cell containing a 2d numeric mat of
%                 stimulus number (row) x repetition (column). mat values
%                 are the sample of the response peak. note: Any traces
%                 removed by the user will receive nan as their markings.
%   starts      - standard cell. same as peaks but with mat values
%                 representing the sample of the response start.
%   peaks_avg   - 2d cell array of channel (row) x intensity (column), with
%                 each cell containing a numeric column vector. vec values
%                 are samples of the response peak, one for each stimulus.
%                 This is similar to "standard cell" (see README) but
%                 inside each cell there is only 1 column for the average
%                 trace instead of a column for every repetition (trace).
%   starts_avg  - same as peaks_avg but with vec values representing
%                 the sample of the response start.
%
% CALL:
%   marking_fig = fepsp_markings("traces", <cell array>, "fs", <numeric
%   scalar>, "protocol_id", <string scalar>, "base_path", <folder path>,
%   "intens", <numeric vector>, "traces_xlim", <numeric 2 elements or
%   empty>, "traces_ylim", <numeric 2 elements, or empty>, "dt", <numeric
%   non-negative scalar>, "max_jitter", <numeric non-negative
%   scalar>, "fast_mark", <logical flag>);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Give Option to enter inputs via struct, instead of name-value
p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('traces',        [],     @(x) validateattributes(x,{'cell'},{'2d'}))
p.addParameter('fs',            [],     @(x) validateattributes(x,{'numeric'},{'scalar'}))
p.addParameter('protocol_id',   [],     @(x) validateattributes(x,{'string','char'},{'scalartext'}))
p.addParameter('intens',        [],     @(x) validateattributes(x,{'numeric'},{'vector'}))
p.addParameter('base_path',     pwd,    @isfolder)
p.addParameter('traces_xlim',   [],     @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('traces_ylim',   [],     @(x) (isnumeric(x) && numel(x)==2) || isempty(x))
p.addParameter('dt',            2,      @(x) validateattributes(x,{'numeric'},{'vector','nonnegative'}))
p.addParameter('max_jitter',    0.5,    @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
p.addParameter('fast_mark',     false,  @(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))
p.addParameter('hold_dt',       false,  @(x) validateattributes(x,{'logical','numeric'},{'binary','scalar'}))

p.parse(varargin{:})

traces          = p.Results.traces;
fs              = p.Results.fs;
protocol_id     = p.Results.protocol_id;
intens          = p.Results.intens;
base_path       = p.Results.base_path;
traces_xlim     = sort(p.Results.traces_xlim);
traces_ylim     = sort(p.Results.traces_ylim);
dt              = p.Results.dt;
max_jitter      = p.Results.max_jitter;
fast_mark       = logical(p.Results.fast_mark);
hold_dt         = logical(p.Results.hold_dt);

% validate intens
if isempty(intens)
    intens = -(size(traces, 2) : -1 : 1);
elseif numel(intens) ~= size(traces, 2)
    error('number of intensities does not fit size of traces')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get protocol info
protocol_info = fepsp_getProtocol("protocol_id",protocol_id,"fs",fs,"dt",dt);
% if user gave traces_xlim, take them over protocol default
if ~isempty(traces_xlim)
    protocol_info.traces_xlim = traces_xlim;
end

% initilialize cell arrays of markings
[mark_peaks, mark_starts] =...
    deal(cellfun(@(x) nan(protocol_info.nStim, size(x, 2)), traces, 'UniformOutput', 0));
[mark_peaks_avg, mark_starts_avg] =...
    deal(repmat({nan(protocol_info.nStim, 1)}, size(traces)));

% convert max_jitter to max index difference from selected time point -
% half in each direction
max_ind_jit = round(((max_jitter/1000)*fs)./2);

% create average trace
for iTraces_group = numel(traces):-1:1
    [Chan,Inten] = ind2sub(size(traces),iTraces_group);
    avg_traces(Chan,Inten,:) = mean(traces{iTraces_group},2,'omitnan'); % electrode-by-intensity-by-sample
end

% find the order of intensities
[intens_sorted,intens_order] = sort(intens,'ascend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

marking_fig = figure();
marking_fig = fepsp_graphics(marking_fig);                        % set graphics
marking_fig.Name = 'fepsp Marking GUI';

% plot the maximum intensity on channel 1. intens_order(end) is the column
% of the highest intensity
traces_plot = plot(protocol_info.Tstamps,traces{1,intens_order(end)});

% organize axis limits
axis tight
xlim(protocol_info.traces_xlim)
if isempty(traces_ylim)
    % use protocol_info.response.win in order to find the minmax values of
    % the Y axis, without including stim artifact 
    yLimit_base = [min(traces{1, intens_order(end)}(protocol_info.response.win,:), [], 'all'),...
        max(traces{1, end}(protocol_info.response.win,:), [], 'all')];
    ylim(yLimit_base*1.1)
else
    % if user defined ylimits, use them
    ylim(traces_ylim)
end

% label & title nicely :)
xlabel('Time [ms]')
ylabel('Amp [mV]')
title({sprintf('Channel %d @ %duA', 1, intens(intens_order(end))),...
    ['\fontsize{18}', 'Move green & red lines to the start & trough of the response'],...
    'Press Enter to save & move on'})

% create Movable Lines.
% their starting location is based on response window, in order to place
% them "around" each stimulation.
lines_start_pos = protocol_info.Tstamps(protocol_info.response.base)';
colors = {'g','r'}; % prepare for loop - switch between green and red lines
line_rule = {'S','P'}; % prepare for loop - switch between start (S) and peak (P) lines
% make lines 20 times bigger then needed to present response. Robust to
% user zoom changes, while maintain label position
lines_Ylim = ylim()'+[-1;1]*(20*max(abs(ylim())));
for iLine = numel(lines_start_pos):-1:1
    rule_ind = mod(iLine,2); % switch between colours (and S/P marks, accordingly)
    rule_ind(rule_ind == 0) = 2; % convert index 0 to 2 in order to complete the switch circle
    mark_lines(iLine) = drawline('Position',[repmat(lines_start_pos(iLine),2,1) lines_Ylim],'InteractionsAllowed','translate','color',colors{rule_ind},...
        'Label',sprintf('%d%s t=%.1f',ceil(iLine/2),line_rule{rule_ind},lines_start_pos(iLine)));
    % 'LabelAlpha',0.3 - Not present in 2019a, however recommended in later
    % versions. Helps preventing label from "obscuring" traces. To add it simply
    % add this name value parameter presented here to the "drawline" call
end
Lis = addlistener(mark_lines,'MovingROI',@(obj,evt) line_move_manage(obj,evt,marking_fig)); % listener for managing lines movement & label

% build GUI Buttons & Interactivity:

% move channel droplist
Ch_strings = split(num2str(1:size(traces,1)));
ch_button = uicontrol(marking_fig,'Style','popupmenu','String',Ch_strings,'Value',1,'Tooltip',sprintf('Channel Selection)'));
% move intensity droplist
intens_button = uicontrol(marking_fig,'Style','popupmenu','String',split(num2str(intens_sorted)),'Value',length(intens_sorted),'Tooltip','Intens Selection');
% callbacks to move intensity & trace droplists - they are dependent on each
% other, so they need to be defined after they were created
intens_button.Callback = @(~,~) switch_traces_group(marking_fig,ch_button,intens_button,'Int');
ch_button.Callback = @(~,~) switch_traces_group(marking_fig,ch_button,intens_button,'Ch');
% export marking to workspace
save_button = uicontrol(marking_fig,'Style','pushbutton','String','Export Marks','Tooltip','Export Markings','Callback',@(~,~) export_marks(marking_fig));
% checkbox for fast analysis option
fast_mark_checkbox = uicontrol(marking_fig,'Style','checkbox','String','Fast Mark','Value',fast_mark);
% checkbox for fixed dt
fixed_dt_checkbox = uicontrol(marking_fig,'Style','checkbox','String','Hold dt','Value',hold_dt,'Callback',@(obj,~) enable_fixed_dt(obj.Value,marking_fig));
% dt 2 hold numeric input box 
fixed_dt_editbox = uicontrol(marking_fig,'Style','edit','String','5','Tooltip','dt to hold [ms]','Callback',@(obj,~) validate_update_editbox(obj,marking_fig),'UserData','5');
% fix buttons with too small width
save_button.Position(3) = 100;
fast_mark_checkbox.Position(3) = 100;
fixed_dt_checkbox.Position(3) = 100;
fixed_dt_editbox.Position(3) = 100;
% align all button to sit next to each other. If you create any more
% uicontrols, add them here to avoid manual placements
align([ch_button, intens_button, save_button, fast_mark_checkbox,fixed_dt_checkbox,fixed_dt_editbox],'Fixed',5,'bottom')

% create context menu to remove or invert traces from it
cm = uicontextmenu('Callback',@CM_show_nTrace);
uimenu(cm,'Label','Trace num -') % using the function trace_trapper, this will mark what plot user right clicked on
uimenu(cm,'Label','Remove Trace','Callback',@(obj,~) remove_trace(obj,marking_fig,ch_button.Value,intens_button.Value));
uimenu(cm,'Label','Invert Trace','Callback',@(obj,~) invert_trace(obj,marking_fig,ch_button.Value,intens_button.Value));
% Note: using 'Callback' property of uicontextmenu is not recommended since
% 2020a. It is recommended to use 'ContextMenuOpeningFcn', its meaning is
% clearer. In order to be compatible to 2019a, we will use
% 'Callback', while remembering this callback run just before uicontextmenu
% is opened.

% add the context menu to the plots
for iTrace = 1:length(traces_plot)
    traces_plot(iTrace).UIContextMenu = cm;
    traces_plot(iTrace).Tag = num2str(iTrace); % tag each trace by its column number in the traces group
    traces_plot(iTrace).ButtonDownFcn = {@trace_trapper,cm}; % improve right click success in finding the right trace
    r = dataTipTextRow('Trace num',ones(size(traces_plot(iTrace).XData))*iTrace); % write trace num in datatips
    traces_plot(iTrace).DataTipTemplate.DataTipRows(end+1) = r;
end
% Note: using 'UIContextMenu' property of line is not recommended since
% 2020a. It is recommended to use 'ContextMenu'. 
% In order to be compatible to 2019a, we will use 'UIContextMenu'.

% add average trace
ax = findobj(marking_fig,'type','axes','-depth',1);
ax = fepsp_graphics(ax);          % set graphics

hold(ax,'on')
traces_plot(end+1) = plot(protocol_info.Tstamps,squeeze(avg_traces(1,intens_order(end),:)),'--b','LineWidth',2);
traces_plot(end).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Avg Trace','');

% make figure callbacks
marking_fig.CloseRequestFcn = @(~,~) close_GUI(marking_fig,Lis); % Export marks on window close
marking_fig.KeyReleaseFcn =  @(~,evt) next_on_enter(evt,marking_fig,ch_button,intens_button); % Create move trace by pressing "enter"/"return"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save between callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save all data related to the fepsp under 1 struct, to keep it organized
fepsp_data.traces               = traces;
fepsp_data.avg_traces           = avg_traces;
fepsp_data.intens               = intens;
fepsp_data.fs                   = fs;
fepsp_data.protocol_info        = protocol_info;
fepsp_data.max_ind_jit          = max_ind_jit;
fepsp_data.mark_starts        = mark_starts;
fepsp_data.mark_peaks         = mark_peaks;
fepsp_data.mark_starts_avg    = mark_starts_avg;
fepsp_data.mark_peaks_avg     = mark_peaks_avg;

% save everything we want to transfer between callbacked under UserData
marking_fig.UserData.fepsp_data                     = fepsp_data;
marking_fig.UserData.traces_handles.traces_plot     = traces_plot;
marking_fig.UserData.traces_handles.cm              = cm;
marking_fig.UserData.mark_lines_handles.Lis         = Lis;
marking_fig.UserData.mark_lines_handles.mark_lines  = mark_lines;
marking_fig.UserData.last_Ch                        = ch_button.Value;
marking_fig.UserData.last_intens                    = intens_button.Value;
marking_fig.UserData.fast_mark_checkbox             = fast_mark_checkbox;
marking_fig.UserData.fixed_dt_checkbox              = fixed_dt_checkbox;
marking_fig.UserData.fixed_dt_editbox               = fixed_dt_editbox;
marking_fig.UserData.intens_order                   = intens_order;
marking_fig.UserData.traces_ylim                  = traces_ylim;
marking_fig.UserData.base_path                    = base_path;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% final touches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update hold dt if needed
enable_fixed_dt(hold_dt,marking_fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_markings(marking_fig,Ch_button_val,intens_button_val)
% saving the current markings according to channel & intensity values.
% before saving, finds the best match within jitter for each trace.
%
% INPUT:
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button_val 
%               - numeric scalar. Index of the channel to save
%                 to.
%   intens_button_val 
%               - numeric scalar. Index of the intensity to save
%                 to.

% extract all needed variables form UserData.fepsp_data
fepsp_data       = marking_fig.UserData.fepsp_data;
protocol_info    = fepsp_data.protocol_info;
max_ind_jit      = fepsp_data.max_ind_jit;
traces           = fepsp_data.traces;
avg_traces       = fepsp_data.avg_traces;
intens_order     = marking_fig.UserData.intens_order;

% intens_button_val is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
traces_col = intens_order(intens_button_val);

% find the X value of each marking line
lines_pos = [marking_fig.UserData.mark_lines_handles.mark_lines.Position];
line_xpos = lines_pos(1,1:2:end);

% calculate the closest Point that exist in data
Tstamp_mat = repmat(protocol_info.Tstamps, 1, protocol_info.nStim);
[~,start_ind_base] = min(abs(Tstamp_mat - line_xpos(1:2:end)), [], 1);
[~,peak_ind_base] = min(abs(Tstamp_mat - line_xpos(2:2:end)), [], 1);

% get mean baseline values for each trace & average trace
baseline_traces = mean(traces{Ch_button_val,traces_col}(protocol_info.baseline,:),1,'omitnan');
baseline_avg_trace = mean(avg_traces(Ch_button_val,traces_col,protocol_info.baseline),'omitnan');

limit_mat(:,1) = [
    min(protocol_info.response.base(1,2),start_ind_base(1)) ; % let user placment override window
    peak_ind_base(1:(end-1))' + max_ind_jit + 1];
limit_mat(:,2) = [
     start_ind_base(2:end)' - max_ind_jit - 1;
     max(protocol_info.response.base(end,2),peak_ind_base(end))]; % let user placment override window
% deal with position jitter (find point in window by max_ind_jit with
% max distance from mean value of baseline)
for iStim = protocol_info.nStim:-1:1
    % create jitter windows for start of the response and for peak:
    jit_win_start = max(start_ind_base(iStim) - max_ind_jit, limit_mat(iStim,1)); % The max makes sure Start tolerance doesn't go before previus stimulus
    jit_win_end = min(start_ind_base(iStim) + max_ind_jit, peak_ind_base(iStim) - max_ind_jit - 1); % The min makes sure Start tolerance doesn't take after peak_ind
    start_jit_win = jit_win_start:jit_win_end;
    start_jit_win = start_jit_win(ismember(start_jit_win,protocol_info.response.win)); % make sure window stays in responce window - avoide artifact (by dt)

    jit_win_start = max(peak_ind_base(iStim) - max_ind_jit, start_ind_base(iStim) + max_ind_jit + 1); % The max makes sure End tolerance doesn't take before start_ind
    jit_win_end = min(peak_ind_base(iStim) + max_ind_jit, limit_mat(iStim,2)); % The min makes sure End tolerance doesn't go before past next stimulus
    peak_jit_win = jit_win_start:jit_win_end;
    peak_jit_win = peak_jit_win(ismember(peak_jit_win,protocol_info.response.win)); % make sure window stays in responce window - avoide artifact (by dt)
    
    % find the maximal distance from mean baseline in all traces:
    % baseline_traces is a row vector, and inside each traces
    % cell there is a mat in which each column is a different trace
    % (repetition). Therefore, the "minus" operation between the mat & the
    % row vector removes the vector from each row of the matrix, resulting
    % in a distance mat
    distance_mat = abs(traces{Ch_button_val,traces_col}(start_jit_win,:) - baseline_traces);
    [start_values,start_ind(iStim,:)] = max(distance_mat,[],1); 
    distance_mat = abs(traces{Ch_button_val,traces_col}(peak_jit_win,:) - baseline_traces);
    [peak_values,peak_ind(iStim,:)] = max(distance_mat,[],1); 
    % right now each index refer to the tolerance window.
    % convert it to full trace
    start_ind(iStim,:) = start_jit_win(1) + start_ind(iStim,:);
    peak_ind(iStim,:) = peak_jit_win(1) + peak_ind(iStim,:);
    % if the minimal value is nan, trace was deleted.
    % however, each index will mistakenly return 1 + Start of tolerance window.
    % fix it to NaN.
    start_ind(iStim,isnan(start_values)) = nan;
    peak_ind(iStim,isnan(peak_values)) = nan;
    
    % find the maximal distance from mean baseline in the average traces:
    distance_vec = abs(avg_traces(Ch_button_val,traces_col,start_jit_win) - baseline_avg_trace);
    [~,start_ind_avg(iStim)] = max(distance_vec); 
    distance_vec = abs(avg_traces(Ch_button_val,traces_col,peak_jit_win) - baseline_avg_trace);
    [~,peak_ind_avg(iStim)] = max(distance_vec); 
    % right now each index refer to the tolerance window.
    % convert it to full trace
    start_ind_avg(iStim) = start_jit_win(1) + start_ind_avg(iStim);
    peak_ind_avg(iStim) = peak_jit_win(1) + peak_ind_avg(iStim);
   
end

% save the markings to the right location in each matching
% UserData.fepsp_data variable
marking_fig.UserData.fepsp_data.mark_starts{Ch_button_val,traces_col}     = start_ind;
marking_fig.UserData.fepsp_data.mark_peaks{Ch_button_val,traces_col}      = peak_ind;
marking_fig.UserData.fepsp_data.mark_peaks_avg{Ch_button_val,traces_col}  = peak_ind_avg';
marking_fig.UserData.fepsp_data.mark_starts_avg{Ch_button_val,traces_col} = start_ind_avg';

end


function switch_traces_group(marking_fig,ch_button,intens_button,Parm)
% manage movement between traces group. Save markings before moving.
%
% INPUT:
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button   - uicontrol handle scalar. The uicontrol in
%                 charge of changing between channels.
%   intens_button  
%               - uicontrol handle scalar. The uicontrol in
%                 charge of changing between intensities.
%   Parm        - char vector, either 'Next', 'Ch' or 'Int'.
%                 What type of transfer to perform:
%                   Next - move to next traces group (same channel with
%                          lower intensity if available, else highest
%                          intensity in the next channel in a rising
%                          circular manner - i.e. 1 to 2 to 3 back to 1).
%                   Ch   - move to user specific channel (according to Ch_button
%                          droplist value), do not change the intensity.
%                   Int  - move to user specific intensity (according to
%                          intens_button droplist value), do not change the channel.

switch Parm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % next
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Next'
        % quick move to the next traces group
        
        % fast analysis case
        if marking_fig.UserData.fast_mark_checkbox.Value
            % save current marks for all intensities
            for iIntens = 1:numel(marking_fig.UserData.fepsp_data.intens)
                save_markings(marking_fig,ch_button.Value,iIntens)
            end
            
            % change intensity to smallest intensity to cause transfer to
            % next channel
            intens_button.Value = 1;
        else
            
            % save the existing markings before moving traces group
            save_markings(marking_fig,ch_button.Value,intens_button.Value)
        end
        
        % move to next group
        if intens_button.Value == 1
            % we are on the smallest intensity. Next group requires a channel switch.
            
            % check if all marks have been already filled - if so, prompt
            % export & close GUI. All mean traces marks are not NaN mean
            % that all marks were filled 
            mark_peaks_avg = [marking_fig.UserData.fepsp_data.mark_peaks_avg{:}];            
            if all(~isnan(mark_peaks_avg),'all')
                Status = questdlg('Finished all channels & intensities. Exit?','Finished Main','Yes','No','Yes');
                switch Status
                    case  "Yes"
                        close_GUI(marking_fig,marking_fig.UserData.mark_lines_handles.Lis)
                        return
                    case []
                        % unclear answer - cancel move
                        msgbox('User cancel, aborting move')
                        
                        % flash droplists colours red, to confirm to user
                        % change was cancelled successfully
                        ch_button.BackgroundColor = 'r';
                        intens_button.BackgroundColor = 'r';
                        pause(0.2)
                        ch_button.BackgroundColor = 'w';
                        intens_button.BackgroundColor = 'w';
                        return
                end
            end
            
            % correct droplists
            if ch_button.Value == size(marking_fig.UserData.fepsp_data.traces,1)
                % if this is the last channel, return to the first 
                ch_button.Value = 1;
            else
                % move to next channel
                ch_button.Value = ch_button.Value+1;
            end
            % move to the highest intensity
            intens_button.Value = length(intens_button.String);
            
            % plot the new traces group - done by getting the new channel &
            % intensity from the droplists
            plot_traces_group(marking_fig,ch_button.Value,intens_button.Value)
            
            % flash droplists colours green, to confirm to user change was
            % done successfully
            ch_button.BackgroundColor     = 'g';
            intens_button.BackgroundColor = 'g';
            pause(0.2)
            ch_button.BackgroundColor     = 'w';
            intens_button.BackgroundColor = 'w';
            
        else
            % we are not on the smallest intensity, move to a smaller one
            intens_button.Value = intens_button.Value-1;
            
            % plot the new traces group - done by getting the new channel &
            % intensity from the droplists
            plot_traces_group(marking_fig,ch_button.Value,intens_button.Value)
            
            % flash droplist colour green, to confirm to user change was
            % done successfully
            intens_button.BackgroundColor = 'g';
            pause(0.2)
            intens_button.BackgroundColor = 'w';
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Ch'
        % change only the channel
        
        % offer the user to save marking before moving to other trace group
        Status = questdlg('Save current markings?','Save current?','Yes','No','Yes');
        switch Status
            % if 'Yes', save. If 'No', continue without saving. If user
            % cancelled (clicked X), it is unclear answer - cancel move.
            case 'Yes'
                save_markings(marking_fig,marking_fig.UserData.last_Ch,intens_button.Value)
            case ''
                % notify user
                msgbox('User cancel, aborting move')
                % restore droplist value to what it was before droplist use
                ch_button.Value = marking_fig.UserData.last_Ch;
                
                % flash droplist colour red, to confirm to user change was
                % cancelled successfully
                ch_button.BackgroundColor = 'r';
                pause(0.2)
                ch_button.BackgroundColor = 'w';
                return
        end
        
        % plot the new traces group - done by getting the new channel &
        % intensity from the droplists
        plot_traces_group(marking_fig,ch_button.Value,intens_button.Value)
        
        % sometime this fix "focus stays on droplist & cause enter key to not work" issue
        figure(marking_fig);
        ch_button.Enable = 'off';
        drawnow
        ch_button.Enable = 'on';
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % intensity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Int'
        % change only the intensity
        
        % offer the user to save marking before moving to other trace group
        Status = questdlg('Save current markings?','Save current?','Yes','No','Yes');
        switch Status
            % if 'Yes', save. If 'No', continue without saving. If user
            % cancelled (clicked X), it is unclear answer - cancel move.
            case 'Yes'
                save_markings(marking_fig,ch_button.Value,marking_fig.UserData.last_intens)
            case ''
                % notify user
                msgbox('User Cancel, Aborting move')
                
                % restore droplist value to what it was before droplist use
                intens_button.Value = marking_fig.UserData.last_intens;
                
                % flash droplist colour red, to confirm to user change was
                % cancelled successfully
                intens_button.BackgroundColor = 'r';
                pause(0.2)
                intens_button.BackgroundColor = 'w';
                return
        end
        
        % plot the new traces group - done by getting the new channel &
        % intensity from the droplists
        plot_traces_group(marking_fig,ch_button.Value,intens_button.Value)
        
        % sometime this fix "focus stays on droplist & cause enter key to not work" issue
        figure(marking_fig);
        intens_button.Enable = 'off';
        drawnow
        intens_button.Enable = 'on';
end

% save current channel & intensity numbers - enable saving marking
% & cancel when changing trace group using droplists, as their
% value is changed before callback
marking_fig.UserData.last_Ch     = ch_button.Value;
marking_fig.UserData.last_intens = intens_button.Value;

end


function plot_traces_group(marking_fig,wanted_ch,wanted_inten)
% change the traces group plotted in marking_fig according to given channel
% & intensity indices
%
% INPUT:
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   wanted_Ch   - numeric scalar. Index of the channel to plot from.
%   wanted_inten 
%               - numeric scalar. Index of the intensity to plot from.

% extract all needed variables form UserData.fepsp_data
fepsp_data      = marking_fig.UserData.fepsp_data;
traces          = fepsp_data.traces;
protocol_info   = fepsp_data.protocol_info;
avg_traces      = fepsp_data.avg_traces;
intens          = fepsp_data.intens;
intens_order    = marking_fig.UserData.intens_order;
traces_ylim   = marking_fig.UserData.traces_ylim;

% wanted_inten is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
traces_col = intens_order(wanted_inten);

% get the axes
ax = findobj(marking_fig,'type','axes','-depth',1);

% return zoom to default
zoom out

% remove old traces group & plot new traces
delete(marking_fig.UserData.traces_handles.traces_plot)
traces_plot = plot(ax,protocol_info.Tstamps,traces{wanted_ch,traces_col});

% add the context menu to the plots
for iTrace = 1:length(traces_plot)
    traces_plot(iTrace).UIContextMenu = marking_fig.UserData.traces_handles.cm;
    traces_plot(iTrace).Tag = num2str(iTrace); % tag each trace by its column number in the trace cell
    traces_plot(iTrace).ButtonDownFcn = {@trace_trapper,marking_fig.UserData.traces_handles.cm}; % improve right click success in finding the right trace
    r = dataTipTextRow('Trace num',ones(size(traces_plot(iTrace).XData))*iTrace); % write trace num in datatips
    traces_plot(iTrace).DataTipTemplate.DataTipRows(end+1) = r;
end

% add average trace
traces_plot(end+1) = plot(protocol_info.Tstamps,squeeze(avg_traces(wanted_ch,traces_col,:)),'--b','LineWidth',2);
traces_plot(end).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Avg Trace','');

% change title to match the new traces group
title(ax, {sprintf('Channel %d @ %duA', wanted_ch, intens(traces_col)),...
    ['\fontsize{18}', 'Move green & red lines to the start & trough of the response'],...
    'Press Enter to save & move on'})

% fix axis limits if necessary:
xlim(protocol_info.traces_xlim)
% get the needed ylimits to present the current traces group (while ignoring the stim artifact)
traces_group_Ylims = [min(traces{wanted_ch, traces_col}(protocol_info.response.win,:), [], 'all'),...
    max(traces{wanted_ch, traces_col}(protocol_info.response.win,:), [], 'all')]*1.1;

if ~isempty(traces_ylim)
    % if user defined ylimits, keep them
    ylim(traces_ylim)
else
    if marking_fig.UserData.last_Ch~= wanted_ch
        % channel was switch - inherited ylimits are irrelevant. Compare
        % ylimits of biggest intensity to current traces group- take the
        % biggest of the two
        
        % get the ylimits matching the biggest intensity (while ignoring the
        % stim artifact). intens_order(end) is the column of the highest intensity
        high_intens_Ylims = [min(traces{wanted_ch, intens_order(end)}(protocol_info.response.win,:), [], 'all'),...
            max(traces{wanted_ch, intens_order(end)}(protocol_info.response.win,:), [], 'all')]*1.1;
        % compare and take the biggest
        ylim([min(traces_group_Ylims(1),high_intens_Ylims(1)) max(traces_group_Ylims(2),high_intens_Ylims(2))]);
    else
        % compare ylimit inherited from last traces group to
        % current - take the biggest of the two.
        
        % get the ylimits inherited
        inherited_Ylims = ylim()*1.1;
        % compare and take the biggest
        ylim([min(traces_group_Ylims(1),inherited_Ylims(1)) max(traces_group_Ylims(2),inherited_Ylims(2))])
    end
end

% make lines 20 times bigger then current ylimits, by adding. Robust to
% user zoom changes, while maintain label position
lines_Ylim = ylim()'+[-1;1]*(20*max(abs(ylim())));
% change markers lines to match new ylimits
for aa = 1:length(marking_fig.UserData.mark_lines_handles.mark_lines)
    marking_fig.UserData.mark_lines_handles.mark_lines(aa).Position(:,2) = lines_Ylim;
end

% save the handle to the new plots
marking_fig.UserData.traces_handles.traces_plot = traces_plot;

end


function remove_trace(obj,marking_fig,Ch_button_val,intens_button_val)
% remove selected trace, and all data connected to it. It uses the context
% menu to find which trace, the window for UserData, and current channel &
% intensity indices in order to remove the data from the correct location
% and redraw the traces group (with the bad trace removed)
%
% INPUT:
%   obj         - uimenu scalar. The caller object for
%                 this callback function.
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button_val  
%               - numeric scalar. Index of the channel to plot from.
%   intens_button_val
%               - numeric scalar. Index of the intensity to plot from.

% get the trace User want to remove, from its plotted line. It is saved
% under the uicontextmenu UserData, the parent of the uimenu that called this
% callback
Trace = obj.Parent.UserData;
if ~isa(Trace,'matlab.graphics.chart.primitive.Line')
    % if no trace was found, give nice error to user, and cancel
    errordlg({'For some reason, no trace have been capture.' 'Please try again'})
    return
end
% get the trace number (its column number in the trace cell) from its tag
nTrace = str2double(Trace.Tag);

% make sure user want to remove trace - it is not a mistake
Status = questdlg({'Are you sure you want to remove trace?' 'This cannot be undone, and will effect the output struct'},'Removal Confirm','Yes','No','No');
if ismember(Status,{'No',''})
    return
end

% intens_button_val is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
intens_order = marking_fig.UserData.intens_order;
traces_col = intens_order(intens_button_val);

% delete the trace data
marking_fig.UserData.fepsp_data.traces{Ch_button_val,traces_col}(:,nTrace)      = nan;
marking_fig.UserData.fepsp_data.mark_peaks{Ch_button_val,traces_col}(:,nTrace)  = nan;
marking_fig.UserData.fepsp_data.mark_starts{Ch_button_val,traces_col}(:,nTrace) = nan;
% recalculate mean trace, for changing display later
marking_fig.UserData.fepsp_data.avg_traces(Ch_button_val,traces_col,:)          = mean(marking_fig.UserData.fepsp_data.traces{Ch_button_val,traces_col},2,'omitnan');

% redraw the traces - the deleted trace won't be present.
plot_traces_group(marking_fig,Ch_button_val,intens_button_val)

if isempty(marking_fig.UserData.traces_ylim)
    % if user defined ylimits, they are fixed, and will not be changed. Else:
    % try to shrink the ylimits, however no smaller than this channel max
    % intens ylimits or the needed to view the traces
    traces = marking_fig.UserData.fepsp_data.traces;
    response_win = marking_fig.UserData.fepsp_data.protocol_info.response.win;
    high_intens_Ylims = [min(traces{Ch_button_val, intens_order(end)}(response_win,:), [], 'all'),...
        max(traces{Ch_button_val, intens_order(end)}(response_win,:), [], 'all')]*1.1;
    traces_group_Ylims = [min(traces{Ch_button_val,traces_col}(response_win,:), [], 'all'),...
        max(traces{Ch_button_val,traces_col}(response_win,:), [], 'all')]*1.1;
    % compare and take the biggest
    ylim([min(traces_group_Ylims(1),high_intens_Ylims(1)) max(traces_group_Ylims(2),high_intens_Ylims(2))]);
end

end


function invert_trace(obj,marking_fig,Ch_button_val,intens_button_val)
% invert a selected trace. It uses the context menu to find which trace, the
% window for UserData, and current channel & intensity indices in order to
% invert the data from the correct location and redraw the traces group
%
% INPUT:
%   obj         - uimenu scalar. The caller object for
%                 this callback function.
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.
%   Ch_button_val 
%               - numeric scalar. Index of the channel to plot from.
%   intens_button_val
%               - numeric scalar. Index of the intensity to plot from.

% get the trace User want to invert, from its plotted line. It is saved
% under the uicontextmenu UserData, the parent of the uimenu that called this
% callback
trace_line_obj = obj.Parent.UserData;
if ~isa(trace_line_obj,'matlab.graphics.chart.primitive.Line')
    errordlg({'For some reason, no trace have been capture.' 'Please try again'})
    return
end
nTrace = str2double(trace_line_obj.Tag);

% intens_button_val is actually indexing intens_order, instead of matching traces
% column number. Convert it back by simple indexing
intens_order = marking_fig.UserData.intens_order;
traces_col = intens_order(intens_button_val);

% invert the trace data
marking_fig.UserData.fepsp_data.traces{Ch_button_val,traces_col}(:,nTrace) = -marking_fig.UserData.fepsp_data.traces{Ch_button_val,traces_col}(:,nTrace);

% recalculate average trace
marking_fig.UserData.fepsp_data.avg_traces(Ch_button_val,traces_col,:) = mean(marking_fig.UserData.fepsp_data.traces{Ch_button_val,traces_col},2,'omitnan');

% redraw trace group
plot_traces_group(marking_fig,Ch_button_val,intens_button_val)

if isempty(marking_fig.UserData.traces_ylim)
    % if user defined ylimits, they are fixed, and will not be changed. Else:
    % try to shrink the ylimits, however no smaller than this channel max
    % intens ylimits or the needed to view the traces
    traces = marking_fig.UserData.fepsp_data.traces;
    response_win = marking_fig.UserData.fepsp_data.protocol_info.response.win;
    high_intens_Ylims = [min(traces{Ch_button_val, intens_order(end)}(response_win,:), [], 'all'),...
        max(traces{Ch_button_val, intens_order(end)}(response_win,:), [], 'all')]*1.1;
    traces_group_Ylims = [min(traces{Ch_button_val,traces_col}(response_win,:), [], 'all'),...
        max(traces{Ch_button_val,traces_col}(response_win,:), [], 'all')]*1.1;
    % compare and take the biggest
    ylim([min(traces_group_Ylims(1),high_intens_Ylims(1)) max(traces_group_Ylims(2),high_intens_Ylims(2))]);
end

end


function export_marks(marking_fig)
% export markings from GUI UserData to base workspace
%
% INPUT:
%   marking_fig - figure scalar, window of fepsp_markings

% extract all outputs. Pack markings for easy export
markings.peaks          = marking_fig.UserData.fepsp_data.mark_peaks;
markings.starts         = marking_fig.UserData.fepsp_data.mark_starts;
markings.peaks_avg      = marking_fig.UserData.fepsp_data.mark_peaks_avg;
markings.starts_avg     = marking_fig.UserData.fepsp_data.mark_starts_avg;
traces                  = marking_fig.UserData.fepsp_data.traces;

% create file name from base_path
base_path = marking_fig.UserData.base_path;
[~,base_name] = fileparts(base_path);
marking_file = [base_path filesep base_name '_fepsp_markings.mat'];

% save
save(marking_file, 'markings', 'traces','-v7.3')


end


function close_GUI(marking_fig,Lis)
% save everything before closing GUI. Delete listeners - they should be
% removed when lineROI are deleted, but deleting them explicitly helps in
% some bug cases.
%
% INPUT:
%   marking_fig - figure scalar, window of fepsp_markings
%   Lis         - listeners vector for lineROI event "MovingROI"

export_marks(marking_fig)
delete(Lis)
delete(marking_fig)

end


function line_move_manage(obj,evt,marking_fig)
% function for managing lines movement.
% change the time written in lineROI label when moving it.
% prevent moving the lines along the Y axis.
% if fixed_dt_checkbox is checked, move the appropriate line to keep the dt.
%
% INPUT:
%   obj         - lineROI object that generated the event "MovingROI".
%   evt         - "MovingROI" event, generated in response to
%                 user moving the lineROI.
%   marking_fig - figure scalar, window of fepsp_markings.Needed
%                 mainly for everything saved under UserData.

% get the new position of the lineROI
XPos = evt.CurrentPosition(1);
% change the label to include it, after the '=' sign
obj.Label = [obj.Label(1:(find(obj.Label == '='))), sprintf('%.1f',XPos)];

% prevent moving on the Y axis, by restoring previous position
obj.Position(:,2) = evt.PreviousPosition(:,2);

% deal with dt holding
if marking_fig.UserData.fixed_dt_checkbox.Value
    update_connected_line(obj,marking_fig)
end

end


function update_connected_line(reference_line,marking_fig)
% move the connected line so the dt from the reference line will be constant.
%
% INPUT:
%   reference_line
%               - lineROI object, the connected line will be position
%                 relatively to its location.
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.

% find the connected line. Note that this only work when there are only
% 2 lines - start (S) and peak (P).
line_tag = regexp(reference_line.Label,'\d*[S,P]','match');
connected_tag = replace(line_tag,["S","P"],["P","S"]);
connected_log = contains({marking_fig.UserData.mark_lines_handles.mark_lines.Label},connected_tag);
connected_line = marking_fig.UserData.mark_lines_handles.mark_lines(connected_log);

% get dt range 2 hold
dt2hold = str2double(marking_fig.UserData.fixed_dt_editbox.String);

% find connected line new distance
if connected_tag{1}(end) == 'S'
    % if I moved the "peak" line, "start" should follow behind it -
    % reverse dt factor
    dt2hold = -dt2hold;
end
connected_pos = reference_line.Position(1,1) + dt2hold;

% update connected line pos & label
connected_line.Position(:,1) = connected_pos;
connected_line.Label = [connected_line.Label(1:(find(connected_line.Label == '='))), sprintf('%.1f',connected_pos)];
end


function enable_fixed_dt(checkbox_Value,marking_fig)
% Simple loop to place all lines inside dt next to their "start", only when
% fixed_dt_checkbox is checked
%
% INPUT:
%   checkbox_Value
%               - logical flag, fixed_dt_checkbox state, only move the lines to position if true
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.

if checkbox_Value
    % if checkbox is checked, make sure editbox is enable, then update all
    % lines for hold dt

    marking_fig.UserData.fixed_dt_editbox.Enable = 'on';
    lines_handles = marking_fig.UserData.mark_lines_handles.mark_lines;
    for iLine = 1:2:numel(lines_handles)
        update_connected_line(lines_handles(iLine),marking_fig)
    end

else
    % just make sure editbox is disabled
    
    marking_fig.UserData.fixed_dt_editbox.Enable = 'off';

end

end

function validate_update_editbox(obj,marking_fig)
% make sure only positive numbers are in the editbox,
% and place the lines accordingly
%
% INPUT:
%   obj         - handle to the hold dt editbox.
%   marking_fig - figure scalar, window of fepsp_markings. Needed
%                 mainly for everything saved under UserData.

dt_value = str2double(obj.String);
if isnan(dt_value) || ~isreal(dt_value) || dt_value <= 0
    % flash red, restore previous value
    obj.BackgroundColor = 'r';
    pause(0.1)
    obj.BackgroundColor = 'w';
    obj.String = obj.UserData;
else
    % flash green, save new value, place lines to match new value
    obj.BackgroundColor = 'g';
    pause(0.1)
    obj.BackgroundColor = 'w';
    obj.UserData = obj.String;
    enable_fixed_dt(marking_fig.UserData.fixed_dt_checkbox.Value,marking_fig)
end

end



function next_on_enter(evt,marking_fig,Ch_button,intens_button)
% cause movement to the next traces group when figure window is in focus,
% and user pressed "enter"/"return" key.
%
% INPUT:
%   evt         - keypress event generated by the figure.
%   marking_fig - figure scalar, window of fepsp_markings.
%   Ch_button   - uicontrol handle scalar. The uicontrol in charge of
%                 changing between channels.
%   intens_button 
%               - uicontrol handle scalar. The uicontrol in charge of
%                 changing between intensities. 

% make sure the key that was pressed was "enter"/"return" - its number is 13
if evt.Character == 13
    % if it was, just call switch_traces_group with parameter 'Next'
    switch_traces_group(marking_fig,Ch_button,intens_button,'Next')
end

end


function trace_trapper(obj,~,cm)
% will help the uicontext miss less often, in the cost of hover datatip
% due to some unknown reason (datacrusermode will still work as usual however)
%
% INPUT:
%   obj         - line obj scalar, of the trace that was clicked
%   ~           - the click event generated by the clicked trace
%   cm          - uicontextmenu scalar, opened when user right-click the trace

% simply place the trace handle under the context menu UserData for easy
% and reliable access
cm.UserData = obj;

end


function CM_show_nTrace(cm,~)
% present in trace number 1 of the uimenues what trace the user had right
% clicked on, when opening the uicontextmenu.
%
% INPUT:
%   cm          - uicontextmenu scalar, opened when user right-click the
%                 trace, object calling this callback
%   ~           - the generated event from opening the uicontextmenu

% extract the tag of the trace from the uicontext menu UserData (it was
% saved there by trace_trapper) and switch in the label of the uimenu what after
% the '-' character
cm.Children(3).Label = [extractBefore(cm.Children(3).Label,'-') '- ' num2str(cm.UserData.Tag)];

end


% EOF