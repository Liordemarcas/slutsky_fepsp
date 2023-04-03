function fepsp_WCP_grp_analysis(basepath,options)
% Analyse a full day of WCP files.
% Assume in each folder there are wcp files named:
% join([out_name,regex_suffix],'')
% where regex_suffix is a suffix that may contain any regular expression (see regexp).
% everything will be saved inside "basepath".
%
% INPUT: (Positional, Required)
%   basepath        - list of folders to be analyse. if empty, defualts to pwd.
%
% INPUT: (Name - Value, Optional)
%   out_names       - Text, list of output names, see fepsp_wcpPipeline.
%                     Default: ["R_IO","R_STP","L_IO","L_STP"] 
%                       (4 types of analysis)
%   regex_suffix    - Text scalar, a single suffix to add to each output
%                     name, to find all relevant files. Assume files are of
%                     format join([out_name,regex_suffix],'').
%                     Default: "_(\d)*uA.wcp"
%   protocol_ids    - Text, list of protocol to match each out_name. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    basepath       (1,:) string {mustBeFolder} = pwd;
    options.out_names      (1,:) string  {mustBeText}   = ["R_IO","R_STP","L_IO","L_STP"];
    options.regex_suffix   (1,:) string  {mustBeText}   = "_(\d)*uA.wcp";
    options.protocol_ids   (1,:) string  {mustBeText}   = ["io_wheel","STP_148ms","io_wheel","STP_148ms"]
    options.fsOut          (1,:) double  {mustBeScalarOrEmpty} = [];
    options.wait_graph     (1,:) logical {mustBeNumericOrLogical} = false;
    options.inverse_all    (1,:) logical {mustBeNumericOrLogical} = false;
    options.overwrite_existing  (1,1) double {mustBeMember(options.overwrite_existing,[-1,0,1])} = 0;
    options.add_fields     (:,2) cell                             = cell.empty([0,2]);
end

out_names    = options.out_names;
regex_suffix = options.regex_suffix;
protocol_ids = options.protocol_ids;
fsOut        = options.fsOut;
wait_graph   = options.wait_graph;
inverse_all  = options.inverse_all;
overwrite_existing = options.overwrite_existing;

% make sure there are the same number of protocol_ids, outnames & regex_suffix
% (all sould match the number of file types to analyse)
if isscalar(protocol_ids)
    protocol_ids = repmat(protocol_ids,size(out_names));
else
    mustBeSameNumel(out_names,protocol_ids)
end

if isscalar(regex_suffix)
    regex_suffix = repmat(regex_suffix,size(out_names));
else
    mustBeSameNumel(out_names,regex_suffix)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all reg_exp, make sure they are not containing each other
reg_exp = join([out_names',regex_suffix'],'');
reg_list = 1:numel(reg_exp);
for iType = reg_list
    if any(contains(reg_exp(reg_list ~= iType),reg_exp(iType)))
        error("Expression based on out_names %d (and prehaps the same index in regex_suffix) is contained within other expression.\n This will cause 'overcatching', so abourting",iType)
    end
end

for iFold = basepath
    % move to the folder to analyse, and get all its files
    cd(iFold)
    [~,fold_name] = fileparts(iFold);
    fprintf('\n***********************\nNow In folder %s\n***********************\n',fold_name)
    files = dir;

    for iType = 1:numel(out_names)
        % find all files to analyse in this type
        out_intens = regexp({files.name},reg_exp(iType),'tokens');
        if all(cellfun(@isempty,out_intens))
            warning("In folder %s, Unable to find files matching %s",fold_name,reg_exp(iType))
            continue
        end
        
        % perform analysis
        try
        fepsp_wcpPipeline("basepath",pwd,"fepsp_protocol",char(protocol_ids(iType)),"fsOut",fsOut,"out_name",out_names(iType),'intens',str2double(string([out_intens{:}])),'wcpfiles',{files(~cellfun(@isempty,out_intens)).name},...
            'overwrite_existing',overwrite_existing,'inverse_all',inverse_all);
        catch err
            % let user skip if they didn't want to overwrite existing file.
            % any other error, rethrow
            if strcmp(err.identifier,'fepsp_wcpPipeline:UserNotOverWrite')
                warning('Skipping %s as user wanted not to overwrite',fullfile(fold_name,out_names(iType)))
            else
                rethrow(err)
            end
        end
        
        % let user inspect the output before moving on
        if wait_graph
            waitfor(gcf)
        else
            close(gcf)
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument Validation Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mustBeSameNumel(arg1,arg2)
% Throw error if number of elements between to variables do not match

    if numel(arg1) ~= numel(arg2)
        errID   = 'mustBeSameNumel:notNumelEqual';
        msgText = 'Input "%s" & input "%s" must have the same number of elements';
        throwAsCaller(MException(errID,msgText,inputname(1),inputname(2)))
    end

end

% EOF