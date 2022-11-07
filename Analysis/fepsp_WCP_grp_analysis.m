function fepsp_WCP_grp_analysis(basepath)
% Analyse a full day of WCP files.
% Assume in each folder there are wcp files named:
% *{Direction}_{Protocol}_(\d)*uA.wcp
% Direction can be R (right) or L (left), protocol can be IO or STP.
% everything will be saved inside "basepath".
%
% INPUT:
%   basepath - list of folders to be analyse. if empty, defualts to pwd.

arguments
    basepath (1,:) string {mustBeFolder} = pwd
end

for iFold = basepath
    cd(iFold)

    files = dir;
    out_intens = regexp({files.name},'R_IO_(\d)*uA.wcp','tokens');
    fepsp_wcpPipeline("basepath",pwd,"fepsp_protocol",'io',"fsOut",5000,"out_name",'R_IO','intens',str2double(string([out_intens{:}])),'wcpfiles',{files(~cellfun(@isempty,out_intens)).name});
    waitfor(gcf)
    files = dir;
    out_intens = regexp({files.name},'L_IO_(\d)*uA.wcp','tokens');
    fepsp_wcpPipeline("basepath",pwd,"fepsp_protocol",'io',"fsOut",5000,"out_name",'L_IO','intens',str2double(string([out_intens{:}])),'wcpfiles',{files(~cellfun(@isempty,out_intens)).name});
    waitfor(gcf)
    files = dir;
    out_intens = regexp({files.name},'L_STP_(\d)*uA.wcp','tokens');
    fepsp_wcpPipeline("basepath",pwd,"fepsp_protocol",'STP_148ms',"fsOut",5000,"out_name",'L_STP','intens',str2double(string([out_intens{:}])),'wcpfiles',{files(~cellfun(@isempty,out_intens)).name});
    waitfor(gcf)
    files = dir;
    out_intens = regexp({files.name},'R_STP_(\d)*uA.wcp','tokens');
    fepsp_wcpPipeline("basepath",pwd,"fepsp_protocol",'STP_148ms',"fsOut",5000,"out_name",'R_STP','intens',str2double(string([out_intens{:}])),'wcpfiles',{files(~cellfun(@isempty,out_intens)).name});
    waitfor(gcf)
end