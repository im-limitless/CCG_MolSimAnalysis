function CP2kInputParameterReplace(Infldr, Outfldr, jobID, param, paramVal)

flname = dir([Infldr jobID '\*-1.restart']);

fid = fopen([Infldr jobID '\' flname.name]);
lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);
lines = lines{1};


relevant =  find(~cellfun(@isempty,strfind(lines,param)));
disp(['Found ' num2str(length(relevant)) ' instances of the parameter ' param ' which will be set to ' paramVal '.']);
for i = 1:length(relevant)
    pIndx = strfind(lines{relevant(i)}, param);
    lines{relevant(i)} = [param ' ' paramVal];
    lines{relevant(i)} = pad(lines{relevant(i)}, length(lines{relevant(i)})+pIndx-1, 'left');
end

fid = fopen([Outfldr jobID '\' jobID '-1.restart'],'w');
for i = 1:length(lines)
    fprintf(fid,'%s\n',lines{i});
end
fclose(fid);


% % % % % % % Input file parameter replacer (standalone version with loop many folders)
% % % % % % % param = 'EXTRAPOLATION_ORDER';paramVal = '0';
% % % % % % % param = 'BACKUP_COPIES';paramVal = '2';
% % % % % % % param = 'mpiexec -n ';paramVal = '768 cp2k.popt ptwater10Na12Cl_CPlike.inp >& out.log';
% % % % % %
% % % % % % Basefldr = 'G:\Imperial\CP-LIKE\10Na12Cl\ASPC0\';
% % % % % % Allfldrs = dir([Basefldr '*_*']);
% % % % % % % flname = '\ptwater10Na12Cl_CPlike.inp';
% % % % % % flname = '\cp2k_parallel.qs';
% % % % % %
% % % % % % for ii = 1:length(Allfldrs)
% % % % % %     fid = fopen([Basefldr Allfldrs(ii).name flname]);
% % % % % %     lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
% % % % % %     fclose(fid);
% % % % % %     lines = lines{1};
% % % % % %     relevant =  find(~cellfun(@isempty,strfind(lines,param)));
% % % % % %     disp(['Found ' num2str(length(relevant)) ' instances of the parameter ' param ' which will be set to ' paramVal '.']);
% % % % % %     for i = 1:length(relevant)
% % % % % %         pIndx = strfind(lines{relevant(i)}, param);
% % % % % % %         lines{relevant(i)} = [param ' ' paramVal];
% % % % % %         lines{relevant(i)} = [param paramVal];
% % % % % %         lines{relevant(i)} = pad(lines{relevant(i)}, length(lines{relevant(i)})+pIndx-1, 'left');
% % % % % %     end
% % % % % %     delete([Basefldr Allfldrs(ii).name flname]);
% % % % % %     fid = fopen([Basefldr Allfldrs(ii).name flname],'w');
% % % % % %     for i = 1:length(lines)
% % % % % %         fprintf(fid,'%s\n',lines{i});
% % % % % %     end
% % % % % %     fclose(fid);
% % % % % % end