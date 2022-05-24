import numpy as np
    
def CP2kInputParameterReplace(Infldr = None,Outfldr = None,jobID = None,param = None,paramVal = None): 
    flname = dir(np.array([Infldr,'*.inp']))
    fid = open(np.array([Infldr,flname.name]))
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    relevant = find(not cellfun(isempty,strfind(lines,param)) )
    print(np.array(['Found ',num2str(len(relevant)),' instances of the parameter ',param,' which will be set to ',paramVal,'.']))
    for i in np.arange(1,len(relevant)+1).reshape(-1):
        pIndx = strfind(lines[relevant(i)],param)
        lines[relevant[i]] = np.array([param,' ',paramVal])
        lines[relevant[i]] = pad(lines[relevant(i)],len(lines[relevant(i)]) + pIndx - 1,'left')
    
    fid = open(np.array([Outfldr,jobID,'.inp']),'w')
    for i in np.arange(1,len(lines)+1).reshape(-1):
        fid.write('%s\n' % (lines[i]))
    
    fid.close()
    # # # # # # # Input file parameter replacer (standalone version with loop many folders)
# # # # # # # param = 'EXTRAPOLATION_ORDER';paramVal = '0';
# # # # # # # param = 'BACKUP_COPIES';paramVal = '2';
# # # # # # # param = 'mpiexec -n ';paramVal = '768 cp2k.popt ptwater10Na12Cl_CPlike.inp >& out.log';
# # # # # #
# # # # # # Basefldr = 'G:\Imperial\CP-LIKE\10Na12Cl\ASPC0\';
# # # # # # Allfldrs = dir([Basefldr '*_*']);
# # # # # # # flname = '\ptwater10Na12Cl_CPlike.inp';
# # # # # # flname = '\cp2k_parallel.qs';
# # # # # #
# # # # # # for ii = 1:length(Allfldrs)
# # # # # #     fid = fopen([Basefldr Allfldrs(ii).name flname]);
# # # # # #     lines = textscan(fid,'#s','delimiter','\n', 'whitespace', '');
fid.close();
# # # # # #     lines = lines{1};
# # # # # #     relevant =  find(~cellfun(@isempty,strfind(lines,param)));
# # # # # #     disp(['Found ' num2str(length(relevant)) ' instances of the parameter ' param ' which will be set to ' paramVal '.']);
# # # # # #     for i = 1:length(relevant)
# # # # # #         pIndx = strfind(lines{relevant(i)}, param);
# # # # # # #         lines{relevant(i)} = [param ' ' paramVal];
# # # # # #         lines{relevant(i)} = [param paramVal];
# # # # # #         lines{relevant(i)} = pad(lines{relevant(i)}, length(lines{relevant(i)})+pIndx-1, 'left');
# # # # # #     end
# # # # # #     delete([Basefldr Allfldrs(ii).name flname]);
# # # # # #     fid = fopen([Basefldr Allfldrs(ii).name flname],'w');
# # # # # #     for i = 1:length(lines)
# # # # # #         fprintf(fid,'#s\n',lines{i});
# # # # # #     end
fid.close();
# # # # # # end