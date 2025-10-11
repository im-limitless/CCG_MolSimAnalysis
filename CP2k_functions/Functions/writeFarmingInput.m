function writeFarmingInput(BaseFldr, Allfldrs, cores, ppn, time)

% if exist([BaseFldr 'Farm'],'dir')
%     warning('Directory already exists! Farming folder not created.');
%     return
% end
mkdir([BaseFldr 'Farm']);

fidout = fopen([BaseFldr 'Farm\Farm-1.restart'], 'w');
fprintf(fidout,['# CP2K farming file created by Matt Darby, Imperial College London' newline]);
fprintf(fidout,[' &GLOBAL' newline]);
fprintf(fidout,['   PROJECT OldMacDonald' newline]);
fprintf(fidout,['   PROGRAM FARMING' newline]);
fprintf(fidout,['   RUN_TYPE  NONE' newline]);
fprintf(fidout,[' &END GLOBAL' newline]);
fprintf(fidout,[' &FARMING' newline]);
fprintf(fidout,['   NGROUPS ' num2str(length(Allfldrs)) newline]);
fprintf(fidout,['   GROUP_SIZE ' num2str(cores*ppn/48) newline]);
for i = 1:length(Allfldrs)
    fprintf(fidout,['   &JOB' newline]);
    fprintf(fidout,['     DIRECTORY ../' Allfldrs(i).name(1:end-4) newline]);
    fprintf(fidout,['     INPUT_FILE_NAME ' Allfldrs(i).name(1:end-4) '-1.restart' newline]);
    fprintf(fidout,['     OUTPUT_FILE_NAME out.log' newline]);
    fprintf(fidout,['   &END JOB' newline]);
end
fprintf(fidout,[' &END FARMING' newline]);

fclose(fidout);

CreateSubmissionFile('LRZ', BaseFldr, 'Farm', time, cores*length(Allfldrs), ppn, 'No');

return