function writeSnaptoxyz(BaseFldr, fldrname, StepNum, XYZ_snap, Atoms, DL, DoubleAnalType)

% if ~exist([BaseFldr fldrname '\1stWL_Trajectories'])
%     mkdir([BaseFldr fldrname '\1stWL_Trajectories']);
% end

% fidout = fopen([BaseFldr fldrname '\1stWL_Trajectories\'  DoubleAnalType '_' num2str(StepNum) '.xyz'], 'w');
fidout = fopen([BaseFldr fldrname '\'  DoubleAnalType '_' num2str(StepNum) '.xyz'], 'w');
fprintf(fidout, [num2str(length(DL)) newline]);
fprintf(fidout, ['i = ' num2str(StepNum) newline]);

for j = 1:length(DL)
    fprintf(fidout,[pad(Atoms(DL(j),:),8) pad(num2str(XYZ_snap(DL(j),1), '%.10g'),20) pad(num2str(XYZ_snap(DL(j),2), '%.10g'),20) pad(num2str(XYZ_snap(DL(j),3), '%.10g'),20) newline]);
end
fclose(fidout);

return
