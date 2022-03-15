% split .xyz file and extract specific (final) geometry, write to .xtd file
% format for materials studio

% Set the location of the calculation output
Basefldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\'; % Base directory containing calculation directory ("\" included at end)
system = 'CP_Pit_20F'; % Name of calculation directory (no "\")
Trajectory = 'CP_Pit_20F-pos-1_7.xyz';

nSampleSteps = 500; % sampling in units of number of steps
InitalStepOverride = []; % Use to set initial step manually. Set to [] to use entire trajectory
% InitalStepOverride = 10000; % Use to set initial step manually. Set to [] to use entire trajectory

fid  = fopen([Basefldr system '\' Trajectory]);
disp(['Reading xyz trajectory file for ' system]);
lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');

fclose(fid);

lines = lines{1};
nAtoms = str2num(lines{1});
relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
nConfigs = length(relevant);

initialI = lines{relevant(1)};
EqI = strfind(initialI, '=');
ComI = strfind(initialI, ',');
indxI = [EqI(1) ComI(1)];
initialI = str2num(initialI(indxI(1)+1:indxI(2)-1));
first = ceil(initialI/nSampleSteps)*nSampleSteps;

finalI = lines{relevant(end)};
EqI = strfind(finalI, '=');
ComI = strfind(finalI, ',');
indxI = [EqI(1) ComI(1)];
finalI = str2num(finalI(indxI(1)+1:indxI(2)-1));
last = floor(finalI/nSampleSteps)*nSampleSteps;

fidout = fopen([Basefldr system '\final.xyz'], 'w');
for i = relevant(end)-1:relevant(end)+nAtoms
    fprintf(fidout, [lines{i} '\n']);
end
fclose(fidout);

if ~isempty(InitalStepOverride)
    first = InitalStepOverride;
end

fidout = fopen([Basefldr system '\' system '_' num2str(first) 'to' num2str(last) '_' num2str(nSampleSteps) 'step.xyz'], 'w'); % AIMD
for i = first-initialI+1:nSampleSteps:last-initialI+1 % AIMD
% fidout = fopen([Basefldr fldrnm '\Sample' num2str(1) '_' num2str(length(relevant)) '.xyz'], 'w'); %GO
%     for i = 1:length(relevant) % GO
        disp(['Writing sample trajectory... ' num2str(100*(i/(last-initialI+1))) ' % complete']);
    for j = relevant(i)-1:relevant(i)+nAtoms
        fprintf(fidout, [lines{j} '\n']);
    end
end
fclose(fidout);

CP2kOptimPathParse(Basefldr,system,[system '_' num2str(first) 'to' num2str(last) '_' num2str(nSampleSteps) 'step.xyz']) % AIMD
% CP2kOptimPathParse(Basefldr,fldrnm, 'CP_Pit_18H22F-1_1000_Metal-pos-1.xyz') % GO