clear all;
close all;

% split .xyz file and extract specific (final) geometry, write to .xtd file
% format for materials studio

% Set the location of the calculation output
BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\'; % Base directory containing calculation directory ("\" included at end)
system = 'CP_Pit_Water'; % Name of calculation directory (no "\")
Trajectory = 'CP_Pit_Water-pos-1_4.xyz';

nSampleSteps = 500; % sampling in units of number of steps
InitalStepOverride = []; % Use to set initial step manually. Set to [] to use entire trajectory
% InitalStepOverride = 16000; % Use to set initial step manually. Set to [] to use entire trajectory

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);

fidout = fopen([BaseFldr system '\final_test.xyz'], 'w');
fprintf(fidout,[num2str(nAtoms) newline]);
fprintf(fidout,['i = ' num2str(StepNum(end)) newline]);
for i = 1:nAtoms
    fprintf(fidout,[Atoms(i,:) pad(extractBefore(num2str(XYZ(end, i, 1), '%#.10g'),10),20) pad(extractBefore(num2str(XYZ(end, i, 2), '%#.10g'),10),20) pad(extractBefore(num2str(XYZ(end, i, 3), '%#.10g'),10),20) newline]);
end
fclose(fidout);

first = ceil(StepNum(1)/nSampleSteps)*nSampleSteps;
last = floor(StepNum(end)/nSampleSteps)*nSampleSteps;

if ~isempty(InitalStepOverride)
    first = InitalStepOverride;
end

fidout = fopen([BaseFldr system '\' system '_' num2str(first) 'to' num2str(last) '_' num2str(nSampleSteps) 'step.xyz'], 'w'); % AIMD
for i = first-StepNum(1)+1:nSampleSteps:last-StepNum(1)+1
    disp(['Writing sample trajectory... ' num2str(100*(i/(last-StepNum(1)+1))) ' % complete']);
    fprintf(fidout,[num2str(nAtoms) newline]);
    fprintf(fidout,['i = ' num2str(StepNum(i)) newline]);

    for j = 1:nAtoms
        fprintf(fidout,[Atoms(j,:) pad(extractBefore(num2str(XYZ(i, j, 1), '%#.10g'),10),20) pad(extractBefore(num2str(XYZ(i, j, 2), '%#.10g'),10),20) pad(extractBefore(num2str(XYZ(i, j, 3), '%#.10g'),10),20) newline]);

    end
end
fclose(fidout);

CP2kOptimPathParse(BaseFldr,system,[system '_' num2str(first) 'to' num2str(last) '_' num2str(nSampleSteps) 'step.xyz']);