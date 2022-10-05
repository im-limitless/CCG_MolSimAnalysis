clear all;  clc;
close all;


BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus2/Al_AlO/';
system = 'Al_water';
Trajectory = 'CP_Pit_20F_43000to73000_500step.xyz';

fldrname = [BaseFldr system '/Bader/'];
ACFfiles = dir([fldrname 'ACF_*.dat']);


% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx, Indxfns] = getAtomNamesFromInputXYZ(BaseFldr, system);

% % find the indices of atoms containing name in first arguement
myAtoms = {'Pt' 'O'};
[Indx, myAtomList, myAtomNums] = detectAtomsOfType(myAtoms, AtomList, Indx, Indxfns);

% % parse coordinates of atoms along trajectory and wrap into cell
[xyz, XYZ, ~, ~, ~, nAtoms, startConfig, nConfigs, StepNum_Traj] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC, [0 0 0]);

% [xyz1, XYZ1, ~, ~, ~, nAtoms1, startConfig1, nConfigs1, StepNum_Traj1] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC, [0 0 0]);
% [xyz2, XYZ2, ~, ~, ~, nAtoms2, startConfig2, nConfigs2, StepNum_Traj2] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0 0 0]);

% % Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
Q = zeros(length(Atoms),length(ACFfiles));
Qnet = zeros(length(Atoms),length(ACFfiles));

for n = 1:length(ACFfiles)
    [Q(:,n), Qnet(:,n), StepNum(n)] = extractBaderCharges(fldrname, ACFfiles(n).name, Atoms, AtomList);
end

MeanCharge = zeros(length(AtomList),length(StepNum));
StdCharge = zeros(length(AtomList),length(StepNum));
SumCharge = zeros(length(AtomList),length(StepNum));

for i = 1:length(StepNum)
    for j = 1:length(AtomList)
        MeanCharge(j,i) = mean(Qnet(Indx.(Indxfns{j}),i));
        StdCharge(j,i) = std(Qnet(Indx.(Indxfns{j}),i));
        SumCharge(j,i) = sum(Qnet(Indx.(Indxfns{j}),i));
    end
end

[StepNum,sortIdx] = sort(StepNum,'ascend');
SumCharge = SumCharge(:,sortIdx);
MeanCharge = MeanCharge(:,sortIdx);
StdCharge = StdCharge(:,sortIdx);

XYZ_snap = reshape(XYZ(1,:,:), [size(XYZ,2) size(XYZ,3)]);

MeanQnet = mean(Qnet(myAtomNums,end-20:end),2);
Bader3DCharge(XYZ_snap(myAtomNums,:), ABC, MeanQnet);
