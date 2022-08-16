clear all;
close all;

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus/Al_AlO/Al_water/Testing/'; 
system = 'Al_water';
Trajectory = 'Al_water_14300to18100_100step.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, ~, ~, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC);
[Atoms, AtomList, AtomIndx] = getAtomNamesFromInputXYZ(BaseFldr, system);

for snap = 39
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
   
    [VecAlAl, DistAlAl] = GetAtomCorrelation(XYZ_snap, AtomIndx.Al1, AtomIndx.Al1, ABC); 

end

for i= 1:length(AtomIndx.Al1)
    DoubleCount= 1:length(AtomIndx.Al1);
    DoubleCount(i)=[];
    AlNN_indx= find(DistAlAl(DoubleCount,i)<3.6)
end

%%%Bader Charge%%%

fldrname = [BaseFldr system '/Bader_Analysis/'];
ACFfiles = dir([fldrname 'ACF_*.dat']);

DoubleAnalType = 'MassDensity';
% DoubleAnalType = 'Radial';

% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx] = getAtomNamesFromInputXYZ(BaseFldr, system);

% % Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
Q = zeros(length(Atoms),length(ACFfiles));
Qnet = zeros(length(Atoms),length(ACFfiles));

for n = 1:length(ACFfiles)
    [Q(:,n), Qnet(:,n), StepNum(n)] = extractBaderCharges(fldrname, ACFfiles(n).name, Atoms, AtomList);
end