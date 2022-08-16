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