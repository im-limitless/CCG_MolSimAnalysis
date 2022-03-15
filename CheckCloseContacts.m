% calculate close contacts
clc
clear all
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/LRZ_superMUC/Al_111_farming_Rashid/Farming/Final/EleventhTimeLucky/';
system = 'AlO_water_0.5ML';
Trajectory = 'AlO_water_0.5ML.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC);
dist = zeros(nAtoms, nAtoms);
for i = 1:nAtoms
    for j = 1:nAtoms
        if j ~= i
        dist(i,j) = sqrt(((XYZ(1,i,1)-XYZ(1,j,1))^2)+((XYZ(1,i,2)-XYZ(1,j,2))^2)+((XYZ(1,i,3)-XYZ(1,j,3))^2));
        else 
            dist(i,j) = 100;
        end
    end
end

min(dist, [], 'all')