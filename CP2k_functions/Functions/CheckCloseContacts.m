% calculate close contacts
clc
clear all
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge_2/Strt_dV/Send_to_ARCHER2_GC2/';
system = 'OH_0.83ML_dV';
Trajectory = 'OH_0.83ML_dV.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);

XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
XYZ_snap(:,:) = XYZ(1,:,:);

[Vec, Dist] = GetAtomCorrelation(XYZ_snap, 1:nAtoms, 1:nAtoms, ABC);

% Dist(boolean(eye(nAtoms)))=nan;
Dist = Dist - eye(size(Dist));
Dist(Dist==-1) = NaN;

minD = min(Dist, [], 'all');

[idx,~] = find(Dist == minD);

disp(['The closest atoms in the system are separated by ' num2str(minD, '%.2f') ' Angstroms']);
disp(['Atom 1 = ' Atoms(idx(1),:) ' ' num2str(idx(1)) ' at coordinates ' num2str(XYZ_snap(idx(1), :), '%.5f  ')]);
disp(['Atom 1 = ' Atoms(idx(2),:) ' ' num2str(idx(2))  ' at coordinates ' num2str(XYZ_snap(idx(2), :), '%.5f  ')]);

[r,c] = find(Dist < 0.1)