clear all;
close all;

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus/Al_AlO/'; 
system = 'Al_water';
Trajectory = 'Al_water_14300to18100_100step.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, ~, ~, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC);
[Atoms, AtomList, AtomIndx] = getAtomNamesFromInputXYZ(BaseFldr, system);

RadFunOH = cell(nConfigs,1);
RadFunFH = cell(nConfigs,1);
RadFunFO = cell(nConfigs,1);
RadFunHH = cell(nConfigs,1);
RadFunOO = cell(nConfigs,1);
RadFunPtO = cell(nConfigs,1);
RadFunPtEO = cell(nConfigs,1);
DistPtO = cell(1,nConfigs);

% for snap = startConfig:startConfig
for snap = startConfig:nConfigs
    disp(['Processing snapshot ' num2str(StepNum(snap)) ' - ' num2str(100*(snap/nConfigs)) ' % complete']);
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
%     AtomIndx.Omid = find(XYZ_snap(AtomIndx.O,3) < ABC(3)/2+5 & XYZ_snap(AtomIndx.O,3) > ABC(3)/2-5);
%    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.Omid, AtomIndx.H, ABC);
%    [VecOO, DistOO] = GetAtomCorrelation(XYZ_snap, AtomIndx.Omid, AtomIndx.O, ABC);
   
    
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.O, AtomIndx.H, ABC);
%     [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, [AtomIndx.O; AtomIndx.OtU; AtomIndx.OtL], AtomIndx.H, ABC);
%     [VecFH, DistFH] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.H, ABC);
%     [VecFO, DistFO] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.O, ABC);
%     [VecHH, DistHH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.H, ABC);
    [VecOO, DistOO] = GetAtomCorrelation(XYZ_snap, AtomIndx.O, AtomIndx.O, ABC);
%     [VecPtO, DistPtO{snap}] = GetAtomCorrelation(XYZ_snap, Indx.Pt, Indx.O, ABC);
%     
    RadFunOH{snap} = reshape(DistOH, [numel(DistOH), 1]);
%     RadFunFH{snap} = reshape(DistFH, [numel(DistFH), 1]);
%     RadFunFO{snap} = reshape(DistFO, [numel(DistFO), 1]);
%     RadFunHH{snap} = reshape(DistHH, [numel(DistHH), 1]);
    RadFunOO{snap} = reshape(DistOO, [numel(DistOO), 1]);
%     RadFunPtO{snap} = reshape(DistPtO{snap}, [numel(DistPtO{snap}), 1]);
    
%         [VecPtEO, DistPtEO] = GetAtomCorrelation(XYZ_snap, AtomIndx.PtE, Indx.O, ABC);
%         RadFunPtEO{snap} = reshape(DistPtEO, [numel(DistPtEO), 1]);
end

RadialDistribution(RadFunOH, ABC, ['O'; 'H'], 1);
% RadialDistribution(RadFunFH, ABC, ['H'; 'F'], 1);
% RadialDistribution(RadFunFO, ABC, ['O'; 'F'], 1);
% RadialDistribution(RadFunHH, ABC, ['H'; 'H'], 1);
RadialDistribution(RadFunOO, ABC, ['O'; 'O'], 1);
% RadialDistribution(RadFunPtO, ABC, ['Pt'; 'O '],1);

% RadialDistribution(RadFunPtEO, ABC, ['PtE'; 'O  '],1)

%%%%%%%%%%%%%%%%% Counting the number of O-H bonds per oxygen atom %%%%%%%%%%%%%%%%%

% the value after DistOH<# is based on the graph of the raidal distribution
% of OH where the value is the extent of the first peak (entered manually
% for now)
OH_bonds= DistOH<1.3;

no_OH_bonds= sum(sum(OH_bonds)==1);
no_H2O_bonds= sum(sum(OH_bonds)==2);
no_H3O_bonds= sum(sum(OH_bonds)==3);