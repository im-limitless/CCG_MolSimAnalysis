clear all;  clc;
close all;
% BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\';
% system = 'CP_Like_1012_Fluoride';
% Trajectory = 'Sample34000_52000.xyz';

% BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\';
% system = 'CP_Like_1010_Fluorine';
% Trajectory = 'Sample50000_58000.xyz';

% BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CorrectVolume\';
% system = 'Pt_12H10F';
% Trajectory = 'Pt_12H10F_0to9500_500step.xyz';

BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
system = 'CP_Pit_20F';
Trajectory = 'Sample41000_52500.xyz';

% BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Vacuum\';
% system = 'Pt_Bulk';
% Trajectory = 'Pt_Bulk.xyz';


fldrname = [BaseFldr system '\Bader\'];
% ACFfiles = dir([fldrname 'ACF_*.dat']);

% DoubleAnalType = 'MassDensity';
DoubleAnalType = 'Radial';

% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx] = getAtomNamesFromInputXYZ(BaseFldr, system);

% % Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
% Q = zeros(length(Atoms),length(ACFfiles));
% Qnet = zeros(length(Atoms),length(ACFfiles));

% for n = 1:length(ACFfiles)
%     [Q(:,n), Qnet(:,n), StepNum(n)] = extractBaderCharges(fldrname, ACFfiles(n).name, Atoms, AtomList);
% end

% % Find indices of all "Al*" atoms
if size(AtomList,1) > 1
    AlList = find(ismember(AtomList, 'Al'));
    AlList = AlList(AlList < length(AtomList)+1);
    Indxfns = fieldnames(Indx);
    Indx.AlAll = [];
    for ii = 1:length(AlList)
        Indx.AlAll = [Indx.AlAll; Indx.(Indxfns{AlList(ii)})];
    end
else
    Indxfns = fieldnames(Indx);
    AlList = 1;
    Indx.AlAll = Indx.(Indxfns{1});
end

% % parse coordinates of atoms along trajectory and wrap into cell
[xyz, XYZ, ~, ~, ~, nAtoms, startConfig, nConfigs, StepNum_Traj] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC);
XYZ = wrapXYZ(XYZ, ABC);

% % compute radial functions for OH, OF and HF
RadFunOH = cell(nConfigs,1);
% RadFunFH = cell(nConfigs,1);
% RadFunFO = cell(nConfigs,1);
RadFunAlO = cell(nConfigs,1);

DistOH = cell(1,nConfigs);
% DistFH = cell(1,nConfigs);
% DistFO = cell(1,nConfigs);
DistAlO = cell(1,nConfigs);

for snap = startConfig:nConfigs
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
    [VecAlO, DistAlO{snap}] = GetAtomCorrelation(XYZ_snap, [Indx.Als; Indx.AlE], Indx.O, ABC);
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.O, Indx.H, ABC);
%     [VecFH, DistFH] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.H, ABC);
%     [VecFO, DistFO] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.O, ABC);
    
    RadFunOH{snap} = reshape(DistOH, [numel(DistOH), 1]);
%     RadFunFH{snap} = reshape(DistFH, [numel(DistFH), 1]);
%     RadFunFO{snap} = reshape(DistFO, [numel(DistFO), 1]);
    RadFunAlO{snap} = reshape(DistAlO{snap}, [numel(DistAlO{snap}), 1]);
end

MinimaOH = RadialDistribution(RadFunOH, ABC, ['O'; 'H'], 0);
% MinimaFH = RadialDistribution(RadFunFH, ABC, ['F'; 'H'], 0);
% MinimaFO = RadialDistribution(RadFunFO, ABC, ['F'; 'O'], 0);
MinimaAlO = RadialDistribution(RadFunAlO, ABC, ['Al'; 'O '], 1);

if strcmp(DoubleAnalType, 'MassDensity')
    
    disp('Determining water layering from mass density profile...');
    
    % % get the O atom distribution and corresponding indices of DL atoms
    [Dens_O, ~, ~, ~, ~, z] = getDensityProfile(xyz, ABC);
    [FirstLayerIndx, SecondLayerIndx] = getWaterLayerIndices(Indx, XYZ, Dens_O, z);
    
    for i = startConfig:nConfigs
        DL1st{i} = [FirstLayerIndx{i}];
        DL2nd{i} = [SecondLayerIndx{i}];
        nonDL{i} = setdiff(Indx.O, [DL1st{i}; DL2nd{i}]);
        
        XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
        XYZ_snap(:,:) = XYZ(i,:,:);
        [~, DistOH1stWL] = GetAtomCorrelation(XYZ_snap, DL1st{i}, Indx.H, ABC);
        [~, DistOH2ndWL] = GetAtomCorrelation(XYZ_snap, DL2nd{i}, Indx.H, ABC);
        [~, DistOHnonDL] = GetAtomCorrelation(XYZ_snap, nonDL{i}, Indx.H, ABC);
        
        for j = 1:length(DL1st{i})
            DL1st{i} = [DL1st{i}; Indx.H(find(DistOH1stWL(:,j)<MinimaOH(1)))];
        end
        
        for j = 1:length(DL2nd{i})
            DL2nd{i} = [DL2nd{i}; Indx.H(find(DistOH2ndWL(:,j)<MinimaOH(1)))];
        end
        
        for j = 1:length(nonDL{i})
            nonDL{i} = [nonDL{i}; Indx.H(find(DistOHnonDL(:,j)<MinimaOH(1)))];
        end
    end
    
elseif strcmp(DoubleAnalType, 'Radial')
    
    disp('Determining water layering from radial distribution...');
    
    for i = startConfig:nConfigs
        [r1st, ~] = find(DistAlO{i} <= MinimaAlO(1));
        DL1st{i} = Indx.O(unique(r1st));
        [r2nd, ~] = find(DistAlO{i} <= MinimaAlO(2) & DistAlO{i} > MinimaAlO(1));
        DL2nd{i} = Indx.O(unique(r2nd));
        nonDL{i} = setdiff(Indx.O, [DL1st{i}; DL2nd{i}]);
        
        XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
        XYZ_snap(:,:) = XYZ(i,:,:);
        
        [~, DistOH1stWL] = GetAtomCorrelation(XYZ_snap, DL1st{i}, Indx.H, ABC);
        [~, DistOH2ndWL] = GetAtomCorrelation(XYZ_snap, DL2nd{i}, Indx.H, ABC);
        [~, DistOHnonDL] = GetAtomCorrelation(XYZ_snap, nonDL{i}, Indx.H, ABC);
        
        for j = 1:length(DL1st{i})
            DL1st{i} = [DL1st{i}; Indx.H(find(DistOH1stWL(:,j)<=MinimaOH(1)))];
        end
        
        for j = 1:length(DL2nd{i})
            
            DL2nd{i} = [DL2nd{i}; Indx.H(find(DistOH2ndWL(:,j)<=MinimaOH(1)))];
        end
        
        for j = 1:length(nonDL{i})
            nonDL{i} = [nonDL{i}; Indx.H(find(DistOHnonDL(:,j)<=MinimaOH(1)))];
        end
        
        nonDL{i} = unique(nonDL{i});
        
    end
end

% MeanCharge = zeros(length(AtomList),length(StepNum));
% StdCharge = zeros(length(AtomList),length(StepNum));
% SumCharge = zeros(length(AtomList),length(StepNum));

% for i = 1:length(StepNum)
%     if any(StepNum(i) == StepNum_Traj)
%         Inter = find(StepNum(i) == StepNum_Traj);
%         DL1st_sum(i) = sum(Qnet(DL1st{Inter},i));
%         DL2nd_sum(i) = sum(Qnet(DL2nd{Inter},i));
%         XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
%         XYZ_snap(:,:) = XYZ(Inter,:,:);
%         writeSnaptoxyz(BaseFldr, system, StepNum(i), XYZ_snap, Atoms, [DL1st{Inter}; Indx.AlAll] , DoubleAnalType)
%     else % redundant?
%         DL1st_sum(i) = 0;
%         DL2nd_sum(i) = 0;
%     end
    
%     for j = 1:length(AtomList)
%         MeanCharge(j,i) = mean(Qnet(Indx.(Indxfns{j}),i));
%         StdCharge(j,i) = std(Qnet(Indx.(Indxfns{j}),i));
%         SumCharge(j,i) = sum(Qnet(Indx.(Indxfns{j}),i));
%     end
% end

[StepNum,sortIdx] = sort(StepNum,'ascend');
% SumCharge = SumCharge(:,sortIdx);
% MeanCharge = MeanCharge(:,sortIdx);
% StdCharge = StdCharge(:,sortIdx);
DL1st_sum = DL1st_sum(sortIdx);
DL2nd_sum = DL2nd_sum(sortIdx);

