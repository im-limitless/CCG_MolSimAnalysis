clear all;
close all;

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/5lyr_systems/Al_AlO_OH/AlO_OH/';
system = 'AlO_1ML_OH';
Trajectory = 'AlO_1ML_OH_0to123000_1000step.xyz';

% % get the names of atoms from original xyz input file
[~, ~, AtomIndx, ~, ~, ~, ~] = getAtomInfoFromInput(BaseFldr, system);

% % get the ABC vectors from CP2K input file
ABC = getABCvectors(BaseFldr, system);

% % Read the xyz data from "Trajectory"
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);


RadFunOH = cell(nConfigs,1);
RadFunFH = cell(nConfigs,1);
RadFunFO = cell(nConfigs,1);
RadFunHH = cell(nConfigs,1);
RadFunOO = cell(nConfigs,1);
RadFunPtO = cell(nConfigs,1);
RadFunPtEO = cell(nConfigs,1);
RadFunPtSO = cell(nConfigs,1);
RadFunAlO = cell(nConfigs,1);
% DistPtO = cell(1,nConfigs);


% for snap = startConfig:startConfig
for snap = startConfig:nConfigs
%     disp(['Processing snapshot ' num2str(StepNum(snap)) ' - ' num2str(100*(snap/nConfigs)) ' % complete']);
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
% % %     AtomIndx.Omid = find(XYZ_snap(AtomIndx.O,3) < ABC(3)/2+5 & XYZ_snap(AtomIndx.O,3) > ABC(3)/2-5);
% % %    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.Omid, AtomIndx.H, ABC);
% % %    [VecOO, DistOO] = GetAtomCorrelation(XYZ_snap, AtomIndx.Omid, AtomIndx.O, ABC);

    
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.O, AtomIndx.H, ABC);
%     [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, [AtomIndx.O; AtomIndx.OtU; AtomIndx.OtL], AtomIndx.H, ABC);
%     [VecFH, DistFH] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.H, ABC);
%     [VecFO, DistFO] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.O, ABC);
%     [VecHH, DistHH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.H, ABC);
    [VecOO, DistOO] = GetAtomCorrelation(XYZ_snap, AtomIndx.O, AtomIndx.O, ABC);
    [VecAlO, DistAlO] = GetAtomCorrelation(XYZ_snap, AtomIndx.Al1, AtomIndx.O, ABC);
%     [VecPtO, DistPtO] = GetAtomCorrelation(XYZ_snap, AtomIndx.Pts, Indx.O, ABC);

% [VecPt, DistPt] = GetAtomCorrelation(XYZ_snap, AtomIndx.Pts, [AtomIndx.Pts; AtomIndx.Ptss], ABC);
%     RadFunOH{snap} = reshape(DistPt, [numel(DistPt), 1]);
%     
    RadFunOH{snap} = reshape(DistOH, [numel(DistOH), 1]);
%     RadFunFH{snap} = reshape(DistFH, [numel(DistFH), 1]);
%     RadFunFO{snap} = reshape(DistFO, [numel(DistFO), 1]);
%     RadFunHH{snap} = reshape(DistHH, [numel(DistHH), 1]);
    RadFunOO{snap} = reshape(DistOO, [numel(DistOO), 1]);
    RadFunAlO{snap} = reshape(DistAlO, [numel(DistAlO), 1]);
%     RadFunPtO{snap} = reshape(DistPtO{snap}, [numel(DistPtO{snap}), 1]);
    
% [VecPtSO, DistPtSO] = GetAtomCorrelation(XYZ_snap, AtomIndx.Pts, Indx.O, ABC);
%         RadFunPtSO{snap} = reshape(DistPtSO, [numel(DistPtSO), 1]);
% 
%         [VecPtEO, DistPtEO] = GetAtomCorrelation(XYZ_snap, AtomIndx.PtE, Indx.O, ABC);
%         RadFunPtEO{snap} = reshape(DistPtEO, [numel(DistPtEO), 1]);
%         [r1all{snap}, ~] = find(DistPtEO <= MinimaPtEO(1)+0.1);
end

RadialDistribution_new(RadFunOH, ABC, ['O'; 'H'], 1);
% RadialDistribution(RadFunFH, ABC, ['H'; 'F'], 1);
% RadialDistribution(RadFunFO, ABC, ['O'; 'F'], 1);
% RadialDistribution(RadFunHH, ABC, ['H'; 'H'], 1);
RadialDistribution_new(RadFunOO, ABC, ['O'; 'O'], 1);
RadialDistribution_new(RadFunAlO, ABC, ['Al'; 'O '], 1);  % "Dimensions of arrays being concatenated are not consistent." Error
% MinimaPtSO = RadialDistribution(RadFunPtSO, ABC, ['Pts'; 'O  '],1);
% [r1stflat, ~] = find(DistPtSO <= MinimaPtSO(1));
% 
% MinimaPtEO = RadialDistribution(RadFunPtEO, ABC, ['PtE'; 'O  '],1);
% [r1st, ~] = find(DistPtEO <= MinimaPtEO(1));
% 
% 
% [~, DistOH1stWL] = GetAtomCorrelation(XYZ_snap, AtomIndx.O(r1st), Indx.H, ABC);
 
 return

        DL1st{i} = Indx.O(unique(r1st));
        [r2nd, ~] = find(DistPtO{i} <= MinimaPtO(2) & DistPtO{i} > MinimaPtO(1));
        DL2nd{i} = Indx.O(unique(r2nd));
        nonDL{i} = setdiff(Indx.O, [DL1st{i}; DL2nd{i}]);
        
        XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
        XYZ_snap(:,:) = XYZ(i,:,:);
        
        [~, DistOH1stWL] = GetAtomCorrelation(XYZ_snap, DL1st{i}, Indx.H, ABC);
        [~, DistOH2ndWL] = GetAtomCorrelation(XYZ_snap, DL2nd{i}, Indx.H, ABC);
        [~, DistOHnonDL] = GetAtomCorrelation(XYZ_snap, nonDL{i}, Indx.H, ABC);
        
        
        AvePos = mean(xyz.O(:, r1st, :),1)
        figure
        hold on
%         plot3(AvePos(1,:,1),AvePos(1,:,2),AvePos(1,:,3),'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 15)
%         text(AvePos(1,:,1),AvePos(1,:,2),AvePos(1,:,3), num2str(r1st))
        for i = 1:length(r1st)
            plot3(xyz.O(:,r1st(i),1), xyz.O(:,r1st(i),2), xyz.O(:,r1st(i),3), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 5)
        end
        %         plot3(xyz.O(:,614,1), xyz.O(:,614,2), xyz.O(:,614,3), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markersize', 15)
        
        r1stflat = setdiff(r1stflat,r1st);
        AvePos1stflat = mean(xyz.O(:, r1stflat, :),1)
%         plot3(AvePos1stflat(1,:,1),AvePos1stflat(1,:,2),AvePos1stflat(1,:,3),'o', 'markerfacecolor', 'm', 'markeredgecolor', 'k', 'markersize', 15)
%        plot3(xyz.O(end,r1stflat,1), xyz.O(end,r1stflat,2), xyz.O(end,r1stflat,3), 'o', 'markerfacecolor', 'm', 'markeredgecolor', 'k', 'markersize', 15)
        
        for i = 1:length(r1stflat)
            plot3(xyz.O(:,r1stflat(i),1), xyz.O(:,r1stflat(i),2), xyz.O(:,r1stflat(i),3), 'o', 'markerfacecolor', 'm', 'markeredgecolor', 'k', 'markersize', 5)
        end
        
                AvePosPtE = mean(xyz.Pt(:, AtomIndx.PtE, :),1)
        plot3(AvePosPtE(1,:,1),AvePosPtE(1,:,2),AvePosPtE(1,:,3),'o', 'markerfacecolor', 'c', 'markeredgecolor', 'k', 'markersize', 30)
        
        AvePosPtS = mean(xyz.Pt(:, AtomIndx.Pts, :),1)
        plot3(AvePosPtS(1,:,1),AvePosPtS(1,:,2),AvePosPtS(1,:,3),'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markersize', 30)

        
        [r2nd, ~] = find(DistPtEO > MinimaPtEO(1) & DistPtEO <= MinimaPtEO(2));
        r2nd = setdiff(r2nd,r1st);
        AvePos2nd = mean(xyz.O(:, r2nd, :),1)
        
        plot3(AvePos2nd(1,:,1),AvePos2nd(1,:,2),AvePos2nd(1,:,3),'o', 'markerfacecolor', [0 0.5 0], 'markeredgecolor', 'k', 'markersize', 15)
        
        

        
