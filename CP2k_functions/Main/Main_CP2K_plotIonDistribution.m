clear all;
close all;

BaseFldr = 'D:\';
system = 'CP_Pit_20F';
Trajectory = 'CP_Pit_20F_43000to73000_500step.xyz';
grFO = [3.28 5.44];
grFH = [1.2 2.4];
grOH = [1.28];

% % get the ABC vectors from CP2K input file
ABC = getABCvectors(BaseFldr, system);

% % Read the xyz data from "Trajectory"
[xyz, XYZ, ~, ~, ~, nAtoms, startConfig, nConfigs, Step] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx, Indxfns, Kinds, Elements, PP] = getAtomInfoFromInput(BaseFldr, system);


% Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
fldrname = [BaseFldr system '\Bader\'];
ACFfiles = dir([fldrname 'ACF_*.dat']);
Q = zeros(length(Atoms),length(Step));
Qnet = zeros(length(Atoms),length(Step));

%% ACF files doesn't match step number - need to matche these in order to extract charges

Tot_SsO_1st = zeros(nConfigs,20);
Tot_SsO_2nd = zeros(nConfigs,20);
Tot_SsH_1st = {};

for snap = startConfig:1
disp(['Extracting charge from ' num2str(Step(snap))]);
    if exist([fldrname 'ACF_' num2str(Step(snap)) '.dat'])
        [Q(:,snap), Qnet(:,snap), StepNum(snap)] = extractBaderCharges(fldrname, ['ACF_' num2str(Step(snap)) '.dat'], Atoms, AtomList, Kinds, PP);
    end
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
    % % get the distances between F and O atoms in each snapshot "snap"
    [VecFO, DistFO] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.O, ABC);
    % % get the distances between F and H atoms in each snapshot "snap"
    [VecFH, DistFH] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.H, ABC);
    % % get the distances between H and O atoms in each snapshot "snap"
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.O, ABC);

    [SsO_1st, SsO_1stF] = find(DistFO <= grFO(1));
    [Hss1Indx, Hss1Indx_F] = find(DistFH <= grFH(1));

    [Hss2Indx, Hss2Indx_F] = find(DistFH > grFH(1) & DistFH <= grFH(2));

    [~, DistOH1st] = GetAtomCorrelation(XYZ_snap, Indx.O, Indx.H, ABC);
    [SsO_1st_H, SsO_1st_HIndx] = find(DistOH1st(:,SsO_1st) <= grOH(1));

    SsO_1st_H_unique = SsO_1st_H(~ismember(SsO_1st_H, Hss1Indx));
    SsO_1st_H_unique = SsO_1st_H_unique(~ismember(SsO_1st_H_unique, Hss2Indx));


    pIndx.H = Indx.H(Hss1Indx); % Closest H to F 
    pIndx.J = Indx.H(Hss2Indx); % 2nd Closest H to F
    pIndx.F = Indx.F([Hss1Indx_F; SsO_1stF]); % F
    pIndx.O = Indx.O(SsO_1st); % Closest H to F
    pIndx.Q = Indx.H(SsO_1st_H_unique); % Closest H to F

    make3dAtomPlot(ABC, XYZ_snap, pIndx, ['H'; 'J'; 'F'; 'O'; 'Q'])
    [CN_O, ~] = hist(SsO_1stF,unique(SsO_1stF));
    text(XYZ_snap(Indx.F, 1), XYZ_snap(Indx.F, 2), XYZ_snap(Indx.F, 3), arrayfun(@num2str, CN_O, 'UniformOutput', 0))

% % % % % % %     % % find the index of O atoms when F-O distances are less than the
% % % % % % %     % first (SsO_1st) and second (SsO_2nd) peaks in the global g(r) for F & O (grFO)
% % % % % % %     [SsO_1st, SsO_1stF] = find(DistFO <= grFO(1));
% % % % % % %     [SsO_2nd, SsO_2ndF] = find(DistFO > grFO(1) & DistFO <= grFO(2));
% % % % % % %     [~, DistOH1st] = GetAtomCorrelation(XYZ_snap, Indx.O, Indx.H, ABC);
% % % % % % %     [SsO_1st_H, SsO_1st_HIndx] = find(DistOH1st(:,SsO_1st) <= grOH(1));
% % % % % % % 
% % % % % % %     % % find the index of H atoms when F-H distances are less than the
% % % % % % %     % first (SsH_1st) peak in the global g(r) for F & H (grFH)
% % % % % % %     [SsH_1st, SsH_1stF] = find(DistFH <= grFH(1));
% % % % % % % 
% % % % % % % %     SsO_1stCharge(snap) = sum(Qnet(SsO_1st, snap));
% % % % % % % %     FCharge(snap) = mean(Qnet(Indx.F, snap));
% % % % % % % 
% % % % % % %     % % here is the charge on the O-F entities (sum, mean and variance)
% % % % % % %     F_OChargeSum(snap) = sum(Qnet(unique(Indx.F(SsO_1stF)), snap));
% % % % % % %     SsO_1stChargeSum(snap) = sum(Qnet(Indx.O(SsO_1st), snap));
% % % % % % %     F_OChargeAve(snap) = mean(Qnet(unique(Indx.F(SsO_1stF)), snap));
% % % % % % %     SsO_1stChargeAve(snap) = mean(Qnet(Indx.O(SsO_1st), snap));
% % % % % % %     F_OChargeVar(snap) = std(Qnet(unique(Indx.F(SsO_1stF)), snap));
% % % % % % %     SsO_1stChargeVar(snap) = std(Qnet(Indx.O(SsO_1st), snap));
% % % % % % % 
% % % % % % %     SsO_1st_HChargeSum(snap) = sum(Qnet(Indx.H(SsO_1st_H), snap));
% % % % % % %     SsO_1st_HChargeAve(snap) = mean(Qnet(Indx.H(SsO_1st_H), snap));
% % % % % % %     SsO_1st_HChargeVar(snap) = std(Qnet(Indx.H(SsO_1st_H), snap));
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % % % % here is the charge on the H-F entities (sum, mean and variance) then
% % % % % % % % not HF
% % % % % % %     F_HChargeSum(snap) = sum(Qnet(unique(Indx.F(SsH_1stF)), snap));
% % % % % % %     SsH_1stChargeSum(snap) = sum(Qnet(Indx.H(SsH_1st), snap));
% % % % % % %     F_HChargeAve(snap) = mean(Qnet(unique(Indx.F(SsH_1stF)), snap));
% % % % % % %     SsH_1stChargeAve(snap) = mean(Qnet(Indx.H(SsH_1st), snap));
% % % % % % %     F_HChargeVar(snap) = std(Qnet(unique(Indx.F(SsH_1stF)), snap));
% % % % % % %     SsH_1stChargeVar(snap) = std(Qnet(Indx.H(SsH_1st), snap));
% % % % % % % 
% % % % % % %     notIndxH = Indx.H;
% % % % % % %     notIndxH(SsH_1st) = [];
% % % % % % %     notIndxF = Indx.F;
% % % % % % %     notIndxF(SsH_1stF) = [];
% % % % % % %     not_F_HChargeSum(snap) = sum(Qnet(unique(notIndxF), snap));
% % % % % % %     not_SsH_1stChargeSum(snap) = sum(Qnet(notIndxH, snap));
% % % % % % %     not_F_HChargeAve(snap) = mean(Qnet(unique(notIndxF), snap));
% % % % % % %     not_SsH_1stChargeAve(snap) = mean(Qnet(notIndxH, snap));
% % % % % % %     not_F_HChargeVar(snap) = std(Qnet(unique(notIndxF), snap));
% % % % % % %     not_SsH_1stChargeVar(snap) = std(Qnet(notIndxH, snap));
% % % % % % % 
% % % % % % % %     F_HChargeSum(snap) = sum(Qnet(SsH_1st, snap));
% % % % % % % 
% % % % % % % 
% % % % % % %     [Tot_SsO_1st(snap,:), ~] = hist(SsO_1stF,unique(SsO_1stF));
% % % % % % %     [Tot_SsO_1st(snap,:), ~] = hist(SsO_1stF,unique(SsO_1stF));
% % % % % % %     [Tot_SsO_2nd(snap,:), ~] = hist(SsO_2ndF,unique(SsO_2ndF));
% % % % % % %     [Tot_SsH_1st{snap}, ~] = hist(SsH_1stF,unique(SsH_1stF));
% % % % % % % Tot_SsH_1st{snap} = nonzeros(Tot_SsH_1st{snap});

end

return
