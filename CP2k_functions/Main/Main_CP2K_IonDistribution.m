clear all;
close all;

BaseFldr = 'D:\';
system = 'CP_Pit_20F';
Trajectory = 'CP_Pit_20F_43000to73000_500step.xyz';
grFO = [3.28 5.44];
grFH = [1.2];
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

for snap = startConfig:nConfigs
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

    % % find the index of O atoms when F-O distances are less than the
    % first (SsO_1st) and second (SsO_2nd) peaks in the global g(r) for F & O (grFO)
    [SsO_1st, SsO_1stF] = find(DistFO <= grFO(1));
    [SsO_2nd, SsO_2ndF] = find(DistFO > grFO(1) & DistFO <= grFO(2));
    [~, DistOH1st] = GetAtomCorrelation(XYZ_snap, Indx.O, Indx.H, ABC);
    [SsO_1st_H, SsO_1st_HIndx] = find(DistOH1st(:,SsO_1st) <= grOH(1));

    % % find the index of H atoms when F-H distances are less than the
    % first (SsH_1st) peak in the global g(r) for F & H (grFH)
    [SsH_1st, SsH_1stF] = find(DistFH <= grFH(1));

%     SsO_1stCharge(snap) = sum(Qnet(SsO_1st, snap));
%     FCharge(snap) = mean(Qnet(Indx.F, snap));

    % % here is the charge on the O-F entities (sum, mean and variance)
    F_OChargeSum(snap) = sum(Qnet(unique(Indx.F(SsO_1stF)), snap));
    SsO_1stChargeSum(snap) = sum(Qnet(Indx.O(SsO_1st), snap));
    F_OChargeAve(snap) = mean(Qnet(unique(Indx.F(SsO_1stF)), snap));
    SsO_1stChargeAve(snap) = mean(Qnet(Indx.O(SsO_1st), snap));
    F_OChargeVar(snap) = std(Qnet(unique(Indx.F(SsO_1stF)), snap));
    SsO_1stChargeVar(snap) = std(Qnet(Indx.O(SsO_1st), snap));

    SsO_1st_HChargeSum(snap) = sum(Qnet(Indx.H(SsO_1st_H), snap));
    SsO_1st_HChargeAve(snap) = mean(Qnet(Indx.H(SsO_1st_H), snap));
    SsO_1st_HChargeVar(snap) = std(Qnet(Indx.H(SsO_1st_H), snap));



% % here is the charge on the H-F entities (sum, mean and variance) then
% not HF
    F_HChargeSum(snap) = sum(Qnet(unique(Indx.F(SsH_1stF)), snap));
    SsH_1stChargeSum(snap) = sum(Qnet(Indx.H(SsH_1st), snap));
    F_HChargeAve(snap) = mean(Qnet(unique(Indx.F(SsH_1stF)), snap));
    SsH_1stChargeAve(snap) = mean(Qnet(Indx.H(SsH_1st), snap));
    F_HChargeVar(snap) = std(Qnet(unique(Indx.F(SsH_1stF)), snap));
    SsH_1stChargeVar(snap) = std(Qnet(Indx.H(SsH_1st), snap));

    notIndxH = Indx.H;
    notIndxH(SsH_1st) = [];
    notIndxF = Indx.F;
    notIndxF(SsH_1stF) = [];
    not_F_HChargeSum(snap) = sum(Qnet(unique(notIndxF), snap));
    not_SsH_1stChargeSum(snap) = sum(Qnet(notIndxH, snap));
    not_F_HChargeAve(snap) = mean(Qnet(unique(notIndxF), snap));
    not_SsH_1stChargeAve(snap) = mean(Qnet(notIndxH, snap));
    not_F_HChargeVar(snap) = std(Qnet(unique(notIndxF), snap));
    not_SsH_1stChargeVar(snap) = std(Qnet(notIndxH, snap));

%     F_HChargeSum(snap) = sum(Qnet(SsH_1st, snap));


    [Tot_SsO_1st(snap,:), ~] = hist(SsO_1stF,unique(SsO_1stF));
    [Tot_SsO_1st(snap,:), ~] = hist(SsO_1stF,unique(SsO_1stF));
    [Tot_SsO_2nd(snap,:), ~] = hist(SsO_2ndF,unique(SsO_2ndF));
    [Tot_SsH_1st{snap}, ~] = hist(SsH_1stF,unique(SsH_1stF));
Tot_SsH_1st{snap} = nonzeros(Tot_SsH_1st{snap});

end

% % plot the charge of F and O in bulk/1st SS vs. time
figure
xlabel('Time (ps)');
ylabel('Bader Charge (e)');
set(gcf, 'position', [940   197   560   420])
hold on
errorbar((1:52)*500/2000 , SsH_1stChargeAve(1:52), SsH_1stChargeVar(1:52), '-ok', 'markerfacecolor', 'b');
errorbar((1:52)*500/2000 , not_SsH_1stChargeAve(1:52), not_SsH_1stChargeVar(1:52), '-ok', 'markerfacecolor', 'c');
hold off
legend('1st Solv H', 'Water H')


figure
xlabel('Time (ps)');
ylabel('Bader Charge (e)');
set(gcf, 'position', [940   197   560   420])
hold on
errorbar((1:52)*500/2000 , F_OChargeAve(1:52), F_OChargeVar(1:52), '-ok', 'markerfacecolor', 'r');
errorbar((1:52)*500/2000 , SsO_1stChargeAve(1:52), SsO_1stChargeVar(1:52), '-ok', 'markerfacecolor', 'b');
errorbar((1:52)*500/2000 , F_OChargeAve(1:52) + SsO_1stChargeAve(1:52), mean([F_OChargeVar(1:52); SsO_1stChargeVar(1:52)]), '-or', 'markerfacecolor', 'b');
errorbar((1:52)*500/2000 , F_OChargeAve(1:52) + SsO_1stChargeAve(1:52) + SsO_1st_HChargeAve(1:52), mean([F_OChargeVar(1:52); SsO_1stChargeVar(1:52); SsO_1st_HChargeVar(1:52)]), '-ok', 'markerfacecolor', 'c');
legend('F (F-OH_x)', 'O (F-OHx)', 'F-O (F-OHx)', 'F-OHx')
hold off


figure
xlabel('Time (ps)');
ylabel('Bader Charge (e)');
set(gcf, 'position', [940   197   560   420])
hold on
errorbar((1:52)*500/2000 , F_HChargeAve(1:52), F_HChargeVar(1:52), '-ok', 'markerfacecolor', 'r');
errorbar((1:52)*500/2000 , not_F_HChargeAve(1:52), not_F_HChargeVar(1:52), '-ok', 'markerfacecolor', 'm');
hold off
legend('F (HF)', 'F (other)')

figure 
xlabel('Time (ps)');
ylabel('Bader Charge (e)');
set(gcf, 'position', [940   197   560   420])
hold on
plot((1:52)*500/2000 , SsH_1stChargeSum(1:52), '-ok', 'markerfacecolor', 'b');
plot((1:52)*500/2000 , F_HChargeSum(1:52), '-ok', 'markerfacecolor', 'r');
plot((1:52)*500/2000 , F_HChargeSum(1:52)+SsH_1stChargeSum(1:52), '-or', 'markerfacecolor', 'b');
legend('H (HF)', 'F (HF)', 'HF total')
hold off

% % plot the charge of F and O in bulk/1st SS vs. time
figure
set(gcf, 'position', [940   197   560   420])
hold on
errorbar((1:52)*500/2000 , mean(Qnet(Indx.O, 1:52)), std(Qnet(Indx.O, 1:52)), '-ok', 'markerfacecolor', 'm');
errorbar((1:52)*500/2000 ,SsO_1stChargeAve(1:52), SsO_1stChargeVar(1:52), '-ok', 'markerfacecolor', 'b');
xlabel('Time (ps)');
ylabel('Bader Charge (e)');
legend('Water O', '1st Solv O')
hold off


% % plot the H 1st SS with time

figure
hold on
xlabel('Time (ps)');
ylabel('Number of H-F');
set(gcf, 'position', [940   197   560   420]);
plot((1:length(cellfun('length', Tot_SsH_1st)))*500/2000, cellfun('length', Tot_SsH_1st), '-ok', 'markerfacecolor', 'r');
plot([0 (length(cellfun('length', Tot_SsH_1st)))*500/2000], [mean(cellfun('length', Tot_SsH_1st)) mean(cellfun('length', Tot_SsH_1st))], '--', 'color', [0.8 0.8 0.8]);


figure
xlabel('Time (ps)');
ylabel('Fraction of F-H')
set(gcf, 'position', [940   197   560   420])
plot((1:length(cellfun('length', Tot_SsH_1st)))*500/2000, cellfun('length', Tot_SsH_1st)/20, '-ok', 'markerfacecolor', 'r')





Ave_SsO_1st =  mean(Tot_SsO_1st, 'all');
Var_SsO_1st = std(Tot_SsO_1st(:));

Ave_SsO_2nd =  mean(Tot_SsO_2nd, 'all');
Var_SsO_2nd = std(Tot_SsO_2nd(:));