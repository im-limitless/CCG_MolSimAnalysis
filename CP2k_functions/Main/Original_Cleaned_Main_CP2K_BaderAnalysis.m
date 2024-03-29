clear all;  clc;
close all;

%%%%%%%%%%%%%% Data collections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge_2/Phase_diagram_sys/';
system = 'AlO_1ML_OH';
Trajectory = 'AlO_1ML_OH_95000to104000_1000step.xyz';


fldrname = [BaseFldr system '/Bader_Analysis/'];
ACFfiles = dir([fldrname 'ACF_*.dat']);

DoubleAnalType = 'MassDensity';
% DoubleAnalType = 'Radial';

% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx, Indxfns, Kinds, Elements, PP] = getAtomInfoFromInput(BaseFldr, system);

% % Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
Q = zeros(length(Atoms),length(ACFfiles));
Qnet = zeros(length(Atoms),length(ACFfiles));

for n = 1:length(ACFfiles)
    [Q(:,n), Qnet(:,n), StepNum(n)] = extractBaderCharges(fldrname, ACFfiles(n).name, Atoms, AtomList,Kinds,PP);
end

% % Find indices of all "Al*" atoms
if size(AtomList,1) > 1
    AlList = find(ismember(AtomList, 'Al'));
    AlList = AlList(AlList < length(AtomList)+1);
    Indxfns = fieldnames(Indx);
    Indx.Al_All = [];
    for ii = 1:length(AlList)
        Indx.Al_All = [Indx.Al_All; Indx.(Indxfns{AlList(ii)})];
    end
else
    Indxfns = fieldnames(Indx);
    AlList = 1;
    Indx.Al_All = Indx.(Indxfns{1});
end

% % parse coordinates of atoms along trajectory and wrap into cell
[xyz, XYZ, ~, ~, ~, nAtoms, startConfig, nConfigs, StepNum_Traj] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0 0 0]);
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

    [VecAlO, DistAlO{snap}] = GetAtomCorrelation(XYZ_snap, [Indx.Al1], Indx.O, ABC);
%     [VecAlO, DistAlO{snap}] = GetAtomCorrelation(XYZ_snap, [Indx.Al11], Indx.O, ABC);
%     [VecAlO11, DistAlO11{snap}] = GetAtomCorrelation(XYZ_snap, [Indx.Al11], Indx.O, ABC);
%     [VecAlO12, DistAlO12{snap}] = GetAtomCorrelation(XYZ_snap, [Indx.Al12], Indx.O, ABC);
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.O, Indx.H, ABC);
%     [VecFH, DistFH] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.H, ABC);
%     [VecFO, DistFO] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.O, ABC);
    
    RadFunOH{snap} = reshape(DistOH, [numel(DistOH), 1]);
% %     RadFunFH{snap} = reshape(DistFH, [numel(DistFH), 1]);
% %     RadFunFO{snap} = reshape(DistFO, [numel(DistFO), 1]);
    RadFunAlO{snap} = reshape(DistAlO{snap}, [numel(DistAlO{snap}), 1]);

end

MinimaOH = RadialDistribution(RadFunOH, ABC, ['O'; 'H'], 1);
% % MinimaFH = RadialDistribution(RadFunFH, ABC, ['F'; 'H'], 0);
% % MinimaFO = RadialDistribution(RadFunFO, ABC, ['F'; 'O'], 0);
MinimaAlO = RadialDistribution(RadFunAlO, ABC, ['Al'; 'O '], 1);


if strcmp(DoubleAnalType, 'MassDensity')
    
    disp('Determining water layering from mass density profile...');
    
    % % get the O atom distribution and corresponding indices of DL atoms
    [Dens_O, ~, ~, ~, z] = getDensityProfile(xyz, ABC);
    %change getwaterlayerIndices to the name of the new function (same name
    %but add persnap at the end) +  [same same minimaZ] = same
    [FirstLayerIndx, SecondLayerIndx] = getWaterLayerIndicesPerSnap_new(Indx, XYZ, Dens_O, z);
    
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

MeanCharge = zeros(length(AtomList),length(StepNum));
StdCharge = zeros(length(AtomList),length(StepNum));
SumCharge = zeros(length(AtomList),length(StepNum));
SumChargeAls= zeros(length(AtomList),length(StepNum));
totalAlsCharge= zeros(1,length(StepNum));

for i = 1:length(StepNum)
    if any(StepNum(i) == StepNum_Traj)
        Inter = find(StepNum(i) == StepNum_Traj);
        DL1st_sum(i) = sum(Qnet(DL1st{Inter},i));
        DL2nd_sum(i) = sum(Qnet(DL2nd{Inter},i));
        nonDL_sum(i) = sum(Qnet(nonDL{Inter},i));
        XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
        XYZ_snap(:,:) = XYZ(Inter,:,:);
        writeSnaptoxyz(BaseFldr, system, StepNum(i), XYZ_snap, Atoms, [DL1st{Inter}; Indx.Al_All] , DoubleAnalType)
    else % redundant?
        DL1st_sum(i) = 0;
        DL2nd_sum(i) = 0;
    end
    
    for j = 1:length(AtomList)
        MeanCharge(j,i) = mean(Qnet(Indx.(Indxfns{j}),i));
        StdCharge(j,i) = std(Qnet(Indx.(Indxfns{j}),i));
        SumCharge(j,i) = sum(Qnet(Indx.(Indxfns{j}),i));
        if contains(Indxfns{j},'Al')
            SumChargeAls(j,i) = sum(Qnet(Indx.(Indxfns{j}),i));
        end
    end
end

for i = 1:length(StepNum)

    totalAlsCharge(:,i)=sum(SumChargeAls(:,i));
end

[StepNum,sortIdx] = sort(StepNum,'ascend');
SumCharge = SumCharge(:,sortIdx);
SumChargeAls = SumChargeAls(:,sortIdx);
MeanCharge = MeanCharge(:,sortIdx);
StdCharge = StdCharge(:,sortIdx);
DL1st_sum = DL1st_sum(sortIdx);
DL2nd_sum = DL2nd_sum(sortIdx);
%nonDL_sum = nonDL_sum(sortIdx);
%%%%%%%%%%%%%% End of data collection, below are graphs and data
%%%%%%%%%%%%%% illustrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ----- Al ----- %%%%
% Total charge per Al type
figure
box on
hold on
xlabel('Time (ps)');
ylabel('Total Charge (e)');
% title(['Total Bader charge for Al Atoms in ' system], 'interpreter', 'none')
C = [218/255 165/255 32/255; 0.25 0.25 0.25; 0 0.5 0; 0 0 1; 1 0 0];
AlC = [];
for i = 1:length(AlList)
    if contains(AtomList(AlList(i), :), 'Al11')
        AlC = C(1,:);
    elseif contains(AtomList(AlList(i), :), 'Al12')
        AlC = C(2,:);
    elseif contains(AtomList(AlList(i), :), 'Al21')
        AlC = C(3,:);
    elseif contains(AtomList(AlList(i), :), 'Al22')
        AlC = C(4,:);
    elseif contains(AtomList(AlList(i), :), 'Alb')
        AlC = C(5,:);
    end
    plot(StepNum/2000, SumCharge(AlList(i),:), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', AlC);
end
if length(AlList) == 2
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), 'interpreter', 'tex')
elseif length(AlList) == 3
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), 'interpreter', 'tex')
elseif length(AlList) == 5
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), AtomList(AlList(4),:), AtomList(AlList(5),:), 'interpreter', 'tex')
end

for i = 1:length(AlList)
    if contains(AtomList(AlList(i), :), 'Al11')
        AlC = C(1,:);
    elseif contains(AtomList(AlList(i), :), 'Al12')
        AlC = C(2,:);
    elseif contains(AtomList(AlList(i), :), 'Al21')
        AlC = C(3,:);
    elseif contains(AtomList(AlList(i), :), 'Al22')
        AlC = C(4,:);
    elseif contains(AtomList(AlList(i), :), 'Alb')
        AlC = C(5,:);
    end
    plot([StepNum(1)/2000 StepNum(end)/2000], [mean(SumCharge(AlList(i),:)) mean(SumCharge(AlList(i),:))], '--', 'color', AlC);
end
if length(AlList) == 2
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), 'interpreter', 'tex')
elseif length(AlList) == 3
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), 'interpreter', 'tex')
elseif length(AlList) == 5
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), AtomList(AlList(4),:), AtomList(AlList(5),:), 'interpreter', 'tex')
end
hold off

% mean charge per Al type
figure
box on
hold on
xlabel('Time (ps)');
ylabel('Ave. Charge per Atom (e)');
% title(['Total Bader charge for Pt Atoms in ' system], 'interpreter', 'none')
AlC = [];
for i = 1:length(AlList)
    if contains(AtomList(AlList(i), :), 'Al11')
        AlC = C(1,:);
    elseif contains(AtomList(AlList(i), :), 'Al12')
        AlC = C(2,:);
    elseif contains(AtomList(AlList(i), :), 'Al21')
        AlC = C(3,:);
    elseif contains(AtomList(AlList(i), :), 'Al22')
        AlC = C(4,:);
    elseif contains(AtomList(AlList(i), :), 'Alb')
        AlC = C(5,:);
    end
    errorbar(StepNum/2000, MeanCharge(AlList(i),:), StdCharge(AlList(i),:), '-o', 'color', AlC, 'markeredgecolor', 'k', 'markerfacecolor', AlC);
end
if length(AlList) == 2
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), 'interpreter', 'tex')
elseif length(AlList) == 3
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), 'interpreter', 'tex')
elseif length(AlList) == 5
    legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), AtomList(AlList(4),:), AtomList(AlList(5),:), 'interpreter', 'tex')
end
hold off

figure
box on
hold on
xlabel('Time (ps)');
ylabel('Total Charge (e)');
title(['Total Bader charge for all electrolyte species in ' system], 'interpreter', 'none')
C = [1 0 0; 1 0 0; 1 0 0; 0 0.5 0; 0 0 1; 218/255 165/255 32/255;219/255 166/255 34/255];
IonList = 1:length(AtomList);
IonList(AlList) = [];
IonSumCharge = sum(SumCharge(IonList,:));
plot(StepNum/2000, IonSumCharge, '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'm');
legend('Total Electrolyte', 'interpreter', 'tex')
hold off

for i = 1:length(AtomList)
    if ~ismember(i, AlList)
        figure
        box on
        hold on
        h1 = plot(StepNum/2000, SumCharge(i,:), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', C(i,:));
        h2 = plot([StepNum(1)/2000 StepNum(end)/2000], [mean(SumCharge(i,2:end)) mean(SumCharge(i,2:end))], '--', 'color', [0.7 0.7 0.7]);
        AveStr = ['Traj. Ave. = ' num2str(mean(SumCharge(i,2:end)), '%.2f') '\pm'  num2str(std(SumCharge(i,2:end)), '%.2f')];
        legend(AtomList(i,:), AveStr, 'interpreter', 'tex')
        xlabel('Time (ps)');
        ylabel('Total Charge (e)');
        title(['Total Bader charge for ' AtomList(i,:) ' in ' system], 'interpreter', 'none')
        uistack(h1, 'top');
        hold off
    end
end

%%%% ----- Water ----- %%%%

% Bader charge and charge density on 1st layer of DL
figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
xlabel('Time (ps)');
plot(StepNum(DL1st_sum~=0)/2000, (DL1st_sum(DL1st_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL1st_sum) mean(DL1st_sum)], '--', 'color', 'r');
ylabel('Total Charge (e)');
legend('1st Water Layer')
hold off

% Bader charge and charge density on 2nd layer of DL
figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(DL2nd_sum~=0)/2000, (DL2nd_sum(DL2nd_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL2nd_sum) mean(DL2nd_sum)], '--', 'color', 'r');
ylabel('Total Charge (e)');
legend('2nd Water Layer')
hold off

% Bader charge and charge density on non-DL "Bulk"
figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(nonDL_sum~=0)/2000, (nonDL_sum(nonDL_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'g');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(nonDL_sum) mean(nonDL_sum)], '--', 'color', 'g');
ylabel('Total Charge (e)');
legend('Bulk Water')
hold off

% Bader charge and charge density on 1st and 2nd layer of DL
figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(DL1st_sum~=0)/2000, (DL1st_sum(DL1st_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL1st_sum) mean(DL1st_sum)], '--', 'color', 'r');
% ylabel('Total Charge (e)');
% legend('1st Water Layer')
plot(StepNum(DL2nd_sum~=0)/2000, (DL2nd_sum(DL2nd_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'b');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL2nd_sum) mean(DL2nd_sum)], '--', 'color', 'b');
ylabel('Total Charge (e)');
% legend('2nd Water Layer')
legend('1st Water Layer', '<1st Water Layer>', '2nd Water Layer', '<2nd Water Layer>', 'interpreter', 'tex')
hold off

% Bader charge and charge density on 1st and 2nd layer of DL and nonDL
figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(DL1st_sum~=0)/2000, (DL1st_sum(DL1st_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL1st_sum) mean(DL1st_sum)], '--', 'color', 'r');

plot(StepNum(DL2nd_sum~=0)/2000, (DL2nd_sum(DL2nd_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'b');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL2nd_sum) mean(DL2nd_sum)], '--', 'color', 'b');

plot(StepNum(nonDL_sum~=0)/2000, (nonDL_sum(nonDL_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'g');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(nonDL_sum) mean(nonDL_sum)], '--', 'color', 'g');
ylabel('Total Charge (e)');

legend('1st Water Layer', '<1st Water Layer>', '2nd Water Layer', '<2nd Water Layer>', 'Bulk Water', '<Bulk Water>', 'interpreter', 'tex')
hold off


%%%% ----- Combo ----- %%%%

% Total charge per all Al types + 1 DL

figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(DL1st_sum~=0)/2000, (DL1st_sum(DL1st_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL1st_sum) mean(DL1st_sum)], '--', 'color', 'r');

plot(StepNum/2000, totalAlsCharge, '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'c');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(totalAlsCharge) mean(totalAlsCharge)], '--', 'color', 'c');
ylabel('Total Charge (e)');

legend('1st Water Layer', '<1st Water Layer>', 'Al Layers', '<Al Layers>', 'interpreter', 'tex')
hold off




% Total charge per Al type + 1&2 DL

figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(DL1st_sum~=0)/2000, (DL1st_sum(DL1st_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL1st_sum) mean(DL1st_sum)], '--', 'color', 'r');

plot(StepNum(DL2nd_sum~=0)/2000, (DL2nd_sum(DL2nd_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'b');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL2nd_sum) mean(DL2nd_sum)], '--', 'color', 'b');

plot(StepNum/2000, totalAlsCharge, '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'c');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(totalAlsCharge) mean(totalAlsCharge)], '--', 'color', 'c');
ylabel('Total Charge (e)');

legend('1st Water Layer', '<1st Water Layer>', '2nd Water Layer', '<2nd Water Layer>', 'Al Layers', '<Al Layers>', 'interpreter', 'tex')
hold off



%%%% ----- Everything ----- %%%%

figure
box on
hold on
set(gca, 'colororder', [0 0 0]);
....
xlabel('Time (ps)');
plot(StepNum(DL1st_sum~=0)/2000, (DL1st_sum(DL1st_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL1st_sum) mean(DL1st_sum)], '--', 'color', 'r');

plot(StepNum(DL2nd_sum~=0)/2000, (DL2nd_sum(DL2nd_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'b');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(DL2nd_sum) mean(DL2nd_sum)], '--', 'color', 'b');

plot(StepNum(nonDL_sum~=0)/2000, (nonDL_sum(nonDL_sum~=0)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'g');
plot([StepNum(1)/2000 StepNum(end)/2000], [mean(nonDL_sum) mean(nonDL_sum)], '--', 'color', 'g');
ylabel('Total Charge (e)');

% legend('1st Water Layer', '<1st Water Layer>', '2nd Water Layer', '<2nd Water Layer>', 'Bulk Water', '<Bulk Water>', 'interpreter', 'tex')

xlabel('Time (ps)');
ylabel('Total Charge (e)');
% title(['Total Bader charge for Al Atoms in ' system], 'interpreter', 'none')
C = [218/255 165/255 32/255; 0.25 0.25 0.25; 0 0.5 0; 0 0 1; 1 0 0];
AlC = [];
for i = 1:length(AlList)
    if contains(AtomList(AlList(i), :), 'Al11')
        AlC = C(1,:);
    elseif contains(AtomList(AlList(i), :), 'Al12')
        AlC = C(2,:);
    elseif contains(AtomList(AlList(i), :), 'Al21')
        AlC = C(3,:);
    elseif contains(AtomList(AlList(i), :), 'Al22')
        AlC = C(4,:);
    elseif contains(AtomList(AlList(i), :), 'Alb')
        AlC = C(5,:);
    end
    plot(StepNum/2000, SumCharge(AlList(i),:), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', AlC);
end

for i = 1:length(AlList)
    if contains(AtomList(AlList(i), :), 'Al11')
        AlC = C(1,:);
    elseif contains(AtomList(AlList(i), :), 'Al12')
        AlC = C(2,:);
    elseif contains(AtomList(AlList(i), :), 'Al21')
        AlC = C(3,:);
    elseif contains(AtomList(AlList(i), :), 'Al22')
        AlC = C(4,:);
    elseif contains(AtomList(AlList(i), :), 'Alb')
        AlC = C(5,:);
    end
    plot([StepNum(1)/2000 StepNum(end)/2000], [mean(SumCharge(AlList(i),:)) mean(SumCharge(AlList(i),:))], '--', 'color', AlC);
end
% if length(AlList) == 2
%     legend(AtomList(AlList(1),:), AtomList(AlList(2),:), 'interpreter', 'tex')
% elseif length(AlList) == 3
%     legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), 'interpreter', 'tex')
% elseif length(AlList) == 5
%     legend(AtomList(AlList(1),:), AtomList(AlList(2),:), AtomList(AlList(3),:), AtomList(AlList(4),:), AtomList(AlList(5),:), 'interpreter', 'tex')
% end

hold off

%%%%%%%%%%%%%% Heat map(s) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Precondition(s)
% <XY coords/all snapshots> --> xyz(nStep, atom, coord) , Bader charge for atoms
% of interest--> <Qnet (atom, nSnap)>

%Output: Heat map for atoms of interest per all snapshots

Al11_xy_tbl= randi(Indx.Al11,n);


