clear all;  clc;
close all;
PathSep = setOSpathSep;

%% Home
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/5lyr_systems/Al_AlO_OH/Al_OH/New-Al_OH/';
system = 'OH_0.185ML';
Trajectory = 'OH_0.185ML_0to71000_1000step.xyz';

% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx, Indxfns, Kinds, Elements, PP] = getAtomInfoFromInput(BaseFldr, system);

% % parse coordinates of atoms along trajectory and wrap into cell
[xyz, XYZ, ~, ~, ~, nAtoms, startConfig, nConfigs, StepNum_Traj] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0 0 0]);

OH_Coverage = zeros(nConfigs, 1);
H2O_Coverage = zeros(nConfigs, 1);
H3O_Coverage = zeros(nConfigs, 1);
O_Coverage = zeros(nConfigs, 1);

for i = 1:nConfigs

    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(i,:,:);

    % get the distances between pairs of atoms
    [~, DistAlO] = GetAtomCorrelation(XYZ_snap, Indx.Al1, Indx.O, ABC);
    [r,c] = find(DistAlO < 2.83); % Rashid to fix this "2" by looking at RDF minimum - r = row aka O atom number, c = column aka Al1 atom number

    [~, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, setdiff(Indx.O,Indx.O(r)), ABC);
%     [~, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.O, ABC);
    [rOH,cOH] = find(DistOH < 1.28);
    [GR, GC] = groupcounts(rOH);
    O_Coverage(i) = sum(GR == 0);
    OH_Coverage(i) = sum(GR == 1);
    H2O_Coverage(i) = sum(GR == 2);

    [gH3O, cH3O]= groupcounts(cOH);
    H3O_Coverage(i) = sum(GR == 3)-sum(gH3O==2); %Potential bug, we are picking up any shared H atom shared between two O atoms that might not be H3O and it might be a H shuttling between OH and H2O!
%One idea is to loop over the shared oxygens and then identify any OH
%molecules and if there are any then they would be excluded from the
%subtraction somehow, i.e. we want to only remove shared O between two
%H2O's and nothing else, those are the suspected double countings.

end

Total_Coverage = OH_Coverage+H2O_Coverage+H3O_Coverage;

%H2O,H3O,OH and total 
figure
hold on
plot((1:nConfigs)*0.5, OH_Coverage, '-ok', 'markerfacecolor', 'r')
plot((1:nConfigs)*0.5, H2O_Coverage, '-ok', 'markerfacecolor', 'b')
plot((1:nConfigs)*0.5, H3O_Coverage, '-ok', 'markerfacecolor', 'auto')
plot((1:nConfigs)*0.5, Total_Coverage, '-ok', 'markerfacecolor', [0 0.8 0])
legend('OH', 'Water','H3O', 'Total')
xlabel('Time (ps)')
ylabel('Number of Molecules')

%H3O and OH only
figure
hold on
plot((1:nConfigs)*0.5, OH_Coverage, '-ok', 'markerfacecolor', 'r')
plot((1:nConfigs)*0.5, H3O_Coverage, '-ok', 'markerfacecolor', 'auto')
legend('OH','H3O')
xlabel('Time (ps)')
ylabel('Number of Molecules')



%Stats
disp(['Ave. coverage of OH = ' num2str(mean(OH_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(OH_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. coverage of H2O = ' num2str(mean(H2O_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(H2O_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. coverage of H3O = ' num2str(mean(H3O_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(H3O_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. Total coverage = ' num2str(mean(Total_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(Total_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. O coverage = ' num2str(mean(O_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(O_Coverage(30:end)/108), '%.2f') ' ML'])