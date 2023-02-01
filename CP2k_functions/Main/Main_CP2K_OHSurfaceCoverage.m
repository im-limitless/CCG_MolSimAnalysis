clear all;  clc;
close all;
PathSep = setOSpathSep;

%% Home
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/5lyr_systems/Al_AlO_OH/Al_OH/New-Al_OH/';
system = 'OH_0.33ML';
Trajectory = 'OH_0.33ML_0to68000_1000step.xyz';

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
    [r,c] = find(DistAlO < 2); % Rashid to fix this "2" by looking at RDF minimum - r = row aka O atom number, c = column aka Al1 atom number

    [~, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.O(r), ABC);
    [rOH,cOH] = find(DistOH < 1.3);
    [GR, GC] = groupcounts(rOH);
    O_Coverage(i) = sum(GR(GR == 0));
    OH_Coverage(i) = sum(GR(GR == 1));
    H2O_Coverage(i) = sum(GR(GR == 2));
    H3O_Coverage(i) = sum(GR(GR == 3));


end

Total_Coverage = OH_Coverage+H2O_Coverage;

figure
hold on
plot((1:nConfigs)*0.5, OH_Coverage, '-ok', 'markerfacecolor', 'r')
plot((1:nConfigs)*0.5, H2O_Coverage, '-ok', 'markerfacecolor', 'b')
plot((1:nConfigs)*0.5, Total_Coverage, '-ok', 'markerfacecolor', [0 0.8 0])
legend('OH', 'Water', 'Total')
xlabel('Time (ps)')
ylabel('Number of Molecules')

disp(['Ave. coverage of OH = ' num2str(mean(OH_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(OH_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. coverage of H2O = ' num2str(mean(H2O_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(H2O_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. coverage of H3O = ' num2str(mean(H3O_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(H3O_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. Total coverage = ' num2str(mean(Total_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(Total_Coverage(30:end)/108), '%.2f') ' ML'])
disp(['Ave. O coverage = ' num2str(mean(O_Coverage(30:end)/108), '%.2f') ' +/- ' num2str(std(O_Coverage(30:end)/108), '%.2f') ' ML'])