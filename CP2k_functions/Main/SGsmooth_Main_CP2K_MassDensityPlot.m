clear all; 
close all; clc;

PathSep =  setOSpathSep;

%% Set the location of the calculation output
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge_2/Phase_diagram_Sys/';
system = 'Al_water';
Trajectory = 'Al_water-pos-1.xyz';

% BaseFldr = '/Users/mtdarby/Dropbox/Mac/Documents/MattProjects/TempAnalysis/';
% system = 'CP_Pit_20F';
% Trajectory = 'CP_Pit_20F_43000to73000_500step.xyz';

% BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\';
% system = 'CP_Like_1012_Fluoride';
% Trajectory = 'Sample34000_52000.xyz';
%%  

% % get the names of atoms from original xyz input file
[~, ~, AtomIndx, ~, ~, ~, ~] = getAtomInfoFromInput(BaseFldr, system);

ABC = getABCvectors(BaseFldr, system);
% [xyz, XYZ, Indx, ~, ~, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC, [pi/2; 1; 1]);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);
% [xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [pi/2; 1; 1]);


% % % % % %% Modify which O atoms go into mass density accoring to explicit naming in
% % % % % % input xyz. AtomIndx is the Index of atoms by name from input.xyz. Modify
% % % % % % % get the names of atoms from original xyz input file
% % % % % % [~, OIndxxyz] = ismember([AtomIndx.O; AtomIndx.OtL; AtomIndx.OtU; AtomIndx.OtS] , Indx.O);
% % % % % [~, OIndxxyz] = ismember([AtomIndx.O] , Indx.O);
% % % % % xyz.O = xyz.O(:,OIndxxyz,:);
%%

[Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
% [Dens_O, Dens_H, TotDen, AveDen, z] = getCylindricalDensity(xyz, ABC);
% [Dens_O, Dens_H, Dens_F, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
% [Dens_O, Dens_H, Dens_Na, Dens_Cl, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);

% [FirstLayerIndx, SecondLayerIndx, ThirdLayerIndx, FourthLayerIndx, MinimaZ] = getWaterLayerIndices(Indx, XYZ, Dens_O, z);

%%%%%%%%%%%%%%%%% Macroscopic ave (uncomment) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prompt = "Macroscopic average using (1)Iterating the average with fixed <vox> [M1], or, (2)Varying <voxels> [M2]? ('1'/'2'/'both'): ";
% oxide = input(prompt);
% if (oxide == '1')
%     M1_getBulkMacroscopicAve(TotDen, z, ABC, BaseFldr, system);
% elseif(oxide=='2')
%     M2_getBulkMacroscopicAve(TotDen, z, ABC, BaseFldr, system);
% elseif(oxide=='both')
%     M1_getBulkMacroscopicAve(TotDen, z, ABC, BaseFldr, system);
%     M2_getBulkMacroscopicAve(TotDen, z, ABC, BaseFldr, system);
% else
%     error ('Incorrect value for the usr prompt');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmax = ABC(3);
bins = length(z) + 1;
% Plot each snapshot
figure
hold on
set(gcf, 'position', [377         423        1123         420]);
xlabel('z (Ang)');
ylabel('Density (kgm^{-3})');
set(gca, 'xlim', [0 ABC(3)], 'ylim', [0 2500]);
% include profile for every snapshot
for i = 1:size(TotDen,2)
    plot(z, TotDen(:,i), '--', 'linewidth', 0.25)
%     plot(z, TotDen(:,i), '--', 'linewidth', 0.25, 'color', 'r')
end
plot(z, sum(TotDen,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'k')
plot([z(1) z(end)], [mean(AveDen) mean(AveDen)], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)-round(3/(zmax/(bins-1))/2)) z((bins/2)-round(3/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)+round(3/(zmax/(bins-1))/2)) z((bins/2)+round(3/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6])
hold off

% figure
% hold on
% set(gcf, 'position', [377         423        1123         420]);
% xlabel('z (Ang)');
% ylabel('Density (kgm^{-3})');
% set(gca, 'xlim', [0 ABC(3)], 'ylim', [0 2500]);

%% no smoothing
% plot(z, sum(TotDen,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'k')
% plot(z, sum(Dens_H,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'b')
% plot(z, sum(Dens_O,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'r')

% %% Savitzky-Golay filtering (smoothing)
% plot(z, sgolayfilt(sum(TotDen,2)/(nConfigs-startConfig+1),2,3), 'linewidth', 1.5, 'color', 'k')
% plot(z, sgolayfilt(sum(Dens_H,2)/(nConfigs-startConfig+1),2,3), 'linewidth', 1.5, 'color', 'b')
% plot(z, sgolayfilt(sum(Dens_O,2)/(nConfigs-startConfig+1),2,3), 'linewidth', 1.5, 'color', 'r')
% 
% % plot(z, sum(Dens_F,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', [34 177 76]/255)
% plot([z(1) z(end)], [mean(AveDen) mean(AveDen)], ':', 'color', [0.6 0.6 0.6])
% plot([z((bins/2)-round(5/(zmax/(bins-1))/2)) z((bins/2)-round(5/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6]) % there might be a bug here since zmax isn't max(z)? 
% plot([z((bins/2)+round(5/(zmax/(bins-1))/2)) z((bins/2)+round(5/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6])
% % legend('Water+Ions', 'H', 'O', 'F', 'Ave. Bulk Density', 'location', 'northeast');
% legend('Water', 'H', 'O', 'Ave. Bulk Density', 'location', 'northeast');
% % legend('Water+Ions', 'H', 'O', 'Na', 'Cl', 'Ave. Bulk Density', 'location', 'northeast');
% hold off

% Sort the data by Z
[Z_sorted, sort_idx] = sort(z);
TotDen_sorted = TotDen(sort_idx);

% Define segments with their Z ranges and filter parameters
segments = [
    struct('Z_start', 0,  'Z_end', 10, 'order', 2, 'frameLength', 3);
    struct('Z_start', 10, 'Z_end', 30, 'order', 3, 'frameLength', 19);
    struct('Z_start', 30, 'Z_end', 50, 'order', 2, 'frameLength', 3);
];

% Initialize arrays for combining results
filtered_density = zeros(size(TotDen_sorted));
count = zeros(size(TotDen_sorted));

% Process each segment
for i = 1:length(segments)
    seg = segments(i);
    
    % Indices within the current Z range
    in_range = Z_sorted >= seg.Z_start & Z_sorted <= seg.Z_end;
    indices = find(in_range);
    if isempty(indices), continue; end
    
    % Segment start and end indices
    idx1 = indices(1);
    idx2 = indices(end);
    
    % Buffer based on frame length
    buffer = (seg.frameLength - 1)/2;
    idx1_ext = max(1, idx1 - buffer);
    idx2_ext = min(length(Z_sorted), idx2 + buffer);
    
    % Extract extended segment
    Z_ext = Z_sorted(idx1_ext:idx2_ext);
    TotDen_ext = TotDen_sorted(idx1_ext:idx2_ext);
    
    % Apply Savitzky-Golay filter
    filtered_ext = sgolayfilt(TotDen_ext, seg.order, seg.frameLength);
    
    % Trim to original segment
    start_offset = idx1 - idx1_ext;
    filtered_segment = filtered_ext(1 + start_offset : end - (idx2_ext - idx2));
    
    % Accumulate results
    filtered_density(idx1:idx2) = filtered_density(idx1:idx2) + filtered_segment;
    count(idx1:idx2) = count(idx1:idx2) + 1;
end

% Average overlapping regions and handle uncovered areas
filtered_density = filtered_density ./ count;
filtered_density(count == 0) = TotDen_sorted(count == 0); % Use original where no filter applied

% Restore original order
[~, unsort_idx] = sort(sort_idx);
filtered_density_unsorted = filtered_density(unsort_idx);

% Plot results
figure;
plot(z, TotDen, 'k', z, filtered_density_unsorted, 'r-');
legend('Original Data', 'Filtered Data');
xlabel('Z');
ylabel('Density');
title('Savitzky-Golay Filter Applied to Different Z Ranges');




%% uncomment to save a jpg of the mass density
% if exist([BaseFldr 'MassDensityProfiles'],'dir')
%     warning('Directory already exists!');
% else
%     mkdir([BaseFldr 'MassDensityProfiles']);
% end

% saveas(gcf, [BaseFldr 'MassDensityProfiles' PathSep system '.jpg']);

 [FirstLayerIndx, SecondLayerIndx, MinimaZ] = getWaterLayerIndicesPerSnapRestrictedRev(Indx, XYZ, Dens_O, z, [100 -100]);
 [FirstLayerIndx_low, SecondLayerIndx_low, MinimaZ] = getWaterLayerIndicesPerSnapRestricted(Indx, XYZ, Dens_O, z, [100 -100]);

