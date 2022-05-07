clear all; 
close all; clc;

%% Set the location of the calculation output
Basefldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus/Al_AlO/'; % Base directory containing calculation directory ("\" included at end)
system = 'Al_water'; % Name of calculation directory (no "\")
Trajectory = 'Al_water_14300to18100_100step.xyz';
%%

ABC = getABCvectors(Basefldr, system);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(Basefldr, system, Trajectory, ABC);
XYZ = wrapXYZ(XYZ, ABC);

%% Modify which O atoms go into mass density accoring to explicit naming in
% input xyz. AtomIndx is the Index of atoms by name from input.xyz. Modify
% ismember argument to choose which atoms are used.
[~, ~, AtomIndx] = getAtomNamesFromInputXYZ(Basefldr, system);
% [~, OIndxxyz] = ismember([AtomIndx.O; AtomIndx.OtL; AtomIndx.OtU; AtomIndx.OtS] , Indx.O);
[~, OIndxxyz] = ismember([AtomIndx.O] , Indx.O);
xyz.O = xyz.O(:,OIndxxyz,:);
%%

% [Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
[Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
% [Dens_O, Dens_H, Dens_Na, Dens_Cl, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);


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
    plot(z, smooth(TotDen(:,i)), '--', 'linewidth', 0.25, 'color', 'r')
end
plot(z, smooth(sum(TotDen,2)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', 'k')
plot([z(1) z(end)], [mean(AveDen) mean(AveDen)], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)-round(3/(zmax/(bins-1))/2)) z((bins/2)-round(3/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)+round(3/(zmax/(bins-1))/2)) z((bins/2)+round(3/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6])
hold off

figure
hold on
set(gcf, 'position', [377         423        1123         420]);
xlabel('z (Ang)');
ylabel('Density (kgm^{-3})');
set(gca, 'xlim', [0 ABC(3)], 'ylim', [0 2500]);
plot(z, smooth(sum(TotDen,2)/(nConfigs-startConfig+1),3), 'linewidth', 1.5, 'color', 'k')
plot(z, smooth(sum(Dens_H,2)/(nConfigs-startConfig+1),3), 'linewidth', 1.5, 'color', 'b')
plot(z, smooth(sum(Dens_O,2)/(nConfigs-startConfig+1),3), 'linewidth', 1.5, 'color', 'r')
%plot(z, smooth(sum(Dens_F,2)/(nConfigs-startConfig+1),3), 'linewidth', 1.5, 'color', [34 177 76]/255)
% plot(z, smooth(sum(Dens_Na,2)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', [128 0 128]/255)
% plot(z, smooth(sum(Dens_Cl,2)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', [34 177 76]/255)
% smoothing off
% plot(z, sum(TotDen,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'k')
% plot(z, sum(Dens_H,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'b')
% plot(z, sum(Dens_O,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'r')
% plot(z, sum(Dens_F,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', [34 177 76]/255)
plot([z(1) z(end)], [mean(AveDen) mean(AveDen)], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)-round(3/(zmax/(bins-1))/2)) z((bins/2)-round(3/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6]) % there might be a bug here since zmax isn't max(z)? 
plot([z((bins/2)+round(3/(zmax/(bins-1))/2)) z((bins/2)+round(3/(zmax/(bins-1))/2))], [0 2500], ':', 'color', [0.6 0.6 0.6])
legend('Water+Ions', 'H', 'O', 'Ave. Bulk Density', 'location', 'northeast');
% legend('Water', 'H', 'O', 'Ave. Bulk Density', 'location', 'northeast');
% legend('Water+Ions', 'H', 'O', 'Na', 'Cl', 'Ave. Bulk Density', 'location', 'northeast');
hold off

%% uncomment to save a jpg of the mass density
% if exist([Basefldr 'MassDensityProfiles'],'dir')
%     warning('Directory already exists!');
% else
%     mkdir([Basefldr 'MassDensityProfiles']);
% end
% 
% saveas(gcf, [Basefldr 'MassDensityProfiles\' system '.jpg']);

