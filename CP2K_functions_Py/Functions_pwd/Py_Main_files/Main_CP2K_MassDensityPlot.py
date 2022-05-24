import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
## Set the location of the calculation output
Basefldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\'

system = 'CP_Pit_20F'

Trajectory = 'CP_Pit_20F_43000to73000_500step.xyz'
##

ABC = getABCvectors(Basefldr,system)
xyz,XYZ,Indx,Atoms,AtomList,nAtoms,startConfig,nConfigs,StepNum = ReadAndParsexyz(Basefldr,system,Trajectory,ABC)
XYZ = wrapXYZ(XYZ,ABC)
## Modify which O atoms go into mass density accoring to explicit naming in
# input xyz. AtomIndx is the Index of atoms by name from input.xyz. Modify
# ismember argument to choose which atoms are used.
__,__,AtomIndx = getAtomNamesFromInputXYZ(Basefldr,system)
# [~, OIndxxyz] = ismember([AtomIndx.O; AtomIndx.OtL; AtomIndx.OtU; AtomIndx.OtS] , Indx.O);
__,OIndxxyz = ismember(np.array([AtomIndx.O]),Indx.O)
xyz.O = xyz.O(:,OIndxxyz,:)
##

# [Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
Dens_O,Dens_H,Dens_F,TotDen,AveDen,z = getDensityProfile(xyz,ABC)
# [Dens_O, Dens_H, Dens_Na, Dens_Cl, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);

zmax = ABC(3)
bins = len(z) + 1
# Plot each snapshot
figure
hold('on')
set(gcf,'position',np.array([377,423,1123,420]))
plt.xlabel('z (Ang)')
plt.ylabel('Density (kgm^{-3})')
set(gca,'xlim',np.array([0,ABC(3)]),'ylim',np.array([0,2500]))
# include profile for every snapshot
for i in np.arange(1,TotDen.shape[2-1]+1).reshape(-1):
    plt.plot(z,smooth(TotDen(:,i)),'--','linewidth',0.25,'color','r')

plt.plot(z,smooth(np.sum(TotDen, 2-1) / (nConfigs - startConfig + 1)),'linewidth',1.5,'color','k')
plt.plot(np.array([z(1),z(end())]),np.array([mean(AveDen),mean(AveDen)]),':','color',np.array([0.6,0.6,0.6]))
plt.plot(np.array([z((bins / 2) - np.round(3 / (zmax / (bins - 1)) / 2)),z((bins / 2) - np.round(3 / (zmax / (bins - 1)) / 2))]),np.array([0,2500]),':','color',np.array([0.6,0.6,0.6]))
plt.plot(np.array([z((bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2)),z((bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2))]),np.array([0,2500]),':','color',np.array([0.6,0.6,0.6]))
hold('off')
figure
hold('on')
set(gcf,'position',np.array([377,423,1123,420]))
plt.xlabel('z (Ang)')
plt.ylabel('Density (kgm^{-3})')
set(gca,'xlim',np.array([0,ABC(3)]),'ylim',np.array([0,2500]))
plt.plot(z,smooth(np.sum(TotDen, 2-1) / (nConfigs - startConfig + 1),3),'linewidth',1.5,'color','k')
plt.plot(z,smooth(np.sum(Dens_H, 2-1) / (nConfigs - startConfig + 1),3),'linewidth',1.5,'color','b')
plt.plot(z,smooth(np.sum(Dens_O, 2-1) / (nConfigs - startConfig + 1),3),'linewidth',1.5,'color','r')
plt.plot(z,smooth(np.sum(Dens_F, 2-1) / (nConfigs - startConfig + 1),3),'linewidth',1.5,'color',np.array([34,177,76]) / 255)
# plot(z, smooth(sum(Dens_Na,2)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', [128 0 128]/255)
# plot(z, smooth(sum(Dens_Cl,2)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', [34 177 76]/255)
# smoothing off
# plot(z, sum(TotDen,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'k')
# plot(z, sum(Dens_H,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'b')
# plot(z, sum(Dens_O,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'r')
# plot(z, sum(Dens_F,2)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', [34 177 76]/255)
plt.plot(np.array([z(1),z(end())]),np.array([mean(AveDen),mean(AveDen)]),':','color',np.array([0.6,0.6,0.6]))
plt.plot(np.array([z((bins / 2) - np.round(3 / (zmax / (bins - 1)) / 2)),z((bins / 2) - np.round(3 / (zmax / (bins - 1)) / 2))]),np.array([0,2500]),':','color',np.array([0.6,0.6,0.6]))
plt.plot(np.array([z((bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2)),z((bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2))]),np.array([0,2500]),':','color',np.array([0.6,0.6,0.6]))
plt.legend('Water+Ions','H','O','F','Ave. Bulk Density','location','northeast')
# legend('Water', 'H', 'O', 'Ave. Bulk Density', 'location', 'northeast');
# legend('Water+Ions', 'H', 'O', 'Na', 'Cl', 'Ave. Bulk Density', 'location', 'northeast');
hold('off')
## uncomment to save a jpg of the mass density
# if exist([Basefldr 'MassDensityProfiles'],'dir')
#     warning('Directory already exists!');
# else
#     mkdir([Basefldr 'MassDensityProfiles']);
# end

# saveas(gcf, [Basefldr 'MassDensityProfiles\' system '.jpg']);