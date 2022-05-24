import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
Basefldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\'

system = 'CP_Like_1010_Fluorine'

Trajectory = 'Sample36000_58000.xyz'
ABC = getABCvectors(Basefldr,system)
xyz,XYZ,Indx,Atoms,AtomList,nAtoms,startConfig,nConfigs,StepNum = ReadAndParsexyz(Basefldr,system,Trajectory,ABC)
## Modify which O atoms go into mass density accoring to explicit naming in
# input xyz. AtomIndx is the Index of atoms by name from input.xyz. Modify
# ismember argument to choose which atoms are used.
__,__,AtomIndx = getAtomNamesFromInputXYZ(Basefldr,system)
# [~, OIndxxyz] = ismember([AtomIndx.O; AtomIndx.OtL; AtomIndx.OtU; AtomIndx.OtS] , Indx.O);
__,OIndxxyz = ismember(AtomIndx.O,Indx.O)
xyz.O = xyz.O(:,OIndxxyz,:)
##

# [Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
Dens_O,Dens_H,Dens_F,TotDen,AveDen,z = getDensityProfile(xyz,ABC)
# [Dens_O, Dens_H, Dens_Na, Dens_Cl, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);

FirstLayerIndx,SecondLayerIndx = getWaterLayerIndices(AtomIndx,XYZ,Dens_O,z)
for snap in np.arange(startConfig,nConfigs+1).reshape(-1):
    print(np.array(['Processing snapshot ',num2str(StepNum(snap)),' - ',num2str(100 * (snap / nConfigs)),' % complete']))
    XYZ_snap = np.zeros((XYZ.shape[2-1],XYZ.shape[3-1]))
    XYZ_snap[:,:] = XYZ(snap,:,:)
    VecOH,DistOH = GetAtomCorrelation(XYZ_snap,AtomIndx.O,AtomIndx.H,ABC)
    for j in np.arange(1,len(AtomIndx.O)+1).reshape(-1):
        BondOH = find(DistOH(:,j) < 1.33)
        # find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
        if len(BondOH) == 2:
            # calculate net vector bisecting H-O-H
            VecH2O = VecOH(BondOH(1),:) + VecOH(BondOH(2),:)
            # # #             # if we want to assume only one direction for both halves if the cell
# # #                     zVec = [0 0 1];
# # #                        Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
            # split the cell into half and determine if O is in upper or lower half
# then take the dot product of resultant water vector with both
# surface normals (assuming the normal is c)
            if XYZ_snap(AtomIndx.O(j),3) < ABC(3) / 2:
                zVec = np.array([0,0,1])
                Theta = acosd(np.dot(VecH2O,zVec) / (norm(VecH2O) * norm(zVec)))
                Phi1 = acosd(np.dot(VecOH(BondOH(1),:),zVec) / (norm(VecOH(BondOH(1),:)) * norm(zVec)))
                Phi2 = acosd(np.dot(VecOH(BondOH(2),:),zVec) / (norm(VecOH(BondOH(2),:)) * norm(zVec)))
            else:
                if XYZ_snap(AtomIndx.O(j),3) >= ABC(3) / 2:
                    zVec = np.array([0,0,- 1])
                    Theta = acosd(np.dot(VecH2O,zVec) / (norm(VecH2O) * norm(zVec)))
                    Phi1 = acosd(np.dot(VecOH(BondOH(1),:),zVec) / (norm(VecOH(BondOH(1),:)) * norm(zVec)))
                    Phi2 = acosd(np.dot(VecOH(BondOH(2),:),zVec) / (norm(VecOH(BondOH(2),:)) * norm(zVec)))
        ThetaSnaps[snap,j] = Theta
        PhiSnaps[snap,j] = np.array([Phi1,Phi2])
    DL_Indx = np.array([[FirstLayerIndx[snap]],[SecondLayerIndx[snap]]])
    #     DL_Indx = SecondLayerIndx{snap};
    DL_Indx_H = []
    #     DL_Indx = FirstLayerIndx{snap};
    for j in np.arange(1,len(DL_Indx)+1).reshape(-1):
        BondOH = find(DistOH(:,find(AtomIndx.O == DL_Indx(j))) < 1.33)
        # find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
        if len(BondOH) == 2:
            # calculate net vector bisecting H-O-H
            VecH2O = VecOH(BondOH(1),:) + VecOH(BondOH(2),:)
            # # #             # if we want to assume only one direction for both halves if the cell
# # #                     zVec = [0 0 1];
# # #                        Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
            # split the cell into half and determine if O is in upper or lower half
# then take the dot product of resultant water vector with both
# surface normals (assuming the normal is c)
            if XYZ_snap(DL_Indx(j),3) < ABC(3) / 2:
                zVec = np.array([0,0,1])
                ThetaDL = acosd(np.dot(VecH2O,zVec) / (norm(VecH2O) * norm(zVec)))
                Phi1DL = acosd(np.dot(VecOH(BondOH(1),:),zVec) / (norm(VecOH(BondOH(1),:)) * norm(zVec)))
                Phi2DL = acosd(np.dot(VecOH(BondOH(2),:),zVec) / (norm(VecOH(BondOH(2),:)) * norm(zVec)))
            else:
                if XYZ_snap(DL_Indx(j),3) >= ABC(3) / 2:
                    zVec = np.array([0,0,- 1])
                    ThetaDL = acosd(np.dot(VecH2O,zVec) / (norm(VecH2O) * norm(zVec)))
                    Phi1DL = acosd(np.dot(VecOH(BondOH(1),:),zVec) / (norm(VecOH(BondOH(1),:)) * norm(zVec)))
                    Phi2DL = acosd(np.dot(VecOH(BondOH(2),:),zVec) / (norm(VecOH(BondOH(2),:)) * norm(zVec)))
        ThetaSnapsDL[snap,j] = ThetaDL
        PhiSnapsDL[snap,j] = np.array([Phi1DL,Phi2DL])
        DL_Indx_H = np.array([[DL_Indx_H],[AtomIndx.H(BondOH)]])
    # write a sample snapshot to .xyz and convert to MS file (only for last
# snapshot - simple bug that number of atoms in DL changes so arrays
# are non-consistent for concatenation when trying to do all snaps together - could modify using cell array)
    if snap == nConfigs:
        writeSnaptoxyz(Basefldr,system,snap,XYZ_snap,Atoms,np.array([[DL_Indx],[DL_Indx_H],[AtomIndx.Pts],[AtomIndx.Ptb],[AtomIndx.Ptss]]),'DL_Density')
        CP2kOptimPathParse(Basefldr,system,np.array(['DL_Density_',num2str(snap),'.xyz']))
    DL_Indx = FirstLayerIndx[snap]
    for j in np.arange(1,len(DL_Indx)+1).reshape(-1):
        BondOH = find(DistOH(:,find(AtomIndx.O == DL_Indx(j))) < 1.33)
        # find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
        if len(BondOH) == 2:
            # calculate net vector bisecting H-O-H
            VecH2O = VecOH(BondOH(1),:) + VecOH(BondOH(2),:)
            # # #             # if we want to assume only one direction for both halves if the cell
# # #                     zVec = [0 0 1];
# # #                        Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
            # split the cell into half and determine if O is in upper or lower half
# then take the dot product of resultant water vector with both
# surface normals (assuming the normal is c)
            if XYZ_snap(DL_Indx(j),3) < ABC(3) / 2:
                zVec = np.array([0,0,1])
                ThetaDL = acosd(np.dot(VecH2O,zVec) / (norm(VecH2O) * norm(zVec)))
                Phi1DL = acosd(np.dot(VecOH(BondOH(1),:),zVec) / (norm(VecOH(BondOH(1),:)) * norm(zVec)))
                Phi2DL = acosd(np.dot(VecOH(BondOH(2),:),zVec) / (norm(VecOH(BondOH(2),:)) * norm(zVec)))
            else:
                if XYZ_snap(DL_Indx(j),3) >= ABC(3) / 2:
                    zVec = np.array([0,0,- 1])
                    ThetaDL = acosd(np.dot(VecH2O,zVec) / (norm(VecH2O) * norm(zVec)))
                    Phi1DL = acosd(np.dot(VecOH(BondOH(1),:),zVec) / (norm(VecOH(BondOH(1),:)) * norm(zVec)))
                    Phi2DL = acosd(np.dot(VecOH(BondOH(2),:),zVec) / (norm(VecOH(BondOH(2),:)) * norm(zVec)))
        ThetaSnaps1stWL[snap,j] = ThetaDL
        PhiSnaps1stWL[snap,j] = np.array([Phi1DL,Phi2DL])
        #                 DL_Indx_H = [DL_Indx_H; AtomIndx.H(BondOH)];

# # # # # # # # Theta is the angle of the water H-O-H bisector with the surface normal/c-vector
Theta = vertcat(ThetaSnaps[:,:])
Theta_DL = vertcat(ThetaSnapsDL[:,:])
Theta_1stWL = vertcat(ThetaSnaps1stWL[:,:])
for j in np.arange(1,ThetaSnapsDL.shape[1-1]+1).reshape(-1):
    nH2O_DL[j] = ThetaSnapsDL.shape[2-1] - nnz(cellfun(isempty,ThetaSnapsDL(j,:)))

for j in np.arange(1,ThetaSnaps1stWL.shape[1-1]+1).reshape(-1):
    nH2O_1stWL[j] = ThetaSnaps1stWL.shape[2-1] - nnz(cellfun(isempty,ThetaSnaps1stWL(j,:)))

print(np.array(['Ave. number of water molecules in 1st + 2nd WL = ',num2str(mean(nH2O_DL) / 2),' \pm ',num2str(std(nH2O_DL) / 2)]))
print(np.array(['Ave. number of water molecules in 1st WL = ',num2str(mean(nH2O_1stWL) / 2),' \pm ',num2str(std(nH2O_1stWL) / 2)]))
figure
hold('on')
plt.title(np.array(['H-O-H Bisector Angle Distribution for ',system]),'interpreter','none')
f,x = hist(Theta,60)
# histogram(Theta,60);
plt.plot(x,smooth(smooth(f)) / len(AtomIndx.O) / nConfigs,'k','linewidth',2)
f,x = hist(Theta_DL,60)
# histogram(Theta_DL,60);
plt.plot(x,smooth(smooth(f)) / len(AtomIndx.O) / nConfigs,'r','linewidth',2)
f,x = hist(Theta_1stWL,60)
# histogram(Theta_DL,60);
plt.plot(x,smooth(smooth(f)) / len(AtomIndx.O) / nConfigs,'b','linewidth',2)
plt.xlabel('\Theta (deg)')
plt.ylabel('Abundance')
plt.legend('Global','1st+2nd Layer','1st Layer')
hold('off')
# # # # # # # # Phi is the angle of the O-H bond vector made with the surface normal
Phi = np.reshape(vertcat(PhiSnaps[:]), tuple(np.array([2 * len(vertcat(PhiSnaps[:])),1])), order="F")
Phi_DL = np.reshape(vertcat(PhiSnapsDL[:]), tuple(np.array([2 * len(vertcat(PhiSnapsDL[:])),1])), order="F")
Phi_1stWL = np.reshape(vertcat(PhiSnaps1stWL[:]), tuple(np.array([2 * len(vertcat(PhiSnaps1stWL[:])),1])), order="F")
figure
hold('on')
plt.title(np.array(['O-H Angle Distribution for ',system]),'interpreter','none')
f,x = hist(Phi,60)
# histogram(Phi,60);
plt.plot(x,smooth(smooth(f)) / len(AtomIndx.O) / 2 / nConfigs,'k','linewidth',2)
f,x = hist(Phi_DL,60)
# histogram(Phi_DL,60);
plt.plot(x,smooth(smooth(f)) / len(AtomIndx.O) / 2 / nConfigs,'r','linewidth',2)
f,x = hist(Phi_1stWL,60)
# histogram(Phi_DL,60);
plt.plot(x,smooth(smooth(f)) / len(AtomIndx.O) / 2 / nConfigs,'b','linewidth',2)
plt.xlabel('\phi (deg)')
plt.ylabel('Abundance')
plt.legend('Global','1st+2nd Layer','1st Layer')
hold('off')
# # # # # # # # plot the average value of theta as a function of z-coordinate
# # # # # # # p=200;
# # # # # # # for i = 1:p
# # # # # # #     Position(i) = (i*max(zCoords)/p)+(((i*max(zCoords)/p)-((i-1)*max(zCoords)/p))/2);
# # # # # # #     AveAng(i) = mean(Theta(find(zCoords < i*max(zCoords)/p & zCoords > (i-1)*max(zCoords)/p)));
# # # # # # #     StdAng(i) = std(Theta(find(zCoords < i*max(zCoords)/p & zCoords > (i-1)*max(zCoords)/p)));
# # # # # # #
# # # # # # # end
# # # # # # # figure
# # # # # # # errorbar(Position, AveAng, StdAng, '-o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'color', 'k');
# # # # # # # xlabel('z (Ang)');
# # # # # # # ylabel('Ave. \theta');
# # # # # # # set(gca, 'xlim', [0 ABC(3)]);
# # # # # # # set(gcf, 'position', [-1646         235        1196         420]);
# # # # # # #
# # # # # # # return
# # # # # # #
# # # # # # #
