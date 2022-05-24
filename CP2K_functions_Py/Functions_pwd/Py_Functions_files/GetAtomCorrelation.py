import numpy as np
    
def GetAtomCorrelation(XYZ = None,Indx1 = None,Indx2 = None,ABC = None): 
    Dist12 = np.zeros((len(Indx2),len(Indx1)))
    for j in np.arange(1,len(Indx1)+1).reshape(-1):
        # scan through each atom 1 and find all atom1-atom2 vectors
        Vec12 = XYZ(Indx2,:) - XYZ(Indx1(j),:)
        # find all atom1-atom2 vector magnitudes
        Dist12[:,j] = np.sqrt((Vec12(:,1) ** 2) + (Vec12(:,2) ** 2) + (Vec12(:,3) ** 2))
        # check for neighbours separated by PBC in x and y only
        Vec12,Dist12[:,j] = searchAcrossPBC(Vec12,Dist12(:,j),XYZ,Indx1(j),Indx2,ABC)
    
    return Vec12,Dist12
    # # # #     # find the O atoms with three O-H bonds within 1.3 Ang, i.e. hydronium
# # # #     if length(Bond12) == 3
# # # #         #print the indices of the O atom and the furthest of the 3 H atoms from the O, i.e. the hydronium H
# # # #         HydroIdx = Bond12(find(max(Dist12(Bond12))));
# # # #         disp(['O atom ' num2str(Indx1(j)) ' is bonded to 3 x H atoms @ XYZ = ' num2str(XYZ(Indx1(j), :)) '. Hydronium H (' num2str(Indx2(HydroIdx)) ') @ XYZ = ' num2str(XYZ(Indx2(HydroIdx),:))])
# # # #         HydroList = [HydroList; Indx2(HydroIdx) Indx1(j) Indx2(Bond12')];
# # # #     end
# # # #
# # # #     # find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
# # # #     if length(Bond12) == 2
# # # #
# # # #         # calculate net vector bisecting H-O-H
# # # #         VecH2O = Vec12(Bond12(1),:)+Vec12(Bond12(2),:);
# # # #
# # # #         # # #             # if we want to assume only one direction for both halves if the cell
# # # #         # # #                     zVec = [0 0 1];
# # # #         # # #                        Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
# # # #
# # # #         # split the cell into half and determine if O is in upper or lower half
# # # #         # then take the dot product of resultant water vector with both
# # # #         # surface normals (assuming the normal is c)
# # # #         if XYZ(Indx1(j),3) < ABC(3)/2
# # # #             zVec = [0 0 1];
# # # #             Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
# # # #             Phi1 = acosd(dot(Vec12(Bond12(1),:),zVec)/(norm(Vec12(Bond12(1),:))*norm(zVec)));
# # # #             Phi2 = acosd(dot(Vec12(Bond12(2),:),zVec)/(norm(Vec12(Bond12(2),:))*norm(zVec)));
# # # #         elseif XYZ(Indx1(j),3) >= ABC(3)/2
# # # #             zVec = [0 0 -1];
# # # #             Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
# # # #             Phi1 = acosd(dot(Vec12(Bond12(1),:),zVec)/(norm(Vec12(Bond12(1),:))*norm(zVec)));
# # # #             Phi2 = acosd(dot(Vec12(Bond12(2),:),zVec)/(norm(Vec12(Bond12(2),:))*norm(zVec)));
# # # #         end
# # # #     end
# # # #
# # # #     HydroListSnaps{snap} = HydroList; # mod: could identify num of hydronium wanted and set if length HydroList > n hydronium, take only first n values according to shortest distance
# # # #     ThetaSnaps{snap,j} = Theta;
# # # #     PhiSnaps{snap,j} = [Phi1 Phi2];
# # # #     RadFunOHSnaps{snap,j} = Dist12;
# # # #     zPosition{snap,j} = XYZ(Indx1(j),3);
# # # #
# # # #     if sum(DistOPt<2.75 & DistOPt > 0) > 0
# # # #         #         if sum(DistOPt<3.75 & DistOPt > 2.75) > 0 & abs(XYZ(Indx1(j),3)-XYZ(IndxPt,3)) > 1.5
# # # #         #         if sum(DistOPt<4.75 & DistOPt > 3.75) > 0
# # # #         #         if sum(DistOPt<10 & DistOPt > 4.75) > 0
# # # #         SurfOSnaps{snap,j} = 1;
# # # #     else
# # # #         SurfOSnaps{snap,j} = 0;
# # # #     end
# # # # end
    return Vec12,Dist12