clear all;
close all;

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/5lyr_systems/Al_AlO_OH/AlO_OH/';
    system = 'AlO_0.33ML_OH';
    Trajectory = 'AlO_0.33ML_OH_50000to134000_1000step.xyz';


% BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\';
%     system = 'CP_Like_1010_Fluorine';
%     Trajectory = 'CP_Like_1010_Fluorine_46000to58000_500step.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);

[~, ~, AtomIndx, ~, ~, ~, ~] = getAtomInfoFromInput(BaseFldr, system);
%[Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC, 120);
[Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);


for snap = startConfig:nConfigs
%     disp(['Processing snapshot ' num2str(StepNum(snap)) ' - ' num2str(100*(snap/nConfigs)) ' % complete']);
    [FirstLayerIndx, SecondLayerIndx, MinimaZ] = getWaterLayerIndicesPerSnap_new(Indx, XYZ, Dens_O, z);

    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.O, AtomIndx.H, ABC);
    
    for j = 1:length(AtomIndx.O)
        BondOH = find(DistOH(:,j) <= 1.28);
        % find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
        if length(BondOH) == 2
            
            % calculate net vector bisecting H-O-H
            VecH2O = VecOH(BondOH(1),:,j)+VecOH(BondOH(2),:,j);
            
            % % %             % if we want to assume only one direction for both halves if the cell
            zVec = [0 0 1];
            Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
            Phi1 = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
            Phi2 = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
            
            % split the cell into half and determine if O is in upper or lower half
            % then take the dot product of resultant water vector with both
            % surface normals (assuming the normal is c)
%             if XYZ_snap(AtomIndx.O(j),3) < ABC(3)/2
%                 zVec = [0 0 1];
%                 Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
%                 Phi1 = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
%                 Phi2 = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
%             elseif XYZ_snap(AtomIndx.O(j),3) >= ABC(3)/2
%                 zVec = [0 0 -1];
%                 Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
%                 Phi1 = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
%                 Phi2 = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
%             end
            ThetaSnaps{snap,j} = Theta;
            PhiSnaps{snap,j} = [Phi1 Phi2];

        else
            ThetaSnaps{snap,j} = [];
            PhiSnaps{snap,j} = [];
        end

    end
    
    DL_Indx = [SecondLayerIndx{snap}];
    %     DL_Indx = SecondLayerIndx{snap};
    DL_Indx_H = [];
    %     DL_Indx = FirstLayerIndx{snap};
    
    for j = 1:length(DL_Indx)
        BondOH = find(DistOH(:,find(AtomIndx.O == DL_Indx(j))) <= 1.28);
        % find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
        if length(BondOH) == 2
            
            % calculate net vector bisecting H-O-H
            VecH2O = VecOH(BondOH(1),:,j)+VecOH(BondOH(2),:,j);
            
            % % %             % if we want to assume only one direction for both halves if the cell
            zVec = [0 0 1];
            ThetaDL = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
            Phi1DL = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
            Phi2DL = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
            
            % split the cell into half and determine if O is in upper or lower half
            % then take the dot product of resultant water vector with both
            % surface normals (assuming the normal is c)
%             if XYZ_snap(DL_Indx(j),3) < ABC(3)/2
%                 zVec = [0 0 1];
%                 ThetaDL = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
%                 Phi1DL = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
%                 Phi2DL = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
%             elseif XYZ_snap(DL_Indx(j),3) >= ABC(3)/2
%                 zVec = [0 0 -1];
%                 ThetaDL = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
%                 Phi1DL = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
%                 Phi2DL = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
%             end
        end
        
        ThetaSnapsDL{snap,j} = ThetaDL;
        PhiSnapsDL{snap,j} = [Phi1DL Phi2DL];
        
        DL_Indx_H = [DL_Indx_H; AtomIndx.H(BondOH)];
    end
    
    % write a sample snapshot to .xyz and convert to MS file (only for last
    % snapshot - simple bug that number of atoms in DL changes so arrays
    % are non-consistent for concatenation when trying to do all snaps together - could modify using cell array)
%     if snap == nConfigs
%         writeSnaptoxyz(BaseFldr, system, snap, XYZ_snap, Atoms, [DL_Indx; DL_Indx_H; AtomIndx.Pts; AtomIndx.Ptb; AtomIndx.Ptss] , 'DL_Density');
%         CP2kOptimPathParse(BaseFldr,system, ['DL_Density_' num2str(snap) '.xyz'])
%     end
    
    DL_Indx = FirstLayerIndx{snap};
    for j = 1:length(DL_Indx)
        BondOH = find(DistOH(:,find(AtomIndx.O == DL_Indx(j))) <= 1.28);
        % find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
        if length(BondOH) == 2
            
            % calculate net vector bisecting H-O-H
            VecH2O = VecOH(BondOH(1),:,j)+VecOH(BondOH(2),:,j);
            
            % % %             % if we want to assume only one direction for both halves if the cell
            zVec = [0 0 1];
            ThetaDL = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
            Phi1DL = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
            Phi2DL = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
            
            % split the cell into half and determine if O is in upper or lower half
            % then take the dot product of resultant water vector with both
            % surface normals (assuming the normal is c)
%             if XYZ_snap(DL_Indx(j),3) < ABC(3)/2
%                 zVec = [0 0 1];
%                 ThetaDL = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
%                 Phi1DL = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
%                 Phi2DL = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
%             elseif XYZ_snap(DL_Indx(j),3) >= ABC(3)/2
%                 zVec = [0 0 -1];
%                 ThetaDL = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
%                 Phi1DL = acosd(dot(VecOH(BondOH(1),:,j),zVec)/(norm(VecOH(BondOH(1),:,j))*norm(zVec)));
%                 Phi2DL = acosd(dot(VecOH(BondOH(2),:,j),zVec)/(norm(VecOH(BondOH(2),:,j))*norm(zVec)));
%             end
        end
        
        ThetaSnaps1stWL{snap,j} = ThetaDL;
        PhiSnaps1stWL{snap,j} = [Phi1DL Phi2DL];
        
        %                 DL_Indx_H = [DL_Indx_H; AtomIndx.H(BondOH)];
    end
    
end

% % % % % % % % Theta is the angle of the water H-O-H bisector with the surface normal/c-vector
Theta=vertcat(ThetaSnaps{:,:});
Theta_DL=vertcat(ThetaSnapsDL{:,:});
Theta_1stWL=vertcat(ThetaSnaps1stWL{:,:});

for j = 1:size(ThetaSnapsDL,1)
    nH2O_DL(j) = size(ThetaSnapsDL,2) - nnz(cellfun(@isempty,ThetaSnapsDL(j,:)));
end

for j = 1:size(ThetaSnaps1stWL,1)
    nH2O_1stWL(j) = size(ThetaSnaps1stWL,2) - nnz(cellfun(@isempty,ThetaSnaps1stWL(j,:)));
end

disp(['Ave. number of water molecules in 2nd WL = ' num2str(mean(nH2O_DL)/2) ' \pm ' num2str(std(nH2O_DL)/2)]);
disp(['Ave. number of water molecules in 1st WL = ' num2str(mean(nH2O_1stWL)/2) ' \pm ' num2str(std(nH2O_1stWL)/2)]);

figure
hold on
title(['H-O-H Bisector Angle Distribution for ' system], 'interpreter', 'none')
[f, x] = hist(Theta,60);
% histogram(Theta,60);
plot(x,smooth(smooth(f))/length(AtomIndx.O)/nConfigs, 'k', 'linewidth', 2)
[f, x] = hist(Theta_DL,60);
% histogram(Theta_DL,60);
plot(x,smooth(smooth(f))/length(AtomIndx.O)/nConfigs, 'r', 'linewidth', 2)
[f, x] = hist(Theta_1stWL,60);
% histogram(Theta_DL,60);
plot(x,smooth(smooth(f))/length(AtomIndx.O)/nConfigs, 'b', 'linewidth', 2)
xlabel('\Theta (deg)');
ylabel('Abundance');
legend('Global', '1st+2nd Layer', '1st Layer');
hold off

% % % % % % % % Phi is the angle of the O-H bond vector made with the surface normal
Phi=reshape(vertcat(PhiSnaps{:}), [2*length(vertcat(PhiSnaps{:})),1]);
Phi_DL=reshape(vertcat(PhiSnapsDL{:}), [2*length(vertcat(PhiSnapsDL{:})),1]);
Phi_1stWL=reshape(vertcat(PhiSnaps1stWL{:}), [2*length(vertcat(PhiSnaps1stWL{:})),1]);

figure
hold on
title(['O-H Angle Distribution for ' system], 'interpreter', 'none')
[f, x] = hist(Phi,60);
% histogram(Phi,60);
plot(x,smooth(smooth(f))/length(AtomIndx.O)/2/nConfigs, 'k', 'linewidth', 2)
[f, x] = hist(Phi_DL,60);
% histogram(Phi_DL,60);
plot(x,smooth(smooth(f))/length(AtomIndx.O)/2/nConfigs, 'r', 'linewidth', 2)
[f, x] = hist(Phi_1stWL,60);
% histogram(Phi_DL,60);
plot(x,smooth(smooth(f))/length(AtomIndx.O)/2/nConfigs, 'b', 'linewidth', 2)
xlabel('\phi (deg)');
ylabel('Abundance');
legend('Global', '1st+2nd Layer', '1st Layer');
hold off

% % % % % % % % plot the average value of theta as a function of z-coordinate
% % % % % % % p=200;
% % % % % % % for i = 1:p
% % % % % % %     Position(i) = (i*max(zCoords)/p)+(((i*max(zCoords)/p)-((i-1)*max(zCoords)/p))/2);
% % % % % % %     AveAng(i) = mean(Theta(find(zCoords < i*max(zCoords)/p & zCoords > (i-1)*max(zCoords)/p)));
% % % % % % %     StdAng(i) = std(Theta(find(zCoords < i*max(zCoords)/p & zCoords > (i-1)*max(zCoords)/p)));
% % % % % % %
% % % % % % % end
% % % % % % % figure
% % % % % % % errorbar(Position, AveAng, StdAng, '-o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'color', 'k');
% % % % % % % xlabel('z (Ang)');
% % % % % % % ylabel('Ave. \theta');
% % % % % % % set(gca, 'xlim', [0 ABC(3)]);
% % % % % % % set(gcf, 'position', [-1646         235        1196         420]);
% % % % % % %
% % % % % % % return
% % % % % % %
% % % % % % %


