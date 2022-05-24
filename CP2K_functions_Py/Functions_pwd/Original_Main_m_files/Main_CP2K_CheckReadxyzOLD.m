clear all;
close all;
newlinechar = char(10);

%     BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
%     fldrname = 'CP_Pit_18H22F';
%     Trajectory = 'Sample11500_15000.xyz';
BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CorrectVolume\Metadynamics\';
fldrname = 'Pt_08H12F';
Trajectory = 'Pt_08H12F.xyz';

ABC = getABCvectors(BaseFldr, fldrname);
[xyz, XYZ, AtomIndx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, fldrname, Trajectory, ABC);
AtomsNew = Atoms;
% % Theta = [];
% % Phi1 = [];
% % Phi2 = [];

RadFunOH = cell(nConfigs,1);
RadFunFH = cell(nConfigs,1);
RadFunHH = cell(nConfigs,1);
RadFunOO = cell(nConfigs,1);

% % HydroListSnaps = cell(length(Sampling),1);
% % ThetaSnaps = cell(length(Sampling),length(IndxO));
% % SurfOSnaps = cell(length(Sampling),length(IndxO));
% % PhiSnaps = cell(length(Sampling),length(IndxO));

for snap = startConfig:nConfigs
%     disp(['Processing snapshot ' num2str(StepNum(snap)) ' - ' num2str(100*(StepNum(snap)/max(StepNum))) ' % complete']);
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
    [VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.O, AtomIndx.H, ABC);
%     [VecFH, DistFH] = GetAtomCorrelation(XYZ_snap, Indx.F, Indx.H, ABC);
%     [VecHH, DistHH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.H, ABC);
%     [VecOO, DistOO] = GetAtomCorrelation(XYZ_snap, Indx.O, Indx.O, ABC);
    
%     RadFunOH{snap} = reshape(DistOH, [numel(DistOH), 1]);
%     RadFunFH{snap} = reshape(DistFH, [numel(DistFH), 1]);
%     RadFunHH{snap} = reshape(DistHH, [numel(DistHH), 1]);
%     RadFunOO{snap} = reshape(DistOO, [numel(DistOO), 1]);
    
end

% RadialDistribution(RadFunOH, ABC, ['O'; 'H'])
% RadialDistribution(RadFunFH, ABC, ['H'; 'F'])
% RadialDistribution(RadFunHH, ABC, ['H'; 'H'])
% RadialDistribution(RadFunOO, ABC, ['O'; 'O'])
HydroList = [];
% return
    for j = 1:length(AtomIndx.O)
[VecOH, DistOH] = GetAtomCorrelation(XYZ_snap, AtomIndx.O(j), AtomIndx.H, ABC);
        % scan through each O atom and find all O-H vectors
% % % % %         VecOH = XYZ_snap(Indx.H,:) - XYZ_snap(Indx.O(j),:);
% % % % %         % find all O-H vector magnitudes
% % % % %         DistOH = sqrt((VecOH(:,1).^2)+(VecOH(:,2).^2)+(VecOH(:,3).^2));
% % % % % % % % %
% % % % % % % % %         % scan through each O atom and find all O-Pt vectors
% % % % % % % % %         VecOPt = XYZ(IndxPt,:) - XYZ(IndxO(j),:);
% % % % % % % % %         % find all O-Pt vector magnitudes
% % % % % % % % %         DistOPt = sqrt((VecOPt(:,1).^2)+(VecOPt(:,2).^2)+(VecOPt(:,3).^2));
% % % % % % % % %
% % % % %         % check for neighbours separated by PBC in x and y only
% % % % %         [VecOH, DistOH] = searchAcrossPBC(PBC_mat, VecOH, DistOH, XYZ_snap, Indx.O(j), Indx.H, ABC); % searchAcrossPBC(PBC_mat, non PBC Vector, non PBC distance, XYZ, looping atom IndxAtom(j), other atom IndxAtom, ABC)
% % % % % % % % %         [VecOPt, DistOPt] = searchAcrossPBC(PBC_mat, VecOPt, DistOPt, XYZ, IndxO(j), IndxPt, ABC);
% % % % % % % % %
        % find index of proximal H atoms to O within 1.3 Angstroms
        BondOH = find(DistOH < 1.4);
        % find the O atoms with three O-H bonds within 1.3 Ang, i.e. hydronium
        if length(BondOH) == 3
            %print the indices of the O atom and the furthest of the 3 H atoms from the O, i.e. the hydronium H
            HydroIdx = BondOH(find(max(DistOH(BondOH))));
%             disp(['O atom ' num2str(AtomIndx.O(j)) ' is bonded to 3 x H atoms @ XYZ = ' num2str(XYZ_snap(AtomIndx.O(j), :)) '. Hydronium H (' num2str(AtomIndx.H(HydroIdx)) ') @ XYZ = ' num2str(XYZ_snap(AtomIndx.H(HydroIdx),:))])
disp(['O atom ' num2str(AtomIndx.O(j)) ' is bonded to 3 x H atoms ' num2str(BondOH(1)) ' ' num2str(BondOH(2)) ' ' num2str(BondOH(3)) ])
AtomsNew(AtomIndx.O(j),:) = pad('Li', size(AtomsNew,2));
AtomsNew(AtomIndx.H(BondOH(1)),:) = pad('S', size(AtomsNew,2));
AtomsNew(AtomIndx.H(BondOH(2)),:) = pad('S', size(AtomsNew,2));
AtomsNew(AtomIndx.H(BondOH(3)),:) = pad('S', size(AtomsNew,2));

            %             HydroList = [HydroList; AtomIndx.H(HydroIdx) AtomIndx.O(j) AtomIndx.H(BondOH')];
        end
    end

    fidout = fopen([BaseFldr fldrname '\New_' Trajectory], 'w');
    fprintf(fidout,[num2str(nAtoms) newlinechar]);
    fprintf(fidout,[' i =   input geom' newlinechar]);
    for kk = 1:size(AtomsNew,1)
       fprintf(fidout,[pad(AtomsNew(kk,:),8) pad(num2str(XYZ_snap(kk,1), '%.10g'),20) pad(num2str(XYZ_snap(kk,2), '%.10g'),20) pad(num2str(XYZ_snap(kk,3), '%.10g'),20) newlinechar]); 
    end
    
    fclose(fidout);
    
        % find the O atoms with ONLY two O-H bonds within 1.3 Ang, i.e. water
% % %         if length(BondOH) == 2
% % % 
% % %             % calculate net vector bisecting H-O-H
% % %             VecH2O = VecOH(BondOH(1),:)+VecOH(BondOH(2),:);
% % % %
% % % %             % % %             % if we want to assume only one direction for both halves if the cell
% % % %             % % %                     zVec = [0 0 1];
% % % %             % % %                        Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
% % % %
% % % %             % split the cell into half and determine if O is in upper or lower half
% % % %             % then take the dot product of resultant water vector with both
% % % %             % surface normals (assuming the normal is c)
% % % %             if XYZ(IndxO(j),3) < ABC(3)/2
% % % %                 zVec = [0 0 1];
% % % %                 Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
% % % %                 Phi1 = acosd(dot(VecOH(BondOH(1),:),zVec)/(norm(VecOH(BondOH(1),:))*norm(zVec)));
% % % %                 Phi2 = acosd(dot(VecOH(BondOH(2),:),zVec)/(norm(VecOH(BondOH(2),:))*norm(zVec)));
% % % %             elseif XYZ(IndxO(j),3) >= ABC(3)/2
% % % %                 zVec = [0 0 -1];
% % % %                 Theta = acosd(dot(VecH2O,zVec)/(norm(VecH2O)*norm(zVec)));
% % % %                 Phi1 = acosd(dot(VecOH(BondOH(1),:),zVec)/(norm(VecOH(BondOH(1),:))*norm(zVec)));
% % % %                 Phi2 = acosd(dot(VecOH(BondOH(2),:),zVec)/(norm(VecOH(BondOH(2),:))*norm(zVec)));
% % % %             end
% % % %         end
% % % %
% % % %         HydroListSnaps{snap} = HydroList; % mod: could identify num of hydronium wanted and set if length HydroList > n hydronium, take only first n values according to shortest distance
% % % %         ThetaSnaps{snap,j} = Theta;
% % % %         PhiSnaps{snap,j} = [Phi1 Phi2];
% % % %         RadFunOHSnaps{snap,j} = DistOH;
% % % %         zPosition{snap,j} = XYZ(IndxO(j),3);
% % % %
% % % %         if sum(DistOPt<2.75 & DistOPt > 0) > 0
% % % %             %         if sum(DistOPt<3.75 & DistOPt > 2.75) > 0 & abs(XYZ(IndxO(j),3)-XYZ(IndxPt,3)) > 1.5
% % % %             %         if sum(DistOPt<4.75 & DistOPt > 3.75) > 0
% % % %             %         if sum(DistOPt<10 & DistOPt > 4.75) > 0
% % % %             SurfOSnaps{snap,j} = 1;
% % % %         else
% % % %             SurfOSnaps{snap,j} = 0;
% % % %         end
% % % %     end
% % % % end
% % % % return
% % % %
% % % % RadialDistribution(RadFunOHSnaps, ABC, ['O'; 'H']);


% % % % % % %
% % % % % % % zCoords=vertcat(zPosition{:,:});
% % % % % % % % Theta is the angle of the water H-O-H bisector with the surface normal/c-vector
% % % % % % % Theta=vertcat(ThetaSnaps{:,:});
% % % % % % % SurfO=vertcat(SurfOSnaps{:,:});
% % % % % % % Theta_DL = nonzeros(Theta.*SurfO);
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % title(['Global H-O-H Bisector Angle Distribution for ' fldrname], 'interpreter', 'none')
% % % % % % % [f, x] = hist(Theta,60);
% % % % % % % histogram(Theta,60);
% % % % % % % plot(x,smooth(f), 'r')
% % % % % % % xlabel('\Theta (deg)');
% % % % % % % ylabel('Abundance');
% % % % % % % hold off
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % title(['Double Layer H-O-H Bisector Angle Distribution for ' fldrname], 'interpreter', 'none')
% % % % % % % [f, x] = hist(Theta_DL,60);
% % % % % % % histogram(Theta_DL,60);
% % % % % % % plot(x,smooth(f), 'r')
% % % % % % % xlabel('\Theta (deg)');
% % % % % % % ylabel('Abundance');
% % % % % % % hold off
% % % % % % %
% % % % % % % % Phi is the angle of the O-H bond vector made with the surface normal
% % % % % % % Phi=reshape(vertcat(PhiSnaps{:}), [2*length(vertcat(PhiSnaps{:})),1]);
% % % % % % % Phi_DL = nonzeros(Phi.*[SurfO; SurfO]);
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % title(['Global O-H Angle Distribution for ' fldrname], 'interpreter', 'none')
% % % % % % % [f, x] = hist(Phi,60);
% % % % % % % histogram(Phi,60);
% % % % % % % plot(x,smooth(f), 'r')
% % % % % % % xlabel('\phi (deg)');
% % % % % % % ylabel('Abundance');
% % % % % % % hold off
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % title(['Double Layer O-H Angle Distribution for ' fldrname], 'interpreter', 'none')
% % % % % % % [f, x] = hist(Phi_DL,60);
% % % % % % % histogram(Phi_DL,60);
% % % % % % % plot(x,smooth(f), 'r')
% % % % % % % xlabel('\phi (deg)');
% % % % % % % ylabel('Abundance');
% % % % % % % hold off
% % % % % % %
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



% % % % % % %
% % % % % % %
% % % % % % %
% % % % % % %
% % % % % % % for j = 1:length(IndxF)
% % % % % % %     VecOF = XYZ(IndxO,:) - XYZ(IndxF(j),:);
% % % % % % %     DistOF = sqrt((VecOF(:,1).^2)+(VecOF(:,2).^2)+(VecOF(:,3).^2));
% % % % % % %
% % % % % % %     % check for neighbours separated by PBC in x and y only
% % % % % % %     for jj = 1:length(PBC_mat)
% % % % % % %         VecOF_PBC = XYZ(IndxO,:) - (XYZ(IndxF(j),:) + PBC_mat(jj,:).*ABC);
% % % % % % %         DistOF_PBC = sqrt((VecOF_PBC(:,1).^2)+(VecOF_PBC(:,2).^2)+(VecOF_PBC(:,3).^2));
% % % % % % %
% % % % % % %         % determine the shortest vector between O atom IndxO(j) and every H after each PBC element is applied
% % % % % % %         if any(DistOF_PBC<DistOF)
% % % % % % %             OFPBCIndx = find(DistOF_PBC<DistOF);
% % % % % % %             VecOF(OFPBCIndx,:) = VecOF_PBC(OFPBCIndx,:);
% % % % % % %         end
% % % % % % %
% % % % % % %         % determine the shortest distance between O atom IndxO(j) and every H after each PBC element is applied
% % % % % % %         DistOF = min(DistOF_PBC, DistOF);
% % % % % % %     end
% % % % % % %
% % % % % % %     RadFunOF = [RadFunOF; DistOF];
% % % % % % %
% % % % % % %     VecHF = XYZ(IndxH,:) - XYZ(IndxF(j),:);
% % % % % % %     DistHF = sqrt((VecHF(:,1).^2)+(VecHF(:,2).^2)+(VecHF(:,3).^2));
% % % % % % %
% % % % % % %     % check for neighbours separated by PBC in x and y only
% % % % % % %     for jj = 1:length(PBC_mat)
% % % % % % %         VecHF_PBC = XYZ(IndxH,:) - (XYZ(IndxF(j),:) + PBC_mat(jj,:).*ABC);
% % % % % % %         DistHF_PBC = sqrt((VecHF_PBC(:,1).^2)+(VecHF_PBC(:,2).^2)+(VecHF_PBC(:,3).^2));
% % % % % % %
% % % % % % %         % determine the shortest vector between O atom IndxO(j) and every H after each PBC element is applied
% % % % % % %         if any(DistHF_PBC<DistHF)
% % % % % % %             HFPBCIndx = find(DistHF_PBC<DistHF);
% % % % % % %             VecHF(HFPBCIndx,:) = VecHF_PBC(HFPBCIndx,:);
% % % % % % %         end
% % % % % % %
% % % % % % %         % determine the shortest distance between O atom IndxO(j) and every H after each PBC element is applied
% % % % % % %         DistHF = min(DistHF_PBC, DistHF);
% % % % % % %     end
% % % % % % %
% % % % % % %
% % % % % % %     RadFunHF = [RadFunHF; DistHF];
% % % % % % %
% % % % % % %
% % % % % % %     %         VecFF = XYZ(IndxF,:) - XYZ(IndxF(j),:);
% % % % % % %     %     DistHF = sqrt((VecHF(:,1).^2)+(VecHF(:,2).^2)+(VecHF(:,3).^2));
% % % % % % %     %
% % % % % % %     %     % check for neighbours separated by PBC in x and y only
% % % % % % %     %     for jj = 1:length(PBC_mat)
% % % % % % %     %         VecHF_PBC = XYZ(IndxF,:) - (XYZ(IndxH(j),:) + PBC_mat(jj,:).*ABC);
% % % % % % %     %         DistHF_PBC = sqrt((VecHF_PBC(:,1).^2)+(VecHF_PBC(:,2).^2)+(VecHF_PBC(:,3).^2));
% % % % % % %     %
% % % % % % %     %         % determine the shortest vector between O atom IndxO(j) and every H after each PBC element is applied
% % % % % % %     %         if any(DistHF_PBC<DistHF)
% % % % % % %     %             HFPBCIndx = find(DistHF_PBC<DistHF);
% % % % % % %     %             VecHF(HFPBCIndx,:) = VecHF_PBC(HFPBCIndx,:);
% % % % % % %     %         end
% % % % % % %     %
% % % % % % %     %         % determine the shortest distance between O atom IndxO(j) and every H after each PBC element is applied
% % % % % % %     %         DistHF = min(DistHF_PBC, DistHF);
% % % % % % %     %     end
% % % % % % %
% % % % % % %
% % % % % % %     %     RadFunHF = [RadFunHF; DistHF];
% % % % % % % end
% % % % % % %
% % % % % % %
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % % RadFunSrtd = sort(RadFunOF);
% % % % % % % [f, x] = hist(RadFunSrtd,500);
% % % % % % % % histogram(RadFunSrtd,500);
% % % % % % % plot(x,smooth(f), 'b')
% % % % % % %
% % % % % % % RadFunSrtd = sort(RadFunHF);
% % % % % % % [f, x] = hist(RadFunSrtd,1000);
% % % % % % % % histogram(RadFunSrtd,500);
% % % % % % % plot(x,smooth(f), 'r')
% % % % % % % xlabel('r (Ang)');
% % % % % % % ylabel('g(r)');
% % % % % % % hold off
% % % % % % %
% % % % % % % % % % % % % % % %% Find H atom indexes for those bound to the surface
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % LowerPt = XYZ(IndxPt,3) < ABC(3)/2; %redundant?
% % % % % % % % % % % % % % % [LowVal LowIndx] = sort(XYZ(IndxPt(LowerPt),3), 'descend');%redundant?
% % % % % % % % % % % % % % % LowIndx = IndxPt(LowIndx(1:(length(IndxPt)/4)));%redundant?
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % UpperPt = IndxPt(XYZ(IndxPt,3) > ABC(3)/2);
% % % % % % % % % % % % % % % [UpperVal UpperIndx] = sort(XYZ(UpperPt,3), 'ascend');
% % % % % % % % % % % % % % % UpperIndx = UpperPt(UpperIndx(1:(length(IndxPt)/4)));
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % LowerPt = IndxPt(XYZ(IndxPt,3) < ABC(3)/2);
% % % % % % % % % % % % % % % [LowerVal LowerIndx] = sort(XYZ(LowerPt,3), 'descend');
% % % % % % % % % % % % % % % LowerIndx = LowerPt(LowerIndx(1:(length(IndxPt)/4)));
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % count =0;Hsurfs = [];
% % % % % % % % % % % % % % % for j = 1:length(IndxH)
% % % % % % % % % % % % % % %     VecPtH = XYZ(IndxH(j),:) - XYZ(UpperIndx,:);
% % % % % % % % % % % % % % %     DistPtH = sqrt((VecPtH(:,1).^2)+(VecPtH(:,2).^2)+(VecPtH(:,3).^2));
% % % % % % % % % % % % % % %     BondPtH = find(DistPtH < 2.29);
% % % % % % % % % % % % % % %     if ~isempty(BondPtH)
% % % % % % % % % % % % % % %        count = count +1;
% % % % % % % % % % % % % % %        disp(['H atom ' num2str(IndxH(j)) ' is ' num2str(mean(DistPtH(BondPtH))) ' from the surface']);
% % % % % % % % % % % % % % %        Hsurfs = [Hsurfs; IndxH(j)];
% % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % count =0;
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % for j = 1:length(IndxH)
% % % % % % % % % % % % % % %     VecPtH = XYZ(IndxH(j),:) - XYZ(LowerIndx,:);
% % % % % % % % % % % % % % %     DistPtH = sqrt((VecPtH(:,1).^2)+(VecPtH(:,2).^2)+(VecPtH(:,3).^2));
% % % % % % % % % % % % % % %     BondPtH = find(DistPtH < 2.6);
% % % % % % % % % % % % % % %     if ~isempty(BondPtH)
% % % % % % % % % % % % % % %        count = count +1;
% % % % % % % % % % % % % % %        disp(['H atom ' num2str(IndxH(j)) ' is ' num2str(mean(DistPtH(BondPtH))) ' from the surface']);
% % % % % % % % % % % % % % %        Hsurfs = [Hsurfs; IndxH(j)];
% % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % end
% % % % % % %
