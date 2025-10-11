clear all;  clc;
% close all;
if 1
    BaseFldr = 'D:\';
    system = 'CP_Pit_20F';
    Trajectory = 'CP_Pit_20F_43000to73000_500step.xyz';
    
    % % Call function to find ABC vectors from .inp file
    ABC = getABCvectors(BaseFldr, system);
    
% % get the names of atoms from original xyz input file
[Atoms, AtomList, AtomIndx, ~, ~, ~, ~] = getAtomInfoFromInput(BaseFldr, system);
    

    % % Read the xyz data from "Trajectory"
    [xyz, XYZ, Indx, ~, ~, nAtoms, startConfig, nConfigs, StepNum_Traj] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0; 0; 0]);

    
%     [Dens_O, Dens_H, Dens_F, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
    [Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC);
    
    [FirstLayerIndx, SecondLayerIndx, ThirdLayerIndx, FourthLayerIndx, Minima] = getWaterLayerIndices(AtomIndx, XYZ, Dens_O, z);
    
    
end

r1st = FirstLayerIndx{end};
r2nd = SecondLayerIndx{end};
r3rd = ThirdLayerIndx{end};
r4th = FourthLayerIndx{end};

%         AvePos = mean(xyz.O(:, r1st, :),1)
figure
hold on
axis equal
set(gcf, 'Position', [-1919          41        1920         963])
% set(gca, 'Position', [0.1300    0.1108    0.7750*ABC(1)/ABC(2) 0.7750])

XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
XYZ_snap(:,:) = XYZ(end,:,:);

[VecOH1stWL, DistOH1stWL] = GetAtomCorrelation(XYZ_snap, r1st, AtomIndx.H, ABC);
[VecOH2ndWL, DistOH2ndWL] = GetAtomCorrelation(XYZ_snap, r2nd, AtomIndx.H, ABC);
[VecOH3rdWL, DistOH3rdWL] = GetAtomCorrelation(XYZ_snap, r3rd, AtomIndx.H, ABC);
[VecOH4thWL, DistOH4thWL] = GetAtomCorrelation(XYZ_snap, r4th, AtomIndx.H, ABC);
H1st = [];
H2nd = [];
H3rd = [];
H4th = [];

WL1 = [r1st];
WL2 = [r2nd];
WL3 = [r3rd];
WL4 = [r4th];

Replica = [0 0 0];
for rx = 0:Replica(1)
    for ry =0:Replica(2)
        % AvePosPtE = mean(xyz.Pt(:, AtomIndx.PtE, :),1)
        % plot3(AvePosPtE(1,:,1),AvePosPtE(1,:,2),AvePosPtE(1,:,3),'o', 'markerfacecolor', 'c', 'markeredgecolor', 'k', 'markersize', 30)
        
% % % %         AvePosPt = mean(XYZ(:, [AtomIndx.Pts; AtomIndx.Ptss; AtomIndx.Ptb], :),1);
% % % %         AvePosPtb = mean(XYZ(:, [AtomIndx.Ptb], :),1);
% % % %         % AvePosPtS = XYZ(end, [AtomIndx.Pts; AtomIndx.Ptss; AtomIndx.Ptb], :);
% % % % 
% % % %         plot3(AvePosPt(1,:,1)+rx*ABC(1),AvePosPt(1,:,2)+ry*ABC(2),AvePosPt(1,:,3),'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markersize', 30)
% % % %         plot3(AvePosPtb(1,:,1)+rx*ABC(1),AvePosPtb(1,:,2)+ry*ABC(2),AvePosPtb(1,:,3)-ABC(3),'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markersize', 30)
    


%  AvePosPt = mean(XYZ(:, [AtomIndx.Al1; AtomIndx.Al2; AtomIndx.Alb], :),1);
%         AvePosPtb = mean(XYZ(:, [AtomIndx.Alb], :),1);
        AvePosPt = XYZ(end, [AtomIndx.Pts; AtomIndx.PtE; AtomIndx.Ptb], :);
        AvePosPtb = mean(XYZ(:, [AtomIndx.Ptb], :),1);

        plot3(AvePosPt(1,:,1)+rx*ABC(1),AvePosPt(1,:,2)+ry*ABC(2),AvePosPt(1,:,3),'o', 'markerfacecolor', [0.9 0.9  0.9], 'markeredgecolor', 'k', 'markersize', 30)
        
%         %%plot an extra z-replica of the metal
%         plot3(AvePosPtb(1,:,1)+rx*ABC(1),AvePosPtb(1,:,2)+ry*ABC(2),AvePosPtb(1,:,3)-ABC(3),'o', 'markerfacecolor', [0.9 0.9  0.9], 'markeredgecolor', 'k', 'markersize', 30)
        
        
        for i = 1:length(r1st)
            plot3(XYZ(end,r1st(i),1)+rx*ABC(1), XYZ(end,r1st(i),2)+ry*ABC(2), XYZ(end,r1st(i),3), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markersize', 8)
        end
        
        for i = 1:length(r2nd)
            plot3(XYZ(end,r2nd(i),1)+rx*ABC(1), XYZ(end,r2nd(i),2)+ry*ABC(2), XYZ(end,r2nd(i),3), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 8)
        end
        
        % for i = 1:length(r3rd)
        %     plot3(XYZ(end,r3rd(i),1), XYZ(end,r3rd(i),2), XYZ(end,r3rd(i),3), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 7)
        % end
        %
        % for i = 1:length(r4th)
        %     plot3(XYZ(end,r4th(i),1), XYZ(end,r4th(i),2), XYZ(end,r4th(i),3), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 7)
        % end
        
        for j = 1:length(r1st)
            H1stIndx = find(DistOH1stWL(:,j)<=1.28);
            H1st = AtomIndx.H(find(DistOH1stWL(:,j)<=1.28));
            
            WL1 = [WL1; H1st];
            
            for k = 1:length(H1st)
                plot3([XYZ(end,r1st(j),1)+rx*ABC(1) XYZ(end,r1st(j),1)+VecOH1stWL(H1stIndx(k),1,j)+rx*ABC(1)], [XYZ(end,r1st(j),2)+ry*ABC(2) XYZ(end,r1st(j),2)+VecOH1stWL(H1stIndx(k),2,j)+ry*ABC(2)], [XYZ(end,r1st(j),3) XYZ(end,r1st(j),3)+VecOH1stWL(H1stIndx(k),3,j)], '-', 'color', 'k', 'linewidth', 1);
            end
            
            HB1stIndx = find(DistOH1stWL(:,j)>1.28 & DistOH1stWL(:,j) <= 2.44);
            HB1st = AtomIndx.H(find(DistOH1stWL(:,j)>1.28 & DistOH1stWL(:,j) <= 2.44));
            for k = 1:length(HB1st)
                plot3([XYZ(end,r1st(j),1)+rx*ABC(1) XYZ(end,r1st(j),1)+VecOH1stWL(HB1stIndx(k),1,j)+rx*ABC(1)], [XYZ(end,r1st(j),2)+ry*ABC(2) XYZ(end,r1st(j),2)+VecOH1stWL(HB1stIndx(k),2,j)+ry*ABC(2)], [XYZ(end,r1st(j),3) XYZ(end,r1st(j),3)+VecOH1stWL(HB1stIndx(k),3,j)], '--m');
            end
            
            plot3(XYZ(end,H1st,1)+rx*ABC(1), XYZ(end,H1st,2), XYZ(end,H1st,3), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', 6)
            
        end
        
        for j = 1:length(r2nd)
            H2ndIndx = find(DistOH2ndWL(:,j)<1.28);
            H2nd = AtomIndx.H(find(DistOH2ndWL(:,j)<1.28));
            WL2 = [WL2; H2nd];
            
            for k = 1:length(H2nd)
                plot3([XYZ(end,r2nd(j),1)+rx*ABC(1) XYZ(end,r2nd(j),1)+VecOH2ndWL(H2ndIndx(k),1,j)+rx*ABC(1)], [XYZ(end,r2nd(j),2)+ry*ABC(2) XYZ(end,r2nd(j),2)+VecOH2ndWL(H2ndIndx(k),2,j)+ry*ABC(2)], [XYZ(end,r2nd(j),3) XYZ(end,r2nd(j),3)+VecOH2ndWL(H2ndIndx(k),3,j)], '-', 'color', 'k', 'linewidth', 1);
            end
            
            HB2ndIndx = find(DistOH2ndWL(:,j)>1.28 & DistOH2ndWL(:,j) <= 2.44);
            HB2nd = AtomIndx.H(find(DistOH2ndWL(:,j)>1.28 & DistOH2ndWL(:,j) <= 2.44));
            for k = 1:length(HB2nd)
                plot3([XYZ(end,r2nd(j),1)+rx*ABC(1) XYZ(end,r2nd(j),1)+VecOH2ndWL(HB2ndIndx(k),1,j)+rx*ABC(1)], [XYZ(end,r2nd(j),2)+ry*ABC(2) XYZ(end,r2nd(j),2)+VecOH2ndWL(HB2ndIndx(k),2,j)+ry*ABC(2)], [XYZ(end,r2nd(j),3) XYZ(end,r2nd(j),3)+VecOH2ndWL(HB2ndIndx(k),3,j)], '--m');
            end
            
            plot3(XYZ(end,H2nd,1)+rx*ABC(1), XYZ(end,H2nd,2)+ry*ABC(2), XYZ(end,H2nd,3), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', 6)
            
        end
        
%         writeSnaptoxyz(BaseFldr, system, 1, XYZ_snap, Atoms, [DL_Indx; DL_Indx_H; AtomIndx.Pts; AtomIndx.Ptb; AtomIndx.Ptss] , 'DL_Density');
%         CP2kOptimPathParse(BaseFldr,system, ['DL_Density_' num2str(snap) '.xyz'])
        
        % for j = 1:length(r3rd)
        %     H3rdIndx = find(DistOH3rdWL(:,j)<1.28);
        %     H3rd = AtomIndx.H(find(DistOH3rdWL(:,j)<1.28));
        %
        %     for k = 1:length(H3rd)
        %         plot3([XYZ(end,r3rd(j),1) XYZ(end,r3rd(j),1)+VecOH3rdWL(H3rdIndx(k),1,j)], [XYZ(end,r3rd(j),2) XYZ(end,r3rd(j),2)+VecOH3rdWL(H3rdIndx(k),2,j)], [XYZ(end,r3rd(j),3) XYZ(end,r3rd(j),3)+VecOH3rdWL(H3rdIndx(k),3,j)], '-k');
        %     end
        %
        %     HB3rdIndx = find(DistOH3rdWL(:,j)>1.28 & DistOH3rdWL(:,j) <= 2.44);
        %     HB3rd = AtomIndx.H(find(DistOH3rdWL(:,j)>1.28 & DistOH3rdWL(:,j) <= 2.44));
        %     for k = 1:length(HB3rd)
        %         plot3([XYZ(end,r3rd(j),1) XYZ(end,r3rd(j),1)+VecOH3rdWL(HB3rdIndx(k),1,j)], [XYZ(end,r3rd(j),2) XYZ(end,r3rd(j),2)+VecOH3rdWL(HB3rdIndx(k),2,j)], [XYZ(end,r3rd(j),3) XYZ(end,r3rd(j),3)+VecOH3rdWL(HB3rdIndx(k),3,j)], '--m');
        %     end
        %
        %     plot3(XYZ(end,H3rd,1), XYZ(end,H3rd,2), XYZ(end,H3rd,3), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', 5)
        %
        % end
        %
        % for j = 1:length(r4th)
        %     H4thIndx = find(DistOH4thWL(:,j)<1.28);
        %     H4th = AtomIndx.H(find(DistOH4thWL(:,j)<1.28));
        %     for k = 1:length(H4th)
        %         plot3([XYZ(end,r4th(j),1) XYZ(end,r4th(j),1)+VecOH4thWL(H4thIndx(k),1,j)], [XYZ(end,r4th(j),2) XYZ(end,r4th(j),2)+VecOH4thWL(H4thIndx(k),2,j)], [XYZ(end,r4th(j),3) XYZ(end,r4th(j),3)+VecOH4thWL(H4thIndx(k),3,j)], '-k');
        %     end
        %
        %     HB4thIndx = find(DistOH4thWL(:,j)>1.28 & DistOH4thWL(:,j) <= 2.44);
        %     HB4th = AtomIndx.H(find(DistOH4thWL(:,j)>1.28 & DistOH4thWL(:,j) <= 2.44));
        %     for k = 1:length(HB4th)
        %         plot3([XYZ(end,r4th(j),1) XYZ(end,r4th(j),1)+VecOH4thWL(HB4thIndx(k),1,j)], [XYZ(end,r4th(j),2) XYZ(end,r4th(j),2)+VecOH4thWL(HB4thIndx(k),2,j)], [XYZ(end,r4th(j),3) XYZ(end,r4th(j),3)+VecOH4thWL(HB4thIndx(k),3,j)], '--m');
        %     end
        %
        %     plot3(XYZ(end,H4th,1), XYZ(end,H4th,2), XYZ(end,H4th,3), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', 5)
        %
        % end
        
    end
end
LwD = 1;
Vec = diag(ABC);
plot3([0 Vec(1,1)],[0 Vec(1,2)],[0 Vec(1,3)],'-k','LineWidth',LwD)
plot3([0 Vec(2,1)],[0 Vec(2,2)],[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3([0 Vec(3,1)],[0 Vec(3,2)],[0 Vec(3,3)],'-k','LineWidth',LwD)

plot3(Vec(1,1)+[0 Vec(2,1)],Vec(1,2)+[0 Vec(2,2)],Vec(1,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(2,1)+[0 Vec(1,1)],Vec(2,2)+[0 Vec(1,2)],Vec(2,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)

plot3(Vec(1,1)+[0 Vec(3,1)],Vec(1,2)+[0 Vec(3,2)],Vec(1,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)
plot3(Vec(2,1)+[0 Vec(3,1)],Vec(2,2)+[0 Vec(3,2)],Vec(2,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)
plot3(Vec(1,1)+Vec(2,1)+[0 Vec(3,1)],Vec(1,2)+Vec(2,2)+[0 Vec(3,2)],Vec(1,3)+Vec(2,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)

plot3(Vec(3,1)+[0 Vec(1,1)],Vec(3,2)+[0 Vec(1,2)],Vec(3,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+[0 Vec(2,1)],Vec(3,2)+[0 Vec(2,2)],Vec(3,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+Vec(1,1)+[0 Vec(2,1)],Vec(3,2)+Vec(1,2)+[0 Vec(2,2)],Vec(3,3)+Vec(1,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+Vec(2,1)+[0 Vec(1,1)],Vec(3,2)+Vec(2,2)+[0 Vec(1,2)],Vec(3,3)+Vec(2,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)



return
%         plot3(AvePos(1,:,1),AvePos(1,:,2),AvePos(1,:,3),'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 15)
%         text(AvePos(1,:,1),AvePos(1,:,2),AvePos(1,:,3), num2str(r1st))
for i = 1:length(r1st)
    plot3(XYZ(end,r1st(i),1), XYZ(end,r1st(i),2), XYZ(end,r1st(i),3), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markersize', 5)
end
%         plot3(xyz.O(:,614,1), xyz.O(:,614,2), xyz.O(:,614,3), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markersize', 15)

r1stflat = setdiff(r1stflat,r1st);
AvePos1stflat = mean(xyz.O(:, r1stflat, :),1)
%         plot3(AvePos1stflat(1,:,1),AvePos1stflat(1,:,2),AvePos1stflat(1,:,3),'o', 'markerfacecolor', 'm', 'markeredgecolor', 'k', 'markersize', 15)
%        plot3(xyz.O(end,r1stflat,1), xyz.O(end,r1stflat,2), xyz.O(end,r1stflat,3), 'o', 'markerfacecolor', 'm', 'markeredgecolor', 'k', 'markersize', 15)

for i = 1:length(r1stflat)
    plot3(xyz.O(:,r1stflat(i),1), xyz.O(:,r1stflat(i),2), xyz.O(:,r1stflat(i),3), 'o', 'markerfacecolor', 'm', 'markeredgecolor', 'k', 'markersize', 5)
end




[r2nd, ~] = find(DistPtEO > MinimaPtEO(1) & DistPtEO <= MinimaPtEO(2));
r2nd = setdiff(r2nd,r1st);
AvePos2nd = mean(xyz.O(:, r2nd, :),1)

plot3(AvePos2nd(1,:,1),AvePos2nd(1,:,2),AvePos2nd(1,:,3),'o', 'markerfacecolor', [0 0.5 0], 'markeredgecolor', 'k', 'markersize', 15)


