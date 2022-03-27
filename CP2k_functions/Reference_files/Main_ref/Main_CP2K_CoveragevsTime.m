clear all;
close all;

BaseFldr = 'Z:\Imperial\MattProjects\OxidesOER\RutheniumOxide\'; 
system = 'RuO2_u1_1ML_u2_1ML';
Trajectory = 'RuO2_u1_1ML_u2_1ML_0to15000_500step.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, ~, ~, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC);
[Atoms, AtomList, AtomIndx] = getAtomNamesFromInputXYZ(BaseFldr, system);

Coverage = zeros(nConfigs,2);

for snap = startConfig:nConfigs
    disp(['Processing snapshot ' num2str(StepNum(snap)) ' - ' num2str(100*(snap/nConfigs)) ' % complete']);
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
    
%     [VecOtUH, DistOtUH] = GetAtomCorrelation(XYZ_snap, AtomIndx.OtU, AtomIndx.H, ABC);
%     [VecOtLH, DistOtLH] = GetAtomCorrelation(XYZ_snap, AtomIndx.OtL, AtomIndx.H, ABC);
    [VecOtUH, DistOtUH] = GetAtomCorrelation(XYZ_snap, AtomIndx.OtU, [AtomIndx.H; AtomIndx.Hsurf], ABC);
    [VecOtLH, DistOtLH] = GetAtomCorrelation(XYZ_snap, AtomIndx.OtL, [AtomIndx.H; AtomIndx.Hsurf], ABC);
    
    Coverage(snap,1) = sum(DistOtUH < 1.25, 'all');
    Coverage(snap,2) = sum(DistOtLH < 1.25, 'all');
end

figure
box on
hold on
xlabel('Time (ps)');
ylabel('\theta_H (ML)');
plot(StepNum/2000, Coverage(:,1)/length(AtomIndx.OtU), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'b');
plot(StepNum/2000, Coverage(:,2)/length(AtomIndx.OtL), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', 'r');
plot(StepNum/2000, (Coverage(:,1)+Coverage(:,2))/(length(AtomIndx.OtL)+length(AtomIndx.OtU)), '-o', 'color', 'k', 'markeredgecolor', 'k', 'markerfacecolor', [0 0.5 0]);
set(gca, 'YTick', 0:0.1:2, 'fontsize', 12)
legend('\mu_1', '\mu_2', 'Total', 'interpreter', 'tex', 'location', 'best')
hold off