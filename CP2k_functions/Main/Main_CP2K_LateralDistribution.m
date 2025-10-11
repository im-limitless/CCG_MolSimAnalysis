% CP2K post-processing script to track lateral O atom diffusion as a
% function of height from the surface. Created by Matt Darby 7th July 2022
% - v1.0
if 1
clear all;
close all;
clc;

%% Set the location of the calculation output
BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
system = 'CP_Pit_Water';
Trajectory = 'CP_Pit_Water_29500to40500_500step.xyz';

%% Call functions to parse the lattice vectors and xyz trajectory files
ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC, [pi/2; 1; 0]);
XYZ = wrapXYZ(XYZ, ABC);
end
%%

% Colors = [ones(length(Qmc),1) linspace(0,1,length(Qmc))' zeros(length(Qmc),1)];
% colormap(Colors);
% 
%     hsurf = surf(xAt,yAt,zAt,'FaceColor', Colors(cIndx,:),'EdgeColor','None');
%     set(hsurf,'AmbientStrength',0.1,...
close all
figure
hold on
xlabel('RMS lateral displacement (Ang)');
ylabel('z (Ang)')

for i = 1:size(xyz.O,2)
% for i = 1:10
    for j = 2:size(xyz.O,1)
        if xyz.O(j,i,2)-xyz.O(j-1,i,2) > ABC(2)/2
            xyz.O(j,i,2) = xyz.O(j,i,2) - ABC(2);
        elseif xyz.O(j,i,2)-xyz.O(j-1,i,2) < -ABC(2)/2
            xyz.O(j,i,2) = xyz.O(j,i,2) + ABC(2);
        end
        
        if xyz.O(j,i,1)-xyz.O(j-1,i,1) > ABC(1)/2
            xyz.O(j,i,1) = xyz.O(j,i,1) - ABC(1);
        elseif xyz.O(j,i,1)-xyz.O(j-1,i,1) < -ABC(1)/2
            xyz.O(j,i,1) = xyz.O(j,i,1) + ABC(1);
        end
    end
   
   plot(smooth(sqrt(((xyz.O(:,i,1)-xyz.O(1,i,1)).^2) + ((xyz.O(:,i,2)-xyz.O(1,i,2)).^2))), smooth(xyz.O(:,i,3))) 
   
   
end
hold off

figure
hold on

for i = 1:size(xyz.O,2)
%     [f, x] = hist(xyz.O(:,i,3), 9);
    histogram(xyz.O(:,i,3), 9);
       
% plot(xx,yy, 'b', 'linewidth', 2)
end
hold off


figure
scatter(mean(xyz.O(:,:,3),1), mean(sqrt((xyz.O(:,:,1)-xyz.O(1,:,1)).^2+(xyz.O(:,:,2)-xyz.O(1,:,2)).^2+(xyz.O(:,:,3)-xyz.O(1,:,3)).^2)), 'markerfacecolor', 'b', 'markeredgecolor', 'k')
xlabel('z (Ang)');
ylabel('RMS Displacement (Ang)');

figure
scatter(mean(xyz.O(:,:,3),1), std(sqrt((xyz.O(:,:,1)-xyz.O(1,:,1)).^2+(xyz.O(:,:,2)-xyz.O(1,:,2)).^2+(xyz.O(:,:,3)-xyz.O(1,:,3)).^2)), 'markerfacecolor', 'r', 'markeredgecolor', 'k')
xlabel('z (Ang)');
ylabel('Variance in RMS Displacement (Ang)');

figure
scatter(mean(xyz.O(:,:,3),1), mean(sqrt((xyz.O(:,:,1)-xyz.O(1,:,1)).^2+(xyz.O(:,:,2)-xyz.O(1,:,2)).^2)), 'markerfacecolor', 'm', 'markeredgecolor', 'k')
xlabel('z (Ang)');
ylabel('Lateral RMS Displacement (Ang)');

figure
scatter(mean(xyz.O(:,:,3),1), std(sqrt((xyz.O(:,:,1)-xyz.O(1,:,1)).^2+(xyz.O(:,:,2)-xyz.O(1,:,2)).^2)), 'markerfacecolor', 'y', 'markeredgecolor', 'k')
xlabel('z (Ang)');
ylabel('Variance in Lateral RMS Displacement (Ang)');


% figure
% muXY = mean(sqrt(((xyz.O(:,:,1)-xyz.O(1,:,1)).^2) + ((xyz.O(:,:,2)-xyz.O(1,:,2)).^2)),1);
% muZ = mean(xyz.O(:,:,3),1);
% sigXY = std(sqrt(((xyz.O(:,:,1)-xyz.O(1,:,1)).^2) + ((xyz.O(:,:,2)-xyz.O(1,:,2)).^2)),1);
% sigZ = std(xyz.O(:,:,3),1);
% scatter(muXY, muZ);
% scatter(sigXY, muZ);
% plot(muXY, smooth(xyz.O(:,:,3)))
% 
% 
% figure
% plot3(xyz.O(:,1:30,1), xyz.O(:,1,2), xyz.O(:,1:30,3))