clear all; 
close all; clc;

fldrname = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Vacuum\CleanSlab\CP_Pit_18H22F-1_1000_Metal\Epot\';

allMicros = dir([fldrname 'aver_z*.dat']);
allMacros = dir([fldrname 'avermacro_z*.dat']);

figure
hold on

for i = 1:length(allMicros)

fid = fopen([fldrname allMicros(i).name]);
Micro{i} = fscanf(fid, '%f %f', [2 inf])';
fclose(fid);

plot(Micro{i}(:,1), Micro{i}(:,2)*27.211386245988, 'color', 'k')

fid = fopen([fldrname allMacros(i).name]);
Macro{i} = fscanf(fid, '%f %f %f', [3 inf])';
fclose(fid);
end

AveMicro = mean(cat(3,Micro{:}),3);
AveMacro = mean(cat(3,Macro{:}),3);

set(gca, 'xlim', [0 max(AveMicro(:,1))]);
set(gcf, 'position', [680   558  1062         420]);
xlabel('z-coordinate (Ang)');
ylabel('Electrostatic Potential (V)');

plot([0 max(AveMicro(:,1))], [0 0], '--', 'color', [0.6 0.6 0.6])
% plot(AveMicro(:,1), AveMicro(:,2)*27.211386245988, 'color', 'b')
plot(AveMacro(:,1), AveMacro(:,3)*27.211386245988, 'color', 'r', 'linewidth', 1.5)
hold off
