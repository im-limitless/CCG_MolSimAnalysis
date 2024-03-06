clear all;
close all; clc;

%% Set the location of the calculation output
BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
system = 'CP_Pit_Water';

Cube = {'CP_Pit_Water-ELECTRON_DENSITY-1_10000.cube'; 'CP_Pit_Water-1_10000_Metal-ELECTRON_DENSITY-1_0.cube'; 'CP_Pit_Water-1_10000_Electrolyte-ELECTRON_DENSITY-1_0.cube'};

ABC = getABCvectors(BaseFldr, system);

for i=1:3

    [GridDens.(['C' num2str(i)]), RawDens, xGrid, yGrid, zGrid] = ReadGaussianCubeFile([BaseFldr system '\ChargeDensityDiff\'], Cube{i});

end

plot(linspace(0, ABC(3), 480), smooth(reshape(mean(GridDens.C1 - GridDens.C2-GridDens.C3, [1 ,2]),480,1)), 'Color', 'r', 'LineWidth', 1.5);