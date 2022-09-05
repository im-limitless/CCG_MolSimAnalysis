clear all;
close all; clc;

%% Set the location of the calculation output
BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
system = 'CP_Pit_Water';

AllFiles = dir([BaseFldr '*' system '\Cubes\*ELEC*.cube']);
Cube = {AllFiles.name};

SumGrid = [];

for i = 1:length(Cube)

    ABC = getABCvectors(BaseFldr, system);

    [GridDens.(['C' num2str(i)]), RawDens, xGrid, yGrid, zGrid] = ReadGaussianCubeFile([BaseFldr system '\Cubes\'], Cube{i});

end

SumGrid = GridDens.C1;
figure
hold on
xlabel('z (Ang)');
ylabel('Charge Density (e\cdotbohr^{-3})')
set(gca, 'xlim', [0 ABC(3)])
set(gcf, 'position', [-1390         270        1105         420])
for j = 1:length(Cube)-1
    SumGrid = SumGrid + GridDens.(['C' num2str(j+1)]);

plot(linspace(0, ABC(3), 480), smooth(reshape(mean(GridDens.(['C' num2str(j)]), [1 ,2]),480,1)), 'Color', 'r', 'LineWidth', 1.5);
end


AveGridDens = SumGrid/length(Cube);
plot(linspace(0, ABC(3), 480), smooth(reshape(mean(AveGridDens, [1 ,2]),480,1)), 'Color', 'k', 'LineWidth', 1.5);

