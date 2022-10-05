function [GridDens, RawDens, xGrid, yGrid, zGrid] = ReadGaussianCubeFile(BaseFldr, Cube)

pathSec =setOSpathSep;


fid = fopen([BaseFldr Cube]);
disp(['Reading cube file ' Cube '...']);
fgetl(fid);
FileType = fgetl(fid);
lines = strsplit(fgetl(fid));
nAtoms = str2num(lines{2});

lines = strsplit(fgetl(fid));
xGrid = str2num(lines{2});
lines = strsplit(fgetl(fid));
yGrid = str2num(lines{2});
lines = strsplit(fgetl(fid));
zGrid = str2num(lines{2});
fclose(fid);

Grid = [xGrid(1) yGrid(1) zGrid(1)];

RawDens = readmatrix([BaseFldr Cube], 'FileType', 'text', 'NumHeaderLines', nAtoms+6, 'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join');

RawDens = RawDens(:,2:end);

GridDens = reshape(RawDens', Grid);

return