function [FirstLayerIndx, SecondLayerIndx] = getWaterLayerIndices(Indx, XYZ, Dens_O, z)

Minima = LocateStationaryPoints(mean(Dens_O,2));

for i = 1:size(XYZ,1)
    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= z(Minima{1}(1)))); intersect(Indx.O,find(XYZ(i,:,3) >= z(Minima{1}(end))))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(1)) & XYZ(i,:,3) < z(Minima{1}(2)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end)) & XYZ(i,:,3) > z(Minima{1}(end-1))))];
end

return