function [FirstLayerIndx, SecondLayerIndx, ThirdLayerIndx, FourthLayerIndx, MinimaZ] = getWaterLayerIndices(Indx, XYZ, Dens_O, z)

Minima = LocateStationaryPoints(mean(Dens_O,2));
% Minima = {[35 46 153 164]} ;disp('Using hard-coded WL override');

MinimaZ = z(Minima{:});


for i = 1:size(XYZ,1)
    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= z(Minima{1}(1)))); intersect(Indx.O,find(XYZ(i,:,3) >= z(Minima{1}(end))))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(1)) & XYZ(i,:,3) <= z(Minima{1}(2)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end)) & XYZ(i,:,3) >= z(Minima{1}(end-1))))];
    ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(2)) & XYZ(i,:,3) <= z(Minima{1}(3)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end-1)) & XYZ(i,:,3) >= z(Minima{1}(end-2))))];
    FourthLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(3)) & XYZ(i,:,3) <= z(Minima{1}(4)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end-2)) & XYZ(i,:,3) >= z(Minima{1}(end-3))))];
end

return