function [FirstLayerIndx, SecondLayerIndx, ThirdLayerIndx, FourthLayerIndx, MinimaZ] = getWaterLayerIndicesPerSnap(Indx, XYZ, Dens_O, z)

Minima = zeros(size(Dens_O,2), 2);




for j = 1:size(Dens_O,2)
    Minima{j} = LocateStationaryPoints(Dens_O(:,j));
end
MinimaZ = z(Minima{:});

for i = 1:size(XYZ,1)
    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= z(Minima{1}(1)))); intersect(Indx.O,find(XYZ(i,:,3) >= z(Minima{1}(end))))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(1)) & XYZ(i,:,3) <= z(Minima{1}(2)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end)) & XYZ(i,:,3) >= z(Minima{1}(end-1))))];
    ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(2)) & XYZ(i,:,3) <= z(Minima{1}(3)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end-1)) & XYZ(i,:,3) >= z(Minima{1}(end-2))))];
    FourthLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > z(Minima{1}(3)) & XYZ(i,:,3) <= z(Minima{1}(4)))); intersect(Indx.O, find(XYZ(i,:,3) < z(Minima{1}(end-2)) & XYZ(i,:,3) >= z(Minima{1}(end-3))))];
end

return