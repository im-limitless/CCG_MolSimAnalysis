function [FirstLayerIndx, SecondLayerIndx, MinimaZ] = getWaterLayerIndicesPerSnap(Indx, XYZ, Dens_O, z)

GlobalMinima = LocateStationaryPoints(mean(Dens_O,2));
Minima = zeros(size(Dens_O,1), size(Dens_O,2));

for i = 1:size(Dens_O,2)
    Minima(:, i) = LocateStationaryPoints(Dens_O(:,i));
    MinimaZ = z(find(Minima(:,i)));
    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= MinimaZ(1))); intersect(Indx.O,find(XYZ(i,:,3) >= MinimaZ(end)))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
end

return