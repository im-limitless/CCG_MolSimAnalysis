function [FirstLayerIndx_H, SecondLayerIndx_H, MinimaZ_H] = H_getWaterLayerIndicesPerSnap_new(Indx, XYZ, Dens_H, z)

GlobalMinima = LocateStationaryPoints(mean(Dens_H,2));
Minima = zeros(size(Dens_H,1), size(Dens_H,2));

for i = 1:size(Dens_H,2)
    Minima(:, i) = LocateStationaryPoints(Dens_H(:,i));
    MinimaZ = z(find(Minima(:,i)));
    FirstLayerIndx_H{i} = [intersect(Indx.H,find(XYZ(i,:,3) <= MinimaZ(1))); intersect(Indx.H,find(XYZ(i,:,3) >= MinimaZ(end)))];
    SecondLayerIndx_H{i} = [intersect(Indx.H, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.H, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
end

return