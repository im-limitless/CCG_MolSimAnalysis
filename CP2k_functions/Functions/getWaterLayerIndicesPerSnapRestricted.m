function [FirstLayerIndx, SecondLayerIndx, MinimaZ, ThirdLayerIndx] = getWaterLayerIndicesPerSnapRestricted(Indx, XYZ, Dens_O, z, restX)

GlobalMinima = LocateStationaryPoints(mean(Dens_O,2));
GlobalZ = z(find(GlobalMinima));
Minima = zeros(size(Dens_O,1), size(Dens_O,2));

for i = 1:size(Dens_O,2)
    Minima(:, i) = LocateStationaryPoints(Dens_O(:,i));
    MinimaIndx = find(Minima(:,i));
    MinimaZ{i} = z(MinimaIndx);

% % % % %     % recursion to move minima per snap to within tolerance of global min
% % % % %     if MinimaZ{i}(1) < GlobalZ(1)
% % % % %         n = 1;
% % % % %         while ~ismembertol(MinimaZ{i}(1), GlobalZ(1), 0.1)
% % % % %             MinimaZ{i}(1) = z(MinimaIndx(1)+n);
% % % % %             n = n + 1;
% % % % %         end
% % % % %     elseif MinimaZ{i}(1) > GlobalZ(1)
% % % % %         n = 1;
% % % % %         while ~ismembertol(MinimaZ{i}(1), GlobalZ(1), 0.03)
% % % % %             MinimaZ{i}(1) = z(MinimaIndx(1)-n);
% % % % %             n = n + 1;
% % % % %         end
% % % % %     else
% % % % % 
% % % % %     end
% 
% % % %     if abs(MinimaZ{i}(1) - GlobalZ(1)) <= abs(MinimaZ{i}(1) - GlobalZ(2))
% % %         MinimaZ{i}(1) = MinimaZ{i}(1);
% % % %         disp('hello')
% % % %     elseif abs(MinimaZ{i}(1) - GlobalZ(1)) > abs(MinimaZ{i}(1) - GlobalZ(2))
% % %         MinimaZ{i}(1) = GlobalZ(1);
% % % %         disp('bye')
% % % %     end

% if MinimaZ{i}(1) <= GlobalZ(1)
% elseif MinimaZ{i}(1) > GlobalZ(1)
%     MinimaZ{i}(1) = GlobalZ(1);
% end

    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= MinimaZ{i}(1) & (XYZ(i,:,1) < restX(1) | XYZ(i,:,1) > restX(2))))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ{i}(1) & XYZ(i,:,3) <= MinimaZ{i}(2) & (XYZ(i,:,1) < restX(1) | XYZ(i,:,1) > restX(2))))];
    ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ{i}(2) & XYZ(i,:,3) <= MinimaZ{i}(3) & (XYZ(i,:,1) < restX(1) | XYZ(i,:,1) > restX(2))))];

end



return