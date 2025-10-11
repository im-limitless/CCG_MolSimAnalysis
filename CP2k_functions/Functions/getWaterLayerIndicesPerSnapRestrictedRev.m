function [FirstLayerIndx, SecondLayerIndx, MinimaZ] = getWaterLayerIndicesPerSnapRestrictedRev(Indx, XYZ, Dens_O, z, restX)

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

% % % %     if abs(MinimaZ{i}(end) - GlobalZ(end)) <= abs(MinimaZ{i}(end) - GlobalZ(end-1))
% % % %         MinimaZ{i}(end) = MinimaZ{i}(end);
% % % %         disp('Hello')
% % % %     elseif abs(MinimaZ{i}(end) - GlobalZ(end)) > abs(MinimaZ{i}(end) - GlobalZ(end-1))
% % %         MinimaZ{i}(end) = GlobalZ(end);
% % % %         disp('Bye')
% % % %     end

% if MinimaZ{i}(end) >= GlobalZ(end)
% elseif MinimaZ{i}(end) < GlobalZ(end)
%     MinimaZ{i}(end) = GlobalZ(end);
% end


    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) >= MinimaZ{i}(end) & (XYZ(i,:,1) < restX(1) | XYZ(i,:,1) > restX(2))))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ{i}(end) & XYZ(i,:,3) >= MinimaZ{i}(end-1) & (XYZ(i,:,1) < restX(1) | XYZ(i,:,1) > restX(2))))];
end

return