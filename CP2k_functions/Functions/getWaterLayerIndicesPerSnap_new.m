function [FirstLayerIndx, SecondLayerIndx, ThirdLayerIndx, MinimaZ] = getWaterLayerIndicesPerSnap_new(Indx, XYZ, Dens_O, z)

GlobalMinima = LocateStationaryPoints(mean(Dens_O,2));
Minima = zeros(size(Dens_O,1), size(Dens_O,2));

for i = 1:size(Dens_O,2)
    Minima(:, i) = LocateStationaryPoints(Dens_O(:,i));
    MinimaZ = z(find(Minima(:,i)));
    FirstLayerIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= MinimaZ(1))); intersect(Indx.O,find(XYZ(i,:,3) >= MinimaZ(end)))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
    ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2)))]; %/This is added for the cases when we have Oxide on the surface then this becomes the second water layer !!

     %% Loop to solve the small peaks issue with second and third layers !
    % if (MinimaZ(2)-MinimaZ(1)) >=1.5
    %     if  (MinimaZ(end) -  MinimaZ(end-1)) >=1.5
    %         SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
    %     else
    %         SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1-1)))];
    %     end
    % else
    %     if  (MinimaZ(end) -  MinimaZ(end-1)) >=1.5
    %         SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2+1))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
    %     else
    %         SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2+1))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1-1)))];
    %     end
    % end
    % 
    % if (MinimaZ(3)-MinimaZ(2)) >=1.5
    %     if  (MinimaZ(end-1) -  MinimaZ(end-2)) >=1.5
    %         ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2)))]; %/This is added for the cases when we have Oxide on the surface then this becomes the second water layer !!
    %     else
    %         ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2-1)))];
    %     end
    % else
    %     if  (MinimaZ(end-1) -  MinimaZ(end-2)) >=1.5
    %         ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3+1))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2)))];
    %     else
    %         ThirdLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3+1))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2-1)))];
    %     end
    % end
        %% end of Loop to solve the small peaks issue with second and third layers !

end

return