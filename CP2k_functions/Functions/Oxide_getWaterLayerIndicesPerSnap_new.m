function [FirstLayerIndx, SecondLayerIndx, OxideIndx, MinimaZ] = Oxide_getWaterLayerIndicesPerSnap_new(Indx, XYZ, Dens_O, z, ABC)

GlobalMinima = LocateStationaryPoints(mean(Dens_O,2));
Minima = zeros(size(Dens_O,1), size(Dens_O,2));

temp_Al_z_ULim=[]; %collect Al_z_Ulim over time stamps 

for i = 1:size(Dens_O,2)
    
    %% This section is to find the extent of Al1 at both interfaces %%
    Al_z_ULim_Indices = intersect(Indx.Al_All,find(XYZ(i,:,3)<=(ABC(3)/2)));
    Al_z_ULim = max(XYZ(i,Al_z_ULim_Indices,3));

    temp_Al_z_ULim=[temp_Al_z_ULim; Al_z_ULim];

    Al_z_LLim_Indices = intersect(Indx.Al_All,find(XYZ(i,:,3)>=(ABC(3)/2)));
    Al_z_LLim = min(XYZ(i,Al_z_LLim_Indices,3));
    %% END %%

    %% Set the minimaz according to the Al1 limits found above %%
    Minima(:, i) = LocateStationaryPoints(sgolayfilt(Dens_O(:,i),2,3)); %using Savitzky-Golay filtering
    % Minima(:, i) = LocateStationaryPoints(Dens_O(:,i)); %no filtering
    MinimaZ = z(find(Minima(:,i)));
    MinimaZ = MinimaZ(find(MinimaZ>Al_z_ULim & MinimaZ<Al_z_LLim));
    %% END %%

    OxideIndx{i} = [intersect(Indx.O,find(XYZ(i,:,3) <= MinimaZ(1))); intersect(Indx.O,find(XYZ(i,:,3) >= MinimaZ(end)))]; %\This includes both the Oxides and the 1WL
    FirstLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
    SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2)))]; %/This is added for the cases when we have Oxide on the surface then this becomes the second water layer !!
end 
   
return


% FirstLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(1) & XYZ(i,:,3) <= MinimaZ(2))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end) & XYZ(i,:,3) >= MinimaZ(end-1)))];
    % SecondLayerIndx{i} = [intersect(Indx.O, find(XYZ(i,:,3) > MinimaZ(2) & XYZ(i,:,3) <= MinimaZ(3))); intersect(Indx.O, find(XYZ(i,:,3) < MinimaZ(end-1) & XYZ(i,:,3) >= MinimaZ(end-2)))];