function [Vec12, Dist12] = GetAtomCorrelation(XYZ, Indx1, Indx2, ABC)

Dist12 = zeros(length(Indx2), length(Indx1));
Vec12 = zeros(length(Indx2), 3, length(Indx1));

for j = 1:length(Indx1)
    
    % scan through each atom 1 and find all atom1-atom2 vectors
    Vec12(:,:,j) = XYZ(Indx2,:) - XYZ(Indx1(j),:);
    % find all atom1-atom2 vector magnitudes
    Dist12(:,j) = sqrt((Vec12(:,1,j).^2)+(Vec12(:,2,j).^2)+(Vec12(:,3,j).^2));
    
    % check for neighbours separated by PBC in x and y only
    [Vec12(:,:,j), Dist12(:,j)] = searchAcrossPBC(Vec12(:,:,j), Dist12(:,j), XYZ, Indx1(j), Indx2, ABC); % searchAcrossPBC(PBC_mat, non PBC Vector, non PBC distance, XYZ, looping atom IndxAtom(j), other atom IndxAtom, ABC)
        
end



return
