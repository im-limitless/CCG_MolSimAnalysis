function [Vec, Dist] = searchAcrossPBC(Vec, Dist, XYZ, IndxA, IndxB, ABC)

PBC_mat = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 1 1 0; -1 1 0; 1 -1 0; -1 -1 0; 1 0 1; -1 0 1; 0 1 1; 0 -1 1; 1 1 1; -1 1 1; 1 -1 1; -1 -1 1; 1 0 -1; -1 0 -1; 0 1 -1; 0 -1 -1; 1 1 -1; -1 1 -1; 1 -1 -1; -1 -1 -1; 0 0 1; 0 0 -1;];


        % check for neighbours separated by PBC in x and y only - would be
        % good to make this a function, the code is a bit busy here
        for jj = 1:length(PBC_mat)
            Vec_PBC = XYZ(IndxB,:) - (XYZ(IndxA,:) + PBC_mat(jj,:).*ABC);
            Dist_PBC = sqrt((Vec_PBC(:,1).^2)+(Vec_PBC(:,2).^2)+(Vec_PBC(:,3).^2));
            
            % determine the shortest vector between O atom IndxO(j) and every H after each PBC element is applied
            if any(Dist_PBC<Dist)
                PBCIndx = find(Dist_PBC<Dist);
                Vec(PBCIndx,:) = Vec_PBC(PBCIndx,:);
            end
            
            % determine the shortest distance between O atom IndxO(j) and every H after each PBC element is applied
            Dist = min(Dist_PBC, Dist);
        end
        
end