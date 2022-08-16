function [MinimaIndx] = LocateStationaryPoints(f)

dydx = zeros(size(f,1), size(f,2));
for i=1:size(f,2)
    %     dydx(:,i) = gradient(smooth(f(:,i))); % smoothed function misses
    %     sharp minima... 
    dydx(:,i) = gradient(f(:,i)); % non-smoothed function picks up noise
    
    count = [];
    for j = 1:length(dydx)-1
        if dydx(j) < 0 & dydx(j+1) > 0
            count = [count j];
        end
    end
    
    MinimaIndx{i} = count;
end

return

