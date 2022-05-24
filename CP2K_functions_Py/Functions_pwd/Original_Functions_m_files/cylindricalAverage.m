clear all
close all

ABC = [1 1 1];

M = [0 1 1;
    1 1 1;
    0 0 1;
    1 0 1;
    1 0 0;
    1 1 0;
    0 1 0;
    0.5 0.5 0.5;
    0 0 0];

plot3(M(:,1),M(:,2),M(:,3), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'k');
xlabel('x');
ylabel('y');
zlabel('z');

% map rectangular coordinates to cylindrical (account for points @ (0 0 z))
rho = [sqrt(M(:,1).^2+M(:,2).^2),  atan(M(:,2)./M(:,1)), M(:,3)];
rho(isnan(rho(:,2)),2) = 0;

% portion z into n "planes" with thickness dz - between 0 and 1 for now...

dz = 0.1;
zrange = 0:dz:ABC(3);

dtheta = 2*pi/180;

% search and find values within pla
for i = 1:length(zrange)
    for j = 
    
end