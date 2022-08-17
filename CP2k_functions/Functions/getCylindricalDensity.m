% function [Dens_O, Dens_H, Dens_F, TotDen, AveDen, z] = getDensityProfile(xyz, ABC)
% function [Dens_O, Dens_H, Dens_Na, Dens_Cl, TotDen, AveDen, z] = getDensityProfile(xyz, ABC)
function [Dens_O, Dens_H, TotDen, AveDen, z] = getCylindricalDensity(xyz, ABC)

%% map from cartesian into cylindrical space
theta = pi/2;

Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% translate origin to edge row
OriginCart = [7.18555517 12.50629757 7.72451052]; % Should modify to be avearage of edge row
OriginCartR = OriginCart*Rx;

xyz.O = reshape(reshape(xyz.O, [size(xyz.O,1)*size(xyz.O,2) size(xyz.O,3)])-OriginCartR, [size(xyz.O,1) size(xyz.O,2) size(xyz.O,3)]);
xyz.H = reshape(reshape(xyz.H, [size(xyz.H,1)*size(xyz.H,2) size(xyz.H,3)])-OriginCartR, [size(xyz.H,1) size(xyz.H,2) size(xyz.H,3)]);

[alpha.O,rho.O,zeta.O] = cart2pol(xyz.O(:,:,1),xyz.O(:,:,2),xyz.O(:,:,3));
[alpha.H,rho.H,zeta.H] = cart2pol(xyz.H(:,:,1),xyz.H(:,:,2),xyz.H(:,:,3));

Vec = diag(ABC);
VecR = Vec*Rx;
zmax = VecR(1,1);
rhoMax = max([rho.O rho.H], [], 'all');


% zmin;
bins = 200;
% zRange = linspace(0, zmax, bins);
rhoRange = linspace(0, rhoMax, bins);
Vol = zeros(length(rhoRange)-1,1);
angle1 = pi;
angle0 = -pi/3;

for i = 1:length(rhoRange)-1
    countO(i,:) = sum(rho.O < rhoRange(i+1) & rho.O > rhoRange(i), 'all'); 
    countH(i,:) = sum(rho.H < rhoRange(i+1) & rho.H > rhoRange(i), 'all'); 
%     countF(i,:) = sum(xyz.F(:,:,3) < zRange(i+1) & xyz.F(:,:,3) > zRange(i) ,2); 
%     countNa(i,:) = sum(xyz.Na(:,:,3) < zRange(i+1) & xyz.Na(:,:,3) > zRange(i) ,2); 
%     countCl(i,:) = sum(xyz.Cl(:,:,3) < zRange(i+1) & xyz.Cl(:,:,3) > zRange(i) ,2); 
    z(i) = (rhoRange(i+1)+rhoRange(i))/2;
    
%     Vol(i) = pi*((rhoRange(i+1)^2)-(rhoRange(i)^2))*Vec(2,2);
    
    Vol(i) = 0.5*Vec(2,2)*(angle1-angle0)*((rhoRange(i+1)^2)-(rhoRange(i)^2));
end

% zmax = ABC(3);
% bins = 200;
% zRange = linspace(0, zmax, bins);
% for i = 1:length(zRange)-1
%      countO(i,:) = sum(xyz.O(:,:,3) < zRange(i+1) & xyz.O(:,:,3) > zRange(i) ,2); 
%     countH(i,:) = sum(xyz.H(:,:,3) < zRange(i+1) & xyz.H(:,:,3) > zRange(i) ,2); 
% %     countF(i,:) = sum(xyz.F(:,:,3) < zRange(i+1) & xyz.F(:,:,3) > zRange(i) ,2); 
% %     countNa(i,:) = sum(xyz.Na(:,:,3) < zRange(i+1) & xyz.Na(:,:,3) > zRange(i) ,2); 
% %     countCl(i,:) = sum(xyz.Cl(:,:,3) < zRange(i+1) & xyz.Cl(:,:,3) > zRange(i) ,2); 
%     z(i) = (zRange(i+1)+zRange(i))/2;
% end


RMM_O = 15.99;
RMM_H = 1.008;
RMM_Na = 22.99;
RMM_Cl = 35.453;
RMM_F = 18.998;

Mass_O = countO*RMM_O/(6.0221409e+23*1000);
Mass_H = countH*RMM_H/(6.0221409e+23*1000);
% Mass_F = countF*RMM_F/(6.0221409e+23*1000);
% Mass_Na = countNa*RMM_Na/(6.0221409e+23*1000);
% Mass_Cl = countCl*RMM_Cl/(6.0221409e+23*1000);

Dens_O = Mass_O./((Vol)*(1e-30));
Dens_H = Mass_H./((Vol)*(1e-30));
% Dens_F = Mass_F/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
% Dens_Na = Mass_Na/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
% Dens_Cl = Mass_Cl/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));

TotDen = (Dens_O+Dens_H);
% TotDen = (Dens_O+Dens_H+Dens_F);
% TotDen = (Dens_O+Dens_H+Dens_Na+Dens_Cl);
AveDen = mean(TotDen((bins/2)-round(5/(zmax/(bins-1))/2):(bins/2)+round(5/(zmax/(bins-1))/2),:));
% AveDen = mean(TotDen((bins/2)-round(10/(zmax/(bins-1))/2):(bins/2)+round(10/(zmax/(bins-1))/2),:));

return