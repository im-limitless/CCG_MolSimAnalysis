clear all; 
close all; clc;

% BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
% fldrname = 'CP_Pit_18H22F';

BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\';
fldrname = 'CP_Like_0812_Fluoride';

% fid  = fopen([BaseFldr fldrname '\' fldrname '.xyz']);
% fid  = fopen([BaseFldr fldrname '\' fldrname '-pos-1.xyz']);
fid  = fopen([BaseFldr fldrname '\Sample0_5500.xyz']);
disp('Reading xyz data...');
lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
nAtoms = str2num(lines{1});
relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
HAtomNum =  find(~cellfun(@isempty,strfind({lines{1:nAtoms+2}}','H ')));
OAtomNum =  find(~cellfun(@isempty,strfind({lines{1:nAtoms+2}}','O ')));
NaAtomNum =  find(~cellfun(@isempty,strfind({lines{1:nAtoms+2}}','Na ')));
ClAtomNum =  find(~cellfun(@isempty,strfind({lines{1:nAtoms+2}}','Cl ')));
FAtomNum =  find(~cellfun(@isempty,strfind({lines{1:nAtoms+2}}','F ')));

startConfig = 1;
nConfigs = length(relevant);

% startConfig = 26;
% nConfigs = 2;
% 

XYZ = zeros(nConfigs, nAtoms, 3);
xyz_O2 = zeros(nConfigs, length(OAtomNum), 3);
xyz_H2 = zeros(nConfigs, length(HAtomNum), 3);
xyz_F2 = zeros(nConfigs, length(FAtomNum), 3);
xyz_O =[];
xyz_H =[];
xyz_Na =[];
xyz_Cl =[];
xyz_F =[];

for j = startConfig:nConfigs
    
    disp(['Processing configuration number ' num2str(j) ' of ' num2str(nConfigs)]);
    for k = 1:nAtoms
        %         for j = 1:3200
        XYZ(j,k, 1:3) =  str2num(lines{relevant(j)+k}(5:end)); % add atom number relevant(i) e.g. relevant(i)+(AtomNumber) to get that atoms coordinates
    end
%     xyz_O = [xyz_O; XYZ(OAtomNum-2,j,1:3)];
%     xyz_H = [xyz_H; XYZ(HAtomNum-2,j,1:3)];
%     xyz_Na = [xyz_Na; XYZ(NaAtomNum-2,j,1:3)];
%     xyz_Cl = [xyz_Cl; XYZ(ClAtomNum-2,j,1:3)];
%     xyz_F = [xyz_F; XYZ(FAtomNum-2,j,1:3)];
    
    xyz_O2(j,:, 1:3) = XYZ(j, OAtomNum-2,1:3);
    xyz_H2(j,:, 1:3) = XYZ(j, HAtomNum-2,1:3);
    xyz_F2(j,:, 1:3) = XYZ(j, FAtomNum-2,1:3);
end

% [~, ~, ABC, ~] = getLastPosAndVel(BaseFldr, fldrname);
ABC = getABCvectors(BaseFldr, fldrname);

%non general volume
zmax = ABC(3);

bins = 300;
% zRange = linspace(0, zmax, 1000);
zRange = linspace(0, zmax, bins);

for i = 1:length(zRange)-1
    
    IndxO = xyz_O(:,3) < zRange(i+1);
    IndxO = xyz_O(IndxO,3) > zRange(i);
    numO(i) = sum(IndxO);
    
    IndxH = xyz_H(:,3) < zRange(i+1);
    IndxH = xyz_H(IndxH,3) > zRange(i);
    numH(i) = sum(IndxH);
    
    IndxNa = xyz_Na(:,3) < zRange(i+1);
    IndxNa = xyz_Na(IndxNa,3) > zRange(i);
    numNa(i) = sum(IndxNa);
    
    IndxCl = xyz_Cl(:,3) < zRange(i+1);
    IndxCl = xyz_Cl(IndxCl,3) > zRange(i);
    numCl(i) = sum(IndxCl);
    
    IndxF = xyz_F(:,3) < zRange(i+1);
    IndxF = xyz_F(IndxF,3) > zRange(i);
    numF(i) = sum(IndxF);
    
    z(i) = (zRange(i+1)+zRange(i))/2;
end

RMM_O = 15.99;
RMM_H = 1.008;
RMM_Na = 22.99;
RMM_Cl = 35.453;
RMM_F = 18.998;

Mass_O = numO*RMM_O/(6.0221409e+23*1000);
Mass_H = numH*RMM_H/(6.0221409e+23*1000);
Mass_Na = numNa*RMM_Na/(6.0221409e+23*1000);
Mass_Cl = numCl*RMM_Cl/(6.0221409e+23*1000);
Mass_F = numF*RMM_F/(6.0221409e+23*1000);

%non general volume
Dens_O = Mass_O/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
Dens_H = Mass_H/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
Dens_Na = Mass_Na/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
Dens_Cl = Mass_Cl/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
Dens_F = Mass_F/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));

TotDen = (Dens_O+Dens_H+Dens_Na+Dens_Cl+Dens_F)/(nConfigs-startConfig+1);
AveDen = mean(TotDen((bins/2)-round(3/(zmax/(bins-1))/2):(bins/2)+round(3/(zmax/(bins-1))/2)));

% % sum(Mass_O+Mass_H+Mass_Na+Mass_Cl)/(ABC(1)*ABC(2)*(z(find((Dens_O+Dens_H+Dens_Na+Dens_Cl), 1, 'last'))-z(find((Dens_O+Dens_H+Dens_Na+Dens_Cl), 1, 'first')))*(1e-30));
% sum(Mass_O+Mass_H+Mass_Na+Mass_Cl)/(ABC(1)*ABC(2)*(zmax)*(1e-30))
% % AvailVolDen = sum(TotDen)*(z(find((Dens_O+Dens_H+Dens_Na+Dens_Cl), 1, 'last'))-z(find((Dens_O+Dens_H+Dens_Na+Dens_Cl), 1, 'first')))/zmax;

figure
hold on
set(gcf, 'position', [377         423        1123         420]);
xlabel('z (Ang)');
ylabel('Density (kgm^{-3})');
set(gca, 'xlim', [0 ABC(3)], 'ylim', [0 2000]);



% % % %% for HF w/out smoothing
% % % 
% % % plot(z, (Dens_O+Dens_H+Dens_F)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'k')
% % % plot(z, Dens_H/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', 'b')
% % % plot(z, Dens_O/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', 'r')
% % % plot(z, Dens_F/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', [34 177 76]/255)
% % % plot([z(1) z(end)], [AveDen AveDen], ':', 'color', [0.6 0.6 0.6])
% % % plot([z((bins/2)-round(3/(zmax/(bins-1))/2)) z((bins/2)-round(3/(zmax/(bins-1))/2))], [0 2000], ':', 'color', [0.6 0.6 0.6])
% % % plot([z((bins/2)+round(3/(zmax/(bins-1))/2)) z((bins/2)+round(3/(zmax/(bins-1))/2))], [0 2000], ':', 'color', [0.6 0.6 0.6])
% % % legend('Water+Ions', 'H', 'O', 'F', 'Ave. Bulk Density', 'location', 'north');


%% for HF w smoothing
plot(z, smooth((Dens_O+Dens_H+Dens_F)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', 'k')
plot(z, smooth(Dens_H/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', 'b')
plot(z, smooth(Dens_O/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', 'r')
plot(z, smooth(Dens_F/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', [34 177 76]/255)
plot([z(1) z(end)], [AveDen AveDen], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)-round(3/(zmax/(bins-1))/2)) z((bins/2)-round(3/(zmax/(bins-1))/2))], [0 2000], ':', 'color', [0.6 0.6 0.6])
plot([z((bins/2)+round(3/(zmax/(bins-1))/2)) z((bins/2)+round(3/(zmax/(bins-1))/2))], [0 2000], ':', 'color', [0.6 0.6 0.6])
legend('Water+Ions', 'H', 'O', 'F', 'Ave. Bulk Density', 'location', 'northeast');

% % % %% for NaCl w/out smoothing
% % % plot(z, (Dens_O+Dens_H)/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', 'k')
% % % plot(z, (Dens_O+Dens_H+Dens_Na+Dens_Cl)/(nConfigs-startConfig+1), 'linewidth', 1.5, 'color', 'k')
% % % plot(z, Dens_H/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', 'b')
% % % plot(z, Dens_O/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', 'r')
% % % plot(z, Dens_Na/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', [128 0 255]/255)
% % % plot(z, Dens_Cl/(nConfigs-startConfig+1), ':', 'linewidth', 1.5, 'color', [34 177 76]/255)
% % % plot([z(1) z(end)], [AveDen AveDen], ':', 'color', [0.6 0.6 0.6])
% % % plot([ABC(3)/2 ABC(3)/2], [0 3000], ':', 'color', [0.6 0.6 0.6])
% % % legend('Water', 'Water+Ions', 'H', 'O', 'Na', 'Cl', 'Ave. Bulk Density', 'location', 'north');
% % % 
% % % 
% % % % % % %% for NaCl w smoothing
% % % plot(z, smooth((Dens_O+Dens_H)/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', 'k')
% % % plot(z, smooth((Dens_O+Dens_H+Dens_Na+Dens_Cl)/(nConfigs-startConfig+1)), 'linewidth', 1.5, 'color', 'k')
% % % plot(z, smooth(Dens_H/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', 'b')
% % % plot(z, smooth(Dens_O/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', 'r')
% % % plot(z, smooth(Dens_Na/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', [128 0 255]/255)
% % % plot(z, smooth(Dens_Cl/(nConfigs-startConfig+1)), ':', 'linewidth', 1.5, 'color', [34 177 76]/255)
% % % plot([z(1) z(end)], [AveDen AveDen], ':', 'color', [0.6 0.6 0.6])
% % % plot([ABC(3)/2 ABC(3)/2], [0 3000], ':', 'color', [0.6 0.6 0.6])
% % % legend('Water', 'Water+Ions', 'H', 'O', 'Na', 'Cl', 'Ave. Bulk Density', 'location', 'north');



hold off
