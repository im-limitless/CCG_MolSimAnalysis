import numpy as np
    
def getDensityProfile(xyz = None,ABC = None): 
    # function [Dens_O, Dens_H, Dens_Na, Dens_Cl, TotDen, AveDen, z] = getDensityProfile(xyz, ABC)
# function [Dens_O, Dens_H, TotDen, AveDen, z] = getDensityProfile(xyz, ABC)
    
    zmax = ABC(3)
    bins = 300
    zRange = np.linspace(0,zmax,bins)
    for i in np.arange(1,len(zRange) - 1+1).reshape(-1):
        countO[i,:] = np.sum(xyz.O(:,:,3) < np.logical_and(zRange(i + 1),xyz.O(:,:,3)) > zRange(i), 2-1)
        countH[i,:] = np.sum(xyz.H(:,:,3) < np.logical_and(zRange(i + 1),xyz.H(:,:,3)) > zRange(i), 2-1)
        countF[i,:] = np.sum(xyz.F(:,:,3) < np.logical_and(zRange(i + 1),xyz.F(:,:,3)) > zRange(i), 2-1)
        #     countNa(i,:) = sum(xyz.Na(:,:,3) < zRange(i+1) & xyz.Na(:,:,3) > zRange(i) ,2);
#     countCl(i,:) = sum(xyz.Cl(:,:,3) < zRange(i+1) & xyz.Cl(:,:,3) > zRange(i) ,2);
        z[i] = (zRange(i + 1) + zRange(i)) / 2
    
    RMM_O = 15.99
    RMM_H = 1.008
    RMM_Na = 22.99
    RMM_Cl = 35.453
    RMM_F = 18.998
    Mass_O = countO * RMM_O / (6.0221409e+23 * 1000)
    Mass_H = countH * RMM_H / (6.0221409e+23 * 1000)
    Mass_F = countF * RMM_F / (6.0221409e+23 * 1000)
    # Mass_Na = countNa*RMM_Na/(6.0221409e+23*1000);
# Mass_Cl = countCl*RMM_Cl/(6.0221409e+23*1000);
    
    Dens_O = Mass_O / (ABC(1) * ABC(2) * (ABC(3) / (bins - 1)) * (1e-30))
    Dens_H = Mass_H / (ABC(1) * ABC(2) * (ABC(3) / (bins - 1)) * (1e-30))
    Dens_F = Mass_F / (ABC(1) * ABC(2) * (ABC(3) / (bins - 1)) * (1e-30))
    # Dens_Na = Mass_Na/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
# Dens_Cl = Mass_Cl/(ABC(1)*ABC(2)*(ABC(3)/(bins-1))*(1e-30));
    
    # TotDen = (Dens_O+Dens_H);
    TotDen = (Dens_O + Dens_H + Dens_F)
    # TotDen = (Dens_O+Dens_H+Dens_Na+Dens_Cl);
    AveDen = mean(TotDen(np.arange((bins / 2) - np.round(5 / (zmax / (bins - 1)) / 2),(bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2)+1),:))
    # AveDen = mean(TotDen((bins/2)-round(10/(zmax/(bins-1))/2):(bins/2)+round(10/(zmax/(bins-1))/2),:));
    
    return Dens_O,Dens_H,Dens_F,TotDen,AveDen,z
    return Dens_O,Dens_H,Dens_F,TotDen,AveDen,z