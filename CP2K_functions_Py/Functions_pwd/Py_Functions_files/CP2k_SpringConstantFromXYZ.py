import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
# split .xyz file and extract specific (final) geometry

fid = open('G:\Imperial\Training\CP-LIKE\10Na10Cl\RemoveH_1ML\RemoveH_1ML_CPlike24\ptwater-pos-1.xyz')
# fid  = fopen('G:\Imperial\CP-LIKE\10Na10Cl\RemoveH_1ML\RemoveH_1ML_CPlike24\ptwater-pos-1.xyz');
# fid  = fopen('G:\Imperial\Federico\10Na10Cl\1ML\__dynamic16\ptwater-pos-1.xyz');
lines = textscan(fid,'%s','delimiter','\n','whitespace','')
fclose(fid)
lines = lines[0]
nAtoms = str2num(lines[0])
relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
nConfigs = len(relevant)
# fidout = fopen('G:\Imperial\CP-LIKE\10Na10Cl\RemoveH_1ML\BOMD_NVT\ptwater-pos-final.xyz', 'w');

# Hsurfs is hard coded from output of CP2k_CheckReadxyz.m - it is the index
# of H atoms attached to the surface
Hsurfs = np.array([184,186,188,190,192,194,196,198,200,204,206,208,210,212,216,218,220,222,224,226,228,230,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,274,276,278,280,282,286,288,290,682,76,80,82,84,86,88,90,92,94,100,102,104,106,108,110,112,114,116,118,120,122,126,128,130,132,134,136,138,140,142,146,148,150,152,156,158,160,162,164,166,168,170,172,176,178,180,182])
# for i = 1:10
for i in np.arange(1,len(Hsurfs)+1).reshape(-1):
    for j in np.arange(1,nConfigs+1).reshape(-1):
        #         for j = 1:3200
        XYZ[i,j,np.arange[1,3+1]] = str2num(lines[relevant(j) + Hsurfs(i)](np.arange(5,end()+1)))
    x = 0.5 * (np.arange(1,nConfigs+1))
    y = XYZ(i,:,3)
    p = polyfit(x,y,19)
    f = polyval(p,x)
    figure
    plt.plot(x,y,'o',x,f,'-','linewidth',2.0,'color','k','markerfacecolor','r','markeredgecolor','r','markersize',2)
    plt.xlabel('Time (fs)')
    plt.ylabel('z-coordinate (Ang)')
    set(gcf,'position',np.array([675,203,560,219]))
    AveZPos[i] = mean(y)

AveZUp = mean(AveZPos(np.arange(1,48+1)))
AveZdown = mean(np.array([AveZPos(np.arange(50,53+1)),AveZPos(np.arange(55,59+1)),AveZPos(np.arange(61,end()+1))]))
fidout.close();