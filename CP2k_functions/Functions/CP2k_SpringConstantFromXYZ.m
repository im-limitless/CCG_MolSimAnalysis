clear all
close all

% split .xyz file and extract specific (final) geometry

fid  = fopen('G:\Imperial\Training\CP-LIKE\10Na10Cl\RemoveH_1ML\RemoveH_1ML_CPlike24\ptwater-pos-1.xyz');
% fid  = fopen('G:\Imperial\CP-LIKE\10Na10Cl\RemoveH_1ML\RemoveH_1ML_CPlike24\ptwater-pos-1.xyz');
% fid  = fopen('G:\Imperial\Federico\10Na10Cl\1ML\__dynamic16\ptwater-pos-1.xyz');
lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
nAtoms = str2num(lines{1});
relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
nConfigs = length(relevant);

% fidout = fopen('G:\Imperial\CP-LIKE\10Na10Cl\RemoveH_1ML\BOMD_NVT\ptwater-pos-final.xyz', 'w');

% Hsurfs is hard coded from output of CP2k_CheckReadxyz.m - it is the index
% of H atoms attached to the surface
Hsurfs =[         184
         186
         188
         190
         192
         194
         196
         198
         200
         204
         206
         208
         210
         212
         216
         218
         220
         222
         224
         226
         228
         230
         236
         238
         240
         242
         244
         246
         248
         250
         252
         254
         256
         258
         260
         262
         264
         266
         268
         270
         274
         276
         278
         280
         282
         286
         288
         290
         682
          76
          80
          82
          84
          86
          88
          90
          92
          94
         100
         102
         104
         106
         108
         110
         112
         114
         116
         118
         120
         122
         126
         128
         130
         132
         134
         136
         138
         140
         142
         146
         148
         150
         152
         156
         158
         160
         162
         164
         166
         168
         170
         172
         176
         178
         180
         182
     ];



 % for i = 1:10
 for i = 1:length(Hsurfs)
     for j = 1:nConfigs
         %         for j = 1:3200
         XYZ(i,j, 1:3) =  str2num(lines{relevant(j)+Hsurfs(i)}(5:end)); % add atom number relevant(i) e.g. relevant(i)+(AtomNumber) to get that atoms coordinates
         
     end
     x=0.5*(1:nConfigs);
     y=XYZ(i,:,3);
     p = polyfit(x,y,19);
     f = polyval(p,x);
     figure
     
     plot(x,y,'o',x,f,'-', 'linewidth', 2.0, 'color', 'k', 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'markersize', 2)
     xlabel('Time (fs)');
     ylabel('z-coordinate (Ang)');
     set(gcf, 'position', [675   203   560   219]);
     
     AveZPos(i) = mean(y);
     
 end

AveZUp = mean(AveZPos(1:48));
AveZdown = mean([AveZPos(50:53),AveZPos(55:59),AveZPos(61:end)]);
% fclose(fidout);
