# check water density
import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
mWater = 18 / 6.02e+23 / 1000
mF = 19 / 6.02e+23 / 1000
mH = 1 / 6.02e+23 / 1000
mO = 16 / 6.02e+23 / 1000
a = 14.5683
b = 25.233
# #0812
# # PtU = 35.76466286; #Original
# PtU = 34.90746286;
# PtL = 4.55143801;
# nF = 12;
# nH = 8;
# c = 39.4800;
# # 40.3372; #Original
# nO = 0;

#1210
# PtU = 34.62549736;
# PtL = 4.43248015;
PtU = 34.63455446
PtL = 4.43248015
nF = 10
nH = 12
# c = 39.2788;
c = 39.2879
nO = 0
# #1010
# # PtU = 34.62549736; #Original
# # PtL = 4.53588015; #Original
# PtU = 34.63455446;
# PtL = 4.43248015;
# nF = 10;
# nH = 10;
# # c = 39.1754; #Original
# c = 39.2879;
# nO = 0;

# #1012
# # PtU = 35.80412939;
# # PtL = 4.54651915;
# PtU = 34.91159939;
# PtL = 4.54651915;
# nF = 12;
# nH = 10;
# # c = 40.3372;
# c = 39.4447;
# nO = 0;

# #1010 + O 1/3 ML
# # PtU = 35.97901096; #O
# # PtL = 5.78674566; #O
# PtU = 37.19098332;
# PtL = 4.58453624;
# nF = 10;
# nH = 10;
# c = 41.7754;
# nO = 108/3;

nPt = 108
vPt = nPt / 2 * ((4 / 3) * pi * (2.804 / 2) ** 3)
vO = nO / 2 * ((4 / 3) * pi * 1.52 ** 3)
vSphere = vPt + vO
# V=a*b*((mean(PtU) - mean(PtL))) - vPt;
V = a * b * (np.array([np.arange(25,35+0.001,0.001)])) - vSphere - vO
# V=a*b*([25:0.001:35]) - vPt;
V = V * 1e-30
# x = 300:400;
# y=(mWater*x+mF*nF+mH*nH)/V;
# plot(x,y);
x = 329
# y=(mWater*x+mF*nF+mH*nH)./V;
# Vsolv = (mWater*x+mF*nF+mH*nH)/998;

y = (mWater * x + mF * nF + mH * nH + mO * nO) / V
Vsolv = (mWater * x + mF * nF + mH * nH + mO * nO) / 998
figure
hold('on')
plt.plot(((V * 1e+30 + vSphere) / (a * b)) + c - (PtU - PtL),y,'k')
l1 = get(gca,'ylim')
l2 = get(gca,'xlim')
plt.plot(np.array([((Vsolv * 1e+30 + vSphere) / (a * b)) + c - (PtU - PtL),((Vsolv * 1e+30 + vSphere) / (a * b)) + c - (PtU - PtL)]),np.array([0,998]),'--','color','r')
plt.plot(np.array([0((Vsolv * 1e+30 + vSphere) / (a * b)) + c - (PtU - PtL)]),np.array([998,998]),'--','color','r')
set(gca,'ylim',l1,'xlim',l2)
plt.xlabel('z length (Ang)')
plt.ylabel('Density (kg\cdotm^{-3})')
print(np.array(['New z length = ',num2str(((Vsolv * 1e+30 + vSphere) / (a * b)) + c - (PtU - PtL)),' changed by ',num2str(((Vsolv * 1e+30 + vSphere) / (a * b)) + c - (PtU - PtL) - c)]))