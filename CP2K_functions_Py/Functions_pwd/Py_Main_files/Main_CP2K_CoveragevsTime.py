import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
BaseFldr = 'Z:\Imperial\MattProjects\OxidesOER\RutheniumOxide\'
system = 'RuO2_u1_1ML_u2_1ML'
Trajectory = 'RuO2_u1_1ML_u2_1ML_0to15000_500step.xyz'
ABC = getABCvectors(BaseFldr,system)
xyz,XYZ,Indx,__,__,nAtoms,startConfig,nConfigs,StepNum = ReadAndParsexyz(BaseFldr,system,Trajectory,ABC)
Atoms,AtomList,AtomIndx = getAtomNamesFromInputXYZ(BaseFldr,system)
Coverage = np.zeros((nConfigs,2))
for snap in np.arange(startConfig,nConfigs+1).reshape(-1):
    print(np.array(['Processing snapshot ',num2str(StepNum(snap)),' - ',num2str(100 * (snap / nConfigs)),' % complete']))
    XYZ_snap = np.zeros((XYZ.shape[2-1],XYZ.shape[3-1]))
    XYZ_snap[:,:] = XYZ(snap,:,:)
    #     [VecOtUH, DistOtUH] = GetAtomCorrelation(XYZ_snap, AtomIndx.OtU, AtomIndx.H, ABC);
#     [VecOtLH, DistOtLH] = GetAtomCorrelation(XYZ_snap, AtomIndx.OtL, AtomIndx.H, ABC);
    VecOtUH,DistOtUH = GetAtomCorrelation(XYZ_snap,AtomIndx.OtU,np.array([[AtomIndx.H],[AtomIndx.Hsurf]]),ABC)
    VecOtLH,DistOtLH = GetAtomCorrelation(XYZ_snap,AtomIndx.OtL,np.array([[AtomIndx.H],[AtomIndx.Hsurf]]),ABC)
    Coverage[snap,1] = np.sum(DistOtUH < 1.25, 'all'-1)
    Coverage[snap,2] = np.sum(DistOtLH < 1.25, 'all'-1)

figure
box('on')
hold('on')
plt.xlabel('Time (ps)')
plt.ylabel('\theta_H (ML)')
plt.plot(StepNum / 2000,Coverage(:,1) / len(AtomIndx.OtU),'-o','color','k','markeredgecolor','k','markerfacecolor','b')
plt.plot(StepNum / 2000,Coverage(:,2) / len(AtomIndx.OtL),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
plt.plot(StepNum / 2000,(Coverage(:,1) + Coverage(:,2)) / (len(AtomIndx.OtL) + len(AtomIndx.OtU)),'-o','color','k','markeredgecolor','k','markerfacecolor',np.array([0,0.5,0]))
set(gca,'YTick',np.arange(0,2+0.1,0.1),'fontsize',12)
plt.legend('\mu_1','\mu_2','Total','interpreter','tex','location','best')
hold('off')