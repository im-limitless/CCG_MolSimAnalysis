import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
fldrname = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Vacuum\CleanSlab\CP_Pit_18H22F-1_1000_Metal\Epot\'
allMicros = dir(np.array([fldrname,'aver_z*.dat']))
allMacros = dir(np.array([fldrname,'avermacro_z*.dat']))
figure
hold('on')
for i in np.arange(1,len(allMicros)+1).reshape(-1):
    fid = open(np.array([fldrname,allMicros(i).name]))
    Micro[i] = np.transpose(fscanf(fid,'%f %f',np.array([2,inf])))
    fid.close()
    plt.plot(Micro[i](:,1),Micro[i](:,2) * 27.211386245988,'color','k')
    fid = open(np.array([fldrname,allMacros(i).name]))
    Macro[i] = np.transpose(fscanf(fid,'%f %f %f',np.array([3,inf])))
    fid.close()

AveMicro = mean(cat(3,Micro[:]),3)
AveMacro = mean(cat(3,Macro[:]),3)
set(gca,'xlim',np.array([0,np.amax(AveMicro(:,1))]))
set(gcf,'position',np.array([680,558,1062,420]))
plt.xlabel('z-coordinate (Ang)')
plt.ylabel('Electrostatic Potential (V)')
plt.plot(np.array([0,np.amax(AveMicro(:,1))]),np.array([0,0]),'--','color',np.array([0.6,0.6,0.6]))
# plot(AveMicro(:,1), AveMicro(:,2)*27.211386245988, 'color', 'b')
plt.plot(AveMacro(:,1),AveMacro(:,3) * 27.211386245988,'color','r','linewidth',1.5)
hold('off')