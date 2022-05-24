import numpy as np
import matplotlib.pyplot as plt
import warnings
clear('all')
close_('all')
# parse the .ener file from cp2k

Basefldr = 'G:\Imperial\MattProjects\Pt_Clean\BOMD\BOMD_NVT_1010_Clean\CP_Like\Noisy\BOMD_NVT_1010_Clean\Noisy\'
Allfldrs = dir(np.array([Basefldr,'*_*']))
figure
hold('on')
set(gcf,'color','w')
plt.ylabel('Total Energy (eV)')
plt.xlabel('Time (ps)')
# set(gca, 'Ylim', [-3.2270e4 -3.2269e4])
# set(gca, 'Xlim', [0 1.5])
C = colormap(lines)
for l in np.arange(1,len(Allfldrs)+1).reshape(-1):
    if not exist(np.array([Basefldr,Allfldrs(l).name,'\out.log'])) :
        warnings.warn(np.array(['Output for ',Allfldrs(l).name,'does not exist!']))
        continue
    else:
        print(np.array(['Processing ',Allfldrs(l).name]))
    ParamScoreIndx = strfind(Allfldrs(l).name,'_')
    Stepsize[l] = str2num(Allfldrs(l).name(np.arange(ParamScoreIndx(1) + 1,ParamScoreIndx(2) - 1+1)))
    EXOR[l] = str2num(Allfldrs(l).name(np.arange(ParamScoreIndx(3) + 1,ParamScoreIndx(4) - 1+1)))
    Gamma[l] = str2num(Allfldrs(l).name(np.arange(ParamScoreIndx(5) + 1,ParamScoreIndx(6) - 1+1)))
    Noisy[l] = str2num(Allfldrs(l).name(np.arange(ParamScoreIndx(7) + 1,end()+1)))
    ConvIndx,NoConvIndx,Perc[l],AveSteps[l],StdSteps[l] = CP2k_parse_log(np.array([Basefldr,Allfldrs(l).name]),0)
    AveConverge[l],StdConverge[l] = CP2k_parse_SCF_Loops(np.array([Basefldr,Allfldrs(l).name]))
    flname = dir(np.array([Basefldr,Allfldrs(l).name,'\*.ener']))
    EnerFile = np.array([Basefldr,Allfldrs(l).name,'\',flname.name])
    fid = open(EnerFile)
    eofstat = False
    # First line is the number of atoms and is the same for all new sections
    HeadLine = fgetl(fid)
    eofstat = feof(fid)
    i = 0
    data = []
    while not eofstat :

        i = i + 1
        textLine = fgetl(fid)
        eofstat = feof(fid)
        data = np.array([[data],[str2num(textLine)]])

    fid.close()
    #         Cc(l) = randi([1 size(C,[1])]);
    Cc[l] = l
    h[l] = plt.plot(data(np.arange(1,100+1),2) / 1000,data(np.arange(1,100+1),6),'-','color',C(Cc(l),:),'markeredgecolor',C(Cc(l),:),'markerfacecolor',C(Cc(l),:),'markersize',3)
    #     text((data(100,2)/1000)+0.01, data(100,6), num2str(Step(l)), 'color', C(Cc(l),:))
    if len(Allfldrs) == 1:
        h[l] = plt.plot(data(np.arange(1,end()+1),2) / 1000,data(np.arange(1,end()+1),6),'-','color',C(Cc(l),:),'markeredgecolor',C(Cc(l),:),'markerfacecolor',C(Cc(l),:),'markersize',3,'linewidth',2)
        plt.legend(h,Allfldrs.name,'interpreter','none','Location','east')
        hold('off')
        figure
        set(gcf,'color','w')
        plt.ylabel('Temperature (K)')
        plt.xlabel('Time (ps)')
        hold('on')
        h2[l] = plt.plot(data(np.arange(1,end()+1),2) / 1000,data(np.arange(1,end()+1),4),'-','color',C(Cc(l),:),'markeredgecolor',C(Cc(l),:),'markerfacecolor',C(Cc(l),:),'markersize',3)
        hold('off')
        Step[l] = str2num(Allfldrs(l).name(np.arange(end() - 7,end()+1)))
        text((data(end(),2) / 1000) + 0.01,data(end() - 100,6),num2str(Step(l)),'color',C(Cc(l),:))
    else:
        Step[l] = str2num(Allfldrs(l).name(np.arange(end() - 7,end()+1)))
        # #     ## Energy
# #     h(l) = plot(data(1:end,2)/1000, data(1:end,6), '-', 'color', C(Cc(l),:), 'markeredgecolor', C(Cc(l),:), 'markerfacecolor', C(Cc(l),:), 'markersize', 3);
# #     text((data(end,2)/1000)+0.01, data(end,6), num2str(Step(l)), 'color', C(Cc(l),:));
#     Temp
        h[l] = plt.plot(data(np.arange(1,end()+1),2) / 1000,data(np.arange(1,end()+1),4),'-','color',C(Cc(l),:),'markeredgecolor',C(Cc(l),:),'markerfacecolor',C(Cc(l),:),'markersize',3,'linewidth',2)
        text((data(end(),2) / 1000) + 0.01,data(end(),4),num2str(Step(l)),'color',C(Cc(l),:))

# # # legend(h, Allfldrs.name, 'interpreter', 'none', 'Location', 'east')
# # # hold off

# [srtd, indx] = sort(Step);

# figure
# hold on
# set(gcf, 'color', 'w')
# ylabel('Mean SCF Steps')
# xlabel('STEPSIZE')
# plot(srtd, AveSteps(indx), 'o', 'color', 'k', 'markerfacecolor', 'b', 'markeredgecolor', 'k')
# # errorbar(srtd, AveSteps(indx), StdSteps(indx), 'o', 'color', 'k', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'color', 'k')
# hold off

figure
c = np.array(['r','g','b'])
for i in np.arange(1,np.amax(EXOR) + 1+1).reshape(-1):
    semilogy(Stepsize(EXOR == i - 1),AveConverge(EXOR == i - 1),'o')
    hold('on')
    plt.ylabel('OT Convergence')
    plt.xlabel('STEPSIZE')
    h[i] = errorbar(Stepsize(EXOR == i - 1),AveConverge(EXOR == i - 1),StdConverge(EXOR == i - 1),'o','markerfacecolor',c[i],'markeredgecolor','k','color','k')

N = np.array([np.arange(np.amin(EXOR),np.amax(EXOR)+1)])
legendCell = cellstr(num2str(np.transpose(N),'EXOR=%-d'))
plt.legend(h,legendCell)
hold('off')