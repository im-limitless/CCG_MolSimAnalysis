import numpy as np
import warnings
import matplotlib.pyplot as plt
    
def CP2K_CalcEffectivePotentialDrop(BaseFldr = None,System = None): 
    # BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
# # System = 'CP_Pit_18H22F';
# # System = 'CP_Pit_20F';
# System = 'CP_Pit_20H22F';
    
    # BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\';
# System = 'CP_Like_1010_Fluorine';
    
    PotenFldr = np.array([BaseFldr,System,'\EffectivePD\'])
    PotenSnaps = dir(np.array([PotenFldr,System,'-1*']))
    # fldrname = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\CP_Pit_18H22F\ElectroPots\';
    
    # allMicros = dir([fldrname 'aver_z*.dat']);
# allMacros = dir([fldrname 'avermacro_z*.dat']);
    
    for i in np.arange(1,len(PotenSnaps)+1).reshape(-1):
        StepIndx = find(PotenSnaps(i).name == '_')
        StepNum[i] = str2num(PotenSnaps(i).name(np.arange(StepIndx(end()) + 1,end()+1)))
        if np.logical_and(np.logical_and(exist(np.array([PotenFldr,PotenSnaps(i).name,'\',PotenSnaps(i).name,'_Metal\Epot\aver_Z_0.dat'])),exist(np.array([PotenFldr,PotenSnaps(i).name,'\',PotenSnaps(i).name,'_Electrolyte\Epot\aver_Z_0.dat']))),exist(np.array([PotenFldr,'aver_Z_',PotenSnaps(i).name(np.arange(StepIndx(end()) + 1,end()+1)),'.dat']))):
            print(np.array(['Processing ',PotenSnaps(i).name,'...']))
            fid = open(np.array([PotenFldr,PotenSnaps(i).name,'\',PotenSnaps(i).name,'_Metal\Epot\aver_Z_0.dat']))
            MicroMetal[i] = np.transpose(fscanf(fid,'%f %f',np.array([2,inf])))
            fid.close()
            fid = open(np.array([PotenFldr,PotenSnaps(i).name,'\',PotenSnaps(i).name,'_Electrolyte\Epot\aver_Z_0.dat']))
            MicroElectrolyte[i] = np.transpose(fscanf(fid,'%f %f',np.array([2,inf])))
            fid.close()
            fid = open(np.array([PotenFldr,'aver_Z_',PotenSnaps(i).name(np.arange(StepIndx(end()) + 1,end()+1)),'.dat']))
            Micro[i] = np.transpose(fscanf(fid,'%f %f',np.array([2,inf])))
            fid.close()
            MicroEffective[i] = Micro[i](:,2) - MicroMetal[i](:,2) - MicroElectrolyte[i](:,2)
            bins = len(Micro[i](:,1))
            zmax = np.amax(Micro[i](:,1))
            BulkAveMicro = mean(MicroEffective[i](np.arange((bins / 2) - np.round(3 / (zmax / (bins - 1)) / 2),(bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2)+1))) * 27.211386245988
            EffectiveMicroPot[i] = BulkAveMicro - mean(np.array([[MicroEffective[i](np.arange(1,20+1))],[MicroEffective[i](np.arange(end() - 19,end()+1))]])) * 27.211386245988
            #         figure
#         hold on
#         set(gca, 'xlim', [0 max(Micro{i}(:,1))]);
#         set(gcf, 'position', [680   558  1062         420]);
#         xlabel('z-coordinate (Ang)');
#         ylabel('Electrostatic Potential (V)');
#         plot([0 max(Micro{i}(:,1))], [0 0], '--', 'color', [0.6 0.6 0.6])
#         plot(Micro{i}(:,1), Micro{i}(:,2)*27.211386245988, 'color', 'b')
#         plot(MicroElectrolyte{i}(:,1), MicroElectrolyte{i}(:,2)*27.211386245988, 'color', 'r')
#         plot(MicroMetal{i}(:,1), MicroMetal{i}(:,2)*27.211386245988, 'color', 'k')
#         hold off
            #         figure
#         hold on
#         set(gca, 'xlim', [0 max(Micro{i}(:,1))]);
#         set(gcf, 'position', [680   558  1062         420]);
#         xlabel('z-coordinate (Ang)');
#         ylabel('Electrostatic Potential (V)');
#         plot(Micro{i}(:,1), (MicroEffective{i}(:,1)-mean([MicroEffective{i}(1:20); MicroEffective{i}(end-19:end)]))*27.211386245988, 'color', [0.7 0.7 0.7])
#         hold off
            fid = open(np.array([PotenFldr,PotenSnaps(i).name,'\',PotenSnaps(i).name,'_Metal\Epot\avermacro_Z_0.dat']))
            MacroMetal[i] = np.transpose(fscanf(fid,'%f %f',np.array([3,inf])))
            fid.close()
            fid = open(np.array([PotenFldr,PotenSnaps(i).name,'\',PotenSnaps(i).name,'_Electrolyte\Epot\avermacro_Z_0.dat']))
            MacroElectrolyte[i] = np.transpose(fscanf(fid,'%f %f',np.array([3,inf])))
            fid.close()
            fid = open(np.array([PotenFldr,'avermacro_Z_',PotenSnaps(i).name(np.arange(StepIndx(end()) + 1,end()+1)),'.dat']))
            Macro[i] = np.transpose(fscanf(fid,'%f %f',np.array([3,inf])))
            fid.close()
            MacroEffective[i] = Macro[i](:,3) - MacroMetal[i](:,3) - MacroElectrolyte[i](:,3)
            # BulkAve = mean(MacroEffective{i}((length(Macro{i})/2)-round(3/(max(Macro{i}(:,1)/(length(Macro{i})/2-1))/2):(length(Macro{i})/2/2)+round(3/(max(Macro{i}(:,1))/(length(Macro{i})/2-1))/2)));
            bins = len(Macro[i](:,1))
            zmax = np.amax(Macro[i](:,1))
            BulkAve = mean(MacroEffective[i](np.arange((bins / 2) - np.round(3 / (zmax / (bins - 1)) / 2),(bins / 2) + np.round(3 / (zmax / (bins - 1)) / 2)+1))) * 27.211386245988
            EffectiveMacroPot[i] = BulkAve - mean(np.array([[MacroEffective[i](np.arange(1,20+1))],[MacroEffective[i](np.arange(end() - 19,end()+1))]])) * 27.211386245988
            # figure
# hold on
# set(gca, 'xlim', [0 max(Macro{i}(:,1))]);
# set(gcf, 'position', [680   558  1062         420]);
# xlabel('z-coordinate (Ang)');
# ylabel('Electrostatic Potential (V)');
# plot([0 max(Macro{i}(:,1))], [0 0], '--', 'color', [0.6 0.6 0.6])
# plot(Macro{i}(:,1), Macro{i}(:,3)*27.211386245988, 'color', 'b')
# plot(MacroElectrolyte{i}(:,1), MacroElectrolyte{i}(:,3)*27.211386245988, 'color', 'r')
# plot(MacroMetal{i}(:,1), MacroMetal{i}(:,3)*27.211386245988, 'color', 'k')
# # hold off
# #
#         figure
#         hold on
#         set(gca, 'xlim', [0 max(Macro{i}(:,1))]);
#         set(gcf, 'position', [680   558  1062         420]);
#         xlabel('z-coordinate (Ang)');
#         ylabel('Electrostatic Potential (V)');
#         plot(Macro{i}(:,1), (MacroEffective{i}(:,1)-mean([MacroEffective{i}(1:20); MacroEffective{i}(end-19:end)]))*27.211386245988, 'color', [0.7 0.7 0.7])
#         hold off
        else:
            EffectiveMacroPot[i] = 0
            EffectiveMicroPot[i] = 0
            warnings.warn(np.array(['Electrostatic potential for ',System,' is missing, continuing with next snapshot...']))
            continue
    
    StepNum = StepNum(EffectiveMacroPot != 0)
    EffectiveMacroPot = EffectiveMacroPot(EffectiveMacroPot != 0)
    EffectiveMicroPot = EffectiveMicroPot(EffectiveMicroPot != 0)
    StepNum,sortIdx = __builtint__.sorted(StepNum,'ascend')
    EffectiveMacroPot = EffectiveMacroPot(sortIdx)
    EffectiveMicroPot = EffectiveMicroPot(sortIdx)
    # figure
# hold on
# set(gca, 'xtick',0:0.5:max(StepNum)/2000, 'xminortick', 'on', 'yminortick', 'on');
# set(gcf, 'position', [680   558  640 480]);
# xlabel('Time (ps)');
# ylabel('Effective Potential Drop (\DeltaV, eV)');
# plot(StepNum/2000, EffectiveMacroPot, '-o', 'color', 'r', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'linewidth', 1)
# plot([0 max(StepNum)/2000], [mean(EffectiveMacroPot) mean(EffectiveMacroPot)], '--', 'color', [0.7 0.7 0.7])
# title(['\DeltaV_{Ave.} = ' num2str(mean(EffectiveMacroPot), '#.2f') ' \pm ' num2str(std(EffectiveMacroPot),'#.2f') ' eV']);
# legend(System, 'Average', 'interpreter', 'none')
# hold off
    
    figure
    hold('on')
    set(gca,'xtick',np.arange(0,np.amax(StepNum) / 2000+1,1),'xminortick','on','yminortick','on')
    set(gcf,'position',np.array([680,300,640,480]))
    plt.xlabel('Time (ps)')
    plt.ylabel('Effective Potential Drop (\DeltaV, eV)')
    plt.plot(StepNum / 2000,EffectiveMicroPot,'-o','color','r','markerfacecolor','r','markeredgecolor','k','linewidth',1)
    plt.plot(np.array([np.amin(StepNum) / 2000,np.amax(StepNum) / 2000]),np.array([mean(EffectiveMicroPot),mean(EffectiveMicroPot)]),'--','color',np.array([0.7,0.7,0.7]))
    plt.title(np.array(['\DeltaV_{Ave.} = ',num2str(mean(EffectiveMicroPot),'%.2f'),' \pm ',num2str(std(EffectiveMicroPot),'%.2f'),' eV']))
    plt.legend(System,'Average','interpreter','none')
    hold('off')
    return StepNum,EffectiveMicroPot