import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
# BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\';
# system = 'CP_Like_1012_Fluoride';
# Trajectory = 'Sample34000_52000.xyz';

# BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CP_Like\';
# system = 'CP_Like_1010_Fluorine';
# Trajectory = 'Sample50000_58000.xyz';

# BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CorrectVolume\';
# system = 'Pt_12H10F';
# Trajectory = 'Pt_12H10F_0to9500_500step.xyz';

BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\'
system = 'CP_Pit_20F'
Trajectory = 'Sample41000_52500.xyz'
# BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Vacuum\';
# system = 'Pt_Bulk';
# Trajectory = 'Pt_Bulk.xyz';

fldrname = np.array([BaseFldr,system,'\Bader\'])
ACFfiles = dir(np.array([fldrname,'ACF_*.dat']))
# DoubleAnalType = 'MassDensity';
DoubleAnalType = 'Radial'
# # Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr,system)
# # get the names of atoms from original xyz input file
Atoms,AtomList,Indx = getAtomNamesFromInputXYZ(BaseFldr,system)
# # Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
Q = np.zeros((len(Atoms),len(ACFfiles)))
Qnet = np.zeros((len(Atoms),len(ACFfiles)))
for n in np.arange(1,len(ACFfiles)+1).reshape(-1):
    Q[:,n],Qnet[:,n],StepNum[n] = extractBaderCharges(fldrname,ACFfiles(n).name,Atoms,AtomList)

# # Find indices of all "Pt*" atoms
if AtomList.shape[1-1] > 1:
    PtList = find(ismember(AtomList,'Pt'))
    PtList = PtList(PtList < len(AtomList) + 1)
    Indxfns = fieldnames(Indx)
    Indx.PtAll = []
    for ii in np.arange(1,len(PtList)+1).reshape(-1):
        Indx.PtAll = np.array([[Indx.PtAll],[getattr(Indx,(Indxfns[PtList(ii)]))]])
else:
    Indxfns = fieldnames(Indx)
    PtList = 1
    Indx.PtAll = getattr(Indx,(Indxfns[0]))

# # parse coordinates of atoms along trajectory and wrap into cell
xyz,XYZ,__,__,__,nAtoms,startConfig,nConfigs,StepNum_Traj = ReadAndParsexyz(BaseFldr,system,Trajectory,ABC)
XYZ = wrapXYZ(XYZ,ABC)
# # compute radial functions for OH, OF and HF
RadFunOH = cell(nConfigs,1)
RadFunFH = cell(nConfigs,1)
RadFunFO = cell(nConfigs,1)
RadFunPtO = cell(nConfigs,1)
DistOH = cell(1,nConfigs)
DistFH = cell(1,nConfigs)
DistFO = cell(1,nConfigs)
DistPtO = cell(1,nConfigs)
for snap in np.arange(startConfig,nConfigs+1).reshape(-1):
    XYZ_snap = np.zeros((XYZ.shape[2-1],XYZ.shape[3-1]))
    XYZ_snap[:,:] = XYZ(snap,:,:)
    VecPtO,DistPtO[snap] = GetAtomCorrelation(XYZ_snap,np.array([[Indx.Pts],[Indx.PtE]]),Indx.O,ABC)
    VecOH,DistOH = GetAtomCorrelation(XYZ_snap,Indx.O,Indx.H,ABC)
    VecFH,DistFH = GetAtomCorrelation(XYZ_snap,Indx.F,Indx.H,ABC)
    VecFO,DistFO = GetAtomCorrelation(XYZ_snap,Indx.F,Indx.O,ABC)
    RadFunOH[snap] = np.reshape(DistOH, tuple(np.array([np.asarray(DistOH).size,1])), order="F")
    RadFunFH[snap] = np.reshape(DistFH, tuple(np.array([np.asarray(DistFH).size,1])), order="F")
    RadFunFO[snap] = np.reshape(DistFO, tuple(np.array([np.asarray(DistFO).size,1])), order="F")
    RadFunPtO[snap] = np.reshape(DistPtO[snap], tuple(np.array([np.asarray(DistPtO[snap]).size,1])), order="F")

MinimaOH = RadialDistribution(RadFunOH,ABC,np.array([['O'],['H']]),0)
MinimaFH = RadialDistribution(RadFunFH,ABC,np.array([['F'],['H']]),0)
MinimaFO = RadialDistribution(RadFunFO,ABC,np.array([['F'],['O']]),0)
MinimaPtO = RadialDistribution(RadFunPtO,ABC,np.array([['Pt'],['O ']]),1)
if str(DoubleAnalType) == str('MassDensity'):
    print('Determining water layering from mass density profile...')
    # # get the O atom distribution and corresponding indices of DL atoms
    Dens_O,__,__,__,__,z = getDensityProfile(xyz,ABC)
    FirstLayerIndx,SecondLayerIndx = getWaterLayerIndices(Indx,XYZ,Dens_O,z)
    for i in np.arange(startConfig,nConfigs+1).reshape(-1):
        DL1st[i] = np.array([FirstLayerIndx[i]])
        DL2nd[i] = np.array([SecondLayerIndx[i]])
        nonDL[i] = setdiff(Indx.O,np.array([[DL1st[i]],[DL2nd[i]]]))
        XYZ_snap = np.zeros((XYZ.shape[2-1],XYZ.shape[3-1]))
        XYZ_snap[:,:] = XYZ(i,:,:)
        __,DistOH1stWL = GetAtomCorrelation(XYZ_snap,DL1st[i],Indx.H,ABC)
        __,DistOH2ndWL = GetAtomCorrelation(XYZ_snap,DL2nd[i],Indx.H,ABC)
        __,DistOHnonDL = GetAtomCorrelation(XYZ_snap,nonDL[i],Indx.H,ABC)
        for j in np.arange(1,len(DL1st[i])+1).reshape(-1):
            DL1st[i] = np.array([[DL1st[i]],[Indx.H(find(DistOH1stWL(:,j) < MinimaOH(1)))]])
        for j in np.arange(1,len(DL2nd[i])+1).reshape(-1):
            DL2nd[i] = np.array([[DL2nd[i]],[Indx.H(find(DistOH2ndWL(:,j) < MinimaOH(1)))]])
        for j in np.arange(1,len(nonDL[i])+1).reshape(-1):
            nonDL[i] = np.array([[nonDL[i]],[Indx.H(find(DistOHnonDL(:,j) < MinimaOH(1)))]])
else:
    if str(DoubleAnalType) == str('Radial'):
        print('Determining water layering from radial distribution...')
        for i in np.arange(startConfig,nConfigs+1).reshape(-1):
            r1st,__ = find(DistPtO[i] <= MinimaPtO(1))
            DL1st[i] = Indx.O(unique(r1st))
            r2nd,__ = find(DistPtO[i] <= np.logical_and(MinimaPtO(2),DistPtO[i]) > MinimaPtO(1))
            DL2nd[i] = Indx.O(unique(r2nd))
            nonDL[i] = setdiff(Indx.O,np.array([[DL1st[i]],[DL2nd[i]]]))
            XYZ_snap = np.zeros((XYZ.shape[2-1],XYZ.shape[3-1]))
            XYZ_snap[:,:] = XYZ(i,:,:)
            __,DistOH1stWL = GetAtomCorrelation(XYZ_snap,DL1st[i],Indx.H,ABC)
            __,DistOH2ndWL = GetAtomCorrelation(XYZ_snap,DL2nd[i],Indx.H,ABC)
            __,DistOHnonDL = GetAtomCorrelation(XYZ_snap,nonDL[i],Indx.H,ABC)
            for j in np.arange(1,len(DL1st[i])+1).reshape(-1):
                DL1st[i] = np.array([[DL1st[i]],[Indx.H(find(DistOH1stWL(:,j) <= MinimaOH(1)))]])
            for j in np.arange(1,len(DL2nd[i])+1).reshape(-1):
                DL2nd[i] = np.array([[DL2nd[i]],[Indx.H(find(DistOH2ndWL(:,j) <= MinimaOH(1)))]])
            for j in np.arange(1,len(nonDL[i])+1).reshape(-1):
                nonDL[i] = np.array([[nonDL[i]],[Indx.H(find(DistOHnonDL(:,j) <= MinimaOH(1)))]])
            nonDL[i] = unique(nonDL[i])

MeanCharge = np.zeros((len(AtomList),len(StepNum)))
StdCharge = np.zeros((len(AtomList),len(StepNum)))
SumCharge = np.zeros((len(AtomList),len(StepNum)))
for i in np.arange(1,len(StepNum)+1).reshape(-1):
    if np.any(StepNum(i) == StepNum_Traj):
        Inter = find(StepNum(i) == StepNum_Traj)
        DL1st_sum[i] = sum(Qnet(DL1st[Inter],i))
        DL2nd_sum[i] = sum(Qnet(DL2nd[Inter],i))
        XYZ_snap = np.zeros((XYZ.shape[2-1],XYZ.shape[3-1]))
        XYZ_snap[:,:] = XYZ(Inter,:,:)
        writeSnaptoxyz(BaseFldr,system,StepNum(i),XYZ_snap,Atoms,np.array([[DL1st[Inter]],[Indx.PtAll]]),DoubleAnalType)
    else:
        DL1st_sum[i] = 0
        DL2nd_sum[i] = 0
    for j in np.arange(1,len(AtomList)+1).reshape(-1):
        MeanCharge[j,i] = mean(Qnet(getattr(Indx,(Indxfns[j])),i))
        StdCharge[j,i] = std(Qnet(getattr(Indx,(Indxfns[j])),i))
        SumCharge[j,i] = sum(Qnet(getattr(Indx,(Indxfns[j])),i))

StepNum,sortIdx = __builtint__.sorted(StepNum,'ascend')
SumCharge = SumCharge(:,sortIdx)
MeanCharge = MeanCharge(:,sortIdx)
StdCharge = StdCharge(:,sortIdx)
DL1st_sum = DL1st_sum(sortIdx)
DL2nd_sum = DL2nd_sum(sortIdx)
# Bader charge and charge density on half-electrode w/out DL
figure
box('on')
set(gca,'colororder',np.array([0,0,0]))
plt.xlabel('Time (ps)')
yyaxis('left')
if len(PtList) == 1:
    HalfElectro = SumCharge(PtList,:) / 2
else:
    if len(PtList) > 1:
        HalfElectro = sum(SumCharge(PtList,:)) / 2

plt.plot(StepNum / 2000,HalfElectro,'-o','color','k','markeredgecolor','k','markerfacecolor','b')
Lyax = gca
etoC = ((1.60218e-19) * (1000000.0)) / (2 * ABC(1) * ABC(2) * (1e-08) * (1e-08))
Lyaxt = get(Lyax,'YTick')
LyaxDegC = np.round(10 * Lyaxt * etoC) / 10
plt.ylabel('Total Charge (e)')
yyaxis('right')
plt.plot(StepNum / 2000,HalfElectro,'-o','color','k','markeredgecolor','k','markerfacecolor','b')
Ryax = gca
set(Ryax,'YTick',Lyaxt,'YTickLabel',LyaxDegC)
plt.ylabel('Charge Density (\muC\cdotcm^{-2})')
plt.legend('Half-Electrode')
# Bader charge and charge density on half-electrode with DL
figure
box('on')
set(gca,'colororder',np.array([0,0,0]))
plt.xlabel('Time (ps)')
yyaxis('left')
plt.plot(StepNum(DL1st_sum != 0) / 2000,HalfElectro(DL1st_sum != 0) + (DL1st_sum(DL1st_sum != 0) / 2) + (DL2nd_sum(DL1st_sum != 0) / 2),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
Lyax = gca
etoC = ((1.60218e-19) * (1000000.0)) / (2 * ABC(1) * ABC(2) * (1e-08) * (1e-08))
Lyaxt = get(Lyax,'YTick')
LyaxDegC = np.round(10 * Lyaxt * etoC) / 10
plt.ylabel('Total Charge (e)')
yyaxis('right')
plt.plot(StepNum(DL1st_sum != 0) / 2000,HalfElectro(DL1st_sum != 0) + (DL1st_sum(DL1st_sum != 0) / 2) + (DL2nd_sum(DL1st_sum != 0) / 2),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
Ryax = gca
set(Ryax,'YTick',Lyaxt,'YTickLabel',LyaxDegC)
plt.ylabel('Charge Density (\muC\cdotcm^{-2})')
plt.legend('Half-Electrode+DL')
# Bader charge and charge density on 1st layer of DL
figure
box('on')
set(gca,'colororder',np.array([0,0,0]))
plt.xlabel('Time (ps)')
yyaxis('left')
plt.plot(StepNum(DL1st_sum != 0) / 2000,(DL1st_sum(DL1st_sum != 0) / 2),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
Lyax = gca
etoC = ((1.60218e-19) * (1000000.0)) / (2 * ABC(1) * ABC(2) * (1e-08) * (1e-08))
Lyaxt = get(Lyax,'YTick')
LyaxDegC = np.round(10 * Lyaxt * etoC) / 10
plt.ylabel('Total Charge (e)')
yyaxis('right')
plt.plot(StepNum(DL1st_sum != 0) / 2000,(DL1st_sum(DL1st_sum != 0) / 2),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
Ryax = gca
set(Ryax,'YTick',Lyaxt,'YTickLabel',LyaxDegC)
plt.ylabel('Charge Density (\muC\cdotcm^{-2})')
plt.legend('1st Water Layer')
# Bader charge and charge density on 2nd layer of DL
figure
box('on')
set(gca,'colororder',np.array([0,0,0]))
plt.xlabel('Time (ps)')
yyaxis('left')
plt.plot(StepNum(DL2nd_sum != 0) / 2000,(DL2nd_sum(DL2nd_sum != 0) / 2),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
Lyax = gca
etoC = ((1.60218e-19) * (1000000.0)) / (2 * ABC(1) * ABC(2) * (1e-08) * (1e-08))
Lyaxt = get(Lyax,'YTick')
LyaxDegC = np.round(10 * Lyaxt * etoC) / 10
plt.ylabel('Total Charge (e)')
yyaxis('right')
plt.plot(StepNum(DL2nd_sum != 0) / 2000,(DL2nd_sum(DL2nd_sum != 0) / 2),'-o','color','k','markeredgecolor','k','markerfacecolor','r')
Ryax = gca
set(Ryax,'YTick',Lyaxt,'YTickLabel',LyaxDegC)
plt.ylabel('Charge Density (\muC\cdotcm^{-2})')
plt.legend('2nd Water Layer')
# Total charge per Pt type
figure
box('on')
hold('on')
plt.xlabel('Time (ps)')
plt.ylabel('Total Charge (e)')
# title(['Total Bader charge for Pt Atoms in ' system], 'interpreter', 'none')
C = np.array([[218 / 255,165 / 255,32 / 255],[0,0.5,0],[0,0,1],[1,0,0]])
PtC = []
for i in np.arange(1,len(PtList)+1).reshape(-1):
    if contains(AtomList(PtList(i),:),'PtE'):
        PtC = C(1,:)
    else:
        if contains(AtomList(PtList(i),:),'Ptb'):
            PtC = C(2,:)
        else:
            if np.logical_and(contains(AtomList(PtList(i),:),'Pts'),not contains(AtomList(PtList(i),:),'Ptss') ):
                PtC = C(3,:)
            else:
                if contains(AtomList(PtList(i),:),'Ptss'):
                    PtC = C(4,:)
    plt.plot(StepNum / 2000,SumCharge(PtList(i),:) / 2,'-o','color','k','markeredgecolor','k','markerfacecolor',PtC)

if len(PtList) == 2:
    plt.legend(AtomList(PtList(1),:),AtomList(PtList(2),:),'interpreter','tex')
else:
    if len(PtList) == 3:
        plt.legend(AtomList(PtList(1),:),AtomList(PtList(2),:),AtomList(PtList(3),:),'interpreter','tex')

hold('off')
# mean charge per Pt type
figure
box('on')
hold('on')
plt.xlabel('Time (ps)')
plt.ylabel('Ave. Charge per Atom (e)')
# title(['Total Bader charge for Pt Atoms in ' system], 'interpreter', 'none')
PtC = []
for i in np.arange(1,len(PtList)+1).reshape(-1):
    if contains(AtomList(PtList(i),:),'PtE'):
        PtC = C(1,:)
    else:
        if contains(AtomList(PtList(i),:),'Ptb'):
            PtC = C(2,:)
        else:
            if np.logical_and(contains(AtomList(PtList(i),:),'Pts'),not contains(AtomList(PtList(i),:),'Ptss') ):
                PtC = C(3,:)
            else:
                if contains(AtomList(PtList(i),:),'Ptss'):
                    PtC = C(4,:)
    errorbar(StepNum / 2000,MeanCharge(PtList(i),:),StdCharge(PtList(i),:),'-o','color',PtC,'markeredgecolor','k','markerfacecolor',PtC)

if len(PtList) == 2:
    plt.legend(AtomList(PtList(1),:),AtomList(PtList(2),:),'interpreter','tex')
else:
    if len(PtList) == 3:
        plt.legend(AtomList(PtList(1),:),AtomList(PtList(2),:),AtomList(PtList(3),:),'interpreter','tex')

hold('off')
figure
box('on')
hold('on')
plt.xlabel('Time (ps)')
plt.ylabel('Total Charge (e)')
plt.title(np.array(['Total Bader charge for all electrolyte species in ',system]),'interpreter','none')
C = np.array([[1,0,0],[0,0.5,0],[0,0,1]])
IonList = np.arange(1,len(AtomList)+1)
IonList[PtList] = []
IonSumCharge = sum(SumCharge(IonList,:))
plt.plot(StepNum / 2000,IonSumCharge,'-o','color','k','markeredgecolor','k','markerfacecolor','m')
plt.legend('Total Electrolyte','interpreter','tex')
hold('off')
for i in np.arange(1,len(AtomList)+1).reshape(-1):
    if not ismember(i,PtList) :
        figure
        box('on')
        hold('on')
        h1 = plt.plot(StepNum / 2000,SumCharge(i,:),'-o','color','k','markeredgecolor','k','markerfacecolor',C(i,:))
        h2 = plt.plot(np.array([StepNum(1) / 2000,StepNum(end()) / 2000]),np.array([mean(SumCharge(i,np.arange(2,end()+1))),mean(SumCharge(i,np.arange(2,end()+1)))]),'--','color',np.array([0.7,0.7,0.7]))
        AveStr = np.array(['Traj. Ave. = ',num2str(mean(SumCharge(i,np.arange(2,end()+1))),'%.2f'),'\pm',num2str(std(SumCharge(i,np.arange(2,end()+1))),'%.2f')])
        plt.legend(AtomList(i,:),AveStr,'interpreter','tex')
        plt.xlabel('Time (ps)')
        plt.ylabel('Total Charge (e)')
        plt.title(np.array(['Total Bader charge for ',AtomList(i,:),' in ',system]),'interpreter','none')
        uistack(h1,'top')
        hold('off')

StepNumPot,EffPotDrop = CP2K_CalcEffectivePotentialDrop(BaseFldr,system)
tf,indx = ismember(StepNumPot,StepNum)
figure
hold('on')
plt.ylabel('Total Charge (e)')
plt.xlabel('Electrostatic Potential (V)')
# plot(EffPotDrop, HalfElectro(indx), 'o', 'color', 'k', 'markerfacecolor', 'r')

WLElectro = HalfElectro(DL1st_sum != 0) + (DL1st_sum(DL1st_sum != 0) / 2) + (DL2nd_sum(DL1st_sum != 0) / 2)
plt.plot(EffPotDrop,etoC * WLElectro(indx),'o','color','k','markerfacecolor','r')
for i in np.arange(1,len(indx)+1).reshape(-1):
    text(EffPotDrop(i),etoC * WLElectro(indx(i)),np.array([num2str(StepNum(indx(i)) / 2000),' s']),'verticalalignment','top','horizontalalignment','right')

hold('off')
PtNums = []
for i in np.arange(1,len(PtList)+1).reshape(-1):
    PtNums = np.array([[PtNums],[getattr(Indx,(Indxfns[PtList(i)]))]])

# XYZ_ave(:,:) = mean(XYZ, 1);

MeanQnet = mean(Qnet(PtNums,np.arange(1,end()+1)),2)
Bader3DCharge(XYZ_snap(PtNums,:),ABC,MeanQnet)
# XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
# XYZ_snap(:,:) = XYZ(1,:,:);
# MeanQnet = mean(Qnet(PtNums,:),2);
# Bader3DCharge(XYZ_snap(PtNums,:), ABC, Qnet(PtNums,1));