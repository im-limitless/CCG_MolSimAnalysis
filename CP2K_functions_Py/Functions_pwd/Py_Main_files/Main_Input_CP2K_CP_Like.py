import numpy as np
fclose('all')
clear('all')
newlinechar = char(10)
# # ########################################################################
# Read Materials Studio file
BaseInFldr = 'G:\Imperial\MattProjects\Pt_O\10Na12Cl\'
fldrname = 'BOMD_NVT_1012_O_0.66ML'
Type = 'DMP'
# Type = 'Noisy';

BaseOutFldr = np.array([BaseInFldr,fldrname,'\',Type,'\'])
## ########################################################################
# ##### Choose a Machine #####
# MyCluster = 'Legion';
# MyCluster = 'Grace'; # Check the cell size and k mesh...
# MyCluster = 'Rhea';
# MyCluster = 'LRZ';
MyCluster = 'LRZ-test'
# MyCluster = 'IBchem';

# # ##### Set the number of cores #####
# ncores = 1;
# ncores = 16;
# ncores = 32;
# ncores = 48;
# ncores = 64;
# ncores = 432;
ncores = 816
# ##### Set the walltime in hours #####
# Walltime = 4;
# Walltime = 12;
# Walltime = 24;
# Walltime = 48;
Walltime = 0.5
VarNames = np.array(['SS','EX','GA','NG'])
# Vars.Params = {[0.0025 0.005 0.0075 0.01 0.0125] [0] [0] [3e-4]};
Vars.Params = np.array([np.array([0.0075]),np.array([0]),np.array([0]),np.array([0.0003])])
Vars.BOMD_NVT_1010_Clean = np.array([np.array([0.01]),np.array([0]),np.array([0]),np.array([0.0005,0.001,0.005,0.01,0.05])])
## ########################################################################

if str(Type) == str('DMP'):
    Vars = Vars.Params
else:
    if str(Type) == str('Noisy'):
        Vars = getattr(Vars,(fldrname))

Coord,Vel,ABC,Fix = getLastPosAndVel(BaseInFldr,fldrname)
for i in np.arange(1,len(Vars[0])+1).reshape(-1):
    for j in np.arange(1,len(Vars[2])+1).reshape(-1):
        for k in np.arange(1,len(Vars[3])+1).reshape(-1):
            for l in np.arange(1,len(Vars[4])+1).reshape(-1):
                DirString = np.array([VarNames[0],'_',num2str(Vars[0](i)),'_',VarNames[2],'_',num2str(Vars[2](j)),'_',VarNames[3],'_',num2str(Vars[3](k)),'_',VarNames[4],'_',num2str(Vars[4](l))])
                CreateCP2KInputfile_CP_Like(BaseInFldr,BaseOutFldr,fldrname,Walltime,Coord,Vel,ABC,Fix,Vars[0](i),Vars[2](j),Vars[3](k),Vars[4](l),DirString)
                CreateSubmissionFile(MyCluster,BaseOutFldr,DirString,Walltime,ncores,'No')
                copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS',np.array([BaseOutFldr,DirString,'\GTH_POTENTIALS']))
                copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS',np.array([BaseOutFldr,DirString,'\GTH_BASIS_SETS']))
                copyfile('G:\Imperial\CP2k_files\dftd3.dat',np.array([BaseOutFldr,DirString,'\dftd3.dat']))
