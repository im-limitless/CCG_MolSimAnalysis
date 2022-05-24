import numpy as np
import os
import warnings
fclose('all')
clear('all')
newlinechar = char(10)
# # ########################################################################
# Read Materials Studio file
BaseInFldr = 'E:\Materials Studio Projects\TempProject Files\Documents\Flat\Adjusted\'

BaseOutFldr = 'G:\Imperial\ScriptTestingSuite\'

FileList = dir(np.array([BaseInFldr,'*_*.xsd']))

## ########################################################################
##### Choose a Machine #####
# MyCluster = 'Legion';
# MyCluster = 'Grace';
# MyCluster = 'Rhea';
MyCluster = 'LRZ'
# MyCluster = 'IBchem';

##### Set the number of cores #####
# ncores = 1;
# ncores = 16;
# ncores = 32;
# ncores = 48;
# ncores = 64;
# ncores = 624;
# ncores = 768;
ncores = 816
# ncores = 1152;
# ncores = 1776;

##### Set the number of processors per node #####
ppn = 40
# ppn = 24;

##### Set the calculation type #####
CalcType = 'BOMD'
# CalcType = 'CP';
# CalcType = 'GEO';

##### Set the walltime in hours #####
Walltime = 8
# Walltime = 12;
# Walltime = 24;
# Walltime = 48;
# Walltime = 120;

##### Choose if the jobs are to be farmed #####
# Farming = 'Yes';
Farming = 'No'
##### Make unique input file to converge wavefunction with no MD #####
# WFN = 'Yes';
WFN = 'No'
## ########################################################################

for i in np.arange(1,len(FileList)+1).reshape(-1):
    flname = FileList(i).name
    print(np.array(['Processing ',FileList(i).name,'...']))
    fldrname = flname(np.arange(1,end() - 4+1))
    fldrname = fldrname.replace(' ','_')
    fldrname = fldrname.replace('(','_')
    fldrname = fldrname.replace(')','_')
    if os.path.exist(str(np.array([BaseOutFldr,fldrname]))):
        warnings.warn('Directory already exists! Continuing with next file.')
        continue
    mkdir(np.array([BaseOutFldr,fldrname]))
    XSDFileName = np.array([BaseInFldr,flname])
    NAtoms,srtdAtomElement,srtdAtomPosition,n,AtomCur,AtomPos,Vec = XSDFileRead(XSDFileName)
    ListOfElems,Fix,Restrain,Vel = convertXSDtoXYZ(BaseOutFldr,fldrname,NAtoms,n,srtdAtomPosition,srtdAtomElement,Vec)
    CreateCP2KInputfile(BaseOutFldr,fldrname,Vec,Walltime,Fix,Restrain,CalcType,Vel,ListOfElems,WFN)
    #     CreateSubmissionFile(MyCluster, BaseOutFldr, fldrname, Walltime, ncores, ppn, 'No');
    # Copy the potential, basis and vdw kernel files from a base location
# into each directory
    copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS',np.array([BaseOutFldr,fldrname,'\GTH_POTENTIALS']))
    copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS',np.array([BaseOutFldr,fldrname,'\GTH_BASIS_SETS']))
    copyfile('G:\Imperial\CP2k_files\dftd3.dat',np.array([BaseOutFldr,fldrname,'\dftd3.dat']))
    #     # Copy the chainjob scripts to each directory
#     copyfile('G:\Imperial\CP2k_files\chainjob.qs', [BaseOutFldr fldrname '\chainjob.qs']);
#     copyfile('G:\Imperial\CP2k_files\AddChain.qs', [BaseOutFldr fldrname '\AddChain.qs']);

if str(Farming) == str('Yes'):
    writeFarmingInput(BaseOutFldr,FileList,ncores,ppn,Walltime)
