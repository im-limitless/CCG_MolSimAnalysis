import numpy as np
import os
import warnings
clear('all')
close_('all')
# script to find coordinates from restart file and make metal only +
# eletrolyte only inputs for Electrostatic Potentials

BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\CP_Pit_18H22F\'
# BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CorrectVolume\Pt_12H10F\';

Restarts = dir(np.array([BaseFldr,'*-1_*.restart']))
for i in np.arange(1,len(Restarts)+1).reshape(-1):
    print(np.array(['Processing restart file ',Restarts(i).name]))
    OutFldr = np.array([BaseFldr,'EffectivePD\',Restarts(i).name(np.arange(1,end() - 8+1)),'\'])
    if os.path.exist(str(OutFldr)):
        warnings.warn('Directory already exists! Continuing with next file.')
        continue
    mkdir(OutFldr)
    fid = open(np.array([BaseFldr,Restarts(i).name]))
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    startCELL = find(not cellfun(isempty,strfind(lines,'&CELL')) )
    endCELL = find(not cellfun(isempty,strfind(lines,'&END CELL')) )
    IndxCell = np.arange(startCELL,endCELL+1)
    startXYZ = find(not cellfun(isempty,strfind(lines,'&COORD')) )
    endXYZ = find(not cellfun(isempty,strfind(lines,'&END COORD')) )
    IndxPt = find(not cellfun(isempty,strfind(lines,'Pt')) )
    IndxPt = IndxPt(IndxPt < np.logical_and(endXYZ,IndxPt) > startXYZ)
    IndxElectrolyte = np.transpose((np.arange(startXYZ + 1,endXYZ - 1+1)))
    IndxElectrolyte = IndxElectrolyte(find(not ismember(IndxElectrolyte,IndxPt) ))
    ## write the file
    MetalOut = np.array([Restarts(i).name(np.arange(1,end() - 8+1)),'_Metal'])
    if os.path.exist(str(np.array([OutFldr,MetalOut]))):
        warnings.warn('Directory already exists! Will not overwrite.')
    else:
        mkdir(np.array([OutFldr,MetalOut]))
        InputFromRestart(OutFldr,MetalOut,IndxPt,IndxCell,lines)
        CreateSubmissionFile('LRZ',OutFldr,MetalOut,1800,768,24,'No')
        copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS',np.array([OutFldr,MetalOut,'\GTH_POTENTIALS']))
        copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS',np.array([OutFldr,MetalOut,'\GTH_BASIS_SETS']))
        copyfile('G:\Imperial\CP2k_files\dftd3.dat',np.array([OutFldr,MetalOut,'\dftd3.dat']))
    ElectrolyteOut = np.array([Restarts(i).name(np.arange(1,end() - 8+1)),'_Electrolyte'])
    if os.path.exist(str(np.array([OutFldr,'\',Restarts(i).name(np.arange(1,end() - 8+1)),'_Electrolyte']))):
        warnings.warn('Directory already exists! Will not overwrite.')
    else:
        mkdir(np.array([OutFldr,ElectrolyteOut]))
        InputFromRestart(OutFldr,ElectrolyteOut,IndxElectrolyte,IndxCell,lines)
        CreateSubmissionFile('LRZ',OutFldr,ElectrolyteOut,1800,768,24,'No')
        copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS',np.array([OutFldr,ElectrolyteOut,'\GTH_POTENTIALS']))
        copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS',np.array([OutFldr,ElectrolyteOut,'\GTH_BASIS_SETS']))
        copyfile('G:\Imperial\CP2k_files\dftd3.dat',np.array([OutFldr,ElectrolyteOut,'\dftd3.dat']))

fclose('all')