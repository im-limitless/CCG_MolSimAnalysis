import numpy as np
    
def CP2kOptimPathParse(basefldr = None,directory = None,flname = None): 
    # Matt Darby. 21-12-2020. Imperial College London.
# Function that parses a CP2K .xyz file output from basefldr
# Version 3.0
# clear all;
# close all;
# clc
    
    # basefldr = 'G:\Imperial\MattProjects\Pt_Clean\';
# directory = 'CP_Like_1012_Fluoride';
# flname = ['Sample' num2str(startConf) '_' num2str(endConf) '.xyz'];
    print('Creating xtd file...')
    foldback = 0.02
    AtomCur = []
    AtomPosImages = []
    ABC = getABCvectors(basefldr,directory)
    Vec = diag(ABC)
    if not len(len(basefldr))==0 :
        if not str(basefldr(end())) == str('\') :
            basefldr = np.array([basefldr,'\'])
    
    fidxdatcar = open(np.array([basefldr,directory,'\',flname]))
    eofstat = False
    # First line is the number of atoms and is the same for all new sections
    textLine = fgetl(fidxdatcar)
    eofstat = feof(fidxdatcar)
    n = str2num(textLine)
    # Second line contains the iteration (i), time, and total energy (E)
    SectionLine = fgetl(fidxdatcar)
    eofstat = feof(fidxdatcar)
    iconf = 1
    while not eofstat :

        # loop over number of atoms (n), extract names of atoms (AtomCur) for
# config 1 only, extract xyz coordinates (xyzpos) and amalgamate for
# all iterations (AtomPosImages)
        for i in np.arange(1,n+1).reshape(-1):
            textLine = fgetl(fidxdatcar)
            eofstat = feof(fidxdatcar)
            # Read in atoms names
            if iconf == 1:
                AtomCur[1,i] = textLine(isstrprop(textLine,'alpha'))
            # Read the xyz coordinates and wrap into unit cell
            xyzpos = str2num(textLine(np.arange(5,end()+1)))
            AtomPosImages[i,:,iconf] = xyzpos
        if not eofstat :
            textLine = fgetl(fidxdatcar)
            eofstat = feof(fidxdatcar)
            textLine = fgetl(fidxdatcar)
            eofstat = feof(fidxdatcar)
        # progress the iteration counter if a line contains " i =" flag
        if contains(textLine,SectionLine(np.arange(1,4+1))):
            iconf = iconf + 1

    
    AtomCur,sortIdx = __builtint__.sorted(AtomCur)
    for ii in np.arange(1,iconf+1).reshape(-1):
        AtomTemp = AtomPosImages(:,:,ii)
        AtomPosImages[:,:,ii] = AtomTemp(sortIdx,:)
    
    Snapshots = np.arange(1,iconf+1)
    AtomPosImages = AtomPosImages - np.multiply(ABC,(AtomPosImages > ABC))
    AtomPosImages = AtomPosImages + np.multiply(ABC,(AtomPosImages < 0))
    fidxdatcar.close()
    if contains(flname,'.xyz'):
        flname = erase(flname,'.xyz')
    
    XTDFileName = np.array([basefldr,directory,'\',flname,'.xtd'])
    # XTDFileName = [basefldr directory '\Sample' num2str(startConf) '_' num2str(endConf) '.xtd'];
    XTDFileWriteCP2k(XTDFileName,n,AtomCur,AtomPosImages,Vec,5,Snapshots)