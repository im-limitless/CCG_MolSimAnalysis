# function [NAtoms,AtomElement,AtomPosition,n,AtomCur,AtomPos,AtomFreeToMove,Vec] = ...
#     CP2kConfigRead(basefldr,flname,foldback)

# Matt Darby - 17-Dec-2020. Imperial College London.
# Function that parses a CP2k .xyz file from basefldr
# flname = '*.xyz'
# Element information is automatically read from .inp if not contained in
# the configuration file
# Version 3.1

import numpy as np
basefldr = 'Z:\Imperial\Federico\10Na10Cl\1ML\__dynamic16'
flname = 'ptwater-pos-1.xyz'
foldback = 0.02
NAtoms = 0
AtomElement = np.array([])
AtomPosition = []
n = []
AtomCur = np.array([])
AtomPos = []
Vec = []
AtomFreeToMove = []
if not len(len(basefldr))==0 :
    if not str(basefldr(end())) == str('\') :
        basefldr = np.array([basefldr,'\'])

# potentialsfile = [basefldr 'POTCAR'];
# if ~exist(potentialsfile,'file')
#     potentialsfile = [basefldr '..\POTCAR'];
# end
configuratfile = np.array([basefldr,flname])
###########################################################################

# Read atom positions and unit cell vectors from .xyz file

fidconf = open(configuratfile)
# First line is a comment
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
# Second line gives the scaling factor for the unit cell
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
ScalF = str2num(textLine)
# Subsequent three lines give the unit vectors
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
va = str2num(textLine)
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
vb = str2num(textLine)
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
vc = str2num(textLine)
Vec = ScalF * np.array([[va],[vb],[vc]])
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
if len(str2num(textLine))==0:
    StrWords = textscan(textLine,'%s')
    AtomCur = np.transpose(StrWords[0])
    textLine = fgetl(fidconf)
    eofstat = feof(fidconf)
else:
    # Read atom types from POTCAR file
    fidpot = open(potentialsfile)
    eofstatpot = False
    while not eofstatpot :

        textLinepot = fgetl(fidpot)
        eofstatpot = feof(fidpot)
        StrWordspot = textscan(textLinepot,'%s')
        if contains(StrWordspot[0][2],'_'):
            newstr = extractBefore(StrWordspot[0][2],'_')
        else:
            newstr = StrWordspot[0][2]
        #         AtomCur = {AtomCur{:} StrWordspot{1}{2}};
        AtomCur = np.array([AtomCur[:],newstr])
        while not eofstatpot :

            textLinepot = fgetl(fidpot)
            eofstatpot = feof(fidpot)
            if eofstatpot or not len(strfind(textLinepot,'End of Dataset'))==0 :
                break


    fidpot.close()

# Next line gives the number of atoms for each element present in the
# system
n = str2num(textLine)
NAtoms = sum(n)
for k in np.arange(1,len(n)+1).reshape(-1):
    for i in np.arange(1,n(k)+1).reshape(-1):
        AtomElement[sum[n[np.arange[1,k - 1+1]]] + i] = AtomCur[k]

# We skip the next two lines
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
textLine = fgetl(fidconf)
eofstat = feof(fidconf)
# We start parsing the atom positions
AtomPosition = np.zeros((NAtoms,3))
AtomFreeToMove = np.zeros((NAtoms,3))
jAtm = 0
while not eofstat :

    textLine = fgetl(fidconf)
    eofstat = feof(fidconf)
    StrWords = textscan(textLine,'%s')
    jAtm = jAtm + 1
    AtomPosition[jAtm,1] = str2num(StrWords[0][0])
    AtomPosition[jAtm,2] = str2num(StrWords[0][2])
    AtomPosition[jAtm,3] = str2num(StrWords[0][3])
    if ('foldback' is not None) and not len(foldback)==0 :
        if len(foldback) == 1:
            foldback = np.array([1,1,1]) * foldback
        for k in np.arange(1,3+1).reshape(-1):
            if foldback(k) > 0:
                if AtomPosition(jAtm,k) > 1 - foldback(k):
                    AtomPosition[jAtm,k] = AtomPosition(jAtm,k) - 1
    for q in np.arange(1,3+1).reshape(-1):
        if str(StrWords[0][q + 3]) == str('T'):
            AtomFreeToMove[jAtm,q] = True
        else:
            if str(StrWords[0][q + 3]) == str('F'):
                AtomFreeToMove[jAtm,q] = False
    if jAtm == NAtoms:
        break


fclose(fidconf)
for k in np.arange(1,len(n)+1).reshape(-1):
    for i in np.arange(1,n(k)+1).reshape(-1):
        AtomPos[k,i,:] = AtomPosition(sum(n(np.arange(1,k - 1+1))) + i,:)

# end