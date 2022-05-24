import numpy as np
    
def extractBaderCharges(fldrname = None,flname = None,Atoms = None,AtomList = None): 
    fid = open(np.array([fldrname,flname]))
    # # find the step number from the Bader file append
    StepIndx = find(flname == '_')
    StepIndx = np.array([StepIndx,find(flname == '.')])
    StepNum = str2num(flname(np.arange(StepIndx(1) + 1,StepIndx(2) - 1+1)))
    eofstat = False
    # # Read the Bader charge ACF file
    HeadLine = fgetl(fid)
    eofstat = feof(fid)
    EndLine = fgetl(fid)
    i = 0
    Bader = []
    while not eofstat :

        i = i + 1
        textLine = fgetl(fid)
        eofstat = feof(fid)
        if str(textLine) == str(EndLine):
            break
        else:
            Bader = np.array([[Bader],[str2num(textLine)]])

    
    fid.close()
    # # Extract the raw Bader charge, Q
    Q = Bader(:,5)
    # # Loop over the atom types/number of each type & convert raw to net charges
    TotalCharge = np.zeros((len(Atoms),1))
    NetCharge = np.zeros((len(Atoms),1))
    StndElec = np.zeros((len(Atoms),1))
    StockAtoms = np.array(['F','H','O','Pt'])
    Elecrons = np.array([[7],[1],[6],[18]])
    nAtoms = np.zeros((len(AtomList),1))
    for i in np.arange(1,AtomList.shape[1-1]+1).reshape(-1):
        for j in np.arange(1,Atoms.shape[1-1]+1).reshape(-1):
            if contains(AtomList(i,:),Atoms(j,:)):
                TotalCharge[j,1] = Bader(j,5)
                Indx[j,i] = 1
                for k in np.arange(1,len(StockAtoms)+1).reshape(-1):
                    if contains(AtomList(i,:),StockAtoms[k]):
                        nAtoms[i,1] = nAtoms(i,1) + 1
                        NetCharge[j,1] = Elecrons(k,1)
    
    PtIndxs = find(ismember(AtomList,'Pt'))
    PtIndxs = PtIndxs(PtIndxs < len(AtomList) + 1)
    NetCharge = np.multiply((NetCharge - TotalCharge),Indx)
    Qnet = np.sum(NetCharge, 2-1)
    return Q,Qnet,StepNum
    return Q,Qnet,StepNum