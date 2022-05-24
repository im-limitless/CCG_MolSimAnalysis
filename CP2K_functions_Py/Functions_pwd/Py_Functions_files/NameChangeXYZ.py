# NameChangeXYZ
import numpy as np
clear('all')
Base = 'G:\Imperial\MattProjects\MD_files\Pit_NVT\'
system = 'CP_Pit2020_H_1.00ML'
BaseFldr = np.array([Base,system,'\'])
AllFldrs = dir(np.array([BaseFldr,'SS1*']))
mkdir(np.array([BaseFldr,'xtd']))
for j in np.arange(1,len(AllFldrs)+1).reshape(-1):
    ABC = getABCvectors(BaseFldr,AllFldrs(j).name)
    fid = open(np.array([BaseFldr,AllFldrs(j).name,'\',AllFldrs(j).name,'.xyz']))
    print('Reading xyz data...')
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    lines[2] = 'i = 0'
    nAtoms = str2num(lines[0])
    relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
    nConfigs = len(relevant)
    Step = []
    Atoms = []
    for i in np.arange(1,nAtoms+1).reshape(-1):
        if contains(lines[i + 2],'PTP'):
            lines[i + 2] = lines[i + 2].replace('PTP','Pt ')
        else:
            if contains(lines[i + 2],'HW1'):
                lines[i + 2] = lines[i + 2].replace('HW1','H  ')
            else:
                if contains(lines[i + 2],'HW2'):
                    lines[i + 2] = lines[i + 2].replace('HW2','H  ')
                else:
                    if contains(lines[i + 2],'OW'):
                        lines[i + 2] = lines[i + 2].replace('OW','O ')
                    else:
                        if contains(lines[i + 2],'H1'):
                            lines[i + 2] = lines[i + 2].replace('H1','H ')
                        else:
                            if contains(lines[i + 2],'H2'):
                                lines[i + 2] = lines[i + 2].replace('H2','H ')
                            else:
                                if contains(lines[i + 2],'H3'):
                                    lines[i + 2] = lines[i + 2].replace('H3','H ')
                                else:
                                    if contains(lines[i + 2],'HS'):
                                        lines[i + 2] = lines[i + 2].replace('HS','B ')
    for i in np.arange(1,nAtoms+1).reshape(-1):
        Atoms = np.array([[Atoms],[lines[(relevant(1) + i)](np.arange(1,5+1))]])
    AtomList = unique(Atoms,'rows')
    XYZ = np.zeros((nConfigs,nAtoms,3))
    for i in np.arange(1,AtomList.shape[1-1]+1).reshape(-1):
        setattr(Indx,strtrim(AtomList(i,:)),find(np.sum(Atoms == AtomList(i,:), 2-1) == len(AtomList(i,:))))
        setattr(xyz,strtrim(AtomList(i,:)),np.zeros((nConfigs,len(getattr(Indx,(strtrim(AtomList(i,:))))),3)))
    fid = open(np.array([BaseFldr,AllFldrs(j).name,'\',AllFldrs(j).name,'_new.xyz']),'w')
    for i in np.arange(1,len(lines)+1).reshape(-1):
        fid.write(np.array([lines[i],newline]) % ())
    fid.close()
    CP2kOptimPathParse(BaseFldr,AllFldrs(j).name,np.array([AllFldrs(j).name,'_new.xyz']))
    copyfile(np.array([BaseFldr,AllFldrs(j).name,'\',AllFldrs(j).name,'_new.xtd']),np.array([BaseFldr,'xtd\',system,'_',num2str(j),'.xtd']))
    copyfile(np.array([BaseFldr,AllFldrs(j).name,'\',AllFldrs(j).name,'_new.arc']),np.array([BaseFldr,'xtd\',system,'_',num2str(j),'.arc']))
