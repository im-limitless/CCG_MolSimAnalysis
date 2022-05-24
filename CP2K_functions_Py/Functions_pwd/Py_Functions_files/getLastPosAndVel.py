import numpy as np
    
def getLastPosAndVel(ParentDir = None,WorkFldr = None,findPos = None,findVel = None): 
    # split .xyz file and extract specific (final) geometry
# ParentDir = 'G:\Imperial\MattProjects\Pt_Clean\';
# WorkFldr = 'BOMD_NVT_1010_Clean';
    if findPos:
        fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'-pos-1.xyz']))
        lines = textscan(fid,'%s','delimiter','\n','whitespace','')
        fid.close()
        lines = lines[0]
        nAtoms = str2num(lines[0])
        relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
        nConfigs = len(relevant)
        PosLines = []
        for i in np.arange(relevant(end()) + 1,relevant(end()) + nAtoms+1).reshape(-1):
            PosLines = np.array([[PosLines],[lines[i]]])
    
    if findVel:
        fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'-vel-1.xyz']))
        lines = textscan(fid,'%s','delimiter','\n','whitespace','')
        fid.close()
        lines = lines[0]
        nAtoms = str2num(lines[0])
        relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
        nConfigs = len(relevant)
        VelLines = []
        for i in np.arange(relevant(end()) + 1,relevant(end()) + nAtoms+1).reshape(-1):
            VelLines = np.array([[VelLines],[lines[i](np.arange(5,end()+1))]])
    
    fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'.inp']))
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    relevant = find(not cellfun(isempty,strfind(lines,'ABC')) )
    VectorLine = lines[relevant]
    fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'.inp']))
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    relevant = find(not cellfun(isempty,strfind(lines,' LIST ')) )
    FixedLine = lines[relevant]
    return PosLines,VelLines,VectorLine,FixedLine
    return PosLines,VelLines,VectorLine,FixedLine