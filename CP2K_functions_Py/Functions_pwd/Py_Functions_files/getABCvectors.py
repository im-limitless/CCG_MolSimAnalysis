import numpy as np
import os
import warnings
    
def getABCvectors(ParentDir = None,WorkFldr = None): 
    if os.path.exist(str(np.array([ParentDir,WorkFldr,'\',WorkFldr,'-1.restart']))):
        fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'-1.restart']))
    else:
        if os.path.exist(str(np.array([ParentDir,WorkFldr,'\',WorkFldr,'.inp']))):
            fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'.inp']))
        else:
            if os.path.exist(str(np.array([ParentDir,WorkFldr,'\',WorkFldr,'.pdb']))):
                fid = open(np.array([ParentDir,WorkFldr,'\',WorkFldr,'.pdb']))
                lines = textscan(fid,'%s','delimiter','\n','whitespace','')
                fid.close()
                lines = lines[0]
                VectorLine = str2num(lines[0](np.arange(7,end() - 15+1)))
                Vectors = VectorLine(np.arange(1,3+1))
                return Vectors
            else:
                warnings.warn('No input or restart file detected...')
                return Vectors
    
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    relevant = find(not cellfun(isempty,strfind(lines,'ABC')) )
    if not len(relevant)==0 :
        VectorLine = lines[relevant]
        vecIndx = find(VectorLine == 'C')
        Vectors = str2num(VectorLine(np.arange(vecIndx + 1,end()+1)))
        return Vectors
    else:
        ALineIndx = find(not cellfun(isempty,strfind(lines,' A ')) )
        BLineIndx = find(not cellfun(isempty,strfind(lines,' B ')) )
        CLineIndx = find(not cellfun(isempty,strfind(lines,' C ')) )
        AVectorLine = lines[ALineIndx]
        vecStartA = find(AVectorLine == 'A')
        A = str2num(AVectorLine(np.arange(vecStartA + 1,end()+1)))
        BVectorLine = lines[BLineIndx]
        vecStartB = find(BVectorLine == 'B')
        B = str2num(BVectorLine(np.arange(vecStartB + 1,end()+1)))
        CVectorLine = lines[CLineIndx]
        vecStartC = find(CVectorLine == 'C')
        C = str2num(CVectorLine(np.arange(vecStartC + 1,end()+1)))
        if norm(A) == np.logical_and(sum(A),norm(B)) == np.logical_and(sum(B),norm(C)) == sum(C):
            Vectors = np.transpose(diag(np.array([[A],[B],[C]])))
        else:
            warnings.warn('Lattice vectors are not orthogonal...')
        return Vectors
    
    return Vectors