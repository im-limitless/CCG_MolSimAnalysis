import numpy as np
import warnings
    
def CP2k_parse_log(Basef = None,printNonConv = None): 
    logfile = np.array([Basef,'\out.log'])
    fid = open(logfile)
    textLines = textscan(fid,'%s','delimiter','\n','whitespace','')
    textLines = textLines[0]
    Conv = find(not cellfun(isempty,strfind(textLines,'*** SCF run converged in ')) )
    NotConv = find(not cellfun(isempty,strfind(textLines,' SCF run NOT ')) )
    nInnerSCF = find(not cellfun(isempty,strfind(textLines,'Leaving inner SCF loop after reaching')) )
    AllStep = __builtint__.sorted(np.array([[Conv],[NotConv]]))
    NotConvStepNum = find(ismember(AllStep,NotConv))
    ConvStepNum = find(ismember(AllStep,Conv))
    if not len(NotConv)==0  and printNonConv == 1:
        for i in np.arange(1,len(NotConvStepNum)+1).reshape(-1):
            warnings.warn(np.array(['SCF steps ',num2str(NotConvStepNum(i)),' failed to converge']))
    
    PercConv = 100 * len(NotConv) / (len(Conv) + len(NotConv))
    print(np.array(['Percentage of SCF steps that did not converge is ',num2str(PercConv),'%']))
    nSteps = np.zeros((len(Conv),1))
    for i in np.arange(1,len(Conv)+1).reshape(-1):
        pIndx = strfind(textLines[Conv(i)],'steps')
        nSteps[i] = str2num(textLines[Conv(i)](np.arange(pIndx - 5,pIndx - 1+1)))
    
    nSCF = np.zeros((len(nInnerSCF),1))
    for i in np.arange(1,len(nInnerSCF)+1).reshape(-1):
        pIndx = strfind(textLines[nInnerSCF(i)],'steps')
        nSCF[i] = str2num(textLines[nInnerSCF(i)](np.arange(pIndx - 5,pIndx - 1+1)))
    
    AveSteps = mean(np.array([[nSCF(np.arange(3,end()+1))],[nSteps(np.arange(3,end()+1))]]))
    StdSteps = std(np.array([[nSCF(np.arange(3,end()+1))],[nSteps(np.arange(3,end()+1))]]))
    print(np.array(['Average number of SCF steps performed is ',num2str(AveSteps)]))
    fid.close()
    return ConvStepNum,NotConvStepNum,PercConv,AveSteps,StdSteps