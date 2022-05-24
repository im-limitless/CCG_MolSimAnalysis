import numpy as np
    
def CP2k_parse_SCF_Loops(Basef = None): 
    # logfile = 'G:\Imperial\CP-LIKE\BulkWater\CPLIKE\Parameterisation\FromShort\STEPSIZE\WaterCP_fromShort_ST_1.00E-01\out.log';
    logfile = np.array([Basef,'\out.log'])
    fid = open(logfile)
    textLines = textscan(fid,'%s','delimiter','\n','whitespace','')
    textLines = textLines[0]
    fid.close()
    OTindx = find(not cellfun(isempty,strfind(textLines,'OT DIIS')) )
    # OTdat =[];
    
    for i in np.arange(1,len(OTindx) - 1+1).reshape(-1):
        if OTindx(i + 1) - OTindx(i) == 1:
            continue
        else:
            if OTindx(i + 1) - OTindx(i) > 1:
                A = strsplit(textLines[OTindx(i)])
                OTdat = str2num(A[7])
                Conv[i] = OTdat
    
    Conv = Conv(Conv != 0)
    AveError = mean(Conv(np.arange(3,end()+1)))
    StdError = std(Conv(np.arange(3,end()+1)))
    return AveError,StdError