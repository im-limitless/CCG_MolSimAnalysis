import numpy as np
    
def writeFarmingInput(BaseFldr = None,Allfldrs = None,cores = None,ppn = None,time = None): 
    # if exist([BaseFldr 'Farm'],'dir')
#     warning('Directory already exists! Farming folder not created.');
#     return
# end
    mkdir(np.array([BaseFldr,'Farm']))
    fidout = open(np.array([BaseFldr,'Farm\Farm-1.restart']),'w')
    fidout.write(np.array(['# CP2K farming file created by Matt Darby, Imperial College London',newline]) % ())
    fidout.write(np.array([' &GLOBAL',newline]) % ())
    fidout.write(np.array(['   PROJECT OldMacDonald',newline]) % ())
    fidout.write(np.array(['   PROGRAM FARMING',newline]) % ())
    fidout.write(np.array(['   RUN_TYPE  NONE',newline]) % ())
    fidout.write(np.array([' &END GLOBAL',newline]) % ())
    fidout.write(np.array([' &FARMING',newline]) % ())
    fidout.write(np.array(['   NGROUPS ',num2str(len(Allfldrs)),newline]) % ())
    fidout.write(np.array(['   GROUP_SIZE ',num2str(cores * ppn / 48),newline]) % ())
    for i in np.arange(1,len(Allfldrs)+1).reshape(-1):
        fidout.write(np.array(['   &JOB',newline]) % ())
        fidout.write(np.array(['     DIRECTORY ../',Allfldrs(i).name(np.arange(1,end() - 4+1)),newline]) % ())
        fidout.write(np.array(['     INPUT_FILE_NAME ',Allfldrs(i).name(np.arange(1,end() - 4+1)),'-1.restart',newline]) % ())
        fidout.write(np.array(['     OUTPUT_FILE_NAME out.log',newline]) % ())
        fidout.write(np.array(['   &END JOB',newline]) % ())
    
    fidout.write(np.array([' &END FARMING',newline]) % ())
    fidout.close()
    CreateSubmissionFile('LRZ',BaseFldr,'Farm',time,cores * len(Allfldrs),ppn,'No')
    return