import numpy as np
    
def writeSnaptoxyz(BaseFldr = None,fldrname = None,StepNum = None,XYZ_snap = None,Atoms = None,DL = None,DoubleAnalType = None): 
    # if ~exist([BaseFldr fldrname '\1stWL_Trajectories'])
#     mkdir([BaseFldr fldrname '\1stWL_Trajectories']);
# end
    
    # fidout = fopen([BaseFldr fldrname '\1stWL_Trajectories\'  DoubleAnalType '_' num2str(StepNum) '.xyz'], 'w');
    fidout = open(np.array([BaseFldr,fldrname,'\',DoubleAnalType,'_',num2str(StepNum),'.xyz']),'w')
    fidout.write(np.array([num2str(len(DL)),newline]) % ())
    fidout.write(np.array(['i = ',num2str(StepNum),newline]) % ())
    for j in np.arange(1,len(DL)+1).reshape(-1):
        fidout.write(np.array([pad(Atoms(DL(j),:),8),pad(num2str(XYZ_snap(DL(j),1),'%.10g'),20),pad(num2str(XYZ_snap(DL(j),2),'%.10g'),20),pad(num2str(XYZ_snap(DL(j),3),'%.10g'),20),newline]) % ())
    
    fidout.close()
    return