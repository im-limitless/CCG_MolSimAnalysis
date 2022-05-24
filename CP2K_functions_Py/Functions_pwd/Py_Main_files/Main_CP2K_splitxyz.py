# split .xyz file and extract specific (final) geometry, write to .xtd file
# format for materials studio

# Set the location of the calculation output
import numpy as np
Basefldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\'

system = 'CP_Pit_20F'

Trajectory = 'CP_Pit_20F-pos-1_7.xyz'
nSampleSteps = 500

InitalStepOverride = []

# InitalStepOverride = 10000; # Use to set initial step manually. Set to [] to use entire trajectory

fid = open(np.array([Basefldr,system,'\',Trajectory]))
print(np.array(['Reading xyz trajectory file for ',system]))
lines = textscan(fid,'%s','delimiter','\n','whitespace','')
fclose(fid)
lines = lines[0]
nAtoms = str2num(lines[0])
relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
nConfigs = len(relevant)
initialI = lines[relevant(1)]
EqI = strfind(initialI,'=')
ComI = strfind(initialI,',')
indxI = np.array([EqI(1),ComI(1)])
initialI = str2num(initialI(np.arange(indxI(1) + 1,indxI(2) - 1+1)))
first = np.ceil(initialI / nSampleSteps) * nSampleSteps
finalI = lines[relevant(end())]
EqI = strfind(finalI,'=')
ComI = strfind(finalI,',')
indxI = np.array([EqI(1),ComI(1)])
finalI = str2num(finalI(np.arange(indxI(1) + 1,indxI(2) - 1+1)))
last = int(np.floor(finalI / nSampleSteps)) * nSampleSteps
fidout = open(np.array([Basefldr,system,'\final.xyz']),'w')
for i in np.arange(relevant(end()) - 1,relevant(end()) + nAtoms+1).reshape(-1):
    fidout.write(np.array([lines[i],'\n']) % ())

fclose(fidout)
if not len(InitalStepOverride)==0 :
    first = InitalStepOverride

fidout = open(np.array([Basefldr,system,'\',system,'_',num2str(first),'to',num2str(last),'_',num2str(nSampleSteps),'step.xyz']),'w')

for i in np.arange(first - initialI + 1,last - initialI + 1+nSampleSteps,nSampleSteps).reshape(-1):
    # fidout = fopen([Basefldr fldrnm '\Sample' num2str(1) '_' num2str(length(relevant)) '.xyz'], 'w'); #GO
#     for i = 1:length(relevant) # GO
    print(np.array(['Writing sample trajectory... ',num2str(100 * (i / (last - initialI + 1))),' % complete']))
    for j in np.arange(relevant(i) - 1,relevant(i) + nAtoms+1).reshape(-1):
        fidout.write(np.array([lines[j],'\n']) % ())

fclose(fidout)
CP2kOptimPathParse(Basefldr,system,np.array([system,'_',num2str(first),'to',num2str(last),'_',num2str(nSampleSteps),'step.xyz']))
# CP2kOptimPathParse(Basefldr,fldrnm, 'CP_Pit_18H22F-1_1000_Metal-pos-1.xyz') # GO