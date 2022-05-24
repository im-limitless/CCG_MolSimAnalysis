import numpy as np
import matplotlib.pyplot as plt
clear('all')
close_('all')
# parse the COLVAR file from cp2k

Basefldr = 'G:\Imperial\MattProjects\Pt_Clean\Metadynamics\\'
Allfldrs = dir(np.array([Basefldr,'*Meta_1010_Fluoride']))
flname = dir(np.array([Basefldr,Allfldrs.name,'\*-COLVAR.metadynLog']))
ColvarFile = np.array([Basefldr,Allfldrs.name,'\',flname.name])
fid = open(ColvarFile)
eofstat = False
# First line is the number of atoms and is the same for all new sections
HeadLine = fgetl(fid)
eofstat = feof(fid)
i = 0
data = []
while not eofstat :

    i = i + 1
    textLine = fgetl(fid)
    eofstat = feof(fid)
    data = np.array([[data],[str2num(textLine)]])


fclose(fid)
figure
hold('on')
plt.plot(data(:,1) - data(1,1),data(:,np.arange(2,11+1)))
# plot(data(:,1)-data(1,1), abs(data(:,2:11)));