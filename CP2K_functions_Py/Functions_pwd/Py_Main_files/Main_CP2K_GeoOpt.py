import numpy as np
import matplotlib.pyplot as plt
import warnings
clear('all')
close_('all')
BaseFldr = 'G:\Imperial\MRes\Kehan\Matt\Clean\'
system = dir(np.array([BaseFldr,'K*']))
mkdir(np.array([BaseFldr,'Energy']))
for j in np.arange(1,len(system)+1).reshape(-1):
    print(np.array(['Processing ',system(j).name,'...']))
    logfile = np.array([BaseFldr,system(j).name,'\out.log'])
    fid = open(logfile)
    textLines = textscan(fid,'%s','delimiter','\n','whitespace','')
    textLines = textLines[0]
    fid.close()
    Conv = find(not cellfun(isempty,strfind(textLines,'Max. gradient ')) )
    TotEng = find(not cellfun(isempty,strfind(textLines,'ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):')) )
    Tol = find(not cellfun(isempty,strfind(textLines,'Conv. limit for gradients')) )
    Tol = (27.2114) * (0.529177 ** 2) * str2num(textLines[Tol(1)](np.arange(33,end()+1)))
    F = np.zeros((len(Conv),1))
    E = np.zeros((len(TotEng),1))
    for i in np.arange(1,len(Conv)+1).reshape(-1):
        F[i,1] = (27.2114) * (0.529177 ** 2) * str2num(textLines[Conv(i)](np.arange(33,end()+1)))
    for i in np.arange(1,len(TotEng)+1).reshape(-1):
        E[i,1] = (27.2114) * str2num(textLines[TotEng(i)](np.arange(50,end()+1)))
    figure
    semilogy(F,'-ok','markerfacecolor','b')
    hold('on')
    semilogy(np.array([0,len(F)]),np.array([Tol,Tol]),'--r')
    plt.ylabel('Force (eV\cdotA^{-2})')
    plt.xlabel('Optimisation Step')
    hold('off')
    drawnow
    figure
    plt.plot(E,'-ok','markerfacecolor','b')
    plt.ylabel('Total Energy (eV)')
    plt.xlabel('Optimisation Step')
    drawnow
    if F(end()) > Tol:
        warnings.warn(np.array(['Forces only converged to ',num2str(F(end())),'...']))
    print(np.array(['Total Energy = ',num2str(E(end())),' eV']))
    dlmwrite(np.array([BaseFldr,system(j).name,'\energy',system(j).name]),E(end()),'precision',15)
    copyfile(np.array([BaseFldr,system(j).name,'\energy',system(j).name]),np.array([BaseFldr,'\Energy\']))
    CP2kOptimPathParse(BaseFldr,system(j).name,np.array([system(j).name,'-pos-1.xyz']))
