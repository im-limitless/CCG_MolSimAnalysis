import numpy as np
clear('all')
close_('all')
fid = open('C:\Users\Matt\Desktop\Debug\PDOS_Al_Bulk_cutoff200\pdos')
data = fscanf(fid,'%f',np.array([12,inf]))
data = np.transpose(data)
sigma = 0.02
de = 0.001
scale = 1
EIGENVALUE_COLUMN = 2
DENSITY_COLUMN = 4
margin = 10 * sigma
emin = np.amin(np.amin(data(:,EIGENVALUE_COLUMN))) - margin
emax = np.amax(np.amax(data(:,EIGENVALUE_COLUMN))) + margin
ncols = data.shape[2-1] - DENSITY_COLUMN + 1
nmesh = int32((emax - emin) / de) + 1
xmesh = np.linspace(emin,emax,nmesh)
ymesh = np.zeros((nmesh,ncols))
fact = de / (sigma * np.sqrt(2.0 * pi))
coloffset = 1
for i in np.arange(1,ncols+1).reshape(-1):
    ncol = data.shape[2-1] - DENSITY_COLUMN
    #     for idx = 1:length(nmesh)
#         func = exp(-(xmesh(idx)-data(:, EIGENVALUE_COLUMN)).^2/(2.0*sigma^2))*fact;
#                     ymesh(idx, coloffset:(coloffset+ncol)) = func.dot(data(:, DENSITY_COLUMN:end));
    g_x = fact * np.exp(- ((np.transpose(xmesh) ** 2) / (2 * sigma ** 2)))
    g_y = fact * np.exp(- ((ymesh ** 2) / (2 * sigma ** 2)))
    g_xy = np.multiply(g_x,g_y)
    ymesh[idx,i] = np.dot(func,data(:,DENSITY_COLUMN + i - 1))
    #     end

zfilt = gaussfilt(data(:,2),np.sum(data(:,np.arange(4,end()+1)), 2-1),0.1)
xmesh = (xmesh - 0.254735) * 27.211384
fclose(fid)