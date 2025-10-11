clear all; close all; clc;

fid = fopen('C:\Users\Matt\Desktop\Debug\PDOS_Al_Bulk_cutoff200\pdos');

data = fscanf(fid, '%f',  [12 inf]);
data=data';

sigma = 0.02;
de = 0.001;
scale = 1;
EIGENVALUE_COLUMN = 2;
DENSITY_COLUMN = 4;



margin = 10*sigma;
emin = min(min(data(:, EIGENVALUE_COLUMN))) - margin;
emax = max(max(data(:, EIGENVALUE_COLUMN)))+ margin;
ncols = size(data,2) - DENSITY_COLUMN + 1;
nmesh = int32((emax-emin)/de)+1;
xmesh = linspace(emin, emax, nmesh);
ymesh = zeros(nmesh, ncols);

fact = de/(sigma*sqrt(2.0*pi));

coloffset = 1;
for i = 1:ncols
ncol = size(data,2) - DENSITY_COLUMN;
%     for idx = 1:length(nmesh)
%         func = exp(-(xmesh(idx)-data(:, EIGENVALUE_COLUMN)).^2/(2.0*sigma^2))*fact;
%                     ymesh(idx, coloffset:(coloffset+ncol)) = func.dot(data(:, DENSITY_COLUMN:end));
        g_x = fact*exp(-((xmesh'.^2)/(2*sigma^2)));
        g_y = fact*exp(-((ymesh.^2)/(2*sigma^2)));
        
        g_xy = g_x.*g_y;
        ymesh(idx, i) = dot(func,data(:, DENSITY_COLUMN+i-1));
%     end
end

[ zfilt ] = gaussfilt(data(:,2),sum(data(:,4:end),2),0.1);

xmesh=(xmesh-0.254735)*27.211384;


fclose(fid);