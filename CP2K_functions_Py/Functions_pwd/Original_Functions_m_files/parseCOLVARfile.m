clear all;
close all;

% parse the COLVAR file from cp2k

Basefldr = 'G:\Imperial\MattProjects\Pt_Clean\Metadynamics\\';
Allfldrs = dir([Basefldr '*Meta_1010_Fluoride']);

flname = dir([Basefldr Allfldrs.name '\*-COLVAR.metadynLog']);
ColvarFile = [Basefldr Allfldrs.name '\' flname.name];

fid = fopen(ColvarFile);

eofstat = false;

% First line is the number of atoms and is the same for all new sections
HeadLine = fgetl(fid); eofstat = feof(fid);
i = 0;
data = [];

while ~eofstat
    i=i+1;
    textLine = fgetl(fid); eofstat = feof(fid);
    data = [data; str2num(textLine)];
end

fclose(fid);

figure
hold on

plot(data(:,1)-data(1,1), data(:,2:11));
% plot(data(:,1)-data(1,1), abs(data(:,2:11)));