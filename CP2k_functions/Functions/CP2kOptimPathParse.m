function CP2kOptimPathParse(basefldr,directory, flname)

% Matt Darby. 21-12-2020. Imperial College London.
% Function that parses a CP2K .xyz file output from basefldr
% Version 3.0
% clear all;
% close all;
% clc

% basefldr = 'G:\Imperial\MattProjects\Pt_Clean\';
% directory = 'CP_Like_1012_Fluoride';
% flname = ['Sample' num2str(startConf) '_' num2str(endConf) '.xyz'];
disp('Creating xtd file...');
foldback = 0.02;

AtomCur = [];
AtomPosImages = [];

ABC = getABCvectors(basefldr, directory);
Vec = diag(ABC);

if ~isempty(length(basefldr))
    if ~strcmp(basefldr(end),'\')
        basefldr = [basefldr '\'];
    end
end

fidxdatcar = fopen([basefldr directory '\' flname]);

eofstat = false;

% First line is the number of atoms and is the same for all new sections
textLine = fgetl(fidxdatcar); eofstat = feof(fidxdatcar);
n = str2num(textLine);

% Second line contains the iteration (i), time, and total energy (E)
SectionLine = fgetl(fidxdatcar); eofstat = feof(fidxdatcar);
iconf = 1;

while ~eofstat
    
    % loop over number of atoms (n), extract names of atoms (AtomCur) for
    % config 1 only, extract xyz coordinates (xyzpos) and amalgamate for
    % all iterations (AtomPosImages)
    for i = 1:n
        textLine = fgetl(fidxdatcar); eofstat = feof(fidxdatcar);
        
        % Read in atoms names
        if iconf == 1
            AtomCur{1,i} = textLine(isstrprop(textLine,'alpha'));
        end
        
        % Read the xyz coordinates and wrap into unit cell
        xyzpos = str2num(textLine(5:end));
        
        
        
        AtomPosImages(i, :, iconf) = xyzpos(:);
        
    end
    
    if ~eofstat
        textLine = fgetl(fidxdatcar); eofstat = feof(fidxdatcar);
        textLine = fgetl(fidxdatcar); eofstat = feof(fidxdatcar);
    end
    
    % progress the iteration counter if a line contains " i =" flag
    if contains(textLine,SectionLine(1:4))
        iconf  = iconf+1;
    end
end

[AtomCur,sortIdx] = sort(AtomCur);

for ii = 1:iconf
    AtomTemp = AtomPosImages(:, :, ii);
    AtomPosImages(:, :, ii) = AtomTemp(sortIdx,:);
end

Snapshots = 1:iconf;
AtomPosImages = AtomPosImages-ABC.*(AtomPosImages>ABC);
AtomPosImages = AtomPosImages+ABC.*(AtomPosImages<0);

fclose(fidxdatcar);

if contains(flname, '.xyz')
    flname = erase(flname, '.xyz');
end

XTDFileName = [basefldr directory '\' flname '.xtd'];
% XTDFileName = [basefldr directory '\Sample' num2str(startConf) '_' num2str(endConf) '.xtd'];
XTDFileWriteCP2k(XTDFileName,n,AtomCur,AtomPosImages,Vec,5, Snapshots)


