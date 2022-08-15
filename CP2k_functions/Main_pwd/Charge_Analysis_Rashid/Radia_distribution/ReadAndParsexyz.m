function [xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, Step] = ReadAndParsexyz(Base, Fldr, Traj, ABC)

% Base = 'G:\Imperial\MattProjects\MD_files\Pit_NVT\CP_Pit1822\';
% Fldr = 'SS1';
% Traj = 'SS1.xyz';

fid  = fopen([Base Fldr '/' Traj]);
disp('Reading xyz data...');
lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
nAtoms = str2num(lines{1});
relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
nConfigs = length(relevant);

Step = [];
Atoms = [];

for i = 1:nAtoms
    Atoms = [Atoms; lines{(relevant(1)+i)}(1:5)];
end

AtomList = unique(Atoms,'rows');
XYZ = zeros(nConfigs, nAtoms, 3);

for i = 1:size(AtomList,1)
    Indx.(strtrim(AtomList(i,:))) =  find(sum(Atoms == AtomList(i,:),2) == length(AtomList(i,:)));
    xyz.(strtrim(AtomList(i,:))) = zeros(nConfigs, length(Indx.(strtrim(AtomList(i,:)))), 3);
end

startConfig = 1;

for j = startConfig:nConfigs
    
    disp(['Processing configuration number ' num2str(j) ' of ' num2str(nConfigs)]);
    for k = 1:nAtoms
        XYZ(j,k, 1:3) =  str2num(lines{relevant(j)+k}(5:end));
    end
    
    XYZ = wrapXYZ(XYZ, ABC);
    
    for k = 1:size(AtomList,1)
        xyz.(strtrim(AtomList(k,:)))(j,:, 1:3) = XYZ(j, Indx.(strtrim(AtomList(k,:))),1:3);
    end
    if nConfigs > 1
        initialI = lines{relevant(j)};
        EqI = strfind(initialI, '=');
        ComI = strfind(initialI, ',');
        posI = [EqI(1) ComI(1)];
        Step = [Step; str2num(initialI(posI(1)+1:posI(2)-1))];
    end
end

return