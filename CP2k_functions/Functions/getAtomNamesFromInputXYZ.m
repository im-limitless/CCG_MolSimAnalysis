function [Atoms, AtomList, AtomIndx] = getAtomNamesFromInputXYZ(BaseFldr, system)

fid  = fopen([BaseFldr system '\' system '.xyz']);
lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
nAtoms = str2num(lines{1});
relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
nConfigs = length(relevant);

PosLines =[];
for i = relevant(end)+1:relevant(end)+nAtoms
    PosLines = [PosLines; lines{i}];
end

Atoms = strtrim(PosLines(:,1:5));
AtomList = unique(Atoms,'rows');
for i = 1:size(AtomList,1)
    AtomIndx.(strtrim(AtomList(i,:))) = find(sum(Atoms == AtomList(i,:),2) == length(AtomList(i,:)));
end


return