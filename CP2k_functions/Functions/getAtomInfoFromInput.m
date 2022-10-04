function [Atoms, AtomList, AtomIndx, AtomIndxfns, Kinds, Elements, PP] = getAtomInfoFromInput(BaseFldr, system)

PathSep = setOSpathSep;

fid  = fopen([BaseFldr system PathSep system '.xyz']);
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

AtomIndxfns = fieldnames(AtomIndx);

% now read the input or restart file to find the elements


if exist([BaseFldr system PathSep system '-1.restart'], 'file')
    fid  = fopen([BaseFldr system PathSep system '-1.restart']);
elseif exist([BaseFldr system PathSep system '.inp'], 'file')
    fid  = fopen([BaseFldr system PathSep system '.inp']);
else
    warning('No ".inp" or "-1.restart" file detected, terminating...');
    return
end

lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
KindLines =  find(~cellfun(@isempty,strfind(lines,'&KIND')));
EndKindLines =  find(~cellfun(@isempty,strfind(lines,'&END KIND')));



if length(KindLines) == length(EndKindLines)
    KindList = {};
    ElementList = {};
    PP = zeros(length(KindLines),1);
    for i = 1:length(KindLines)
        KindTemp = strsplit(lines{KindLines(i)}, 'KIND');
        KindList{i} = pad(strtrim(KindTemp{end}),size(AtomList,2));
        ElementList{i} = [];
        for j = 1:(EndKindLines(i)-KindLines(i))
            if contains(lines{KindLines(i)+j}, 'ELEMENT')
                ElemTemp = strsplit(lines{KindLines(i)+j}, 'ELEMENT');
                ElementList{i} = pad(strtrim(ElemTemp{end}),size(AtomList,2));
            end
            if contains(lines{KindLines(i)+j}, 'POTENTIAL') && ~contains(lines{KindLines(i)+j}, '&') && ~contains(lines{KindLines(i)+j}, 'END')
                PPTemp = strsplit(lines{KindLines(i)+j}, '-q');
                PP(i) = str2num(PPTemp{end});
            end
        end

        if isempty(ElementList{i})
            ElementList{i} = KindList{i};
        end
    end


else
    warning('Kind section in input file does not have corresponding end...')
end

idx=[];
for ii = 1:length(AtomList)
    if any(strcmp(KindList, AtomList(ii,:)))
        idx = [idx; find(strcmp(KindList, AtomList(ii,:)))];
    end
end

Kinds = cell2mat({KindList{idx}}'); Elements = cell2mat({ElementList{idx}}'); PP = PP(idx);

return