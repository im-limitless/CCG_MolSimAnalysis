% NameChangeXYZ
clear all

Base = 'G:\Imperial\MattProjects\MD_files\Pit_NVT\';
system = 'CP_Pit2020_H_1.00ML';

BaseFldr = [Base system '\'];
AllFldrs = dir([BaseFldr 'SS1*']);
mkdir([BaseFldr 'xtd']);

for j = 1:length(AllFldrs)
    [ABC] = getABCvectors(BaseFldr, AllFldrs(j).name);
    fid  = fopen([BaseFldr AllFldrs(j).name '\' AllFldrs(j).name '.xyz']);
    disp('Reading xyz data...');
    lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
    fclose(fid);
    
    lines = lines{1};
    lines{2} = 'i = 0';
    nAtoms = str2num(lines{1});
    relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
    nConfigs = length(relevant);
    
    Step = [];
    Atoms = [];
    
    for i = 1:nAtoms
        if contains(lines{i+2}, 'PTP')
            lines{i+2} =  strrep(lines{i+2}, 'PTP', 'Pt ');
        elseif contains(lines{i+2}, 'HW1')
            lines{i+2} =  strrep(lines{i+2}, 'HW1', 'H  ');
        elseif contains(lines{i+2}, 'HW2')
            lines{i+2} =  strrep(lines{i+2}, 'HW2', 'H  ');
        elseif contains(lines{i+2}, 'OW')
            lines{i+2} =  strrep(lines{i+2}, 'OW', 'O ');
        elseif contains(lines{i+2}, 'H1')
            lines{i+2} =  strrep(lines{i+2}, 'H1', 'H ');
        elseif contains(lines{i+2}, 'H2')
            lines{i+2} =  strrep(lines{i+2}, 'H2', 'H ');
        elseif contains(lines{i+2}, 'H3')
            lines{i+2} =  strrep(lines{i+2}, 'H3', 'H ');
        elseif contains(lines{i+2}, 'HS')
            lines{i+2} =  strrep(lines{i+2}, 'HS', 'B ');
        end
    end
    
    
    for i = 1:nAtoms
        Atoms = [Atoms; lines{(relevant(1)+i)}(1:5)];
    end
    
    AtomList = unique(Atoms,'rows');
    XYZ = zeros(nConfigs, nAtoms, 3);
    
    for i = 1:size(AtomList,1)
        Indx.(strtrim(AtomList(i,:))) =  find(sum(Atoms == AtomList(i,:),2) == length(AtomList(i,:)));
        xyz.(strtrim(AtomList(i,:))) = zeros(nConfigs, length(Indx.(strtrim(AtomList(i,:)))), 3);
    end
    
    fid = fopen([BaseFldr AllFldrs(j).name '\' AllFldrs(j).name '_new.xyz'], 'w');
    for i =  1:length(lines)
        fprintf(fid, [lines{i} newline]);
    end
    fclose(fid);
    
    CP2kOptimPathParse(BaseFldr,AllFldrs(j).name,[AllFldrs(j).name '_new.xyz'])
    copyfile([BaseFldr AllFldrs(j).name '\' AllFldrs(j).name '_new.xtd'], [BaseFldr 'xtd\' system '_' num2str(j) '.xtd']);
    copyfile([BaseFldr AllFldrs(j).name '\' AllFldrs(j).name '_new.arc'], [BaseFldr 'xtd\' system '_' num2str(j) '.arc']);
end
