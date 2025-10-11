function [PosLines, VelLines, VectorLine, FixedLine] = getLastPosAndVel(ParentDir, WorkFldr, findPos, findVel)

% split .xyz file and extract specific (final) geometry
% ParentDir = 'G:\Imperial\MattProjects\Pt_Clean\';
% WorkFldr = 'BOMD_NVT_1010_Clean';
if findPos
    fid  = fopen([ParentDir WorkFldr '\' WorkFldr '-pos-1.xyz']);
    
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
    
end

if findVel
    fid  = fopen([ParentDir WorkFldr '\' WorkFldr '-vel-1.xyz']);
    
    lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
    fclose(fid);
    
    lines = lines{1};
    nAtoms = str2num(lines{1});
    relevant =  find(~cellfun(@isempty,strfind(lines,'i =')));
    nConfigs = length(relevant);
    
    VelLines =[];
    for i = relevant(end)+1:relevant(end)+nAtoms
        VelLines = [VelLines; lines{i}(5:end)];
    end
end
fid  = fopen([ParentDir WorkFldr '\' WorkFldr '.inp']);

lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
relevant =  find(~cellfun(@isempty,strfind(lines,'ABC')));

VectorLine = lines{relevant};

fid  = fopen([ParentDir WorkFldr '\' WorkFldr '.inp']);

lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);

lines = lines{1};
relevant =  find(~cellfun(@isempty,strfind(lines,' LIST ')));

FixedLine = lines{relevant};
return