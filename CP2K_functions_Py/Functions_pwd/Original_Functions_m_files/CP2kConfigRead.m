% function [NAtoms,AtomElement,AtomPosition,n,AtomCur,AtomPos,AtomFreeToMove,Vec] = ...
%     CP2kConfigRead(basefldr,flname,foldback)

% Matt Darby - 17-Dec-2020. Imperial College London.
% Function that parses a CP2k .xyz file from basefldr
% flname = '*.xyz'
% Element information is automatically read from .inp if not contained in
% the configuration file
% Version 3.1

basefldr = 'Z:\Imperial\Federico\10Na10Cl\1ML\__dynamic16';
flname = 'ptwater-pos-1.xyz';
foldback = 0.02;

NAtoms = 0;
AtomElement = {};
AtomPosition = [];
n = [];
AtomCur = {};
AtomPos = [];
Vec = [];
AtomFreeToMove = [];

if ~isempty(length(basefldr)) 
    if ~strcmp(basefldr(end),'\')
        basefldr = [basefldr '\'];
    end
end

% potentialsfile = [basefldr 'POTCAR'];
% if ~exist(potentialsfile,'file')
%     potentialsfile = [basefldr '..\POTCAR'];
% end
configuratfile = [basefldr flname];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read atom positions and unit cell vectors from .xyz file

fidconf = fopen(configuratfile);

% First line is a comment
textLine = fgetl(fidconf); eofstat = feof(fidconf); 

% Second line gives the scaling factor for the unit cell
textLine = fgetl(fidconf); eofstat = feof(fidconf); 
ScalF = str2num(textLine);

% Subsequent three lines give the unit vectors
textLine = fgetl(fidconf); eofstat = feof(fidconf); 
va = str2num(textLine);
textLine = fgetl(fidconf); eofstat = feof(fidconf); 
vb = str2num(textLine);
textLine = fgetl(fidconf); eofstat = feof(fidconf); 
vc = str2num(textLine);

Vec = ScalF*[va; vb; vc];

textLine = fgetl(fidconf); eofstat = feof(fidconf);

if isempty(str2num(textLine)) % CONTCAR file format

    StrWords = textscan(textLine,'%s');
    AtomCur=StrWords{1}.';
    textLine = fgetl(fidconf); eofstat = feof(fidconf);

else
    
    % Read atom types from POTCAR file
    fidpot = fopen(potentialsfile);
    
    eofstatpot = false;
    
    while ~eofstatpot
        textLinepot = fgetl(fidpot); eofstatpot = feof(fidpot);
        StrWordspot = textscan(textLinepot,'%s');
        
        if contains(StrWordspot{1}{2}, '_')
            newstr = extractBefore(StrWordspot{1}{2}, "_");
        else
            newstr = StrWordspot{1}{2};
        end
            
%         AtomCur = {AtomCur{:} StrWordspot{1}{2}};
        AtomCur = {AtomCur{:} newstr};
        
        while ~eofstatpot
            textLinepot = fgetl(fidpot); eofstatpot = feof(fidpot);
            if eofstatpot || ~isempty(strfind(textLinepot,'End of Dataset'))
                break
            end
        end
    end
    
    fclose(fidpot);
    
end

% Next line gives the number of atoms for each element present in the
% system
n = str2num(textLine);
NAtoms = sum(n);
for k = 1:length(n)
    for i = 1:n(k)
        AtomElement{sum(n(1:k-1))+i} = AtomCur{k};
    end
end

% We skip the next two lines
textLine = fgetl(fidconf); eofstat = feof(fidconf); 
textLine = fgetl(fidconf); eofstat = feof(fidconf); 

% We start parsing the atom positions
AtomPosition = zeros(NAtoms,3);
AtomFreeToMove = zeros(NAtoms,3);
jAtm = 0;
while ~eofstat
    
    textLine = fgetl(fidconf); eofstat = feof(fidconf);

    StrWords = textscan(textLine,'%s');
    
    jAtm = jAtm + 1;
    AtomPosition(jAtm,1) = str2num(StrWords{1}{1});
    AtomPosition(jAtm,2) = str2num(StrWords{1}{2});
    AtomPosition(jAtm,3) = str2num(StrWords{1}{3});
    
    if exist('foldback','var') && ~isempty(foldback)
        if length(foldback) == 1;
            foldback = [1 1 1]*foldback;
        end
        for k = 1:3
            if foldback(k) > 0
                if AtomPosition(jAtm,k) > 1-foldback(k)
                    AtomPosition(jAtm,k) = AtomPosition(jAtm,k) - 1;
                end
            end
        end
    end
    
    for q = 1:3
        if strcmp(StrWords{1}{q+3},'T')
            AtomFreeToMove(jAtm,q) = true;
        elseif strcmp(StrWords{1}{q+3},'F')
            AtomFreeToMove(jAtm,q) = false;
        end
    end
    
    if jAtm == NAtoms
        break
    end
    
end

fclose(fidconf);

for k = 1:length(n)
    for i = 1:n(k)
        AtomPos(k,i,:) = AtomPosition(sum(n(1:k-1))+i,:);
    end
end

% end