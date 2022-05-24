function [Vectors] = getABCvectors(ParentDir, WorkFldr)

if exist([ParentDir WorkFldr '\' WorkFldr '-1.restart'], 'file')
    fid  = fopen([ParentDir WorkFldr '\' WorkFldr '-1.restart']);
elseif exist([ParentDir WorkFldr '\' WorkFldr '.inp'], 'file')
    fid  = fopen([ParentDir WorkFldr '\' WorkFldr '.inp']);
elseif exist([ParentDir WorkFldr '\' WorkFldr '.pdb'], 'file')
    fid  = fopen([ParentDir WorkFldr '\' WorkFldr '.pdb']);
    lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
    fclose(fid);
    lines = lines{1};
    VectorLine = str2num(lines{1}(7:end-15));
    Vectors = VectorLine(1:3);
    return
else
    warning('No input or restart file detected...')
    return
end

lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
fclose(fid);
lines = lines{1};
relevant =  find(~cellfun(@isempty,strfind(lines,'ABC')));

if ~isempty(relevant)
    VectorLine = lines{relevant};
    vecIndx = find(VectorLine == 'C');
    Vectors = str2num(VectorLine(vecIndx+1:end));
    return
else
    ALineIndx = find(~cellfun(@isempty,strfind(lines,' A ')));
    BLineIndx = find(~cellfun(@isempty,strfind(lines,' B ')));
    CLineIndx = find(~cellfun(@isempty,strfind(lines,' C ')));
    
    
    AVectorLine = lines{ALineIndx};
    vecStartA = find(AVectorLine == 'A');
    A = str2num(AVectorLine(vecStartA+1:end));
    
    
    BVectorLine = lines{BLineIndx};
    vecStartB = find(BVectorLine == 'B');
    B = str2num(BVectorLine(vecStartB+1:end));
    
    CVectorLine = lines{CLineIndx};
    vecStartC = find(CVectorLine == 'C');
    C = str2num(CVectorLine(vecStartC+1:end));
    
    if norm(A) == sum(A) & norm(B) == sum(B) & norm(C) == sum(C)
        Vectors = diag([A; B; C])';
    else
        warning('Lattice vectors are not orthogonal...')
    end
    return
end