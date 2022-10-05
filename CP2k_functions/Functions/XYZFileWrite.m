function XYZFileWrite(basefldr,n,AtomCur,AtomPos,Vec,name)

pathSec =setOSpathSep;

% Michail Stamatakis 16-Nov-2011. University of Delaware.
% Function that writes a VASP configuration to basefldr
% Version 3.0
%
% Usage:
% VASPConfigWrite(basefldr,potenfldr,n,AtomCur,AtomPos,Vec,fixedatoms,writepotcar,commentstr)
% - basefldr is the directory where the VASP input will be written
% - potenfldr is the directory containing the pseudopotentials
% - n(:) is a vector containing the number of atoms per element
% - AtomCur{:} a list of the elements (must match the order in n)
% - AtomPos(i1,i2,i3) is the x,y,z (for i3 = 1:3) positions of the atom i2 of element i1
% - Vec(:,i1) are the a,b,c (for i1 = 1:3) vectors defining the unit cell
% - fixedatoms(:) is a vector containing the numbers of the fixed atoms
% - writepotcar is a logical variable controlling whether a POTCAR file is written
% - commentstr is a string describing the system that will be written in the first line of POSCAR

% newlinechar = [char(13) char(10)];
newlinechar = char(10);

if ~isempty(length(basefldr)) 
    if strcmp(basefldr(end),'\')
        basefldr = [basefldr '\'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write XYZ

fidout = fopen([basefldr 'input.xyz'],'w');

fprintf(fidout,[num2str(n) newlinechar]);
fprintf(fidout,['ABC = ' num2str(Vec(1,1)) ' ' num2str(Vec(2,2)) ' ' num2str(Vec(3,3)) newlinechar]);
for k = 1:length(n)
    [sortAtomPos,indxs] = sort(AtomPos(k,:,3));
    
    for j = 1:length(sortAtomPos)
        s1 = num2str(AtomPos(k,indxs(j),1)*Vec(1,1),'%20.16f'); n1 = 20-length(s1);
        s2 = num2str(AtomPos(k,indxs(j),2)*Vec(2,2),'%20.16f'); n2 = 20-length(s2);
        s3 = num2str(AtomPos(k,indxs(j),3)*Vec(3,3),'%20.16f'); n3 = 20-length(s3);
        fprintf(fidout,[AtomCur{k} ' ' repmat(' ',1,n1) s1 repmat(' ',1,n2) s2 repmat(' ',1,n3) s3 newlinechar]);
    end
end
fclose(fidout);

end