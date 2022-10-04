function [Indx, myAtomList, myAtomIndxs, myAtomCore] = detectAtomsOfType(myAtom, AtomList, Indx, Indxfns, Kinds, Elements, PP)

myAtomList = [];
myAtomIndxs = [];
Indx.myAtom = [];
myAtomCore = zeros(size(myAtom, 2),1);
% % Find indices of all "myAtom*" atoms
if strcmp(myAtom, 'all')
    myAtomList = [1:length(AtomList)]';
else
    if size(AtomList,1) > 1
        for atoms = 1:size(myAtom, 2)
            for i = 1:length(Elements)
                if strtrim(Elements(i,:)) == myAtom{atoms}
                    for j = 1:length(AtomList)
                        if strcmp(strtrim(AtomList(j,:)), strtrim(Kinds(i,:)))
                            myAtomList = [myAtomList; j];
                            myAtomCore(atoms,1) = PP(i);
                        end
                    end
                end
            end
        end
    else
        warning('myAtoms not found in element list, selecting all atoms');
        myAtomList = 1:length(AtomList);
    end
end

for ii = 1:length(myAtomList)
    Indx.myAtom = [Indx.myAtom; Indx.(Indxfns{myAtomList(ii)})];
    myAtomIndxs = [myAtomIndxs; Indx.(Indxfns{myAtomList(ii)})];
end

return