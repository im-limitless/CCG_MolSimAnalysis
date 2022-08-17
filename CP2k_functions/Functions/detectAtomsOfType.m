function [Indx, myAtomList, myAtomIndxs] = detectAtomsOfType(myAtom, AtomList, Indx, Indxfns)

myAtomList = [];
myAtomIndxs = [];
Indx.myAtom = [];

% % Find indices of all "myAtom*" atoms
if size(AtomList,1) > 1
    for atoms = 1:size(myAtom, 2)
        for i = 1:size(AtomList,1)
            if contains(AtomList(i,:), myAtom{atoms})
                myAtomList = [myAtomList; i];
            end
        end
    end
    
    for ii = 1:length(myAtomList)
        Indx.myAtom = [Indx.myAtom; Indx.(Indxfns{myAtomList(ii)})];
        myAtomIndxs = [myAtomIndxs; Indx.(Indxfns{myAtomList(ii)})];
    end
else
    myAtomList = 1;
    Indx.myAtom = Indx.(Indxfns{1});
end

return