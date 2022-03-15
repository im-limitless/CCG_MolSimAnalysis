function [NAtoms,srtdAtomElement,srtdAtomPosition,n,AtomCur,AtomPos,Vec] = ...
    XSDFileRead(XSDFileName)

% Michail Stamatakis 16-Nov-2011. University of Delaware.
% Function that reads a Materials Studio xsd file
% Version 3.0

XSDdoc = xmlread(XSDFileName);

for k = 0:XSDdoc.getChildNodes.getLength-1
    if strcmpi(XSDdoc.getChildNodes.item(k).getNodeName,'xsd')
        kxsdroot = k; % don't break here, since <!DOCTYPE XSD []> is also a node with name XSD
    end
end

XSDRoot = XSDdoc.getChildNodes.item(kxsdroot);

for k = 0:XSDRoot.getChildNodes.getLength-1
    if strcmpi(XSDRoot.getChildNodes.item(k).getNodeName,'AtomisticTreeRoot')
        katomisticrreetoot = k;
        break
    end
end

AtomisticTreeRoot = XSDRoot.getChildNodes.item(katomisticrreetoot);
for i = 0:AtomisticTreeRoot.getChildNodes.getLength-1
    if ~isempty(AtomisticTreeRoot.getChildNodes.item(i))
        if strcmpi(char(AtomisticTreeRoot.getChildNodes.item(i).getNodeName),...
                'SymmetrySystem')
            IdentityMapping = AtomisticTreeRoot. ...
                getChildNodes.item(i).item(1).item(1).item(1);
        end
    end
end

NAtoms = 0;
AtomPosition = [];
AtomID = [];
AtomElement = {};
Vec = [];
for i = 0:IdentityMapping.getChildNodes.getLength-1
    
    if ~isempty(IdentityMapping.getChildNodes.item(i))
        
        if strcmpi(IdentityMapping.getChildNodes.item(i).getNodeName,'Atom3D')
            
            NAtoms = NAtoms + 1;
            AtomID(NAtoms) = 0;
%             disp([num2str(i) ' ' char(IdentityMapping.getChildNodes.item(i).getNodeName)]);
            
            for j = 0:IdentityMapping.getChildNodes.item(i).getAttributes.getLength-1
                
                if ~isempty(IdentityMapping.getChildNodes.item(i).getAttributes.item(j))
                    
                    if strcmpi(char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getName),'XYZ')
                        AtomPosition(NAtoms,1:3) = str2num(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getValue);
                    end
                    if strcmpi(char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getName),'UserID')
                        AtomID(NAtoms) = str2num(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getValue);
                    end
                    if strcmpi(char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getName),'Components')
                        AtomElement{NAtoms} = char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getValue);
                    end
                    
                end
            end
            
        end
        
        if strcmpi(IdentityMapping.getChildNodes.item(i).getNodeName,'SpaceGroup')
            
            for j = 0:IdentityMapping.getChildNodes.item(i).getAttributes.getLength-1
                
                if ~isempty(IdentityMapping.getChildNodes.item(i).getAttributes.item(j))
                    
                    if strcmpi(char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getName),'AVector')
                        Vec(1,:) = str2num(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getValue);
                    elseif strcmpi(char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getName),'BVector')
                        Vec(2,:) = str2num(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getValue);
                    elseif strcmpi(char(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getName),'CVector')
                        Vec(3,:) = str2num(IdentityMapping.getChildNodes.item(i).getAttributes.item(j).getValue);
                    end
                    
                end
            end
            
        end
        
    end
    
end

[srtdAtomID,indxs] = sort(AtomID);
srtdAtomPosition = AtomPosition(indxs,:);
srtdAtomElement = AtomElement(indxs);

k = 1;
AtomCur{k} = srtdAtomElement{1};
AtomPos(1,1,:) = srtdAtomPosition(1,:);
n(k) = 1;
for i = 2:length(srtdAtomElement)
    
    for j = 1:k
        if strcmp(srtdAtomElement{i},AtomCur{j});
            n(j) = n(j) + 1;
            AtomPos(j,n(j),:) = srtdAtomPosition(i,:);
            break;
        end
    end
    
    if ~strcmp(srtdAtomElement{i},AtomCur{j});
        k = k + 1;
        AtomCur{k} = srtdAtomElement{i};
        n(k) = 1;
        AtomPos(k,n(k),:) = srtdAtomPosition(i,:);
    end
    
end

srtdAtomPosition = [];
srtdAtomElement = {};
for k = 1:length(n)
    for j = 1:n(k)
        s1 = AtomPos(k,j,1);
        s2 = AtomPos(k,j,2);
        s3 = AtomPos(k,j,3);
        
        srtdAtomPosition(sum(n(1:k-1))+j,:) = [s1 s2 s3];
        srtdAtomElement{sum(n(1:k-1))+j} = AtomCur{k};
    end
end
