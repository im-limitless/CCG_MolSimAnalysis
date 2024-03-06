function XSDFileWriteCP2K(XSDFileName,n,AtomCur,AtomPos,Vec)

% Michail Stamatakis 16-Nov-2011. University of Delaware.
% Function that writes a Materials Studio xsd file
% Version 3.0

NAtoms = sum(n);

% ID 1 -> reserved for AtomisticTreeRoot element
% ID 2 -> reserved for SymmetrySystem element
% ID 3 -> reserved for InfiniteMapping element
% IDs 4~NAtoms+3 -> used for each atom present in the system
% IDs NAtoms+4 -> reserved for MappingSet element
% IDs NAtoms+5 -> reserved for SpaceGroup element
% IDs NAtoms+6 -> reserved for MappingFamily element
% IDs NAtoms+7 -> reserved for IdentityMapping element
% IDs NAtoms+8 -> reserved for ReciprocalLattice3D element

docNode = com.mathworks.xml.XMLUtils.createDocument('XSD');

docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('Version','4.0');

AtomisticTreeRootElement = docNode.createElement('AtomisticTreeRoot'); 
AtomisticTreeRootElement.setAttribute('ID','1');
AtomisticTreeRootElement.setAttribute('NumProperties','40');
AtomisticTreeRootElement.setAttribute('NumChildren','1');

Property1 = docNode.createElement('Property');
Property1.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property1.setAttribute('Name','AngleEnergy');
Property1.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property1);
 
Property2 = docNode.createElement('Property');
Property2.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property2.setAttribute('Name','BendBendEnergy');
Property2.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property2);
 
Property3 = docNode.createElement('Property');
Property3.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property3.setAttribute('Name','BendTorsionBendEnergy');
Property3.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property3);
 
Property4 = docNode.createElement('Property');
Property4.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property4.setAttribute('Name','BondEnergy');
Property4.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property4);
 
Property5 = docNode.createElement('Property');
Property5.setAttribute('DefinedOn','Atom');
Property5.setAttribute('Name','EFGAsymmetry');
Property5.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property5);
 
Property6 = docNode.createElement('Property');
Property6.setAttribute('DefinedOn','Atom');
Property6.setAttribute('Name','EFGQuadrupolarCoupling');
Property6.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property6);
 
Property7 = docNode.createElement('Property');
Property7.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property7.setAttribute('Name','ElectrostaticEnergy');
Property7.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property7);
 
Property8 = docNode.createElement('Property');
Property8.setAttribute('DefinedOn','GrowthFace');
Property8.setAttribute('Name','FaceMillerIndex');
Property8.setAttribute('Type','MillerIndex');
AtomisticTreeRootElement.appendChild(Property8);
 
Property9 = docNode.createElement('Property');
Property9.setAttribute('DefinedOn','GrowthFace');
Property9.setAttribute('Name','FacetTransparency');
Property9.setAttribute('Type','Float');
AtomisticTreeRootElement.appendChild(Property9);
 
Property10 = docNode.createElement('Property');
Property10.setAttribute('DefinedOn','Bondable');
Property10.setAttribute('Name','Force');
Property10.setAttribute('Type','CoDirection');
AtomisticTreeRootElement.appendChild(Property10);
 
Property11 = docNode.createElement('Property');
Property11.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property11.setAttribute('Name','HydrogenBondEnergy');
Property11.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property11);
 
Property12 = docNode.createElement('Property');
Property12.setAttribute('DefinedOn','Bondable');
Property12.setAttribute('Name','ImportOrder');
Property12.setAttribute('Type','UnsignedInteger');
AtomisticTreeRootElement.appendChild(Property12);
 
Property13 = docNode.createElement('Property');
Property13.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property13.setAttribute('Name','InversionEnergy');
Property13.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property13);
 
Property14 = docNode.createElement('Property');
Property14.setAttribute('DefinedOn','Atom');
Property14.setAttribute('Name','IsBackboneAtom');
Property14.setAttribute('Type','Boolean');
AtomisticTreeRootElement.appendChild(Property14);
 
Property15 = docNode.createElement('Property');
Property15.setAttribute('DefinedOn','Atom');
Property15.setAttribute('Name','IsChiralCenter');
Property15.setAttribute('Type','Boolean');
AtomisticTreeRootElement.appendChild(Property15);
 
Property16 = docNode.createElement('Property');
Property16.setAttribute('DefinedOn','Atom');
Property16.setAttribute('Name','IsOutOfPlane');
Property16.setAttribute('Type','Boolean');
AtomisticTreeRootElement.appendChild(Property16);
 
Property17 = docNode.createElement('Property');
Property17.setAttribute('DefinedOn','BestFitLineMonitor');
Property17.setAttribute('Name','LineExtentPadding');
Property17.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property17);
 
Property18 = docNode.createElement('Property');
Property18.setAttribute('DefinedOn','Linkage');
Property18.setAttribute('Name','LinkageGroupName');
Property18.setAttribute('Type','String');
AtomisticTreeRootElement.appendChild(Property18);
 
Property19 = docNode.createElement('Property');
Property19.setAttribute('DefinedOn','PropertyList');
Property19.setAttribute('Name','ListIdentifier');
Property19.setAttribute('Type','String');
AtomisticTreeRootElement.appendChild(Property19);
 
Property20 = docNode.createElement('Property');
Property20.setAttribute('DefinedOn','Atom');
Property20.setAttribute('Name','NMRShielding');
Property20.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property20);
 
Property21 = docNode.createElement('Property');
Property21.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property21.setAttribute('Name','NonBondEnergy');
Property21.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property21);
 
Property22 = docNode.createElement('Property');
Property22.setAttribute('DefinedOn','Bondable');
Property22.setAttribute('Name','NormalMode');
Property22.setAttribute('Type','Direction');
AtomisticTreeRootElement.appendChild(Property22);
 
Property23 = docNode.createElement('Property');
Property23.setAttribute('DefinedOn','Bondable');
Property23.setAttribute('Name','NormalModeFrequency');
Property23.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property23);
 
Property24 = docNode.createElement('Property');
Property24.setAttribute('DefinedOn','Bondable');
Property24.setAttribute('Name','OrbitalCutoffRadius');
Property24.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property24);
 
Property25 = docNode.createElement('Property');
Property25.setAttribute('DefinedOn','BestFitPlaneMonitor');
Property25.setAttribute('Name','PlaneExtentPadding');
Property25.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property25);
 
Property26 = docNode.createElement('Property');
Property26.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property26.setAttribute('Name','PotentialEnergy');
Property26.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property26);
 
Property27 = docNode.createElement('Property');
Property27.setAttribute('DefinedOn','ScalarFieldBase');
Property27.setAttribute('Name','QuantizationValue');
Property27.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property27);
 
Property28 = docNode.createElement('Property');
Property28.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property28.setAttribute('Name','RestraintEnergy');
Property28.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property28);
 
Property29 = docNode.createElement('Property');
Property29.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property29.setAttribute('Name','SeparatedStretchStretchEnergy');
Property29.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property29);
 
Property30 = docNode.createElement('Property');
Property30.setAttribute('DefinedOn','Trajectory');
Property30.setAttribute('Name','SimulationStep');
Property30.setAttribute('Type','Integer');
AtomisticTreeRootElement.appendChild(Property30);
 
Property31 = docNode.createElement('Property');
Property31.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property31.setAttribute('Name','StretchBendStretchEnergy');
Property31.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property31);
 
Property32 = docNode.createElement('Property');
Property32.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property32.setAttribute('Name','StretchStretchEnergy');
Property32.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property32);
 
Property33 = docNode.createElement('Property');
Property33.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property33.setAttribute('Name','StretchTorsionStretchEnergy');
Property33.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property33);
 
Property34 = docNode.createElement('Property');
Property34.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property34.setAttribute('Name','TorsionBendBendEnergy');
Property34.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property34);
 
Property35 = docNode.createElement('Property');
Property35.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property35.setAttribute('Name','TorsionEnergy');
Property35.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property35);
 
Property36 = docNode.createElement('Property');
Property36.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property36.setAttribute('Name','TorsionStretchEnergy');
Property36.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property36);
 
Property37 = docNode.createElement('Property');
Property37.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property37.setAttribute('Name','ValenceCrossTermEnergy');
Property37.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property37);
 
Property38 = docNode.createElement('Property');
Property38.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property38.setAttribute('Name','ValenceDiagonalEnergy');
Property38.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property38);
 
Property39 = docNode.createElement('Property');
Property39.setAttribute('DefinedOn','ClassicalEnergyHolder');
Property39.setAttribute('Name','VanDerWaalsEnergy');
Property39.setAttribute('Type','Double');
AtomisticTreeRootElement.appendChild(Property39);
 
Property40 = docNode.createElement('Property');
Property40.setAttribute('DefinedOn','SymmetrySystem');
Property40.setAttribute('Name','_Stress');
Property40.setAttribute('Type','Matrix');
AtomisticTreeRootElement.appendChild(Property40);

SymmSys = docNode.createElement('SymmetrySystem');
SymmSys.setAttribute('ID','2');
SymmSys.setAttribute('Mapping','3');
tmpstr = sprintf('%1.0f,',[4:NAtoms+3 NAtoms+5]);
SymmSys.setAttribute('Children',tmpstr(1:end-1));
SymmSys.setAttribute('Normalized','1');
SymmSys.setAttribute('Name','SymmSys');
SymmSys.setAttribute('UserID',num2str(NAtoms+18));
SymmSys.setAttribute('XYZ','0.00000000000000,0.00000000000000,0.000000000000000');
SymmSys.setAttribute('OverspecificationTolerance','0.05');
SymmSys.setAttribute('PeriodicDisplayType','Original');

MappngSet = docNode.createElement('MappingSet');
MappngSet.setAttribute('ID',num2str(NAtoms+4));
MappngSet.setAttribute('SymmetryDefinition',num2str(NAtoms+5));
MappngSet.setAttribute('ActiveSystem','2');
MappngSet.setAttribute('NumFamilies','1');
MappngSet.setAttribute('OwnsTotalConstraintMapping','1');
MappngSet.setAttribute('TotalConstraintMapping','3');

MappngFamily = docNode.createElement('MappingFamily');
MappngFamily.setAttribute('ID',num2str(NAtoms+6));
MappngFamily.setAttribute('NumImageMappings','0');

IdentMappng = docNode.createElement('IdentityMapping');
IdentMappng.setAttribute('ID',num2str(NAtoms+7));
IdentMappng.setAttribute('Element','1,0,0,0,0,1,0,0,0,0,1,0');
IdentMappng.setAttribute('Constraint','1,0,0,0,0,1,0,0,0,0,1,0');
tmpstr = sprintf('%1.0f,',[4:NAtoms+3]);
IdentMappng.setAttribute('MappedObjects',tmpstr(1:end-1));
tmpstr = sprintf('%1.0f,',[NAtoms+5 NAtoms+8]);
IdentMappng.setAttribute('DefectObjects',tmpstr(1:end-1));
IdentMappng.setAttribute('NumImages',num2str(NAtoms));
IdentMappng.setAttribute('NumDefects','2');

MappngRepairs = docNode.createElement('MappingRepairs');
MappngRepairs.setAttribute('NumRepairs','0');

MappngFamily.appendChild(IdentMappng);
MappngFamily.appendChild(MappngRepairs);
MappngSet.appendChild(MappngFamily);
SymmSys.appendChild(MappngSet);
AtomisticTreeRootElement.appendChild(SymmSys);


for k = 1:length(n)
    for i = 1:n(k)
        
        jAtom = sum(n(1:k-1))+i;
        NewAtom = docNode.createElement('Atom3d');
        NewAtom.setAttribute('ID',num2str(jAtom+3));
        NewAtom.setAttribute('Mapping',num2str(NAtoms+7));
        NewAtom.setAttribute('Parent','2');
        NewAtom.setAttribute('Name',[AtomCur{i} num2str(jAtom)]);
        NewAtom.setAttribute('UserID',num2str(jAtom));
        NewAtom.setAttribute('DisplayStyle',CPK_or_BnS(AtomCur{i}));
        tmpstr = sprintf('%10.16f,',AtomPos(k,i,:));
        NewAtom.setAttribute('XYZ',tmpstr(1:end-1));
        NewAtom.setAttribute('Components',AtomCur{i});
%         NewAtom.setAttribute('FormalCharge','0/1');
        IdentMappng.appendChild(NewAtom);
        
    end
end


SpaceGrp = docNode.createElement('SpaceGroup');
SpaceGrp.setAttribute('ID',num2str(NAtoms+5));
SpaceGrp.setAttribute('Parent','2');
SpaceGrp.setAttribute('Children',num2str(NAtoms+8));
SpaceGrp.setAttribute('DisplayStyle','Solid');
SpaceGrp.setAttribute('XYZ','0.00,0.00,0.00');
SpaceGrp.setAttribute('Color','0,0,0,255');
tmpstr = sprintf('%10.16f,',Vec(1,:));
SpaceGrp.setAttribute('AVector',tmpstr(1:end-1));
tmpstr = sprintf('%10.16f,',Vec(2,:));
SpaceGrp.setAttribute('BVector',tmpstr(1:end-1));
tmpstr = sprintf('%10.16f,',Vec(3,:));
SpaceGrp.setAttribute('CVector',tmpstr(1:end-1));
SpaceGrp.setAttribute('OrientationBase','C along Z, B in YZ plane');
SpaceGrp.setAttribute('Centering','3D Primitive-Centered');
SpaceGrp.setAttribute('Lattice','3D Triclinic');
SpaceGrp.setAttribute('GroupName','GroupName');
SpaceGrp.setAttribute('Operators','1,0,0,0,0,1,0,0,0,0,1,0');
SpaceGrp.setAttribute('DisplayRange','0,1,0,1,0,1');
SpaceGrp.setAttribute('LineThickness','2');
SpaceGrp.setAttribute('CylinderRadius','0.2');
SpaceGrp.setAttribute('LabelAxes','1');
SpaceGrp.setAttribute('ActiveSystem','2');
SpaceGrp.setAttribute('ITNumber','1');
SpaceGrp.setAttribute('LongName','P 1');
SpaceGrp.setAttribute('Qualifier','Origin-1');
SpaceGrp.setAttribute('SchoenfliesName','C1-1');
SpaceGrp.setAttribute('System','Triclinic');
SpaceGrp.setAttribute('Class','1');
IdentMappng.appendChild(SpaceGrp);

RecLattc = docNode.createElement('ReciprocalLattice3D');
RecLattc.setAttribute('ID',num2str(NAtoms+8));
RecLattc.setAttribute('Parent',num2str(NAtoms+5));
IdentMappng.appendChild(RecLattc);

InfiniteMappng = docNode.createElement('InfiniteMapping');
InfiniteMappng.setAttribute('ID','3');
InfiniteMappng.setAttribute('Element','1,0,0,0,0,1,0,0,0,0,1,0');
InfiniteMappng.setAttribute('MappedObjects','2');

MappngSet.appendChild(InfiniteMappng);

docRootNode.appendChild(AtomisticTreeRootElement);

xmlwrite(XSDFileName,docNode);

end

function out1 = CPK_or_BnS(in1)

% out1 = 'CPK';

if strcmpi(in1,'Cu') || ...
        strcmpi(in1,'Ag') || ...
        strcmpi(in1,'Au') || ...
        strcmpi(in1,'Rh') || ...
        strcmpi(in1,'Ni') || ...
        strcmpi(in1,'Co') || ...
        strcmpi(in1,'Pd') || ...
        strcmpi(in1,'Pt')
    out1 = 'CPK';
else
%     strcmpi(in1,'C') || ...
%             strcmpi(in1,'H') || ...
%             strcmpi(in1,'O') || ...
%             strcmpi(in1,'S') || ...
%             strcmpi(in1,'F') || ...
%             strcmpi(in1,'Mo') || ...
%             strcmpi(in1,'N');
    out1 = 'Ball and Stick';
end

end
