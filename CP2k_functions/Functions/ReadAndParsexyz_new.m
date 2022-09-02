function [xyz, XYZ, Indx, Atoms, AtomList, nAtoms, startConfig, nConfigs, Step] = ReadAndParsexyz_new(Base, Fldr, Traj, ABC, Rot)
% Function to read .xyz files as outputted by CP2K. Written by Matt Darby
% version 3
% The function requires the base directory (Base), 
% calculation folder name (Fldr), xyz file name (Traj),
% lattice vectors from "getABCvectors" (ABC)
% and an optional parameter 
% Rot = [angle of rotation, axis of rotation (x=1, y=2, z=3), true (1) or false (0)] see "RotateCoordinateAxes"
% returns "xyz" data structure with fields of "AtomList" (aka atom names from xyz file) e.g. "xyz.Pt". xyz.Pt is a 3D array of (configuration, atom number, coordinates) and size = nConfigs x nAtoms (of AtomList type) x 3. 
% XYZ is a non-atomised array of (configuration, atom number, coordinates) and size = nConfigs x nAtoms x 3
% Indx is a data structure containing the atom numbers (as ordered in the .xyz file) of each atom type, e.g. "Indx.Pt" retunrs the position of all Pt atoms in the xyz file.
% Indx.Pt can be used to access elements of XYZ corresponding to Pt only, i.e. XYZ(1, Indx.Pt, :) will return the x, y and z (3) coordinates for configuration 1 of all Pt atoms.


fid  = fopen([Base Fldr '\' Traj]);
disp('Reading xyz data...');
nAtoms = str2num(fgetl(fid));
InitialI = fgetl(fid);

Atoms = [];

for i = 1:nAtoms
    AtomLine = fgetl(fid);
    Atoms = [Atoms; AtomLine(1:5)];
end

fgetl(fid);
NextI = fgetl(fid);

fclose(fid);

A = readmatrix([Base Fldr '\' Traj], 'headerlines', 0, 'FileType', 'Text', 'CommentStyle', 'i =');
nConfigs = size(A,1)/(nAtoms+1);
removeNAtoms = (0:nConfigs-1)*(nAtoms+1)+1;
A(removeNAtoms,:) = [];
A = A(:,2:4);

XYZ = permute(reshape(A',[3,nAtoms, nConfigs]),[3,2,1]);

XYZ = wrapXYZ(XYZ, ABC);
XYZ = RotateCoordinateAxes(XYZ, ABC, Rot(1), Rot(2), Rot(3));


AtomList = unique(Atoms,'rows');

for i = 1:size(AtomList,1)
    Indx.(strtrim(AtomList(i,:))) =  find(sum(Atoms == AtomList(i,:),2) == length(AtomList(i,:)));
    xyz.(strtrim(AtomList(i,:))) = zeros(nConfigs, length(Indx.(strtrim(AtomList(i,:)))), 3);
end

startConfig = 1;
if nConfigs > 1
    %     EqI = strfind(InitialI, '=');
    %     ComI = strfind(InitialI, ',');
    %     posI = [EqI(1) ComI(1)];
    %     Step1 = str2num(InitialI(posI(1)+1:posI(2)-1));
    Step1str = strsplit(InitialI);
    Step1 = str2num(Step1str{3});
    %     EqI = strfind(NextI, '=');
    %     ComI = strfind(NextI, ',');
    %     posI = [EqI(1) ComI(1)];
    %     Step2 = str2num(NextI(posI(1)+1:posI(2)-1));

    Step2str = strsplit(NextI);
    Step2 = str2num(Step2str{3});

    Step = Step1:(Step2-Step1):(Step1+(Step2-Step1)*(nConfigs-1));

else
    Step = 0;
end

for j = startConfig:nConfigs

    disp(['Processing configuration number ' num2str(j) ' of ' num2str(nConfigs)]);

    for k = 1:size(AtomList,1)
        xyz.(strtrim(AtomList(k,:)))(j,:, 1:3) = XYZ(j, Indx.(strtrim(AtomList(k,:))),1:3);
    end
end

return