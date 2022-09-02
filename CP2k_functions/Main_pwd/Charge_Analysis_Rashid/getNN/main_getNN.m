clear all;
close all;

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus/Al_AlO/Al_water/Testing/'; 
system = 'Al_water';
Trajectory = 'Al_water_14300to18100_100step.xyz';

ABC = getABCvectors(BaseFldr, system);
[xyz, XYZ, Indx, ~, ~, nAtoms, startConfig, nConfigs, StepNum] = ReadAndParsexyz(BaseFldr, system, Trajectory, ABC);
[Atoms, AtomList, AtomIndx] = getAtomNamesFromInputXYZ(BaseFldr, system);

for snap = 39
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
   
    [VecAlAl, DistAlAl] = GetAtomCorrelation(XYZ_snap, AtomIndx.Al1, AtomIndx.Al1, ABC); 

end

% for i= 1:length(AtomIndx.Al1)
%     DoubleCount= 1:length(AtomIndx.Al1);
%     DoubleCount(i)=[];
%     AlNN_indx= find(DistAlAl(DoubleCount,i)<3.6)
% end

%%%Bader Charge%%%

fldrname = [BaseFldr system '/Bader_Analysis/'];
ACFfiles = dir([fldrname 'ACF_*.dat']);

DoubleAnalType = 'MassDensity';
% DoubleAnalType = 'Radial';

% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx] = getAtomNamesFromInputXYZ(BaseFldr, system);

% % Read the Bader charge "ACF" files and extract the raw charge Q/net charge Qnet
Q = zeros(length(Atoms),length(ACFfiles));
Qnet = zeros(length(Atoms),length(ACFfiles));

for n = 1:length(ACFfiles)
    [Q(:,n), Qnet(:,n), StepNum(n)] = extractBaderCharges(fldrname, ACFfiles(n).name, Atoms, AtomList);
end

%%%%%%%%%%%%%%%%%%% 1st NN Bader Charge %%%%%%%%%%%%%%%%%%%
Atom=[]

Atom_1Q=[]
% Qnet_Atom=[]

NN_1=[]

NN_1A=[]
NN_1Q=[]
% Qnet_NN_1=[]

for i= 1:length(AtomIndx.Al1)
    DoubleCount= 1:length(AtomIndx.Al1);
    DoubleCount(i)=[];
    i;
    R_ind_Atom=AtomIndx.Al1(i);
%     Atom=[Atom;R_ind_Atom R_Q_Atom]
    R_Q_Atom=Qnet(AtomIndx.Al1(i),1);
%     Qnet_Atom=[Qnet_Atom;R_Q_Atom]
    Atom=[Atom;R_ind_Atom R_Q_Atom];
    Atom_1Q=[Atom_1Q; R_Q_Atom];


    AlNN_indx= find(DistAlAl(DoubleCount,i)<3.6);
    R_ind_NN=AtomIndx.Al1(AlNN_indx)
    R_Q_NN=Qnet(AtomIndx.Al1(AlNN_indx),1)
%     Qnet(AtomIndx.Al1(AlNN_indx),1)

    NN_1=[NN_1;R_ind_NN R_Q_NN];
    NN_1A=[NN_1A; R_ind_NN];
    NN_1Q=[NN_1Q; R_Q_NN];

end

%%%%%%%%%%%%%%%%%%% 2nd NN Bader Charge %%%%%%%%%%%%%%%%%%%
for snap = 39
    
    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(snap,:,:);
   
    [VecAlAl2, DistAlAl2] = GetAtomCorrelation(XYZ_snap, AtomIndx.Al1, AtomIndx.Al2, ABC); 

end

NN_2=[]

NN_2A=[]
NN_2Q=[]
% Qnet_NN_2=[]

for i= 1:length(AtomIndx.Al2)
    DoubleCount= 1:length(AtomIndx.Al2);
    DoubleCount(i)=[];
    i


    AlNN_indx2= find(DistAlAl2(DoubleCount,i)>3.6 & DistAlAl2(DoubleCount,i)<=4.52);
    AlNN_indx2

    R_ind_NN=AtomIndx.Al2(AlNN_indx2)
    R_Q_NN=Qnet(AtomIndx.Al2(AlNN_indx2),1)
%     Qnet(AtomIndx.Al1(AlNN_indx),1)

    NN_2=[NN_2;R_ind_NN R_Q_NN];
    NN_2A=[NN_2A;i;i;i ;R_ind_NN];
    NN_2Q=[NN_2Q; R_Q_NN];

end