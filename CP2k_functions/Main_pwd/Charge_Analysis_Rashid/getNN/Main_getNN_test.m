clear all;  clc;
close all;

BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus/Al_AlO/';
system = 'Al_water';
Trajectory = 'Al_water_14300to18100_100step.xyz';


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

