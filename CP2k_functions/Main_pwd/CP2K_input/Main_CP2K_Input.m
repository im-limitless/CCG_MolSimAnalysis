fclose('all'); clear all; clc;

newlinechar = char(10);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Materials Studio file
BaseInFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/Systems_prep/Original_structures/'; % Add the path to the folder where xsd files are (end in "\")
BaseOutFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/Systems_prep/'; % Add the path where the output folders will be written

FileList = dir([BaseInFldr '*_*.xsd']); % set the name of the xsd files using "*" wildcard; "*_*.xsd" reads any file within BaseInFldr with and underscore and .xsd file type. Any number of files can be read in one go.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Choose a Machine %%%%%
% MyCluster = 'Legion';
% MyCluster = 'Grace';
% MyCluster = 'Rhea';
MyCluster = 'ARCHER';
% MyCluster = 'IBchem';

%%%%% Set the number of cores %%%%%
% ncores = 1;
% ncores = 16;
% ncores = 32;
% ncores = 48;
% ncores = 64;
% ncores = 624;
% ncores = 768;
ncores = 1280;
% ncores = 1152;
% ncores = 1776;

%%%%% Set the number of processors per node %%%%%
ppn = 40;
% ppn = 24;

%%%%% Set the calculation type %%%%%
CalcType = 'BOMD';
% CalcType = 'CP';
% CalcType = 'GEO';

%%%%% Set the walltime in hours %%%%%
Walltime = 48;
% Walltime = 12;
% Walltime = 24;
% Walltime = 48;
% Walltime = 120;

%%%%% Choose if the jobs are to be farmed %%%%%
% Farming = 'Yes';
Farming = 'No';

%%%%% Make unique input file to converge wavefunction with no MD %%%%%
% WFN = 'Yes';
WFN = 'No';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(FileList)
    
    flname = FileList(i).name;
    disp(['Processing ' FileList(i).name '...']);
    fldrname = flname(1:end-4);
    fldrname = strrep(fldrname,' ','_');
    fldrname = strrep(fldrname,'(','_');
    fldrname = strrep(fldrname,')','_');
    
    if exist([BaseOutFldr fldrname],'dir')
        warning('Directory already exists! Continuing with next file.');
        continue
    end
    mkdir([BaseOutFldr fldrname]);
    
    XSDFileName = [BaseInFldr flname];
    
    [NAtoms,srtdAtomElement,srtdAtomPosition,n,AtomCur,AtomPos,Vec] = XSDFileRead(XSDFileName);
    
    [ListOfElems, Fix, Restrain, Vel] = convertXSDtoXYZ(BaseOutFldr, fldrname,  NAtoms, n, srtdAtomPosition, srtdAtomElement, Vec);
       
    CreateCP2KInputfile(BaseOutFldr, fldrname, Vec, Walltime, Fix, Restrain, CalcType, Vel, ListOfElems, WFN);
    
    CreateSubmissionFile(MyCluster, BaseOutFldr, fldrname, Walltime, ncores, ppn, 'No');
    
    % Copy the potential, basis and vdw kernel files from a base location
    % into each directory
    copyfile('/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/Systems_prep/CP2k_files/GTH_POTENTIALS', [BaseOutFldr fldrname '/GTH_POTENTIALS']);
    copyfile('/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/Systems_prep/CP2k_files/GTH_BASIS_SETS', [BaseOutFldr fldrname '/GTH_BASIS_SETS']);
    copyfile('/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/Systems_prep/CP2k_files/dftd3.dat', [BaseOutFldr fldrname '/dftd3.dat']);
    
%     % Copy the chainjob scripts to each directory
    %     copyfile('G:\Imperial\CP2k_files\chainjob.qs', [BaseOutFldr fldrname '\chainjob.qs']);
    %     copyfile('G:\Imperial\CP2k_files\AddChain.qs', [BaseOutFldr fldrname '\AddChain.qs']);
    
end

if strcmp(Farming, 'Yes')
    writeFarmingInput(BaseOutFldr, FileList, ncores, ppn, Walltime);
end