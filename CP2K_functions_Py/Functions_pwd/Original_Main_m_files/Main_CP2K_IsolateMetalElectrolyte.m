clear all; close all; clc;

% script to find coordinates from restart file and make metal only +
% eletrolyte only inputs for Electrostatic Potentials

BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\CP_Pit_18H22F\';
% BaseFldr = 'G:\Imperial\MattProjects\Pt_Clean\CorrectVolume\Pt_12H10F\';

Restarts = dir([BaseFldr '*-1_*.restart']);

for i = 1:length(Restarts)
    
    disp(['Processing restart file ' Restarts(i).name]); 
    
    OutFldr = [BaseFldr 'EffectivePD\' Restarts(i).name(1:end-8) '\'];
    
    if exist(OutFldr,'dir')
        warning('Directory already exists! Continuing with next file.');
        continue
    end
    
    mkdir(OutFldr);
    
    fid = fopen([BaseFldr Restarts(i).name]);
    lines = textscan(fid,'%s','delimiter','\n', 'whitespace', '');
    fclose(fid);
    
    lines = lines{1};
    startCELL = find(~cellfun(@isempty,strfind(lines,'&CELL')));
    endCELL = find(~cellfun(@isempty,strfind(lines,'&END CELL')));
    IndxCell = startCELL:endCELL;
    
    startXYZ =  find(~cellfun(@isempty,strfind(lines,'&COORD')));
    endXYZ = find(~cellfun(@isempty,strfind(lines,'&END COORD')));
    
    IndxPt = find(~cellfun(@isempty,strfind(lines,'Pt')));
    IndxPt = IndxPt(IndxPt < endXYZ & IndxPt > startXYZ);
    IndxElectrolyte = (startXYZ+1:endXYZ-1)';
    IndxElectrolyte = IndxElectrolyte(find(~ismember(IndxElectrolyte,IndxPt)));
    
    %% write the file
    MetalOut = [Restarts(i).name(1:end-8) '_Metal'];
    if exist([OutFldr MetalOut],'dir')
        warning('Directory already exists! Will not overwrite.');
    else
        
        mkdir([OutFldr MetalOut]);
        InputFromRestart(OutFldr, MetalOut, IndxPt, IndxCell, lines);
        CreateSubmissionFile('LRZ', OutFldr, MetalOut, 1800, 768, 24, 'No');
        copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS', [OutFldr MetalOut '\GTH_POTENTIALS']);
        copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS', [OutFldr MetalOut '\GTH_BASIS_SETS']);
        copyfile('G:\Imperial\CP2k_files\dftd3.dat', [OutFldr MetalOut '\dftd3.dat']);
    end
    
    ElectrolyteOut = [Restarts(i).name(1:end-8) '_Electrolyte'];
    if exist([OutFldr '\' Restarts(i).name(1:end-8) '_Electrolyte'],'dir')
        warning('Directory already exists! Will not overwrite.');
    else
        mkdir([OutFldr ElectrolyteOut]);
        InputFromRestart(OutFldr, ElectrolyteOut, IndxElectrolyte, IndxCell, lines);
        CreateSubmissionFile('LRZ', OutFldr, ElectrolyteOut, 1800, 768, 24, 'No');
        copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS', [OutFldr ElectrolyteOut '\GTH_POTENTIALS']);
        copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS', [OutFldr ElectrolyteOut '\GTH_BASIS_SETS']);
        copyfile('G:\Imperial\CP2k_files\dftd3.dat', [OutFldr ElectrolyteOut '\dftd3.dat']);
    end
    
    
end
fclose('all');
