fclose('all'); clear all; clc;

newlinechar = char(10);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Materials Studio file
BaseInFldr = 'G:\Imperial\MattProjects\Pt_O\10Na12Cl\';
fldrname = 'BOMD_NVT_1012_O_0.66ML';



Type = 'DMP';
% Type = 'Noisy';

BaseOutFldr = [BaseInFldr fldrname '\' Type '\'];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Choose a Machine %%%%%
% MyCluster = 'Legion';
% MyCluster = 'Grace'; % Check the cell size and k mesh...
% MyCluster = 'Rhea';
% MyCluster = 'LRZ';
MyCluster = 'LRZ-test';
% MyCluster = 'IBchem';

% % %%%%% Set the number of cores %%%%%
% ncores = 1;
% ncores = 16;
% ncores = 32;
% ncores = 48;
% ncores = 64;
% ncores = 432;
ncores = 816;

% %%%%% Set the walltime in hours %%%%%
% Walltime = 4;
% Walltime = 12;
% Walltime = 24;
% Walltime = 48;
Walltime = 0.5;

VarNames = {'SS', 'EX', 'GA', 'NG'};
% Vars.Params = {[0.0025 0.005 0.0075 0.01 0.0125] [0] [0] [3e-4]};
Vars.Params = {[0.0075] [0] [0] [3e-4]};

Vars.BOMD_NVT_1010_Clean = {[0.01] [0] [0] [5e-4 1e-3 5e-3 1e-2 5e-2]};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Type, 'DMP')
    Vars = Vars.Params;
elseif strcmp(Type, 'Noisy')
    Vars = Vars.(fldrname);
end

[Coord, Vel, ABC, Fix] = getLastPosAndVel(BaseInFldr, fldrname);

for i = 1:length(Vars{1})
    for j = 1:length(Vars{2})
        for k = 1:length(Vars{3})
            for l = 1:length(Vars{4})
                DirString = [VarNames{1} '_' num2str(Vars{1}(i)) '_' VarNames{2} '_' num2str(Vars{2}(j)) '_' VarNames{3} '_' num2str(Vars{3}(k)) '_' VarNames{4} '_' num2str(Vars{4}(l))];
                CreateCP2KInputfile_CP_Like(BaseInFldr, BaseOutFldr, fldrname, Walltime, Coord, Vel, ABC, Fix, Vars{1}(i), Vars{2}(j), Vars{3}(k), Vars{4}(l), DirString);
                CreateSubmissionFile(MyCluster, BaseOutFldr, DirString, Walltime, ncores, 'No');
                
                copyfile('G:\Imperial\CP2k_files\GTH_POTENTIALS', [BaseOutFldr DirString '\GTH_POTENTIALS']);
                copyfile('G:\Imperial\CP2k_files\GTH_BASIS_SETS', [BaseOutFldr DirString '\GTH_BASIS_SETS']);
                copyfile('G:\Imperial\CP2k_files\dftd3.dat', [BaseOutFldr DirString '\dftd3.dat']);
            end
        end
    end
end
