clear all;  clc;
close all;
PathSep = setOSpathSep;

%% Home
% BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/5lyr_systems/Al_AlO/Al_water/Fix_volume/BOMD/Temperature_fix/';
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge_2/The_rest/';
% BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge/5lyr_systems/Al_AlO/';
system = 'OH_0.2ML';
Trajectory = 'OH_0.2ML_172000to220000_1000step.xyz';

% % Call function to find ABC vectors from .inp file
ABC = getABCvectors(BaseFldr, system);

% % get the names of atoms from original xyz input file
[Atoms, AtomList, Indx, Indxfns, Kinds, Elements, PP] = getAtomInfoFromInput(BaseFldr, system);

% % parse coordinates of atoms along trajectory and wrap into cell
[xyz, XYZ, ~, ~, ~, nAtoms, startConfig, nConfigs, StepNum_Traj] = ReadAndParsexyz_new(BaseFldr, system, Trajectory, ABC, [0 0 0]);

OH_Coverage = zeros(nConfigs, 1);
H2O_Coverage = zeros(nConfigs, 1);
H3O_Coverage = zeros(nConfigs, 1);
O_Coverage = zeros(nConfigs, 1);


for i = 1:nConfigs

    XYZ_snap = zeros(size(XYZ,2), size(XYZ,3));
    XYZ_snap(:,:) = XYZ(i,:,:);

    % get the distances between pairs of atoms
    [~, DistAlO] = GetAtomCorrelation(XYZ_snap, Indx.Al1, Indx.O, ABC);
    [r,c] = find(DistAlO < 2.5); % Rashid to fix this "2" by looking at RDF minimum - r = row aka O atom number, c = column aka Al1 atom number

    [C,~]=unique(r); %How many O close to Al1
    [num,~]=size(C);
    rlen=(1:length(DistAlO(:,1)))';
    O_Bulk=setdiff(rlen,C);
    

    % [~, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, setdiff(Indx.O,Indx.O(C)), ABC);
    [~, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.O(O_Bulk), ABC);
%     [~, DistOH] = GetAtomCorrelation(XYZ_snap, Indx.H, Indx.O, ABC);
    [rOH,cOH] = find(DistOH < 1.28);

    OH=[];
    H2O=[];
    H3O=[];

    for j =1:length(rOH)
        if size(find(rOH==rOH(j)),1)==1
            OH=[OH,find(rOH==rOH(j))];

        elseif size(find(rOH==rOH(j)),1)==2
            H2O=[H2O,find(rOH==rOH(j))];

        elseif size(find(rOH==rOH(j)),1)==3
            H3O=[H3O,find(rOH==rOH(j))];
        
        end
                
    end

    % OH_unique=OH;
    H2O_unique=unique(H2O(1,:))';

    if not(isempty(H3O))
        H3O_unique=unique(H3O(1,:))';
        H3O_Coverage(i) = length(unique(H3O(1,:)));
        H3O_indicies{i}=Indx.O(O_Bulk(rOH(H3O_unique))); %save the indicies for H3Os 
    else
        H3O_Coverage(i) = 0;
    end


     if not(isempty(OH))
        OH_unique=unique(OH(1,:))';
        OH_Coverage(i) = length(unique(OH(1,:)));
        OH_indicies{i}=Indx.O(O_Bulk(rOH(OH_unique))); %save the indicies for OHs 
    else
        OH_Coverage(i) = 0;
     end

    % OH_Coverage(i) =length(unique(OH(1,:)));
    H2O_Coverage(i) = length(unique(H2O(1,:))); 
    
    %save the indicies for H2O
    % OH_indicies{i}=Indx.O(C(rOH(OH_unique)));
    H2O_indicies{i}=Indx.O(O_Bulk(rOH(H2O_unique)));
    


end
%-------------------------------------------------------

Total_Coverage = OH_Coverage+H2O_Coverage+H3O_Coverage;

%H2O,H3O,OH and total 
figure
hold on
plot(StepNum_Traj/2000, OH_Coverage, '-ok', 'markerfacecolor', 'r')
plot(StepNum_Traj/2000, H2O_Coverage, '-ok', 'markerfacecolor', 'b')
plot(StepNum_Traj/2000, H3O_Coverage, '-ok', 'markerfacecolor', 'k')
plot(StepNum_Traj/2000, Total_Coverage, '-ok', 'markerfacecolor', [0 0.8 0])
% plot((1:nConfigs)/200, OH_Coverage, '-ok', 'markerfacecolor', 'r')
% plot((1:nConfigs)/200, H2O_Coverage, '-ok', 'markerfacecolor', 'b')
% plot((1:nConfigs)/200, H3O_Coverage, '-ok', 'markerfacecolor', 'auto')
% plot((1:nConfigs)/200, Total_Coverage, '-ok', 'markerfacecolor', [0 0.8 0])
legend('OH', 'Water','H3O', 'Total')
xlabel('Time (ps)')
ylabel('Number of Molecules')

%H3O and OH only
figure
hold on
plot(StepNum_Traj/2000, OH_Coverage, '-ok', 'markerfacecolor', 'r')
plot(StepNum_Traj/2000, H3O_Coverage, '-ok', 'markerfacecolor', 'k')
% plot((1:nConfigs)/200, OH_Coverage, '-ok', 'markerfacecolor', 'r')
% plot((1:nConfigs)/200, H3O_Coverage, '-ok', 'markerfacecolor', 'auto')
legend('OH','H3O')
xlabel('Time (ps)')
ylabel('Number of Molecules')



%Stats
disp(['Ave. coverage of OH = ' num2str(mean(OH_Coverage(1:end)/108), '%.2f') ' +/- ' num2str(std(OH_Coverage(1:end)/108), '%.2f') ' ML'])
disp(['Ave. coverage of H2O = ' num2str(mean(H2O_Coverage(1:end)/108), '%.2f') ' +/- ' num2str(std(H2O_Coverage(1:end)/108), '%.2f') ' ML'])
disp(['Ave. coverage of H3O = ' num2str(mean(H3O_Coverage(1:end)/108), '%.2f') ' +/- ' num2str(std(H3O_Coverage(1:end)/108), '%.2f') ' ML'])
disp(['Ave. Total coverage = ' num2str(mean(Total_Coverage(1:end)/108), '%.2f') ' +/- ' num2str(std(Total_Coverage(1:end)/108), '%.2f') ' ML'])
disp(['Ave. O coverage = ' num2str(mean(O_Coverage(1:end)/108), '%.2f') ' +/- ' num2str(std(O_Coverage(1:end)/108), '%.2f') ' ML'])