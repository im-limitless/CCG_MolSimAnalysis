fclose('all'); clear all; clc;

newlinechar = char(10);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Materials Studio file
BaseInFldr = 'G:\Imperial\MRes\Kehan\Matt\H\'; % Add the path to the folder where xsd files are (end in "\")
BaseOutFldr = 'G:\Imperial\MRes\Kehan\Matt\H2\';
FileList = dir([BaseInFldr '*_*']); % set the name of the xsd files using "*" wildcard; "*_*.xsd" reads any file within BaseInFldr with and underscore and .xsd file type. Any number of files can be read in one go.

for i = 1:length(FileList)
    
mkdir([BaseOutFldr FileList(i).name]);

CP2kInputParameterReplace(BaseInFldr, BaseOutFldr, FileList(i).name, 'EPS_SCF', '5e-7')

end