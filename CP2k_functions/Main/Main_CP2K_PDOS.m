clear all;
close all;

% plot the CP2K pdos file

BaseFldr = 'G:\Imperial\MattProjects\Edges\PostEquilibration\Pit\HF\';
system = 'CP_Pit_20F';


pDOSflnames = dir([BaseFldr system '\PDOS_16500\' '*_pp.pdos']);
if ~isempty(pDOSflnames)
for i = 1:length(pDOSflnames)
       
    disp(['Parsing PDOS file ' num2str(i) ' of ' num2str(length(pDOSflnames)) ' for ' system '...']);
    
    fid = fopen([BaseFldr system '\PDOS_16500\' pDOSflnames(i).name]);
    eofstat = false;
    
    % First line is the comment with the atom type and E(fermi)
%     SystemLine = fgetl(fid); eofstat = feof(fid);
%     KindIndx = find(ismember(SystemLine, ' '));
%     Kind = SystemLine(KindIndx(6)+1:KindIndx(7)-1);
Kind = 'Pts';
    HeadingLine = fgetl(fid); eofstat = feof(fid);
    j = 0;
    data.(Kind) = [];
    
    while ~eofstat
        j=j+1;
        textLine = fgetl(fid); eofstat = feof(fid);
        data.(Kind) = [data.(Kind); str2num(textLine)];
    end
    
    fclose(fid);

end
else
    warning('PDOS file does not exist, exiting...')
    return
end

