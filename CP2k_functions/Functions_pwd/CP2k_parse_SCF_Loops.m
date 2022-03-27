function [AveError, StdError] = CP2k_parse_SCF_Loops(Basef)
% logfile = 'G:\Imperial\CP-LIKE\BulkWater\CPLIKE\Parameterisation\FromShort\STEPSIZE\WaterCP_fromShort_ST_1.00E-01\out.log';
logfile = [Basef '\out.log'];
fid = fopen(logfile);

textLines = textscan(fid, '%s','delimiter','\n', 'whitespace', '');
textLines = textLines{1};
fclose(fid);

OTindx =  find(~cellfun(@isempty,strfind(textLines,'OT DIIS')));
% OTdat =[];

for i = 1:length(OTindx)-1
    
    if OTindx(i+1)-OTindx(i) == 1
        continue
    elseif OTindx(i+1)-OTindx(i) > 1
        A=strsplit(textLines{OTindx(i)});
        OTdat = str2num(A{7});
        Conv(i) = OTdat;
    end
end

Conv = Conv(Conv ~= 0);
AveError = mean(Conv(3:end));
StdError = std(Conv(3:end));