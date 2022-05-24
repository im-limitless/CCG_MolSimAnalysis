function [ConvStepNum NotConvStepNum PercConv AveSteps StdSteps] = CP2k_parse_log(Basef, printNonConv)

logfile = [Basef '\out.log'];
fid = fopen(logfile);

textLines = textscan(fid, '%s','delimiter','\n', 'whitespace', '');
textLines = textLines{1};

Conv =  find(~cellfun(@isempty,strfind(textLines,'*** SCF run converged in ')));
NotConv = find(~cellfun(@isempty,strfind(textLines,' SCF run NOT ')));
nInnerSCF = find(~cellfun(@isempty,strfind(textLines,'Leaving inner SCF loop after reaching')));
AllStep = sort([Conv; NotConv]);
NotConvStepNum = find(ismember(AllStep, NotConv));
ConvStepNum = find(ismember(AllStep, Conv));

if ~isempty(NotConv) && printNonConv==1
    for i = 1:length(NotConvStepNum)
        warning(['SCF steps ' num2str(NotConvStepNum(i)) ' failed to converge'])
    end
end
PercConv = 100*length(NotConv)/(length(Conv)+length(NotConv));
disp(['Percentage of SCF steps that did not converge is ' num2str(PercConv) '%']);

nSteps = zeros(length(Conv),1);
    for i = 1:length(Conv)
        pIndx = strfind(textLines{Conv(i)}, 'steps');
        nSteps(i) = str2num(textLines{Conv(i)}(pIndx-5:pIndx-1));
    end
    
    nSCF = zeros(length(nInnerSCF),1);
    for i = 1:length(nInnerSCF)
        pIndx = strfind(textLines{nInnerSCF(i)}, 'steps');
        nSCF(i) = str2num(textLines{nInnerSCF(i)}(pIndx-5:pIndx-1));
    end
    
AveSteps = mean([nSCF(3:end); nSteps(3:end)]);
StdSteps = std([nSCF(3:end); nSteps(3:end)]);


disp(['Average number of SCF steps performed is ' num2str(AveSteps)]);

fclose(fid);