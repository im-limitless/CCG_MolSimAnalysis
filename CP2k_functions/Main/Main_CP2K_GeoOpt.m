clear all;
close all;
clc;

BaseFldr = 'G:\Imperial\NaIon\hard_carbon\pore\pristine\sodium\Matt\';
system = dir([BaseFldr 'p*']);

mkdir([BaseFldr 'Energy']);
mkdir([BaseFldr 'xsd']);

for j = 1:length(system)

    disp(['Processing ' system(j).name '...']);
    
logfile = [BaseFldr system(j).name '\out.log'];
fid = fopen(logfile);

textLines = textscan(fid, '%s','delimiter','\n', 'whitespace', '');
textLines = textLines{1};
fclose(fid);

Conv =  find(~cellfun(@isempty,strfind(textLines,'Max. gradient ')));
TotEng =  find(~cellfun(@isempty,strfind(textLines,'ENERGY| Total FORCE_EVAL')));
Tol =  find(~cellfun(@isempty,strfind(textLines,'Conv. limit for gradients')));
Tol = (27.2114)*(0.529177^2)*str2num(textLines{Tol(1)}(33:end));

F = zeros(length(Conv),1);
E = zeros(length(TotEng),1);

for i = 1:length(Conv)
    F(i,1) = (27.2114)*(0.529177^2)*str2num(textLines{Conv(i)}(33:end));
end

for i = 1:length(TotEng)
    E(i,1) = (27.2114)*str2num(textLines{TotEng(i)}(50:end));
end

plot(E, '-ok', 'markerfacecolor', 'r')
ylabel('Total Energy (eV)')
xlabel('Optimisation Step')
drawnow
set(gcf, 'position', [872   349   560   420])

semilogy(F, '-ok', 'markerfacecolor', 'b')
hold on
semilogy([0 length(F)],[Tol Tol], '--r')
ylabel('Force (eV\cdotA^{-2})')
xlabel('Optimisation Step')
hold off
drawnow
set(gcf, 'position', [222   349   560   420])

if F(end) > Tol
    warning(['Forces only converged to ' num2str(F(end)) '...'])
end

disp(['Total Energy = ' num2str(E(end)) ' eV']);

dlmwrite([BaseFldr system(j).name '\energy' system(j).name],E(end),'precision',15);
        copyfile([BaseFldr system(j).name '\energy' system(j).name], [BaseFldr '\Energy\']);

CP2kOptimPathParse(BaseFldr,system(j).name,[system(j).name '-pos-1.xyz'])
copyfile([BaseFldr system(j).name '\' system(j).name '-pos-1.xsd'], [BaseFldr '\xsd\']);

end



