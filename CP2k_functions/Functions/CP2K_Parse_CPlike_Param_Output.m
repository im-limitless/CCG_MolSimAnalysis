clear all
close all

% parse the .ener file from cp2k

Basefldr = 'G:\Imperial\MattProjects\Pt_Clean\BOMD\BOMD_NVT_1010_Clean\CP_Like\Noisy\BOMD_NVT_1010_Clean\Noisy\';
Allfldrs = dir([Basefldr '*_*']);

figure
hold on
set(gcf, 'color', 'w')
ylabel('Total Energy (eV)')
xlabel('Time (ps)')
% set(gca, 'Ylim', [-3.2270e4 -3.2269e4])
% set(gca, 'Xlim', [0 1.5])
C = colormap(lines);

for l=1:length(Allfldrs)
    

    if ~exist([Basefldr Allfldrs(l).name '\out.log'])
       warning(['Output for ' Allfldrs(l).name 'does not exist!'])
        continue
    else
        disp(['Processing ' Allfldrs(l).name]);
    end
    
        ParamScoreIndx = strfind(Allfldrs(l).name, '_');
    
    Stepsize(l) = str2num(Allfldrs(l).name(ParamScoreIndx(1)+1:ParamScoreIndx(2)-1));
    EXOR(l) = str2num(Allfldrs(l).name(ParamScoreIndx(3)+1:ParamScoreIndx(4)-1));
    Gamma(l) = str2num(Allfldrs(l).name(ParamScoreIndx(5)+1:ParamScoreIndx(6)-1));
    Noisy(l) = str2num(Allfldrs(l).name(ParamScoreIndx(7)+1:end));
    
    [ConvIndx NoConvIndx Perc(l) AveSteps(l) StdSteps(l)] = CP2k_parse_log([Basefldr Allfldrs(l).name], 0);
    [AveConverge(l) StdConverge(l)] = CP2k_parse_SCF_Loops([Basefldr Allfldrs(l).name]);

    
    flname = dir([Basefldr Allfldrs(l).name '\*.ener']);
    EnerFile = [Basefldr Allfldrs(l).name '\' flname.name];
    
    fid = fopen(EnerFile);
    
    eofstat = false;
    
    % First line is the number of atoms and is the same for all new sections
    HeadLine = fgetl(fid); eofstat = feof(fid);
    i = 0;
    data = [];
    
    while ~eofstat
        i=i+1;
        textLine = fgetl(fid); eofstat = feof(fid);
        data = [data; str2num(textLine)];
    end
    
    fclose(fid);
    
    
%         Cc(l) = randi([1 size(C,[1])]);
    Cc(l) = l;

    
    h(l) = plot(data(1:100,2)/1000, data(1:100,6), '-', 'color', C(Cc(l),:), 'markeredgecolor', C(Cc(l),:), 'markerfacecolor', C(Cc(l),:), 'markersize', 3);
%     text((data(100,2)/1000)+0.01, data(100,6), num2str(Step(l)), 'color', C(Cc(l),:))
    
    if length(Allfldrs) == 1
        h(l) = plot(data(1:end,2)/1000, data(1:end,6), '-', 'color', C(Cc(l),:), 'markeredgecolor', C(Cc(l),:), 'markerfacecolor', C(Cc(l),:), 'markersize', 3, 'linewidth', 2);
        legend(h, Allfldrs.name, 'interpreter', 'none', 'Location', 'east')
        hold off
        figure
        set(gcf, 'color', 'w')
        ylabel('Temperature (K)')
        xlabel('Time (ps)')
        hold on
        h2(l) = plot(data(1:end,2)/1000, data(1:end,4), '-', 'color', C(Cc(l),:), 'markeredgecolor', C(Cc(l),:), 'markerfacecolor', C(Cc(l),:), 'markersize', 3);
        hold off
        Step(l) = str2num(Allfldrs(l).name(end-7:end));
        text((data(end,2)/1000)+0.01, data(end-100,6), num2str(Step(l)), 'color', C(Cc(l),:))

    else
            Step(l) = str2num(Allfldrs(l).name(end-7:end));
            
% %     %% Energy
% %     h(l) = plot(data(1:end,2)/1000, data(1:end,6), '-', 'color', C(Cc(l),:), 'markeredgecolor', C(Cc(l),:), 'markerfacecolor', C(Cc(l),:), 'markersize', 3);
% %     text((data(end,2)/1000)+0.01, data(end,6), num2str(Step(l)), 'color', C(Cc(l),:));
%     Temp
    h(l) = plot(data(1:end,2)/1000, data(1:end,4), '-', 'color', C(Cc(l),:), 'markeredgecolor', C(Cc(l),:), 'markerfacecolor', C(Cc(l),:), 'markersize', 3, 'linewidth', 2);
    text((data(end,2)/1000)+0.01, data(end,4), num2str(Step(l)), 'color', C(Cc(l),:))
        
    end
    
end
% % % legend(h, Allfldrs.name, 'interpreter', 'none', 'Location', 'east')
% % % hold off

% [srtd, indx] = sort(Step);
% 
% figure
% hold on
% set(gcf, 'color', 'w')
% ylabel('Mean SCF Steps')
% xlabel('STEPSIZE')
% plot(srtd, AveSteps(indx), 'o', 'color', 'k', 'markerfacecolor', 'b', 'markeredgecolor', 'k')
% % errorbar(srtd, AveSteps(indx), StdSteps(indx), 'o', 'color', 'k', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'color', 'k')
% hold off
% 
figure
c = {'r' 'g' 'b'};
for i=1:max(EXOR)+1
semilogy(Stepsize(EXOR==i-1),AveConverge(EXOR==i-1), 'o')
hold on
ylabel('OT Convergence')
xlabel('STEPSIZE')
h(i)=errorbar(Stepsize(EXOR==i-1),AveConverge(EXOR==i-1),  StdConverge(EXOR==i-1), 'o', 'markerfacecolor', c{i}, 'markeredgecolor', 'k', 'color', 'k');
end

N = [min(EXOR):max(EXOR)];
legendCell = cellstr(num2str(N', 'EXOR=%-d'));
legend(h,legendCell)
hold off



