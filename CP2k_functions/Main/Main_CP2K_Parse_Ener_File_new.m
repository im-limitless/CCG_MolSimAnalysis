% parse the .ener file from cp2k
clear all;
close all;
PathSep = setOSpathSep;

% Set the location of the calculation output
BaseFldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/Grand_Challenge_2/The_rest/'; % Base directory containing calculation directory ("\" included at end)
system = 'Al'; % Name of calculation directory (no "\")

% Set the number of steps before the end that are used to compute averages
% SampleRange = 6000;
SampleRange = 20000;

% Choose whether to plot an inset
% Inset = 'Yes';
Inset = 'No';

Allfldrs = dir([BaseFldr system '*']);

for i = 1:length(Allfldrs)

    flname = dir([BaseFldr Allfldrs(i).name setOSpathSep '*.ener']);
    EnerFile = [BaseFldr Allfldrs(i).name setOSpathSep flname.name];

    if ~isempty(flname)
        disp(['Parsing energy file for ' Allfldrs(i).name '...']);

        fid = fopen(EnerFile);
        data = readmatrix(EnerFile, 'NumHeaderLines',1, 'FileType', 'Text');
    else
        warning('Energy file does not exist, exiting...')
        return
    end

    fig1 = figure;
    set(gcf, 'color', 'w', 'position', [388   417   592   420])
    ax1 = axes;
    h1 = plot(ax1, data(1:end,2)/1000, 27.211399*data(1:end,5), '-', 'linewidth', 1, 'color', 'r');
    title(['Ave. PE for last ' num2str(SampleRange/2000) ' ps = ' num2str(mean(27.211399*data(end-SampleRange:end,5)),'%.2f') ' \pm ' num2str(std(27.211399*data(end-SampleRange:end,5)),'%.2f') ' eV']);
    ylabel(ax1, 'Potential Energy (eV)')
    xlabel(ax1, 'Time (ps)')
    set(ax1, 'fontsize', 12)
    if size(data,1) > SampleRange && strcmp(Inset, 'Yes')
        ax2 = axes('Position',[0.64 0.6 0.25 0.2]);
        h2 = plot(ax2, data(end-SampleRange:end,2)/1000, 27.211399*data(end-SampleRange:end,5), '-', 'linewidth',0.5, 'color', 'r');
        ylabel(ax2, 'P. E. (eV)')
        xlabel(ax2, 'Time (ps)')
        set(ax2, 'fontsize', 8, 'xlim', [data(end-SampleRange,2)/1000 data(end,2)/1000] )

    end

    fig2 = figure;
    set(gcf, 'color', 'w', 'position', [388   417   592   420])
    ax1 = axes;
    h1 = plot(ax1, data(1:end,2)/1000, 27.211399*data(1:end,3), '-', 'linewidth', 1, 'color', [0 0.5 0]);
    title(['Ave. KE for last ' num2str(SampleRange/2000) ' ps = ' num2str(mean(27.211399*data(end-SampleRange:end,3)),'%.2f') ' \pm ' num2str(std(27.211399*data(end-SampleRange:end,3)),'%.2f') ' eV']);
    ylabel(ax1, 'Kinetic Energy (eV)')
    xlabel(ax1, 'Time (ps)')
    set(ax1, 'fontsize', 12)
    if size(data,1) > SampleRange && strcmp(Inset, 'Yes')
        ax2 = axes('Position',[0.64 0.6 0.25 0.2]);
        h2 = plot(ax2, data(end-SampleRange:end,2)/1000, 27.211399*data(end-SampleRange:end,3), '-', 'linewidth',0.5, 'color', [0 0.5 0]);
        ylabel(ax2, 'K. E. (eV)')
        xlabel(ax2, 'Time (ps)')
        set(ax2, 'fontsize', 8, 'xlim', [data(end-SampleRange,2)/1000 data(end,2)/1000] )
    end

    fig3 = figure;
    set(gcf, 'color', 'w', 'position', [388   417   592   420])
    ax1 = axes;
    h1 = plot(ax1, data(1:end,2)/1000, 27.211399*data(1:end,6), '-', 'linewidth', 1, 'color', 'b');
    % h2 = plot(ax1, data(1:end,2)/1000, 27.211399*(data(1:end,3)+data(1:end,5)), '-', 'linewidth', 1, 'color', 'm');
    title(['Ave. CQ for last ' num2str(SampleRange/2000) ' ps = ' num2str(mean(27.211399*data(end-SampleRange:end,6)),'%.2f') ' \pm ' num2str(std(27.211399*data(end-SampleRange:end,6)),'%.2f') ' eV']);
    ylabel(ax1, 'Conserved Quantity (eV)')
    xlabel(ax1, 'Time (ps)')
    set(ax1, 'fontsize', 12)
    if size(data,1) > SampleRange && strcmp(Inset, 'Yes')
        ax2 = axes('Position',[0.64 0.6 0.25 0.2]);
        h2 = plot(ax2, data(end-SampleRange:end,2)/1000, 27.211399*data(end-SampleRange:end,6), '-', 'linewidth',0.5, 'color', 'b');
        ylabel(ax2, 'C. Q. (eV)')
        xlabel(ax2, 'Time (ps)')
        set(ax2, 'fontsize', 8, 'xlim', [data(end-SampleRange,2)/1000 data(end,2)/1000] )
    end
    CQDrift = 27.211399*data(end-SampleRange:end,6)-detrend(27.211399*data(end-SampleRange:end,6));
    dCQDriftdt = mean(diff(CQDrift))*2000;

    fig4 = figure;
    set(gcf, 'color', 'w', 'position', [388   417   592   420])
    ax1 = axes;
    h1 = plot(ax1, data(1:end,2)/1000, data(1:end,4), '-', 'linewidth', 1, 'color', 'k');
    title(['Ave. Temp. for last ' num2str(SampleRange/2000) ' ps = ' num2str(mean(data(end-SampleRange:end,4)),'%.2f') ' \pm ' num2str(std(data(end-SampleRange:end,4)),'%.2f') ' K']);
    ylabel(ax1, 'Temperature (K)')
    xlabel(ax1, 'Time (ps)')
    set(ax1, 'fontsize', 12)
    if size(data,1) > SampleRange && strcmp(Inset, 'Yes')
        ax2 = axes('Position',[0.64 0.6 0.25 0.2]);
        h2 = plot(ax2, data(end-SampleRange:end,2)/1000, data(end-SampleRange:end,4), '-', 'linewidth',0.5, 'color', 'k');
        ylabel(ax2, 'Temp. (K)')
        xlabel(ax2, 'Time (ps)')
        set(ax2, 'fontsize', 8, 'xlim', [data(end-SampleRange,2)/1000 data(end,2)/1000] )
    end
    
    TDrift = data(end-SampleRange:end,4)-detrend(data(end-SampleRange:end,4));
    dTDriftdt = mean(diff(TDrift))*2000;

    % exportgraphics(fig1,[BaseFldr Allfldrs(i).name setOSpathSep Allfldrs(i).name '_PE.jpg'],'Resolution',300)
    % exportgraphics(fig2,[BaseFldr Allfldrs(i).name setOSpathSep Allfldrs(i).name '_KE.jpg'],'Resolution',300)
    % exportgraphics(fig3,[BaseFldr Allfldrs(i).name setOSpathSep Allfldrs(i).name '_CQ.jpg'],'Resolution',300)
    % exportgraphics(fig4,[BaseFldr Allfldrs(i).name setOSpathSep Allfldrs(i).name '_Temp.jpg'],'Resolution',300)

    tLast(i) = data(end,2)/1000;
    Range(i) = SampleRange/2000;

    Temp(i) = mean(data(end-SampleRange:end,4));
    TempStd(i) = std(data(end-SampleRange:end,4));
    TempDrift(i) = dTDriftdt;

    TotEng(i) = mean(27.211399*data(end-SampleRange:end,3))+mean(27.211399*data(end-SampleRange:end,5));
    TotEngStd(i) = std(27.211399*data(end-SampleRange:end,6));
    TotEngDrift(i) = dCQDriftdt;

    disp(['Ave. Total Energy for last ' num2str(SampleRange/2000) ' ps = ' num2str(TotEng(i)) ' eV']);
    disp(['Std. Dev. = ' num2str(TotEngStd(i)) ' eV']);
    disp(['Drift = ' num2str(TotEngDrift(i)) ' eV/ps']);

    disp(['Ave. Temperature for last ' num2str(SampleRange/2000) ' ps = ' num2str(mean(data(end-SampleRange:end,4))) ' K']);
    disp(['Std. Dev. = ' num2str(std(data(end-SampleRange:end,4))) ' K']);
    disp(['Drift = ' num2str(dTDriftdt) ' K/ps']);
    % disp(['Ave. Total Energy for last 1 ps = ' num2str(mean(27.211399*data(end-SampleRange:end,6)))]);

    [Atoms, AtomList, Indx, Indxfns, Kinds, Elements, PP] = getAtomInfoFromInput(BaseFldr, Allfldrs(i).name);
    AtomCount(i) = size(Atoms,1);
    OAtomCount(i) = size(Indx.O,1);
    HAtomCount(i) = size(Indx.H,1);
end

if length(Allfldrs) > 1
    T = table({Allfldrs.name}', AtomCount', OAtomCount', HAtomCount', tLast', Range', TotEng', TotEngStd', TotEngDrift', Temp', TempStd', TempDrift');
    T.Properties.VariableNames = {'System', 'Atom Count', 'O Atom Count', 'H Atom Count',  'Time (ps)', 'Sampled (ps)', 'Ave. Total Energy (eV)', 'Std. Dev. (eV)', 'Drift (eV)', 'Ave. Temperature (K)', 'Std. Dev. (K)', 'Drift (K)'};
    writetable(T, [BaseFldr 'Data.xlsx'])
end
