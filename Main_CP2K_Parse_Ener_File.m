% parse the .ener file from cp2k
clear all;
close all;

% Set the location of the calculation output
Basefldr = '/Users/rashidal-heidous/Google Drive (local)/Academic Career (Current:local)/UK Postgrad Journey (ICL)/PhD/PhD/cp2k jobs/Jobs/ARCHER2/AIMD/EleventhTimeLucky_Plus/Al_AlO_OH/'; % Base directory containing calculation directory ("\" included at end)
system = 'AlO_0.5ML_OH'; % Name of calculation directory (no "\")


Allfldrs = dir([Basefldr '*' system]);
flname = dir([Basefldr Allfldrs.name '/*.ener']);
EnerFile = [Basefldr Allfldrs.name '/' flname.name];

if ~isempty(flname)
    disp(['Parsing energy file for ' system '...']);
    
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
else
    warning('Energy file does not exist, exiting...')
    return
end

fig1 = figure;
set(gcf, 'color', 'w', 'position', [388   417   592   420])
ax1 = axes;
h1 = plot(ax1, data(1:end,2)/1000, 27.211399*data(1:end,5), '-', 'linewidth', 1, 'color', 'r');
ylabel(ax1, 'Potential Energy (eV)')
xlabel(ax1, 'Time (ps)')
set(ax1, 'fontsize', 12)
if size(data,1) > 2000
    ax2 = axes('Position',[0.64 0.68 0.25 0.2]);
    h2 = plot(ax2, data(end-2000:end,2)/1000, 27.211399*data(end-2000:end,5), '-', 'linewidth',0.5, 'color', 'r');
    ylabel(ax2, 'P. E. (eV)')
    xlabel(ax2, 'Time (ps)')
    set(ax2, 'fontsize', 8, 'xlim', [data(end-2000,2)/1000 data(end,2)/1000] )
end


fig2 = figure;
set(gcf, 'color', 'w', 'position', [388   417   592   420])
ax1 = axes;
h1 = plot(ax1, data(1:end,2)/1000, 27.211399*data(1:end,3), '-', 'linewidth', 1, 'color', [0 0.5 0]);
ylabel(ax1, 'Kinetic Energy (eV)')
xlabel(ax1, 'Time (ps)')
set(ax1, 'fontsize', 12)
if size(data,1) > 2000
    ax2 = axes('Position',[0.64 0.68 0.25 0.2]);
    h2 = plot(ax2, data(end-2000:end,2)/1000, 27.211399*data(end-2000:end,3), '-', 'linewidth',0.5, 'color', [0 0.5 0]);
    ylabel(ax2, 'K. E. (eV)')
    xlabel(ax2, 'Time (ps)')
    set(ax2, 'fontsize', 8, 'xlim', [data(end-2000,2)/1000 data(end,2)/1000] )
end

fig3 = figure;
set(gcf, 'color', 'w', 'position', [388   417   592   420])
ax1 = axes;
h1 = plot(ax1, data(1:end,2)/1000, 27.211399*data(1:end,6), '-', 'linewidth', 1, 'color', 'b');
ylabel(ax1, 'Conserved Quantity (eV)')
xlabel(ax1, 'Time (ps)')
set(ax1, 'fontsize', 12)
if size(data,1) > 2000
    ax2 = axes('Position',[0.64 0.68 0.25 0.2]);
    h2 = plot(ax2, data(end-2000:end,2)/1000, 27.211399*data(end-2000:end,6), '-', 'linewidth',0.5, 'color', 'b');
    ylabel(ax2, 'C. Q. (eV)')
    xlabel(ax2, 'Time (ps)')
    set(ax2, 'fontsize', 8, 'xlim', [data(end-2000,2)/1000 data(end,2)/1000] )
end

fig4 = figure;
set(gcf, 'color', 'w', 'position', [388   417   592   420])
ax1 = axes;
h1 = plot(ax1, data(1:end,2)/1000, data(1:end,4), '-', 'linewidth', 1, 'color', 'k');
ylabel(ax1, 'Temperature (K)')
xlabel(ax1, 'Time (ps)')
set(ax1, 'fontsize', 12)
if size(data,1) > 2000
    ax2 = axes('Position',[0.64 0.68 0.25 0.2]);
    h2 = plot(ax2, data(end-2000:end,2)/1000, data(end-2000:end,4), '-', 'linewidth',0.5, 'color', 'k');
    ylabel(ax2, 'Temp. (K)')
    xlabel(ax2, 'Time (ps)')
    set(ax2, 'fontsize', 8, 'xlim', [data(end-2000,2)/1000 data(end,2)/1000] )
end

exportgraphics(fig1,[Basefldr Allfldrs.name '/' system '_PE.jpg'],'Resolution',300)
exportgraphics(fig2,[Basefldr Allfldrs.name '/' system '_KE.jpg'],'Resolution',300)
exportgraphics(fig3,[Basefldr Allfldrs.name '/' system '_CQ.jpg'],'Resolution',300)
exportgraphics(fig4,[Basefldr Allfldrs.name '/' system '_Temp.jpg'],'Resolution',300)
