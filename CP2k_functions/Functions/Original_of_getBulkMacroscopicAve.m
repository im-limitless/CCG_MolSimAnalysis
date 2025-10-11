function [AveTotDen]= Original_of_getBulkMacroscopicAve(TotDen, z, ABC)

AveTotDen= mean(TotDen,2);

limits= [ABC(3)/2-7.5 ABC(3)/2+7.5];
% limits= [9.9 30.9];

z_indx= find(z>=limits(1) & z<=limits(2));

iterations=20;

figure
hold on
set(gcf, 'position', [377         423        1123         420]);
xlabel('z (Ang)');
ylabel('Density (kgm^{-3})');
set(gca, 'xlim', [0 ABC(3)], 'ylim', [0 2500],'FontSize',14);

plot(z,mean(TotDen,2),  'linewidth', 1.5, 'color', 'k')


for i=1:iterations

    AveTotDen(z_indx)= smooth(AveTotDen(z_indx),3);
    plot(z,AveTotDen,  'linewidth', 1.5)
    if i==iterations
       text(z(end)/2, 1500, ['Bulk Density= ' num2str(median(AveTotDen(z_indx)),'%.0f') ' kgm^{-3}'],'HorizontalAlignment','center','FontSize',14)
    end
end



return