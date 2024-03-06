function [AveTotDen]= Original_of_getBulkMacroscopicAve(TotDen, z, ABC, BaseFldr, system)

PathSep =  setOSpathSep;

AveTotDen= mean(TotDen,2);
counter=0;
Macro_AveDens=[];

for k=[3:1.5:10]
    AveTotDen= mean(TotDen,2);
    limits= [ABC(3)/2-k ABC(3)/2+k];
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
    title('M1 Box limits ',k);

%     Temp_den= AveTotDen;

    j=3; %this is the value for the #voxels to be smoothed
    for i=1:iterations

       
%         Temp_den= smooth(Temp_den(z_indx),j);
        AveTotDen(z_indx)= smooth(AveTotDen(z_indx),j);
        plot(z,AveTotDen,  'linewidth', 1.5)
        if i==iterations
            text(z(end)/2, 1500, ['Bulk Density(last)= ' num2str(median(AveTotDen(z_indx)),'%.0f') ' kgm^{-3}'],'HorizontalAlignment','center','FontSize',14)
        end
        Macro_AveDens_T=[k,j,median(AveTotDen(z_indx)), i];
        Macro_AveDens=[Macro_AveDens; Macro_AveDens_T];
        Box_Range=Macro_AveDens(:,1);
        Ave_vox=Macro_AveDens(:,2);
        Bulk_dens=Macro_AveDens(:,3);
        Iterations=Macro_AveDens(:,4);


    end
end

Allfldrs = dir([BaseFldr '*' system]);

T = table( Box_Range, Ave_vox, Iterations, Bulk_dens);
T.Properties.VariableNames = {'Box_range(Å)', 'Averaged Voxels', 'Iteration', 'Bulk Density (kgm^-3)'};
writetable(T, [BaseFldr Allfldrs.name PathSep system '_Macro_M1.xlsx'])


return