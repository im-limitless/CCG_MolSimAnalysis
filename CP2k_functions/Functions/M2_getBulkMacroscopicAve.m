function [AveTotDen]= getBulkMacroscopicAve(TotDen, z, ABC, BaseFldr, system)

PathSep =  setOSpathSep;


counter=0;
Macro_AveDens=[];
for k=[3:1.5:10]
    AveTotDen= mean(TotDen,2);
    counter=counter+1;
    limits= [ABC(3)/2-k ABC(3)/2+k];
    % limits= [ABC(3)/2-7.5 ABC(3)/2+7.5];
    % limits= [8.81 34.35];

    z_indx= find(z>=limits(1) & z<=limits(2));

    iterations=1;

    figure
    hold on
%     my_legend = legend();
    set(gcf, 'position', [377         423        1123         420]);
    xlabel('z (Ang)');
    ylabel('Density (kgm^{-3})');
    set(gca, 'xlim', [0 ABC(3)], 'ylim', [0 2500],'FontSize',14);

    plot(z,mean(TotDen,2),  'linewidth', 1.5, 'color', 'k')
    title('M2 Box limits ',k);
    Temp_den= AveTotDen;
    
   
    for i=1:iterations
        n=0;
        for j=1:20
            n=n+1;
            Macro_AveDens_T=[k,j,median(AveTotDen(z_indx))];
            Temp_den= smooth(Temp_den,j);
            AveTotDen(z_indx)= Temp_den(z_indx);
            plot(z,AveTotDen,  'linewidth', 1.5)
            if i==iterations & j==20
                text(z(end)/2, 1500, ['Bulk Density (last)= ' num2str(median(AveTotDen(z_indx)),'%.0f') ' kgm^{-3}'],'HorizontalAlignment','center','FontSize',14)
            end
            Macro_AveDens=[Macro_AveDens; Macro_AveDens_T];
            Box_Range=Macro_AveDens(:,1);
            Ave_vox=Macro_AveDens(:,2);
            Bulk_dens=Macro_AveDens(:,3);
        end
    end
end

Allfldrs = dir([BaseFldr '*' system]);

T = table( Box_Range, Ave_vox, Bulk_dens);
T.Properties.VariableNames = {'Box_range(Å)', 'Averaged Voxels', 'Bulk Density (kgm^-3)'};
writetable(T, [BaseFldr Allfldrs.name PathSep system '_Macro_M2.xlsx'])

return