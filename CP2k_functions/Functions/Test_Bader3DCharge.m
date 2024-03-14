function Bader3DCharge(XYZ, ABC, Qmc) 

disp('Creating 3D charge distribution map...');

[xs,ys,zs] = sphere(10,100);   %//sphere(n) will give an sphere with n X n faces (cross sectional and vertical lines)
Radius = 1;

figure('Color','w')
hold on
axis equal

[Qmsort, QmIndx] = sort(Qmc, 'Ascend');

%% %%% Quick Color refs. %%%%%%%%%%%%%%%%%%%
% cmap(1,:) = [1 0 0];   %// red
% cmap(2,:) = [255 255 0]/255;   %// yellow
% cmap(2,:) = [0 0.8 0];   %// green
% cmap(2,:) = [0.8 0.8 0.8];   %// white
% cmap(2,:) = [0.5 0 0.5];   %// purple
% cmap(3,:) = [0 0 1];   %// blue
% cmap(1,:) = [0.7216    0.9098    0.9882];   %// even lighter blue
% cmap(2,:) = [0    0.4471    0.7412];   %// lighter blue
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% We want to split Qmc based on +ve and -ve values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
PQmc=[];    %store +ve 
NQmc=[];    %store -ve

for j=1:length(Qmc)
    if (Qmc(j) >=0)
        PQmc=[PQmc;Qmc(j)];
    else
        NQmc=[NQmc;Qmc(j)];
    end
end

%%% +ve colors
n = length(PQmc);                %// number of colors

cmap(1,:) = [0.8 0.8 0.8];   %// white
cmap(2,:) = [0.7216    0.9098    0.9882];   %// even lighter blue
cmap(3,:) = [0 0 1];   %// blue

[X,Y] = meshgrid([1:3],[1:n]);  %// mesh of indices

Colors = interp2(X([1,floor(n/2),n],:),Y([1,floor(n/2),n],:),cmap,X,Y); %// interpolate colormap
%// on "interp2(X,Y,V,Xq,Yq)": X and Y contains specific locations to be
%//replaced by V and then interpolation happens between these points and are
%//placed in Xq and Yq

% colormap(cmap) %// set color map

% Colors = [ones(length(Qmc),1) linspace(0,1,length(Qmc))' zeros(length(Qmc),1)];
colormap(Colors);


%%% -ve colors
n = length(NQmc);                %// number of colors

cmap(1,:) = [1 0 0];   %// red
cmap(2,:) = [ 1.0000    0.4314    0.4314];   %// lighter red
cmap(3,:) = [1.0000    0.7412    0.7412];   %// very light red

[X,Y] = meshgrid([1:3],[1:n]);  %// mesh of indices

Colors2 = interp2(X([1,floor(n/2),n],:),Y([1,floor(n/2),n],:),cmap,X,Y); %// interpolate colormap
%// on "interp2(X,Y,V,Xq,Yq)": X and Y contains specific locations to be
%//replaced by V and then interpolation happens between these points and are
%//placed in Xq and Yq


% colormap(cmap) %// set color map

% Colors = [ones(length(Qmc),1) linspace(0,1,length(Qmc))' zeros(length(Qmc),1)];
colormap(Colors2);


%%% Combo
Colors_all=[Colors;Colors2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(XYZ, [1])
    
    xAt = xs*Radius + XYZ(i,1);
    yAt = ys*Radius + XYZ(i,2);
    zAt = zs*Radius + XYZ(i,3);

    if (Qmc(i) >=0)
    
        %     cIndx = find(QmIndx == i);
        cIndx = length(Colors)*(Qmc(i)-min(Qmc))/(max(Qmc)-min(Qmc));
        C=round(cIndx);
        if C == 0
            C=1;
        end
        hsurf = surf(xAt,yAt,zAt,'FaceColor', Colors(C,:),'EdgeColor','None');
        set(hsurf,'AmbientStrength',0.1,...
            'LineStyle','none',...
            'FaceLighting','phong',...
            'FaceAlpha',1.0,...
            'AmbientStrength',0.2,...
            'DiffuseStrength',0.9,...
            'SpecularStrength',0.7);
    else
        %     cIndx = find(QmIndx == i);
        cIndx = length(Colors2)*(Qmc(i)-min(Qmc))/(max(Qmc)-min(Qmc));
        C=round(cIndx);
        if C == 0
            C=1;
        end
        hsurf = surf(xAt,yAt,zAt,'FaceColor', Colors2(C,:),'EdgeColor','None');
        set(hsurf,'AmbientStrength',0.1,...
            'LineStyle','none',...
            'FaceLighting','phong',...
            'FaceAlpha',1.0,...
            'AmbientStrength',0.2,...
            'DiffuseStrength',0.9,...
            'SpecularStrength',0.7);
    end


% xAt = xs*Radius + XYZ(i,1) + ABC(1);
% 
%     
%     cIndx = find(QmIndx == i);
%  
%     hsurf = surf(xAt,yAt,zAt,'FaceColor', Colors(cIndx,:),'EdgeColor','None');
%     set(hsurf,'AmbientStrength',0.1,...
%     'LineStyle','none',...
%     'FaceLighting','phong',...
%     'FaceAlpha',1.0,...
%     'AmbientStrength',0.2,...
%     'DiffuseStrength',0.9,...
%     'SpecularStrength',0.7);
% 
% 
% xAt = xs*Radius + XYZ(i,1) - ABC(1);
% 
%     
%     cIndx = find(QmIndx == i);
%  
%     hsurf = surf(xAt,yAt,zAt,'FaceColor', Colors(cIndx,:),'EdgeColor','None');
%     set(hsurf,'AmbientStrength',0.1,...
%     'LineStyle','none',...
%     'FaceLighting','phong',...
%     'FaceAlpha',1.0,...
%     'AmbientStrength',0.2,...
%     'DiffuseStrength',0.9,...
%     'SpecularStrength',0.7);

%     text(XYZ(i,1), XYZ(i,2)-1, XYZ(i,3), num2str(PtIndx(i)), 'fontsize', 8)
end

caxis([min(Qmc) max(Qmc)])
hcb = colorbar;


colorTitleHandle = get(hcb,'Title');
titleString = 'Bader Charge (e)';
set(colorTitleHandle ,'String',titleString);

Vec = diag(ABC);
LwD = 2;

plot3([0 Vec(1,1)],[0 Vec(1,2)],[0 Vec(1,3)],'-k','LineWidth',LwD)
plot3([0 Vec(2,1)],[0 Vec(2,2)],[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3([0 Vec(3,1)],[0 Vec(3,2)],[0 Vec(3,3)],'-k','LineWidth',LwD)

plot3(Vec(1,1)+[0 Vec(2,1)],Vec(1,2)+[0 Vec(2,2)],Vec(1,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(2,1)+[0 Vec(1,1)],Vec(2,2)+[0 Vec(1,2)],Vec(2,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)

plot3(Vec(1,1)+[0 Vec(3,1)],Vec(1,2)+[0 Vec(3,2)],Vec(1,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)
plot3(Vec(2,1)+[0 Vec(3,1)],Vec(2,2)+[0 Vec(3,2)],Vec(2,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)
plot3(Vec(1,1)+Vec(2,1)+[0 Vec(3,1)],Vec(1,2)+Vec(2,2)+[0 Vec(3,2)],Vec(1,3)+Vec(2,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)

plot3(Vec(3,1)+[0 Vec(1,1)],Vec(3,2)+[0 Vec(1,2)],Vec(3,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+[0 Vec(2,1)],Vec(3,2)+[0 Vec(2,2)],Vec(3,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+Vec(1,1)+[0 Vec(2,1)],Vec(3,2)+Vec(1,2)+[0 Vec(2,2)],Vec(3,3)+Vec(1,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+Vec(2,1)+[0 Vec(1,1)],Vec(3,2)+Vec(2,2)+[0 Vec(1,2)],Vec(3,3)+Vec(2,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)
% light

view(11.2545,8.5505)
% set(gca,'CameraUpVector',[2 1 1],'CameraViewAngle',[4.]);
% set(gca,'CameraUpVector',[0 0 1],'CameraViewAngle',[6.5],'CameraPosition',[-5.3 -77.0 81.0]);
set(gcf,'Position',[599   393   803   600]);
set(gca, 'fontsize', 14);

xlabel('x (Ang)')
ylabel('y (Ang)')
zlabel('z (Ang)')
hold off


axis equal




