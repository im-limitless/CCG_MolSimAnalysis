function Bader3DCharge_fixV2(XYZ, ABC, Qmc, SelectedIndices) 
% Bader3DCharge(XYZ, ABC, Qmc) plots all atoms with charge-based coloring.
% Bader3DCharge(XYZ, ABC, Qmc, SelectedIndices) plots all atoms, but only colors
% the ones specified in SelectedIndices. Other atoms are semi-transparent with borders.

disp('Creating 3D charge distribution map...');

% [xs,ys,zs] = sphere(10,100);   % Sphere with 10x10 faces (version 2023b)
[xs,ys,zs] = sphere(10);   % Sphere with 10x10 faces (post version 2023b)
Radius = 1;

% Handle optional SelectedIndices parameter
if nargin < 4 || isempty(SelectedIndices)
    SelectedIndices = 1:size(XYZ,1); % Default: all atoms
end

% Validate SelectedIndices and create logical mask
isSelected = false(size(XYZ,1), 1);
isSelected(SelectedIndices) = true;
Qmc_selected = Qmc(isSelected);

figure('Color','w')
hold on
axis equal

% Generate colormap only if there are selected atoms
if ~isempty(Qmc_selected)
    min_Q = min(Qmc_selected);
    max_Q = max(Qmc_selected);
    
    % Handle uniform charges
    if max_Q == min_Q
        if max_Q > 0
            cmap = [0 0 1]; % Blue
        elseif max_Q < 0
            cmap = [1 0 0]; % Red
        else
            cmap = [0 0.8 0]; % Green
        end
        cmap_length = 1;
    else
        has_neg = any(Qmc_selected < 0);
        has_pos = any(Qmc_selected > 0);
        
        % Define colormap anchors based on charge signs
        if has_neg && has_pos
            x_points = [min_Q, 0, max_Q];
            color_points = [1 0 0; 0 0.8 0; 0 0 1]; % Red-Green-Blue
        elseif has_pos
            x_points = [0, max_Q];
            color_points = [0 0.8 0; 0 0 1]; % Green-Blue
        else
            x_points = [min_Q, 0];
            color_points = [1 0 0; 0 0.8 0]; % Red-Green
        end
        
        cmap_length = 256;
        xi = linspace(min_Q, max_Q, cmap_length);
        cmap = interp1(x_points, color_points, xi, 'linear', 'extrap');
    end
    
    colormap(cmap);
    caxis([min_Q, max_Q]);
    
    % Add colorbar
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Bader Charge (e)';
    set(colorTitleHandle ,'String',titleString);
else
    warning('No atoms selected for heat map. All atoms will be gray.');
end

% First plot unselected atoms with transparency
for i = find(~isSelected)'
    xAt = xs*Radius + XYZ(i,1);
    yAt = ys*Radius + XYZ(i,2);
    zAt = zs*Radius + XYZ(i,3);
    
    hsurf = surf(xAt, yAt, zAt,...
        'FaceColor', [0.7 0.7 0.7],...  % Light gray
        'EdgeColor', [0.4 0.4 0.4],...   % Darker gray edges
        'FaceAlpha', 0.3);               % Semi-transparent
    
    set(hsurf,...
        'LineWidth', 0.3,...            % Thin borders
        'AmbientStrength', 0.3,...       % Reduced lighting
        'DiffuseStrength', 0.4,...       % Softer shading
        'SpecularStrength', 0.1);        % Minimal highlights
end

% Then plot selected atoms with full color
for i = find(isSelected)'
    xAt = xs*Radius + XYZ(i,1);
    yAt = ys*Radius + XYZ(i,2);
    zAt = zs*Radius + XYZ(i,3);
    
    if ~isempty(Qmc_selected)
        if max_Q == min_Q
            C = 1;
        else
            normalized = (Qmc(i) - min_Q) / (max_Q - min_Q);
            C = round(normalized * (cmap_length - 1)) + 1;
            C = max(1, min(C, cmap_length));
        end
        atomColor = cmap(C,:);
    else
        atomColor = [0.7 0.7 0.7]; % Fallback gray
    end
    
    hsurf = surf(xAt, yAt, zAt,...
        'FaceColor', atomColor,...
        'EdgeColor', 'none');           % No edges for highlighted atoms
    
    set(hsurf,...
        'AmbientStrength', 0.1,...
        'DiffuseStrength', 0.9,...      % Stronger shading
        'SpecularStrength', 0.5,...     % More reflective
        'FaceAlpha', 1.0);              % Fully opaque
end

% Plot cell vectors (unchanged)
Vec = diag(ABC);
LwD = 2;

plot3([0 Vec(1,1)],[0 Vec(1,2)],[0 Vec(1,3)],'-k','LineWidth',LwD)
plot3([0 Vec(2,1)],[0 Vec(2,2)],[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3([0 Vec(3,1)],[0 Vec(3,2)],[0 Vec(3,3)],'-k','LineWidth',LwD)

% Additional cell vector plotting code (unchanged)
% ...
plot3(Vec(1,1)+[0 Vec(2,1)],Vec(1,2)+[0 Vec(2,2)],Vec(1,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(2,1)+[0 Vec(1,1)],Vec(2,2)+[0 Vec(1,2)],Vec(2,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)

plot3(Vec(1,1)+[0 Vec(3,1)],Vec(1,2)+[0 Vec(3,2)],Vec(1,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)
plot3(Vec(2,1)+[0 Vec(3,1)],Vec(2,2)+[0 Vec(3,2)],Vec(2,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)
plot3(Vec(1,1)+Vec(2,1)+[0 Vec(3,1)],Vec(1,2)+Vec(2,2)+[0 Vec(3,2)],Vec(1,3)+Vec(2,3)+[0 Vec(3,3)],'-k','LineWidth',LwD)

plot3(Vec(3,1)+[0 Vec(1,1)],Vec(3,2)+[0 Vec(1,2)],Vec(3,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+[0 Vec(2,1)],Vec(3,2)+[0 Vec(2,2)],Vec(3,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+Vec(1,1)+[0 Vec(2,1)],Vec(3,2)+Vec(1,2)+[0 Vec(2,2)],Vec(3,3)+Vec(1,3)+[0 Vec(2,3)],'-k','LineWidth',LwD)
plot3(Vec(3,1)+Vec(2,1)+[0 Vec(1,1)],Vec(3,2)+Vec(2,2)+[0 Vec(1,2)],Vec(3,3)+Vec(2,3)+[0 Vec(1,3)],'-k','LineWidth',LwD)

% Final formatting
view(11.2545,8.5505)
set(gcf,'Position',[599   393   803   600]);
set(gca, 'fontsize', 14);

xlabel('x (Ang)')
ylabel('y (Ang)')
zlabel('z (Ang)')
hold off
axis equal