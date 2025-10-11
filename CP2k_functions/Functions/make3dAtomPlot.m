function make3dAtomPlot(ABC, XYZ_snap, pIndx, AtomType)

Replica = [0 0 0];




figure
hold on
axis equal
set(gcf, 'Position', [200         30        600         600])

for rx = 0:Replica(1)
    for ry =0:Replica(2)
        for i = 1:size(AtomType,1)
            % set the colour of atoms
            if strcmp(AtomType(i,:), 'F')
                C = [0 1 1]; ms = 9;
            elseif contains(AtomType(i,:), 'H')
                C = [1 1 1]; ms = 5;
            elseif strcmp(AtomType(i,:), 'O')
                C = [1 0 0]; ms = 7;
            elseif strcmp(AtomType(i,:), 'Q')
                C = [1 1 0]; ms = 5;
            else
                C = [0.5 0.5 0.5]; ms = 3;
            end
            
            plot3(XYZ_snap(pIndx.(AtomType(i)), 1)+rx*ABC(1),XYZ_snap(pIndx.(AtomType(i)),2)+ry*ABC(2),XYZ_snap(pIndx.(AtomType(i)),3),'o', 'markerfacecolor', C, 'markeredgecolor', 'k', 'markersize', ms)

        end
    end
end

LwD = 1;
Vec = diag(ABC);
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

return