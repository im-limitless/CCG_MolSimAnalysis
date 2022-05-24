import matplotlib.pyplot as plt
import numpy as np
    
def Bader3DCharge(XYZ = None,ABC = None,Qmc = None): 
    xs,ys,zs = sphere(10,100)
    Radius = 1
    plt.figure('Color','w')
    hold('on')
    plt.axis('equal')
    Qmsort,QmIndx = __builtint__.sorted(Qmc,'Ascend')
    Colors = np.array([np.ones((len(Qmc),1)),np.transpose(np.linspace(0,1,len(Qmc))),np.zeros((len(Qmc),1))])
    colormap(Colors)
    for i in np.arange(1,XYZ.shape[np.array([1])-1]+1).reshape(-1):
        xAt = xs * Radius + XYZ(i,1)
        yAt = ys * Radius + XYZ(i,2)
        zAt = zs * Radius + XYZ(i,3)
        cIndx = find(QmIndx == i)
        hsurf = surf(xAt,yAt,zAt,'FaceColor',Colors(cIndx,:),'EdgeColor','None')
        set(hsurf,'AmbientStrength',0.1,'LineStyle','none','FaceLighting','phong','FaceAlpha',1.0,'AmbientStrength',0.2,'DiffuseStrength',0.9,'SpecularStrength',0.7)
        #     text(XYZ(i,1), XYZ(i,2)-1, XYZ(i,3), num2str(PtIndx(i)), 'fontsize', 8)
    
    caxis(np.array([np.amin(Qmc),np.amax(Qmc)]))
    hcb = colorbar
    colorTitleHandle = get(hcb,'Title')
    titleString = 'Bader Charge (e)'
    set(colorTitleHandle,'String',titleString)
    Vec = diag(ABC)
    LwD = 2
    plot3(np.array([0,Vec(1,1)]),np.array([0,Vec(1,2)]),np.array([0,Vec(1,3)]),'-k','LineWidth',LwD)
    plot3(np.array([0,Vec(2,1)]),np.array([0,Vec(2,2)]),np.array([0,Vec(2,3)]),'-k','LineWidth',LwD)
    plot3(np.array([0,Vec(3,1)]),np.array([0,Vec(3,2)]),np.array([0,Vec(3,3)]),'-k','LineWidth',LwD)
    plot3(Vec(1,1) + np.array([0,Vec(2,1)]),Vec(1,2) + np.array([0,Vec(2,2)]),Vec(1,3) + np.array([0,Vec(2,3)]),'-k','LineWidth',LwD)
    plot3(Vec(2,1) + np.array([0,Vec(1,1)]),Vec(2,2) + np.array([0,Vec(1,2)]),Vec(2,3) + np.array([0,Vec(1,3)]),'-k','LineWidth',LwD)
    plot3(Vec(1,1) + np.array([0,Vec(3,1)]),Vec(1,2) + np.array([0,Vec(3,2)]),Vec(1,3) + np.array([0,Vec(3,3)]),'-k','LineWidth',LwD)
    plot3(Vec(2,1) + np.array([0,Vec(3,1)]),Vec(2,2) + np.array([0,Vec(3,2)]),Vec(2,3) + np.array([0,Vec(3,3)]),'-k','LineWidth',LwD)
    plot3(Vec(1,1) + Vec(2,1) + np.array([0,Vec(3,1)]),Vec(1,2) + Vec(2,2) + np.array([0,Vec(3,2)]),Vec(1,3) + Vec(2,3) + np.array([0,Vec(3,3)]),'-k','LineWidth',LwD)
    plot3(Vec(3,1) + np.array([0,Vec(1,1)]),Vec(3,2) + np.array([0,Vec(1,2)]),Vec(3,3) + np.array([0,Vec(1,3)]),'-k','LineWidth',LwD)
    plot3(Vec(3,1) + np.array([0,Vec(2,1)]),Vec(3,2) + np.array([0,Vec(2,2)]),Vec(3,3) + np.array([0,Vec(2,3)]),'-k','LineWidth',LwD)
    plot3(Vec(3,1) + Vec(1,1) + np.array([0,Vec(2,1)]),Vec(3,2) + Vec(1,2) + np.array([0,Vec(2,2)]),Vec(3,3) + Vec(1,3) + np.array([0,Vec(2,3)]),'-k','LineWidth',LwD)
    plot3(Vec(3,1) + Vec(2,1) + np.array([0,Vec(1,1)]),Vec(3,2) + Vec(2,2) + np.array([0,Vec(1,2)]),Vec(3,3) + Vec(2,3) + np.array([0,Vec(1,3)]),'-k','LineWidth',LwD)
    # light
    
    view(0,0)
    # set(gca,'CameraUpVector',[2 1 1],'CameraViewAngle',[4.]);
# set(gca,'CameraUpVector',[0 0 1],'CameraViewAngle',[6.5],'CameraPosition',[-5.3 -77.0 81.0]);
    set(gcf,'Position',np.array([599,393,803,600]))
    set(gca,'fontsize',14)
    plt.xlabel('x (Ang)')
    plt.ylabel('y (Ang)')
    plt.zlabel('z (Ang)')
    hold('off')
    plt.axis('equal')