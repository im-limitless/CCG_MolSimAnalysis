import numpy as np
import matplotlib.pyplot as plt
    
def RadialDistribution(RadFunSnaps = None,ABC = None,ElemNames = None,plotTF = None): 
    print(np.array(['Computing radial distribution function for ',ElemNames(1,:),' and ',ElemNames(2,:)]))
    # concatenate all bond distances for all scanned atoms from all snapshots
# and sort
    dR = __builtint__.sorted(vertcat(RadFunSnaps[:,:]))
    # set bin-width
    L = 0.04
    r = np.arange(0,30+L,L)
    rho = len(dR) / (ABC(1) * ABC(2) * (ABC(3) - 9))
    
    for b in np.arange(1,len(r) - 1+1).reshape(-1):
        count[b] = sum(dR > np.logical_and(r(b),dR) < r(b + 1))
        V[b] = (4 / 3) * pi * (r(b + 1) ** 3) - (4 / 3) * pi * (r(b) ** 3)
        g[b] = (count(b) / V(b)) * (1 / rho)
    
    if np.logical_and(np.any(ismember(ElemNames,'O')),np.any(ismember(ElemNames,'H'))):
        C = np.array([0,0,0])
        yScale = np.array([0,2.5])
    else:
        if np.logical_and(np.any(ismember(ElemNames,'F')),np.any(ismember(ElemNames,'H'))):
            C = np.array([1,0,1])
            yScale = np.array([0,3])
        else:
            if sum(ismember(ElemNames,'O')) == 2:
                C = np.array([0,0,1])
                yScale = np.array([0,3.5])
            else:
                if sum(ismember(ElemNames,'H')) == 2:
                    C = np.array([0,0.5,0])
                    yScale = np.array([0,2.5])
                else:
                    if sum(ismember(ElemNames,'F')) == 2:
                        C = np.array([1,0,0])
                        yScale = np.array([0,2])
                    else:
                        C = np.array([218 / 255,165 / 255,32 / 255])
                        yScale = np.array([0,2])
    
    Minima = LocateStationaryPoints(smooth(g))
    # Minima = LocateStationaryPoints(g);
    RadialMinima = r(Minima[0] + 1)
    if plotTF == 1:
        figure
        set(gcf,'position',np.array([600,600,459,288]))
        hold('on')
        box('on')
        # plot(r(2:end), g, 'color', C, 'linewidth', 1.5);
        plt.plot(r(np.arange(2,end()+1)),smooth(g),'color',C,'linewidth',2)
        set(gca,'xlim',np.array([0,7]),'ylim',yScale,'ytick',np.array([0,0.5,1,1.5,2,2.5,3,3.5,4]),'yticklabel',np.array(['0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0']),'fontsize',14)
        plt.legend(np.array(['g_{',ElemNames(1,:),'}_{',ElemNames(2,:),'}(r)']))
        plt.xlabel('r (Ang)')
        plt.ylabel(np.array(['g_{',ElemNames(1,:),'}_{',ElemNames(2,:),'}(r)']))
        hold('off')
    
    return RadialMinima
    return RadialMinima