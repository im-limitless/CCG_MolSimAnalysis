import numpy as np
    
def convertXSDtoXYZ(BaseOutFldr = None,fldrname = None,NAtoms = None,n = None,srtdAtomPosition = None,srtdAtomElement = None,Vec = None): 
    newlinechar = char(10)
    while srtdAtomPosition(srtdAtomPosition < - 0.01):

        srtdAtomPosition[srtdAtomPosition < - 0.01] = srtdAtomPosition(srtdAtomPosition < - 0.01) + 1

    
    while srtdAtomPosition(srtdAtomPosition > 1.01):

        srtdAtomPosition[srtdAtomPosition > 1.01] = srtdAtomPosition(srtdAtomPosition > 1.01) - 1

    
    fidout = open(np.array([BaseOutFldr,fldrname,'\',fldrname,'.xyz']),'w')
    AscAtomPosition = []
    for jj in np.arange(1,len(n)+1).reshape(-1):
        AscAtomPosition = np.array([[AscAtomPosition],[sortrows(srtdAtomPosition(np.arange(sum(n(np.arange(1,jj+1))) - n(jj) + 1,sum(n(np.arange(1,jj+1)))+1),:),3)]])
    
    H,O,F,Pt = determineVelocityDistribution
    fidout.write(np.array([num2str(NAtoms),newlinechar]) % ())
    fidout.write(np.array([' i =   input geom',newlinechar]) % ())
    Fix = []
    Restrain = []
    Vel = np.zeros((NAtoms,3))
    StockElem = cell(NAtoms,1)
    for kk in np.arange(1,NAtoms+1).reshape(-1):
        if cellfun(strcmp,srtdAtomElement(kk),np.array(['Au'])):
            srtdAtomElement[kk] = np.array(['Pts'])
            for v in np.arange(1,3+1).reshape(-1):
                Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
            StockElem[kk] = 'Pt'
            # # # # #             # Margherita
# # # # #             srtdAtomElement(kk) = {'Au'};
# # # # #             for v=1:3
# # # # #                 Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
# # # # #             end
# # # # #             StockElem{kk} = 'Au';
        else:
            if cellfun(strcmp,srtdAtomElement(kk),np.array(['Ag'])):
                srtdAtomElement[kk] = np.array(['Ptb'])
                Vel[kk,np.arange[1,3+1]] = 0
                Fix = np.array([[Fix],[kk]])
                StockElem[kk] = 'Pt'
                # # # # #             # Margherita
# # # # #             srtdAtomElement(kk) = {'Ag'};
# # # # #             for v=1:3
# # # # #                 Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
# # # # #             end
# # # # #             StockElem{kk} = 'Ag';
# # # # #
# # # # #             elseif cellfun(@strcmp, srtdAtomElement(kk), {'Cr'})
# # # # #                         srtdAtomElement(kk) = {'Agb'};
# # # # #                         Vel(kk,1:3) = 0;
# # # # #                         Fix = [Fix; kk];
# # # # #                         StockElem{kk} = 'Ag';
# # # # #
# # # # #                         elseif cellfun(@strcmp, srtdAtomElement(kk), {'W'})
# # # # #                         srtdAtomElement(kk) = {'Agb'};
# # # # #                         Vel(kk,1:3) = 0;
# # # # #                         Fix = [Fix; kk];
# # # # #                         StockElem{kk} = 'Ag';
            else:
                if cellfun(strcmp,srtdAtomElement(kk),np.array(['Cu'])):
                    srtdAtomElement[kk] = np.array(['Ptss'])
                    for v in np.arange(1,3+1).reshape(-1):
                        Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                    StockElem[kk] = 'Pt'
                else:
                    if cellfun(strcmp,srtdAtomElement(kk),np.array(['Rg'])):
                        srtdAtomElement[kk] = np.array(['PtE'])
                        for v in np.arange(1,3+1).reshape(-1):
                            Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                        StockElem[kk] = 'Pt'
                    else:
                        if cellfun(strcmp,srtdAtomElement(kk),np.array(['B'])):
                            srtdAtomElement[kk] = np.array(['Hsurf'])
                            for v in np.arange(1,3+1).reshape(-1):
                                Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                            Restrain = np.array([[Restrain],[kk]])
                            #             Fix = [Fix; kk];
                            StockElem[kk] = 'H'
                        else:
                            if cellfun(strcmp,srtdAtomElement(kk),np.array(['H'])):
                                for v in np.arange(1,3+1).reshape(-1):
                                    Vel[kk,v] = H.mid(find(np.random.rand(1) * sum(H.N) < cumsum(H.N) == 1,1)) + (H.edges(2) - H.edges(1)) * (2 * np.random.rand(1) - 1)
                                StockElem[kk] = 'H'
                            else:
                                if cellfun(strcmp,srtdAtomElement(kk),np.array(['S'])):
                                    srtdAtomElement[kk] = np.array(['Ohy'])
                                    for v in np.arange(1,3+1).reshape(-1):
                                        Vel[kk,v] = O.mid(find(np.random.rand(1) * sum(O.N) < cumsum(O.N) == 1,1)) + (O.edges(2) - O.edges(1)) * (2 * np.random.rand(1) - 1)
                                    StockElem[kk] = 'O'
                                else:
                                    if cellfun(strcmp,srtdAtomElement(kk),np.array(['O'])):
                                        for v in np.arange(1,3+1).reshape(-1):
                                            Vel[kk,v] = O.mid(find(np.random.rand(1) * sum(O.N) < cumsum(O.N) == 1,1)) + (O.edges(2) - O.edges(1)) * (2 * np.random.rand(1) - 1)
                                        StockElem[kk] = 'O'
                                        #             srtdAtomElement(kk) = {'Ohy'};
                                    else:
                                        if cellfun(strcmp,srtdAtomElement(kk),np.array(['Li'])):
                                            srtdAtomElement[kk] = np.array(['Hhy'])
                                            for v in np.arange(1,3+1).reshape(-1):
                                                Vel[kk,v] = H.mid(find(np.random.rand(1) * sum(H.N) < cumsum(H.N) == 1,1)) + (H.edges(2) - H.edges(1)) * (2 * np.random.rand(1) - 1)
                                            StockElem[kk] = 'H'
                                            #             srtdAtomElement(kk) = {'Hhy'};
                                        else:
                                            if cellfun(strcmp,srtdAtomElement(kk),np.array(['F'])):
                                                srtdAtomElement[kk] = np.array(['F'])
                                                for v in np.arange(1,3+1).reshape(-1):
                                                    Vel[kk,v] = F.mid(find(np.random.rand(1) * sum(F.N) < cumsum(F.N) == 1,1)) + (F.edges(2) - F.edges(1)) * (2 * np.random.rand(1) - 1)
                                                StockElem[kk] = 'F'
                                            else:
                                                if cellfun(strcmp,srtdAtomElement(kk),np.array(['Cl'])):
                                                    srtdAtomElement[kk] = np.array(['Cl'])
                                                    for v in np.arange(1,3+1).reshape(-1):
                                                        Vel[kk,v] = F.mid(find(np.random.rand(1) * sum(F.N) < cumsum(F.N) == 1,1)) + (F.edges(2) - F.edges(1)) * (2 * np.random.rand(1) - 1)
                                                    StockElem[kk] = 'Cl'
                                                else:
                                                    if cellfun(strcmp,srtdAtomElement(kk),np.array(['Se'])):
                                                        srtdAtomElement[kk] = np.array(['OtU'])
                                                        for v in np.arange(1,3+1).reshape(-1):
                                                            Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                        StockElem[kk] = 'O'
                                                    else:
                                                        if cellfun(strcmp,srtdAtomElement(kk),np.array(['Te'])):
                                                            srtdAtomElement[kk] = np.array(['OtL'])
                                                            for v in np.arange(1,3+1).reshape(-1):
                                                                Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                            StockElem[kk] = 'O'
                                                        else:
                                                            if cellfun(strcmp,srtdAtomElement(kk),np.array(['Po'])):
                                                                srtdAtomElement[kk] = np.array(['OtS'])
                                                                for v in np.arange(1,3+1).reshape(-1):
                                                                    Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                                StockElem[kk] = 'O'
                                                            else:
                                                                if cellfun(strcmp,srtdAtomElement(kk),np.array(['Sb'])):
                                                                    srtdAtomElement[kk] = np.array(['OtSS'])
                                                                    for v in np.arange(1,3+1).reshape(-1):
                                                                        Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                                    StockElem[kk] = 'O'
                                                                else:
                                                                    if cellfun(strcmp,srtdAtomElement(kk),np.array(['Pb'])):
                                                                        srtdAtomElement[kk] = np.array(['Om'])
                                                                        for v in np.arange(1,3+1).reshape(-1):
                                                                            Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                                        StockElem[kk] = 'O'
                                                                    else:
                                                                        if cellfun(strcmp,srtdAtomElement(kk),np.array(['Bi'])):
                                                                            srtdAtomElement[kk] = np.array(['Ob'])
                                                                            Vel[kk,np.arange[1,3+1]] = 0
                                                                            Fix = np.array([[Fix],[kk]])
                                                                            StockElem[kk] = 'O'
                                                                        else:
                                                                            if cellfun(strcmp,srtdAtomElement(kk),np.array(['Mt'])):
                                                                                srtdAtomElement[kk] = np.array(['Rus'])
                                                                                for v in np.arange(1,3+1).reshape(-1):
                                                                                    Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                                                StockElem[kk] = 'Ru'
                                                                            else:
                                                                                if cellfun(strcmp,srtdAtomElement(kk),np.array(['Ds'])):
                                                                                    srtdAtomElement[kk] = np.array(['Rub'])
                                                                                    Vel[kk,np.arange[1,3+1]] = 0
                                                                                    Fix = np.array([[Fix],[kk]])
                                                                                    StockElem[kk] = 'Ru'
                                                                                else:
                                                                                    if cellfun(strcmp,srtdAtomElement(kk),np.array(['Hs'])):
                                                                                        srtdAtomElement[kk] = np.array(['Rum'])
                                                                                        Vel[kk,np.arange[1,3+1]] = 0
                                                                                        Fix = np.array([[Fix],[kk]])
                                                                                        StockElem[kk] = 'Ru'
                                                                                        # Rashid
                                                                                    else:
                                                                                        if cellfun(strcmp,srtdAtomElement(kk),np.array(['Ti'])):
                                                                                            srtdAtomElement[kk] = np.array(['Alb'])
                                                                                            Vel[kk,np.arange[1,3+1]] = 0
                                                                                            Fix = np.array([[Fix],[kk]])
                                                                                            StockElem[kk] = 'Al'
                                                                                        else:
                                                                                            if cellfun(strcmp,srtdAtomElement(kk),np.array(['Zr'])):
                                                                                                srtdAtomElement[kk] = np.array(['Al2'])
                                                                                                for v in np.arange(1,3+1).reshape(-1):
                                                                                                    Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                                                                StockElem[kk] = 'Al'
                                                                                            else:
                                                                                                if cellfun(strcmp,srtdAtomElement(kk),np.array(['Hf'])):
                                                                                                    srtdAtomElement[kk] = np.array(['Al1'])
                                                                                                    for v in np.arange(1,3+1).reshape(-1):
                                                                                                        Vel[kk,v] = Pt.mid(find(np.random.rand(1) * sum(Pt.N) < cumsum(Pt.N) == 1,1)) + (Pt.edges(2) - Pt.edges(1)) * (2 * np.random.rand(1) - 1)
                                                                                                    StockElem[kk] = 'Al'
                                                                                                else:
                                                                                                    if cellfun(strcmp,srtdAtomElement(kk),np.array(['Mo'])):
                                                                                                        srtdAtomElement[kk] = np.array(['Mo'])
                                                                                                        Vel[kk,np.arange[1,3+1]] = 0
                                                                                                        Fix = np.array([[Fix],[kk]])
                                                                                                        StockElem[kk] = 'Mo'
        CartPos = AscAtomPosition(kk,:) * Vec
        fidout.write(np.array([pad(srtdAtomElement[kk],8),pad(num2str(CartPos(1),'%.10g'),20),pad(num2str(CartPos(2),'%.10g'),20),pad(num2str(CartPos(3),'%.10g'),20),newlinechar]) % ())
    
    fidout.close()
    ListOfElems = np.array([StockElem(cumsum(n)),np.transpose(srtdAtomElement(cumsum(n)))])
    return ListOfElems,Fix,Restrain,Vel