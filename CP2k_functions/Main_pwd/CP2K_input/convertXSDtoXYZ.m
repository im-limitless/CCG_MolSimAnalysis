function [ListOfElems, Fix, Restrain, Vel] = convertXSDtoXYZ(BaseOutFldr, fldrname,  NAtoms, n, srtdAtomPosition, srtdAtomElement, Vec)

newlinechar = char(10);

while srtdAtomPosition(srtdAtomPosition < -0.01)
    srtdAtomPosition(srtdAtomPosition < -0.01) =  srtdAtomPosition(srtdAtomPosition < -0.01)+1;
end
while srtdAtomPosition(srtdAtomPosition > 1.01)
    srtdAtomPosition(srtdAtomPosition > 1.01) =  srtdAtomPosition(srtdAtomPosition > 1.01)-1;
end

fidout = fopen([BaseOutFldr fldrname '/' fldrname '.xyz'],'w');

AscAtomPosition = [];
for jj = 1:length(n)
    AscAtomPosition = [AscAtomPosition; sortrows(srtdAtomPosition(sum(n(1:jj))-n(jj)+1:sum(n(1:jj)),:),3)];
end

[H, O, F, Pt] = determineVelocityDistribution;

fprintf(fidout,[num2str(NAtoms) newlinechar]);
fprintf(fidout,[' i =   input geom' newlinechar]);
Fix = [];
Restrain = [];
Vel = zeros(NAtoms,3);
StockElem = cell(NAtoms,1);

for kk = 1:NAtoms
    
    if cellfun(@strcmp, srtdAtomElement(kk), {'Au'})
        srtdAtomElement(kk) = {'Pts'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Pt';
        % % % % %             % Margherita
        % % % % %             srtdAtomElement(kk) = {'Au'};
        % % % % %             for v=1:3
        % % % % %                 Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        % % % % %             end
        % % % % %             StockElem{kk} = 'Au';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Ag'})
        srtdAtomElement(kk) = {'Ptb'};
        Vel(kk,1:3) = 0;
        Fix = [Fix; kk];
        StockElem{kk} = 'Pt';
        
        % % % % %             % Margherita
        % % % % %             srtdAtomElement(kk) = {'Ag'};
        % % % % %             for v=1:3
        % % % % %                 Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        % % % % %             end
        % % % % %             StockElem{kk} = 'Ag';
        % % % % %
        % % % % %             elseif cellfun(@strcmp, srtdAtomElement(kk), {'Cr'})
        % % % % %                         srtdAtomElement(kk) = {'Agb'};
        % % % % %                         Vel(kk,1:3) = 0;
        % % % % %                         Fix = [Fix; kk];
        % % % % %                         StockElem{kk} = 'Ag';
        % % % % %
        % % % % %                         elseif cellfun(@strcmp, srtdAtomElement(kk), {'W'})
        % % % % %                         srtdAtomElement(kk) = {'Agb'};
        % % % % %                         Vel(kk,1:3) = 0;
        % % % % %                         Fix = [Fix; kk];
        % % % % %                         StockElem{kk} = 'Ag';
        
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Cu'})
        srtdAtomElement(kk) = {'Ptss'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Pt';
        
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Rg'})
        srtdAtomElement(kk) = {'PtE'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Pt';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'B'})
        srtdAtomElement(kk) = {'Hsurf'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1); % need to come up with Hsurf velocities, set to Pt for now, H presumably to big in liquid state
        end
        Restrain = [Restrain; kk];
        %             Fix = [Fix; kk];
        StockElem{kk} = 'H';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'H'})
        for v=1:3
            Vel(kk,v)=H.mid(find(rand(1)*sum(H.N)<cumsum(H.N)==1,1))+(H.edges(2)-H.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'H';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'S'})
        srtdAtomElement(kk) = {'Ohy'};
        for v=1:3
            Vel(kk,v)=O.mid(find(rand(1)*sum(O.N)<cumsum(O.N)==1,1))+(O.edges(2)-O.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'O'})
        for v=1:3
            Vel(kk,v)=O.mid(find(rand(1)*sum(O.N)<cumsum(O.N)==1,1))+(O.edges(2)-O.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
        %             srtdAtomElement(kk) = {'Ohy'};
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Li'})
        srtdAtomElement(kk) = {'Hhy'};
        for v=1:3
            Vel(kk,v)=H.mid(find(rand(1)*sum(H.N)<cumsum(H.N)==1,1))+(H.edges(2)-H.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'H';
        %             srtdAtomElement(kk) = {'Hhy'};
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'F'})
        srtdAtomElement(kk) = {'F'};
        for v=1:3
            Vel(kk,v)=F.mid(find(rand(1)*sum(F.N)<cumsum(F.N)==1,1))+(F.edges(2)-F.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'F';
        
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Cl'})
        srtdAtomElement(kk) = {'Cl'};
        for v=1:3
            Vel(kk,v)=F.mid(find(rand(1)*sum(F.N)<cumsum(F.N)==1,1))+(F.edges(2)-F.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Cl';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Se'})
        srtdAtomElement(kk) = {'OtU'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Te'})
        srtdAtomElement(kk) = {'OtL'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Po'})
        srtdAtomElement(kk) = {'OtS'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Sb'})
        srtdAtomElement(kk) = {'OtSS'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Pb'})
        srtdAtomElement(kk) = {'Om'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Bi'})
        srtdAtomElement(kk) = {'Ob'};
        Vel(kk,1:3) = 0;
        Fix = [Fix; kk];
        StockElem{kk} = 'O';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Mt'})
        srtdAtomElement(kk) = {'Rus'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Ru';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Ds'})
        srtdAtomElement(kk) = {'Rub'};
        Vel(kk,1:3) = 0;
        Fix = [Fix; kk];
        StockElem{kk} = 'Ru';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Hs'})
        srtdAtomElement(kk) = {'Rum'};
        Vel(kk,1:3) = 0;
        Fix = [Fix; kk];
        StockElem{kk} = 'Ru';
        % Rashid
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Ti'})
        srtdAtomElement(kk) = {'Alb'};
        Vel(kk,1:3) = 0;
        Fix = [Fix; kk];
        StockElem{kk} = 'Al';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Zr'})
        srtdAtomElement(kk) = {'Al2'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Al';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Hf'})
        srtdAtomElement(kk) = {'Al1'};
        for v=1:3
            Vel(kk,v)=Pt.mid(find(rand(1)*sum(Pt.N)<cumsum(Pt.N)==1,1))+(Pt.edges(2)-Pt.edges(1))*(2*rand(1)-1);
        end
        StockElem{kk} = 'Al';
    elseif cellfun(@strcmp, srtdAtomElement(kk), {'Mo'})
        srtdAtomElement(kk) = {'Mo'};
        Vel(kk,1:3) = 0;
        Fix = [Fix; kk];
        StockElem{kk} = 'Mo';
        
    end
    CartPos = AscAtomPosition(kk,:)*Vec;
    fprintf(fidout,[pad(srtdAtomElement{kk},8) pad(num2str(CartPos(1), '%.10g'),20) pad(num2str(CartPos(2), '%.10g'),20) pad(num2str(CartPos(3), '%.10g'),20) newlinechar]);
end
fclose(fidout);

ListOfElems = {StockElem(cumsum(n)) srtdAtomElement(cumsum(n))'};