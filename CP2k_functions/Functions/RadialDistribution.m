function [RadialMinima] = RadialDistribution(RadFunSnaps, ABC, ElemNames, plotTF)

disp(['Computing radial distribution function for ' ElemNames(1,:) ' and ' ElemNames(2,:)])

% concatenate all bond distances for all scanned atoms from all snapshots
% and sort
dR = sort(vertcat(RadFunSnaps{:,:}));

% set bin-width
L=0.04;
r = 0:L:30;
rho = length(dR)/(ABC(1)*ABC(2)*(ABC(3)-9)); % assumes volume available is calculated ignoring 9 Angstroms of z where metal is...
for b=1:length(r)-1
    count(b) = sum(dR>r(b) & dR<r(b+1));
    V(b) = (4/3)*pi*(r(b+1)^3) - (4/3)*pi*(r(b)^3);
    g(b) = (count(b)/V(b))*(1/rho);
end

if any(ismember(ElemNames, 'O')) & any(ismember(ElemNames, 'H'))
    C = [0 0 0];
    yScale = [0 2.5];
elseif any(ismember(ElemNames, 'F')) & any(ismember(ElemNames, 'H'))
    C = [1 0 1];
    yScale = [0 3];
elseif sum(ismember(ElemNames, 'O')) == 2
    C = [0 0 1];
    yScale = [0 3.5];
elseif sum(ismember(ElemNames, 'H')) == 2
    C = [0 0.5 0];
    yScale = [0 2.5];
elseif sum(ismember(ElemNames, 'F')) == 2
    C = [1 0 0];
    yScale = [0 2];
else
    C = [218/255 165/255 32/255];
    yScale = [0 2];
end

Minima = LocateStationaryPoints(smooth(g));
% Minima = LocateStationaryPoints(g);
RadialMinima = r(Minima{1}+1);

if plotTF == 1
    figure
    set(gcf, 'position', [600 600 459   288]);
    hold on
    box on
    plot(r(2:end), g, 'color', C, 'linewidth', 2);
%     plot(r(2:end),smooth(g), 'color', C, 'linewidth', 2);
    set(gca, 'xlim',[0 7], 'ylim', yScale, 'ytick', [0 0.5 1 1.5 2 2.5 3 3.5 4], 'yticklabel', {'0.0' '0.5' '1.0' '1.5' '2.0' '2.5' '3.0' '3.5' '4.0'},'fontsize', 14)
    legend(['g_{' ElemNames(1,:) '}_{' ElemNames(2,:) '}(r)'])
    xlabel('r (Ang)');
    ylabel(['g_{' ElemNames(1,:) '}_{' ElemNames(2,:) '}(r)']);
    hold off
end

return
