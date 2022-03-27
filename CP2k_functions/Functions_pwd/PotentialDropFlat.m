scatter([-1.59 -1.27 -0.89 -0.6], [-1.01 -0.7 1.5 3.1], 'markerfacecolor', 'r', 'markeredgecolor', 'k')
plot([-1.3 pzc], [0 0], '--', 'color', [0.8 0.8 0.8])
plot([pzc pzc], [-1 0], '--', 'color', [0.8 0.8 0.8])
fitlm([-1.27 -0.89 -0.6], [-0.7 1.5 3.1])
pzc = -6.5234/5.6776;