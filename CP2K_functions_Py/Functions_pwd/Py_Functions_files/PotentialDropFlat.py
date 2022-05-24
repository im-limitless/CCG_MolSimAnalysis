import numpy as np
import matplotlib.pyplot as plt
scatter(np.array([- 1.59,- 1.27,- 0.89,- 0.6]),np.array([- 1.01,- 0.7,1.5,3.1]),'markerfacecolor','r','markeredgecolor','k')
plt.plot(np.array([- 1.3,pzc]),np.array([0,0]),'--','color',np.array([0.8,0.8,0.8]))
plt.plot(np.array([pzc,pzc]),np.array([- 1,0]),'--','color',np.array([0.8,0.8,0.8]))
fitlm(np.array([- 1.27,- 0.89,- 0.6]),np.array([- 0.7,1.5,3.1]))
pzc = - 6.5234 / 5.6776