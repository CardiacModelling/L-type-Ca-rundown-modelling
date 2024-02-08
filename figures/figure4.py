import matplotlib.pyplot as plt
import matplotlib
import numpy as np


fig = plt.figure(figsize=(6.6, 1))
ax = fig.add_subplot(111)

# Set font
matplotlib.rc('font', family='arial', size = 8)

x = np.linspace(1, 422)
ax.plot(x, np.ones(shape=np.shape(x)), color = 'black')
ax.vlines(1, 0.9, 1.1, color = 'black')
ax.vlines(422, 0.9, 1.1, color = 'black')

y = np.linspace(5, 300)
ax.plot(y, np.ones(shape=np.shape(y)), color = 'red')
ax.vlines(5, 0.9, 1.1, color = 'red')
ax.vlines(300, 0.9, 1.1, color = 'red')

z = np.linspace(10, 40)
ax.plot(z, np.ones(shape=np.shape(z)), color = 'blue', ls='--')
ax.vlines(10, 0.9, 1.1, color = 'blue')
ax.vlines(20, 0.9, 1.1, color = 'blue')
ax.vlines(40, 0.9, 1.1, color = 'blue')

ax.text(0.9, 0.6, '0 s')
ax.text(240, 0.6, '282 s')
ax.text(2.1, 1.2, '$t_{diff}$', color = 'black')

ax.text(8, 0.6, '10 s')
ax.text(35, 0.6, '40 s')
ax.text(18, 1.2, '$t_{hold}$', color = 'blue')

ax.text(4, 0.6, '5 s')
ax.text(370, 0.6, '300 s')
ax.text(100, 1.2, '$t_{0}$', color = 'red')

ax.set_ylim(0, 2)
ax.set_yticks([])
ax.axis('off')
ax.set_xscale('log')
plt.savefig('figures/figure4.pdf')
plt.close()