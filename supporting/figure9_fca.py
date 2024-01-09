import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import helpers
import parameters

pca = parameters.pca_default
n_sweep = 20
dr = parameters.dR # cm
R = parameters.R_default # cm
N_div = int((R-parameters.Rh)/dr)
t0_range = [10, 20, 40, 60, 80, 100, 150, 200, 250, 300] # seconds 
t0_range = [x*1000 for x in t0_range] # milli seconds

# Set font
matplotlib.rc('font', family='arial', size = 8)

fig = plt.figure(figsize=(6.6, 3))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)


ax1.set_ylabel('Peak-[$Ca^{2+}]_s$ (normalised)')
ax1.set_xlabel('Time (s)')
ax2.set_xlabel('Time (s)')
ax3.set_xlabel('Time (s)')


colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, len(t0_range))]

for i in range(len(t0_range)):
    print(t0_range[i]/1000)

    print('t1')
    thold = 10000
    n_sweep, mult = parameters.cal_n_sweep(thold)
    d1, _ = helpers.simulate_ical(n_sweep, N_div, R, t0_range[i], pca, thold = thold, peak = True, block_fca=False)
    x_range = np.linspace(0, n_sweep -1, n_sweep)
    ax1.plot(x_range*mult, d1, color = colors[i], label = f'{int(t0_range[i]/1000)}s')

    print('t2')
    thold = 20000
    n_sweep, mult = parameters.cal_n_sweep(thold)
    d2, _ = helpers.simulate_ical(n_sweep, N_div, R, t0_range[i], pca, thold = thold, peak = True, block_fca=False)
    x_range = np.linspace(0, n_sweep -1, n_sweep)
    ax2.plot(x_range*mult, d2, color = colors[i], label = f'{int(t0_range[i]/1000)}s')

    print('t3')
    thold = 40000
    n_sweep, mult = parameters.cal_n_sweep(thold)
    d3, _ = helpers.simulate_ical(n_sweep, N_div, R, t0_range[i], pca, thold = thold, peak = True, block_fca=False)
    x_range = np.linspace(0, n_sweep -1, n_sweep)
    ax3.plot(x_range*mult, d3, color = colors[i], label = f'{int(t0_range[i]/1000)}s')

ax1ymin, ax1ymax = ax1.get_ylim()
ax2ymin, ax2ymax = ax2.get_ylim()
ax3ymin, ax3ymax = ax3.get_ylim()

ymin = min([ax1ymin, ax2ymin, ax3ymin])
ymax = max([ax1ymax, ax2ymax, ax3ymax])

ax1.set_ylim(ymin, ymax)
ax2.set_ylim(ymin, ymax)
ax3.set_ylim(ymin, ymax)

ax1.set_xticks([0, 50, 100, 150], [0, 50, 100, 150])
ax2.set_xticks([0, 50, 100, 150], [0, 50, 100, 150])
ax3.set_xticks([0, 50, 100, 150], [0, 50, 100, 150])

ax2.legend(bbox_to_anchor = (1.6, 0.95, 0.8, 0.2) , ncol = 10, \
    prop={'size': 6}, bbox_transform = ax2.transAxes)

ax1.text(0.05, 0.97, '$t_{hold} = 10s$', transform=ax1.transAxes,
                 verticalalignment='top', horizontalalignment='left', fontsize=8)
ax2.text(0.05, 0.97, '$t_{hold} = 20s$', transform=ax2.transAxes,
                 verticalalignment='top', horizontalalignment='left', fontsize=8)
ax3.text(0.05, 0.97, '$t_{hold} = 40s$', transform=ax3.transAxes,
                 verticalalignment='top', horizontalalignment='left', fontsize=8)

plt.subplots_adjust(bottom=0.17, wspace = 0.35)
plt.savefig('supporting/figure9_fca.png')
plt.close()
