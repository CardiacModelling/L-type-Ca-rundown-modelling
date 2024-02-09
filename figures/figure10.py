
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.optimize import fsolve

import helpers
import parameters

pca = parameters.pca_default # nL/s default
dR = parameters.dR # cm

# Set font
matplotlib.rc('font', family='arial', size = 8)

fig = plt.figure(figsize=(6.6, 2.6))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

ax1.set_ylabel('Rundown')
ax1.set_xlabel('Time (s)')
ax2.set_xlabel('Time (s)')
ax3.set_xlabel('Time (s)')


n_draw_ax1 = 10 # number of draws for each
n_draw_ax2 = 5

# plots of no change with zeta less
count_t = 0
while count_t < n_draw_ax2:
    found_p = False
    while found_p == False:
        thold, t0, tdiff = helpers.draw_t() # milli seconds
        t_const = [thold, t0, tdiff]
        if min(t_const) == tdiff:
            a = t0/thold
            b = t0/tdiff
            if a*b > 1.5:
                found_p = True
    
    count_t += 1
    n_sweep, mult = parameters.cal_n_sweep(thold)
    peak_ca_a = helpers.peak_ical(thold, t0, tdiff, n_sweep, pca, dR)
    x_range = np.linspace(0, n_sweep -1, n_sweep)*mult

    ax1.plot(x_range, peak_ca_a, color  = '#1f77b4')

# plots of no change with zeta greater
count_t = 0
while count_t < n_draw_ax2:
    found_p = False
    while found_p == False:
        thold, t0, tdiff = helpers.draw_t() # milli seconds
        t_const = [thold, t0, tdiff]
        if min(t_const) == tdiff:
            a = t0/thold
            b = t0/tdiff
            if a*b < 0.5:
                found_p = True
    
    count_t += 1
    n_sweep, mult = parameters.cal_n_sweep(thold)
    peak_ca_a = helpers.peak_ical(thold, t0, tdiff, n_sweep, pca, dR)
    x_range = np.linspace(0, n_sweep -1, n_sweep)*mult

    # TODO: Remove
    if peak_ca_a[1] > 0:
        count_t = count_t -1
        continue
    
    ax1.plot(x_range, peak_ca_a, color  = '#ff7f0e')

# plots of inverse
count_t = 0
while count_t < n_draw_ax1:
    found_p = False
    while found_p == False:
        thold, t0, tdiff = helpers.draw_t() # milli seconds
        t_const = [thold, t0, tdiff]
        if min(t_const) == t0:
            found_p = True

    count_t += 1
    n_sweep, mult = parameters.cal_n_sweep(thold)
    peak_ca_a = helpers.peak_ical(thold, t0, tdiff, n_sweep, pca, dR)
    x_range = np.linspace(0, n_sweep -1, n_sweep)*mult

    ax2.plot(x_range, peak_ca_a, color  = '#ff7f0e')


# plots of log with zeta > 1
count_t = 0
while count_t < n_draw_ax2:
    found_p = False
    while found_p == False:
        thold, t0, tdiff = helpers.draw_t() # milli seconds
        t_const = [thold, t0, tdiff]
        if min(t_const) == thold:
            a = t0/thold
            b = t0/tdiff
            if a*b > 1.5:
                found_p = True

    count_t += 1
    n_sweep, mult = parameters.cal_n_sweep(thold)
    peak_ca_a = helpers.peak_ical(thold, t0, tdiff, n_sweep, pca, dR)
    x_range = np.linspace(0, n_sweep -1, n_sweep)*mult
        
    ax3.plot(x_range, peak_ca_a, color  = '#1f77b4')


# plots of log with zeta < 1
count_t = 0
while count_t < n_draw_ax2:
    found_p = False
    while found_p == False:
        thold, t0, tdiff = helpers.draw_t() # milli seconds
        t_const = [thold, t0, tdiff]
        if min(t_const) == thold:
            a = t0/thold
            b = t0/tdiff
            if a*b < 0.5:
                found_p = True

    count_t += 1
    n_sweep, mult = parameters.cal_n_sweep(thold)
    peak_ca_a = helpers.peak_ical(thold, t0, tdiff, n_sweep, pca, dR)
    x_range = np.linspace(0, n_sweep -1, n_sweep)*mult

    ax3.plot(x_range, peak_ca_a, color  = '#ff7f0e')


ax3.plot([], [], color  = '#ff7f0e', label = '$\zeta$ < 1.5')
ax3.plot([], [], color  = '#1f77b4', label = '$\zeta$ > 4.5')
ax3.legend(bbox_to_anchor = (-1.12, 1.05, 0.6, 0.15), loc = 'center', ncol =2,\
    prop={'size': 8})

# ax1.text(0.05, 0.95, f"N = {n_draw_ax1}", transform=ax1.transAxes, verticalalignment='top', horizontalalignment='left')

ax1.set_xticks([0, 50, 100, 150], [0, 50, 100, 150])
ax2.set_xticks([0, 50, 100, 150], [0, 50, 100, 150])
ax3.set_xticks([0, 50, 100, 150], [0, 50, 100, 150])

ax1.set_xlim(0, 160)
ax2.set_xlim(0, 160)
ax3.set_xlim(0, 160)

plt.subplots_adjust(bottom=0.17, wspace = 0.35, top = 0.85)
plt.savefig('figures/figure10.pdf')
plt.close()
