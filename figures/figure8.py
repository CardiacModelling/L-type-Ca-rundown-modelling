import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import helpers
import parameters

pca_range = np.linspace(parameters.pca_min, parameters.pca_max, 10)
thold = parameters.thold_default
n_sweep, mult = parameters.cal_n_sweep(thold)
R = parameters.R_default # cm
dR = parameters.dR # cm
t0_arr = [parameters.t0_min, parameters.t0_max] # milli seconds

N_div = int((R-parameters.Rh)/dR)

# Set font
matplotlib.rc('font', family='arial', size = 8)

fig = plt.figure(figsize=(6.6, 3))
fig_cont = []
for i in range(2):
    fig_cont.append(fig.add_subplot(int(f'12{i+1}')))
fig_cont[0].set_ylabel('Rundown')
fig_cont[0].set_xlabel('Time (s)')
fig_cont[1].set_xlabel('Time (s)')
fig_cont[0].set_xticks([0, 50, 100, 150, 200], [0, 50, 100, 150, 200])
fig_cont[1].set_xticks([0, 50, 100, 150, 200], [0, 50, 100, 150, 200])

colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, len(pca_range))]

for i in range(len(pca_range)):
    #d = helpers.simulate_ical(n_sweep, N_div, R, t0, pca_range[i])
    x_range = np.linspace(0, n_sweep-1, n_sweep)*mult

    d, _ = helpers.simulate_ical_curr(n_sweep, N_div, R, t0_arr[0], pca_range[i])
    fig_cont[0].plot(x_range, d, color = colors[i])

    d, _ = helpers.simulate_ical_curr(n_sweep, N_div, R, t0_arr[1], pca_range[i])
    fig_cont[1].plot(x_range, d, color = colors[i])
    
# Add colorbar between the subplots on the top
cax = fig.add_axes([0.85, 0.18, 0.02, 0.7])  # Adjust position and size as needed
norm = matplotlib.colors.Normalize(vmin=np.min(pca_range), vmax=np.max(pca_range))
sm = plt.cm.ScalarMappable(cmap=matplotlib.cm.viridis, norm=norm)
sm.set_array([])  # You need to set an empty array to the ScalarMappable

# tlabel
fig_cont[0].set_title("Minimum $t_0$")
fig_cont[1].set_title("Maximum $t_0$")

# Create colorbar
cbar = plt.colorbar(sm, cax=cax, orientation='vertical')

# Manually set the colorbar label at the bottom
label_text = '$\overline{P}_{Ca}$ (nL/s)'
cbar.set_label('')
cbar.ax.text(1.8, 1.1, label_text, va='top', ha='center', \
             transform=cbar.ax.transAxes)


plt.subplots_adjust(top=0.92, right=0.82, bottom= 0.15, wspace=0.30)
plt.savefig('figures/figure8.pdf')
plt.close()