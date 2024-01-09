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
t0 = parameters.t0_max # milli seconds

N_div = int((R-parameters.Rh)/dR)

# Set font
matplotlib.rc('font', family='arial', size = 8)

fig = plt.figure(figsize=(6.6, 3))
fig_cont = []
for i in range(2):
    fig_cont.append(fig.add_subplot(int(f'12{i+1}')))
fig_cont[0].set_ylabel('[[$Ca^{2+}]_s$ (mM)')
fig_cont[1].set_ylabel('[B] (mM)')
fig_cont[1].set_xlabel('Time (s)')
fig_cont[0].set_xlabel('Time (s)')
fig_cont[0].set_xticks([0, 50, 100, 150], [0, 50, 100, 150])
fig_cont[1].set_xticks([0, 50, 100, 150], [0, 50, 100, 150])


colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, len(pca_range))]

for i in range(len(pca_range)):
    d = helpers.simulate_ical(n_sweep, N_div, R, t0, pca_range[i], block_fca=False)
    len_x = len(d['calcium_dynamics.Ca0'])
    x_range = np.linspace(0, n_sweep, len_x)*mult
    fig_cont[0].plot(x_range, d['calcium_dynamics.Ca0'], color = colors[i])
    fig_cont[1].plot(x_range, d['calcium_dynamics.B0'], color = colors[i])
    

for i in range(n_sweep):
    fig_cont[0].axvspan(i*mult, (i+1)*mult, color = 'grey', \
    alpha = (i%2) * 0.3 + 0.1, transform=fig_cont[0].get_xaxis_transform())
    fig_cont[1].axvspan(i*mult, (i+1)*mult, color = 'grey', \
    alpha = (i%2) * 0.3 + 0.1, transform=fig_cont[1].get_xaxis_transform())

# Add colorbar between the subplots on the top
cax = fig.add_axes([0.15, 0.92, 0.7, 0.02])  # Adjust position and size as needed
norm = matplotlib.colors.Normalize(vmin=np.min(pca_range), vmax=np.max(pca_range))
sm = plt.cm.ScalarMappable(cmap=matplotlib.cm.viridis, norm=norm)
sm.set_array([])  # You need to set an empty array to the ScalarMappable

# Create colorbar
cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
cbar.set_label('$\overline{P}_{Ca}$ (nL/s)', labelpad = -30)

plt.subplots_adjust(top=0.82, right=0.88)
plt.savefig('supporting/figure8_fca.png')
plt.close()