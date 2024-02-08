import matplotlib.pyplot as plt
import matplotlib
import string
import numpy as np

import helpers
import parameters

pca_range = [parameters.pca_min, parameters.pca_max] # nL/s
r_range = [parameters.R_min +3e-4, parameters.R_max] # cm, increased the min slightly to avoid convergence issues
t0_range = [parameters.t0_min, parameters.t0_max] # ms
thold = parameters.thold_default
n_sweep, mult = parameters.cal_n_sweep(thold)
dR = parameters.dR # cm
Rh = parameters.Rh

# Set font
matplotlib.rc('font', family='arial', size = 8)

fig, axs = plt.subplots(
    nrows=7, # add 3 new rows
    ncols=5,  # add 1 new column ...
    figsize=(6.6, 8),                      
    gridspec_kw={'width_ratios': [1, 1, 0.01, 1, 1], 
                 'height_ratios': [1, 0.01, 1, 0.01, 1, 0.01, 1]}
)

for ax in axs[:, 2]: 
    ax.axis("off")  # Turn off axis lines and labels for this column

for ax in axs[1, :]:
    ax.axis("off")

for ax in axs[3, :]:
    ax.axis("off")

for ax in axs[5, :]:
    ax.axis("off")

fig_cont = []

for i in range(7):
    if i == 1 or i ==3 or i==5:
        continue

    for j in range(5):
        if j == 2:
            continue
        fig_cont.append(fig.add_subplot(axs[i, j]))

for i in range(16):
    if i % 2 == 0:
        pass
    else:
        if i == 9 or i ==11:
            fig_cont[i].set_ylim(0, 10)
        else:
            fig_cont[i].set_ylim(9.7, 10)

fig.text(0.17, 0.97, 'Minimum permeability ($\overline{P}_{Ca}$)', \
         fontsize = 9, family = 'arial')
fig.text(0.68, 0.97, 'Maximum permeability ($\overline{P}_{Ca}$)', \
         fontsize = 9, family = 'arial')


fig_cont[0].set_title('[$Ca^{2+}]_s$ (mM)')
fig_cont[1].set_title('[B] (mM)')
fig_cont[2].set_title('[$Ca^{2+}]_s$ (mM)')
fig_cont[3].set_title('[B] (mM)')
fig_cont[0].set_ylabel('Min $R_0$, Min $t_0$')
fig_cont[4].set_ylabel('Min $R_0$, Max $t_0$')
fig_cont[8].set_ylabel('Max $R_0$, Min $t_0$')
fig_cont[12].set_ylabel('Max $R_0$, Max $t_0$')
fig_cont[12].set_xlabel('Time (s)')
fig_cont[13].set_xlabel('Time (s)')
fig_cont[14].set_xlabel('Time (s)')
fig_cont[15].set_xlabel('Time (s)')

for f in fig_cont:
    f.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0, 2))

gs0 = matplotlib.gridspec.GridSpec(1, 1)

gs0.update(left = 0, right = 0, wspace = 5)

props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
for i in range(8):

    fig_cont[2*i].text(0.05, 0.9, string.ascii_uppercase[i] + '.1', transform = fig_cont[2*i].transAxes, weight = 'bold', bbox=props)
    fig_cont[2*i+1].text(0.05, 0.9, string.ascii_uppercase[i] + '.2', transform = fig_cont[2*i +1].transAxes, weight = 'bold', bbox=props)

i = 0
for r in r_range:
    N_div = int((r-Rh)/dR)
    print('Radius:(microm)', r*10000, 'Number of div:', N_div)

    for t0 in t0_range:
        for pca in pca_range:
            d = helpers.simulate_ical(n_sweep, N_div, r, t0, pca)
            len_x = len(d['calcium_dynamics.Ca0'])
            x_range = np.linspace(0, n_sweep, len_x)
            x_range = mult*x_range
            fig_cont[i].plot(x_range, d['calcium_dynamics.Ca0'])
            fig_cont[i+1].plot(x_range, d['calcium_dynamics.B0'])
            i+= 2

for i in range(n_sweep):
    for j in range(16):
        fig_cont[j].axvspan(i*mult, (i+1)*mult, color = 'grey', \
            alpha = (i%2) * 0.3 + 0.1, transform=fig_cont[j].get_xaxis_transform())      

plt.tight_layout()
plt.subplots_adjust(left = 0.08, bottom = 0.05, right = 0.98, top= 0.925)
plt.savefig('figures/figure7.pdf')
plt.close()