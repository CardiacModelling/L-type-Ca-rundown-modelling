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

fig = plt.figure(figsize=(6.6, 8))
gs = matplotlib.gridspec.GridSpec(4, 4, figure = fig)
fig_cont = []
for i in range(16):
    fig_cont.append(fig.add_subplot(gs[i]))
    if i % 2 == 0:
        pass

fig_cont[0].set_title('$\overline{P}_{Ca}$ Min\n[$\mathrm{Ca}^{2+}$] (mM)')
fig_cont[1].set_title('$\overline{P}_{Ca}$ Min\n[B] (mM)')
fig_cont[2].set_title('$\overline{P}_{Ca}$ Max\n[$\mathrm{Ca}^{2+}$] (mM)')
fig_cont[3].set_title('$\overline{P}_{Ca}$ Max\n[B] (mM)')
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
plt.subplots_adjust(left = 0.08, bottom = 0.05, right = 0.98)
plt.savefig('figures/figure7.pdf')
plt.close()