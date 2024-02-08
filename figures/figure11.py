import matplotlib
import matplotlib.pyplot as plt
import myokit
import copy
import numpy as np
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# import project files
from helpers import multi_compartment_model
import parameters

# Set font
matplotlib.rc('font', family='arial', size = 8)

# Constants
bmax_range = [0.1] #np.linspace(0.1, 20, 10)
bmax_range = bmax_range + [2*(i+1) for i in range(10)]
R = parameters.R_default # cm
dR = parameters.dR # cm
N_div = int((R-parameters.Rh)/dR) #70
t0 = parameters.t0_max # milli seconds
thold = 40000 #parameters.thold_default # milli seconds
n_sweep, mult = parameters.cal_n_sweep(thold) 
colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, len(bmax_range))]

fig = plt.figure(figsize=(6.6, 3))
#run_1 = fig.add_subplot(223)
#ical_1 = fig.add_subplot(221)
run_2 = fig.add_subplot(122)
ical_2 = fig.add_subplot(121)

# array for store
log_b = [f'calcium_dynamics.B{i}' for i in range(N_div)]

model = myokit.load_model('resources/zeng-ical-template.mmt')
model.get('calcium_dynamics.R').set_rhs(R)

m0 = multi_compartment_model(model, N_div, [0 for _ in range(N_div)])

# define B initial
m1 = copy.deepcopy(m0)
m1.get('calcium_dynamics.buffering').set_rhs(0)
m1.get('calcium_dynamics.current').set_rhs(0)

# strip ion-channel model
m1.get('L_type_Ca_current.d').demote
m1.get('L_type_Ca_current.d').set_rhs(0)
m1.get('L_type_Ca_current.f').demote
m1.get('L_type_Ca_current.f').set_rhs(1)

for i in range(N_div):
    m1.get(f'calcium_dynamics.Ca{i}').demote
    m1.get(f'calcium_dynamics.Ca{i}').set_rhs(0)
    m1.get(f'calcium_dynamics.CaB{i}').demote
    m1.get(f'calcium_dynamics.CaB{i}').set_rhs(0)

for i in range(N_div):
    m1.get(f'calcium_dynamics.B{i}').set_state_value(0)

s0 = myokit.Simulation(m1)
s0.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)

s = copy.deepcopy(s0)
m = copy.deepcopy(m0)
d = s.run(t0 + 1, log = log_b, log_times = np.array([t0]))

for b in log_b:
    m.get(b).set_state_value(d[b][0])

# ion-channel model ss
d_inf = m.get('L_type_Ca_current.d_infinity').pyfunc()(-90)
m.get('L_type_Ca_current.d').set_state_value(d_inf)
f_inf = m.get('L_type_Ca_current.f_infinity').pyfunc()(-90)
m.get('L_type_Ca_current.f').set_state_value(f_inf)
volt = m.get('membrane.V')
volt.set_binding(None)

sim = myokit.Simulation(m)
sim.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)


times = np.arange(0, 120, 0.1)
t_plot = times/1000

sim_store = copy.deepcopy(sim)

for j in range(len(bmax_range)):

    if bmax_range[j] == 10:
        color = 'red'
    else:
        color = colors[j]
        
    sim = copy.deepcopy(sim_store)

    sim.set_constant('calcium_dynamics.B_max', bmax_range[j])

    x = []
    y = []
    for i in range(n_sweep):

        sim.set_constant('membrane.V', -90)
        sim.run(thold, log = [])
        sim.set_constant('membrane.V', -0.001) 
        time = sim.time()
        d = sim.run(120, log = ['L_type_Ca_current.i_CaL', 'calcium_dynamics.Ca0', \
                                'L_type_Ca_current.f_Ca'], log_times=np.arange(time, time+ 120, 0.1))
        
        if i == 0:
            ical_2.plot(t_plot, d['L_type_Ca_current.i_CaL'], color = color)

        imin = min(d['L_type_Ca_current.i_CaL'])
        if i == 0:
            i_sweep0 = imin 
        x.append((i+1)*mult)
        y.append(1 - imin/i_sweep0)
    run_2.plot(x, y, color = color)

# Add colorbar between the subplots on the top
cax = fig.add_axes([0.90, 0.18, 0.02, 0.7])  # Adjust position and size as needed
norm = matplotlib.colors.Normalize(vmin=np.min(bmax_range), vmax=np.max(bmax_range))
sm = plt.cm.ScalarMappable(cmap=matplotlib.cm.viridis, norm=norm)
sm.set_array([])  # You need to set an empty array to the ScalarMappable

# Create colorbar
cbar = plt.colorbar(sm, cax=cax, orientation='vertical')

# Manually set the colorbar label at the bottom
label_text = '${B}_{max}$ (mM)'
cbar.set_label('')
cbar.ax.text(1.8, 1.11, label_text, va='top', ha='center', \
             transform=cbar.ax.transAxes)

    
# run_1.set_xlabel('Time (s)')
run_2.set_xlabel('Time (s)')

ical_2.set_ylabel('Current (pA)')
run_2.set_ylabel('Rundown')

plt.subplots_adjust(top=0.95, right=0.85, bottom= 0.15, wspace=0.30)
plt.savefig('figures/figure11.pdf')
plt.close()














