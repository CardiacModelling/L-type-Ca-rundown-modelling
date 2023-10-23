import matplotlib
import matplotlib.pyplot as plt
import myokit
import copy
import numpy as np
from matplotlib import cm

# import project files
from helpers import multi_compartment_model
import parameters

# Set font
matplotlib.rc('font', family='arial', size = 8)

# Constants
R = parameters.R_default # cm
dR = parameters.dR # cm
N_div = int((R-parameters.Rh)/dR) #70
t0 = parameters.t0_default # milli seconds
thold = parameters.thold_default # milli seconds
n_sweep, mult = parameters.cal_n_sweep(thold) 

# array for store
log_b = [f'calcium_dynamics.B{i}' for i in range(N_div)]

model = myokit.load_model('resources/zeng-ical-template.mmt')
model.get('calcium_dynamics.R').set_rhs(R)

m = multi_compartment_model(model, N_div, [0 for _ in range(N_div)])

# define B initial
m1 = copy.deepcopy(m)
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

s = myokit.Simulation(m1)
s.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)

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

# can set initial state here

colors = [cm.viridis(x) for x in np.linspace(0, 1, n_sweep)]

fig = plt.figure(figsize=(6.6, 5.5))
run = fig.add_subplot(224)
ical = fig.add_subplot(223)
cal = fig.add_subplot(221)
fca = fig.add_subplot(222)

ax = run.twinx()


times = np.arange(0, 120, 0.1)
t_plot = times/1000

for i in range(n_sweep):

    sim.set_constant('membrane.V', -90)
    sim.run(thold, log = [])
    sim.set_constant('membrane.V', -0.001) 
    time = sim.time()
    d = sim.run(120, log = ['L_type_Ca_current.i_CaL', 'calcium_dynamics.Ca0', 'L_type_Ca_current.f_Ca'], \
        log_times=np.arange(time, time+ 120, 0.1))

    fca.plot(t_plot, d['L_type_Ca_current.f_Ca'], color = colors[i])
    ical.plot(t_plot, d['L_type_Ca_current.i_CaL'], color = colors[i])
    cal.plot(t_plot, d['calcium_dynamics.Ca0'], color = colors[i])
    imin = min(d['L_type_Ca_current.i_CaL'])
    if i == 0:
        i_sweep0 = imin 
    run.scatter((i+1)*mult, imin, color = colors[i])
    ax.scatter((i+1)*mult, 1 - imin/i_sweep0, color = colors[i])

ax1_ylow, ax1_yhigh = ical.get_ylim()

    
ical.set_xlabel('Time (s)')
run.set_xlabel('Time (s)')
ical.set_ylabel('Current (pA)')
run.set_ylabel('Peak current (pA)')
fca.set_ylabel('[$\mathrm{Ca}^{2+}$]-dependent gate ($f_{Ca}$)')
cal.set_ylabel('[$\mathrm{Ca}^{2+}$] in the outermost shell (mM)')
fca.set_ylim(0,1)
run.set_ylim(ax1_ylow, ax1_yhigh)
ax.set_ylim(1 - ax1_ylow/i_sweep0, 1 - ax1_yhigh/i_sweep0)
ax.set_ylabel('Rundown')

plt.tight_layout()
plt.savefig('figures/figure6.pdf')
plt.close()














