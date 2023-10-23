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
m_fca = copy.deepcopy(m)
d_inf = m_fca.get('L_type_Ca_current.d_infinity').pyfunc()(-90)
m_fca.get('L_type_Ca_current.d').set_state_value(d_inf)
f_inf = m_fca.get('L_type_Ca_current.f_infinity').pyfunc()(-90)
m_fca.get('L_type_Ca_current.f').set_state_value(f_inf)
volt = m_fca.get('membrane.V')
volt.set_binding(None)

d_inf = m.get('L_type_Ca_current.d_infinity').pyfunc()(-90)
m.get('L_type_Ca_current.d').set_state_value(d_inf)
f_inf = m.get('L_type_Ca_current.f_infinity').pyfunc()(-90)
m.get('L_type_Ca_current.f').set_state_value(f_inf)
volt = m.get('membrane.V')
volt.set_binding(None)
# block fca
m.get('L_type_Ca_current.open_prob').set_rhs('d*f')

sim = myokit.Simulation(m)
sim.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)

sim_fca = myokit.Simulation(m_fca)
sim_fca.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)

# can set initial state here

colors = [cm.viridis(x) for x in np.linspace(0, 1, n_sweep)]

fig = plt.figure(figsize=(6.6, 5.5))

ax_cal = fig.add_subplot(221)
ax_run = fig.add_subplot(223)
ax_cal_fca = fig.add_subplot(222)
ax_run_fca = fig.add_subplot(224)


times = np.arange(0, 120, 0.1)
t_plot = times/1000

for i in range(n_sweep):

    sim.set_constant('membrane.V', -90)
    sim.run(thold, log = [])
    sim.set_constant('membrane.V', -0.001) 
    time = sim.time()
    d = sim.run(120, log = ['calcium_dynamics.Ca0'], \
        log_times=np.arange(time, time+ 120, 0.1))

    ax_cal.plot(t_plot, d['calcium_dynamics.Ca0'], color = colors[i])
    camax = max(d['calcium_dynamics.Ca0'])
    # if i == 0:
    #     i_sweep0 = imin 
    ax_run.scatter((i+1)*mult, camax, color = colors[i])
    #ax.scatter((i+1)*mult, 1 - imin/i_sweep0, color = colors[i])

    sim_fca.set_constant('membrane.V', -90)
    sim_fca.run(thold, log = [])
    sim_fca.set_constant('membrane.V', -0.001) 
    time = sim_fca.time()
    d_fca = sim_fca.run(120, log = ['calcium_dynamics.Ca0'], \
        log_times=np.arange(time, time+ 120, 0.1))

    ax_cal_fca.plot(t_plot, d_fca['calcium_dynamics.Ca0'], color = colors[i])
    camax = max(d_fca['calcium_dynamics.Ca0'])
    # if i == 0:
    #     i_sweep0 = imin 
    ax_run_fca.scatter((i+1)*mult, camax, color = colors[i])
    #ax.scatter((i+1)*mult, 1 - imin/i_sweep0, color = colors[i])


ax1_ylow, ax1_yhigh = ax_cal.get_ylim()

    
ax_run_fca.set_xlabel('Time (s)')
ax_run.set_xlabel('Time (s)')
ax_cal.set_ylabel('[$\mathrm{Ca}^{2+}$] in the outermost shell (mM)')
ax_run.set_ylabel('peak-[$\mathrm{Ca}^{2+}$] in the outermost shell (mM)')
ax_cal_fca.set_title('$f_{Ca}$ not blocked')
ax_cal.set_title('$f_{Ca}$ blocked')

plt.tight_layout()
plt.savefig('supporting/figure6_no_fca.png')
plt.close()














