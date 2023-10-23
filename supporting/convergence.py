import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import optimize
import myokit

import sys
sys.path.append('supporting')
import analytical_solution
import parameters
import helpers


DB = parameters.DB #cm2/ms
R = parameters.R_max #cm
rh = parameters.Rh 
B_max = parameters.Bmax
t_mult = DB/pow(R-rh, 2)
t = parameters.t0_max # ms time at which we're investigating the concentration profile
tau = t*t_mult 
x_range = np.linspace(0, 1, 200)
h = 1 - (rh/R)
roots = analytical_solution.check_roots(h)


matplotlib.rc('font', family='arial', size = 8)
fig = plt.figure(figsize=(6.6, 3))
ax = fig.add_subplot(121) 
ax2 = fig.add_subplot(122) 

ax.set_ylim(0, B_max)
ax.set_xlim(0, R*10000)
ax.set_ylabel('Free buffer (mM)')
ax.set_xlabel('Distance from the centre ($\mu$m)')


arr = []
rad = []

for x in x_range:
    term_a = analytical_solution.term1(x, h)
    term_b = analytical_solution.convergence(x, tau, roots, h)
    u = term_a - term_b
    r = R + x*(rh - R)
    arr.append(B_max*u*rh/r)
    rad.append(r*10000)

ax.plot(rad, arr, label = 'Analytical', color = 'grey')

# simulation bit
mod_path = 'resources/zeng-ical-template.mmt'

mod = myokit.load_model(mod_path)
mod.get('calcium_dynamics.buffering').set_rhs(0)
mod.get('calcium_dynamics.current').set_rhs(0)
mod.get('membrane.V').set_binding(None)

# Test for the biggest possible cell
mod.get('calcium_dynamics.R').set_rhs(R)

# remove ion-channel model
mod.get('L_type_Ca_current.d').demote
mod.get('L_type_Ca_current.d').set_rhs(0)
mod.get('L_type_Ca_current.f').demote
mod.get('L_type_Ca_current.f').set_rhs(1)

def compute_mass(N_div):

    log_b = ['calcium_dynamics.B' + str(i) for i in range(N_div)]
    B_arr = [0 for _ in range(N_div)]

    m = helpers.multi_compartment_model(mod, N_div, B_arr)

    log_v = []
    for i in range(N_div):
        log_v.append(m.get(f'calcium_dynamics.V{i}').pyfunc()())

    s = myokit.Simulation(m)
    s.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)
    d = s.run(t + 1, log = log_b, log_times = np.array([t]))

    # the result d is stored for each state-variable, time-wise
    ans_len = []

    for i in range(N_div):
        ans_len.append(d[f'calcium_dynamics.B{i}'][0])

    return ans_len


N_range = [400, 600, 800, 1000] 
x_arr = []
y_arr = []


colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

for j in range(len(N_range)):
    N_div = N_range[j]
    ans_len = compute_mass(N_div)
    x_plot = np.linspace(R*10000, rh*10000, N_div)

    ax.plot(x_plot, ans_len, ls = '--', label = f'N = {N_div}')
    dR = 1e4*(R - rh)/N_div

    error = abs((ans_len[0] - arr[0])/arr[0])
    error = np.log(error)
    
    dR = np.log(dR)
    
    ax2.scatter(dR, error, color = colors[j])
    x_arr.append(dR)
    y_arr.append(error)
    


def fun(x, m, c):
    """
    Linear fit according to inital data
    """
    return m*x + c


popt_n, _ = optimize.curve_fit(fun, x_arr, y_arr)
ax2.plot(x_arr, fun(np.array(x_arr), *popt_n), ls = '--', \
    label = 'Slope = 1.0', color = 'grey')
print(popt_n)

ax2.set_ylabel('ln(normalised error)')
ax2.set_xlabel('ln($L$)')

ax.legend()
ax2.legend()


plt.tight_layout()
plt.savefig('supporting/convergence.png')
plt.close()