import copy
import myokit
import numpy as np
from scipy.optimize import fsolve

import parameters

def multi_compartment_model(m, no_of_div, B_arr):

    """
    no_of_div: number of equally spaced divisions along the radius (shells) [radius - shell width]
    at r= r0: Cai = 0, B = 10
    B_arr: array of the values of B, where i = 0 corrresponds to near the membrane
    Here, i = 0 starts from nearest the membrane.
    """
    # ensure that a new object is created rather than reference to previous
    model = copy.deepcopy(m)  
    calcium_dynamics = model.get('calcium_dynamics')
    N = calcium_dynamics.get('N')
    N.set_rhs(no_of_div)

    # First include all the variables needed as literals
    for i in range(0, no_of_div+1):
        # shell 0 variables (and some shell1) are already defined
        try: 
            calcium_dynamics.add_variable('V' + str(i))
        except:
            pass
        try:
            calcium_dynamics.add_variable('A' + str(i))
        except:
            pass
        try:
            if i < no_of_div:
                calcium_dynamics.add_variable('Ca' + str(i))
                calcium_dynamics.add_variable('B' + str(i))
                calcium_dynamics.add_variable('CaB' + str(i))
        except:
            pass

        vol = calcium_dynamics.get('V%s' %(i))
        vol.set_rhs(f'2*3.14*((R-{i}*L)^3 - (R-{i+1}*L)^3)/3')
        area = calcium_dynamics.get('A%s' %(i))
        area.set_rhs(f'2*3.14*(R-({i}*L))^2') # 

    # create same state variables within this file as created in .mmt above
    Ca = ['Ca' + str(i) for i in range(no_of_div)]
    B = ['B' + str(i) for i in range(no_of_div)]
    CaB = ['CaB' + str(i) for i in range(no_of_div)]
    
    # Promote the variables needed as states (shell 0 defined in file already)
    for i in range(1, no_of_div):

        c = calcium_dynamics.get(Ca[i])
        b = calcium_dynamics.get(B[i])
        cab = calcium_dynamics.get(CaB[i])

        if i + 1 < no_of_div:
            c.set_rhs(f'(D_Ca *(Ca{i-1} - Ca{i}) * A{i}/(L*V{i}))\
                + buffering * (-k_on * Ca{i} * B{i} + k_off * CaB{i}) \
                    - (D_Ca *(Ca{i} - Ca{i+1}) * A{i+1}/(L*V{i}))')
            b.set_rhs(f'D_B *(B{i-1} - B{i})* A{i}/(L*V{i}) + buffering * (- k_on * Ca{i} * B{i} \
                + k_off * CaB{i}) - D_B *(B{i} - B{i+1})* A{i+1}/(L*V{i})')
            cab.set_rhs(f'D_B *(CaB{i-1} - CaB{i})* A{i}/(L*V{i}) - buffering * (- k_on * Ca{i} * B{i} \
                + k_off * CaB{i}) - D_B *(CaB{i} - CaB{i+1})* A{i+1}/(L*V{i})')

            c.promote(0)
            b.promote(B_arr[i])
            cab.promote(0)
        else:
            c.set_rhs(f'(D_Ca *(Ca{i-1} - Ca{i}) * A{i}/(L*V{i}))\
                + buffering * (-k_on * Ca{i} * B{i} + k_off * CaB{i}) \
                    - (D_Ca *(Ca{i} - 0) * A{i+1}/(L*V{i}))')
            b.set_rhs(f'D_B *(B{i-1} - B{i})* A{i}/(L*V{i}) + buffering * (- k_on * Ca{i} * B{i} \
                + k_off * CaB{i}) - D_B *(B{i} - B_max)* A{i+1}/(L*V{i})')
            cab.set_rhs(f'D_B *(CaB{i-1} - CaB{i})* A{i}/(L*V{i}) - buffering * (- k_on * Ca{i} * B{i} \
                + k_off * CaB{i}) - D_B *(CaB{i} - 0)* A{i+1}/(L*V{i})')
            c.promote(0)
            b.promote(B_arr[i])
            cab.promote(0)                             
        
    return model


def simulate_ical(n_sweep, N_div, R, t0, PCa, thold = 10000, peak = False):
    """
    Everything in the "hollow" hemispherical model
    """

    # array for store
    log_b = [f'calcium_dynamics.B{i}' for i in range(N_div)]

    mod = myokit.load_model('resources/zeng-ical-template.mmt')
    mod.get('calcium_dynamics.R').set_rhs(R) # cm
    mod.get('L_type_Ca_current.P_Ca').set_rhs(PCa)
    m = multi_compartment_model(mod, N_div, [0 for _ in range(N_div)])
    
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
    s.set_tolerance(abs_tol=1e-9, rel_tol=1e-09)
    
    d = s.run(t0 + 1, log = log_b, log_times = np.array([t0]))

    for b in log_b:
        m.get(b).set_state_value(d[b][0])


    # ion-channel model ss
    d_inf = m.get('L_type_Ca_current.d_infinity').pyfunc()(-90)
    m.get('L_type_Ca_current.d').set_state_value(d_inf)
    f_inf = m.get('L_type_Ca_current.f_infinity').pyfunc()(-90)
    m.get('L_type_Ca_current.f').set_state_value(f_inf)
    # if block_fca == True:
    #     m.get('L_type_Ca_current.open_prob').set_rhs('d*f')
    m.get('membrane.V').set_binding(None)

    sim = myokit.Simulation(m)
    sim.set_tolerance(abs_tol=1e-9, rel_tol=1e-09)

    sim.set_constant('membrane.V', -90)
    sim.run(thold, log=[]) 
    sim.set_constant('membrane.V', -0.001)
    t = sim.time()
    d = sim.run(120, log = ['calcium_dynamics.Ca0', 'calcium_dynamics.B0'],\
         log_times = np.arange(t, t+120, 1))
    
    if peak == True:
        peak_cal = [max(d['calcium_dynamics.Ca0'])]
        # peak_b = [min(d['calcium_dynamics.B0'])]

    for i in range(n_sweep -1):
        sim.set_constant('membrane.V', -90)
        sim.run(thold, log = [])
        sim.set_constant('membrane.V', -0.001)
        t = sim.time()
        b = sim.run(120, log = d, log_times = np.arange(t, t+120, 1))

        if peak == True:
            peak_cal.append(max(b['calcium_dynamics.Ca0'][-120:]))
            # peak_b.append(min(b['calcium_dynamics.B0']))

    if peak == True:
        return [x/peak_cal[0] - 1 for x in peak_cal], 0 
    
    else:
        return d
    
def simulate_ical_curr(n_sweep, N_div, R, t0, PCa, thold = 10000):
    """
    Everything in the "hollow" hemispherical model
    """

    # array for store
    log_b = [f'calcium_dynamics.B{i}' for i in range(N_div)]

    mod = myokit.load_model('resources/zeng-ical-template.mmt')
    mod.get('calcium_dynamics.R').set_rhs(R) # cm
    mod.get('L_type_Ca_current.P_Ca').set_rhs(PCa)
    m = multi_compartment_model(mod, N_div, [0 for _ in range(N_div)])
    
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
    s.set_tolerance(abs_tol=1e-9, rel_tol=1e-09)
    
    d = s.run(t0 + 1, log = log_b, log_times = np.array([t0]))

    for b in log_b:
        m.get(b).set_state_value(d[b][0])


    # ion-channel model ss
    d_inf = m.get('L_type_Ca_current.d_infinity').pyfunc()(-90)
    m.get('L_type_Ca_current.d').set_state_value(d_inf)
    f_inf = m.get('L_type_Ca_current.f_infinity').pyfunc()(-90)
    m.get('L_type_Ca_current.f').set_state_value(f_inf)
    # if block_fca == True:
    #     m.get('L_type_Ca_current.open_prob').set_rhs('d*f')
    m.get('membrane.V').set_binding(None)

    sim = myokit.Simulation(m)
    sim.set_tolerance(abs_tol=1e-9, rel_tol=1e-09)

    sim.set_constant('membrane.V', -90)
    sim.run(thold, log=[]) 
    sim.set_constant('membrane.V', -0.001)
    t = sim.time()
    d = sim.run(120, log = ['L_type_Ca_current.i_CaL'],\
         log_times = np.arange(t, t+120, 1))
    
    peak_ical = [min(d['L_type_Ca_current.i_CaL'])]

    for i in range(n_sweep -1):
        sim.set_constant('membrane.V', -90)
        sim.run(thold, log = [])
        sim.set_constant('membrane.V', -0.001)
        t = sim.time()
        b = sim.run(120, log = d, log_times = np.arange(t, t+120, 1))

        peak_ical.append(min(b['L_type_Ca_current.i_CaL'][-120:]))
    
    return [1 - x/peak_ical[0] for x in peak_ical], 0


def draw_t():
    # draw all three time variables
    thold_range = [10000, 20000, 40000] # milli seconds
    thold = np.random.choice(thold_range)
    t0 = int(np.random.uniform(5000, 300000, size = 1)[0])
    tdiff = int(np.random.uniform(144, 422601, size = 1)[0])
    
    return thold, t0, tdiff

def peak_ca(thold, t0, tdiff, n_sweep, pca, dR, block_fca = True):

    # compute r from t_diff
    f = lambda x: 0.0163 * x ** 3 - 0.019 * x ** 2 + 0.008 * x - 0.016 - tdiff/1000
    r = fsolve(f, [30])[0] * 1e-4 # convert to cm

    N_div = int((r-parameters.Rh)/dR)
    print('radius (um)', r*1e4, 'Ndiv', N_div)

    peak_ca, _ = simulate_ical(n_sweep, N_div, r, t0, pca, thold = thold, peak = True)

    return peak_ca

def peak_ical(thold, t0, tdiff, n_sweep, pca, dR):

    # compute r from t_diff
    f = lambda x: 0.0163 * x ** 3 - 0.019 * x ** 2 + 0.008 * x - 0.016 - tdiff/1000
    r = fsolve(f, [30])[0] * 1e-4 # convert to cm

    N_div = int((r-parameters.Rh)/dR)
    print('radius (um)', r*1e4, 'Ndiv', N_div)

    peak_ical, _ = simulate_ical_curr(n_sweep, N_div, r, t0, pca, thold = thold)

    return peak_ical

def simple_beeswarm(y, nbins=None):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        nbins = len(y) // 6

    # Get upper bounds of bins
    x = np.zeros(len(y))
    ylo = np.min(y)
    yhi = np.max(y)
    dy = (yhi - ylo) / nbins
    ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

    # Divide indices into bins
    i = np.arange(len(y))
    ibs = [0] * nbins
    ybs = [0] * nbins
    nmax = 0
    for j, ybin in enumerate(ybins):
        f = y <= ybin
        ibs[j], ybs[j] = i[f], y[f]
        nmax = max(nmax, len(ibs[j]))
        f = ~f
        i, y = i[f], y[f]
    ibs[-1], ybs[-1] = i, y
    nmax = max(nmax, len(ibs[-1]))

    # Assign x indices
    dx = 1 / (nmax // 2)
    for i, y in zip(ibs, ybs):
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(y)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x


def lin(t, min = False, max = False):
    """
    c varies from 0.1 t0 0
    """
    c_min = 0
    c_max = 0.008
    
    if min:
        return -1 + c_min*t
    elif max:
        return -1 + c_max*t
    
    c = np.random.uniform(low = c_min, high = c_max, size = 1)
    return -1 + c[0]*t

def logar(t, min = False, max = False):
    """
    c between -1 and 0
    """
    c_min = -0.7
    c_max = -0.1
    
    if min:
        return -1 + t/(t+1) + c_min
    elif max:
        return -1 + t/(t+1) + c_max

    c = np.random.uniform(low = c_min, high = c_max, size = 1)
    return -1 + t/(t+1) + c

def hyper(t, min = False, max = False):
    """
    alpha varies from 0 to 1
    """
    alpha = np.random.uniform(low = 0.01, high = 0.3, size = 1)
    return -1*(1-alpha) - alpha*t/(t+1)

def calcium(x, min = False, max = False):
    ic50_min = 0.0000005
    ic50_max = 0.5

    if min:
        return -ic50_min*(1 + 1/x) 
    elif max:
        return -ic50_max*(1 + 1/x)

    ic50 = np.random.uniform(low = ic50_min, high = ic50_max, size = 1)
    return -ic50[0]*(1 + 1/x) 

