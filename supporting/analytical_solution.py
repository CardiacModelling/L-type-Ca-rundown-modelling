import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import optimize

import parameters

B_max = parameters.Bmax # mM
r_h = parameters.Rh #cm
DB = parameters.DB #cm2/ms


def roots_den(h):
    def f(x):
        lam = x
        return np.sin(lam) * h - np.cos(lam)*lam # roots are given by tan(x) = x/h

    x0 = np.array([3.1*i for i in range(1, 50)])
    sol = optimize.root(f, x0)
    p_range = sol.x
    return p_range

def term1(x, h):
    return (h*x - 1)/(h-1)

def term2(lam, x, tau, h):
    def alpha():
        y = (-1 + (1/h) + (h/pow(lam,2)))*np.sin(lam) - \
            (h/lam)*np.cos(lam)
        return y/(h-1)

    def beta():
        y = pow(lam, 2) + pow(h, 2) + \
            (pow(lam, 2) - pow(h,2))*np.sin(2*lam)/(2*lam) - \
                2*h*pow(np.sin(lam), 2)
        return y/(2*pow(h, 2))

    def func_x():
        return np.sin(lam*x) - (lam/h)*np.cos(lam*x)

    def func_t():
        return np.exp(-pow(lam,2)*tau)

    term = alpha() * func_x() * func_t()/(beta())

    return term


def check_roots(h):
    # The below checks that the roots found are sensible
    roots = roots_den(h)
    roots = np.sort(roots)
    index = []
    for i in range(len(roots)):
        if roots[i] < 0:
            index.append(i)
        elif i == 0:
            continue
        elif roots[i] - roots[i-1] < 3: # root will be in intervals of pi
            index.append(i)

    roots = np.delete(roots, index)

    for i in range(1, len(roots)):
        if roots[i] < roots[i-1]:
            print(roots[i])
            raise ValueError ('Calculation of roots is not in increasing order')
        elif roots[i-1] <= 0:
            print(roots[i-1])
            raise ValueError ('Found a non-positive root')
    return roots


def convergence(x, tau, roots, h):
    answer = 0
    term_store = 0
    
    for root in roots:
        term = term2(root, x, tau, h)
        answer +=term

        if root == roots[0]:
            term_store = term
            continue
        elif term_store == 0:
            break
        elif abs(term)/ abs(answer) < 0.0001: 
            break
        elif root == roots[-1]:
            if abs(term)/ abs(answer) < 0.005:
                return answer
            else:
                print('answer: ', answer)
                raise ValueError(f'Series did not converge for x={x}, tau = {tau}')
        else:
            term_store = term 
            continue
    return answer   

def cal_tdiff(r):
    h = 1 - (r_h/r)

    roots = roots_den(h)
    roots.sort() # roots in ascending order

    arr_ind = []
    for i in range(len(roots)):
        if roots[i] < 0:
            arr_ind.append(i)
        elif i == len(roots) - 1:
            continue
        elif roots[i+1] - roots[i] < 2:
            arr_ind.append(i)
    roots = np.delete(roots, arr_ind) # remove negative roots

    tau_range = np.linspace(0.5, 5000, 10000)
    for tau in tau_range:
        # only need solution at the boundar, which is x = 0
        answer = convergence(0, tau, roots, h) 
        # check that term2 is much smaller than term1, i.e, less than 1%
        if abs(answer) < 0.01 * abs(1/(1 - h)):
            tdiff = (tau * (r - r_h)**2)/DB
            return tdiff
        elif tau == tau_range[-1]:
            raise ValueError('tau out of range')

def fun(x, a, b, c, d):
    """
    Curve fit
    """
    return a*x**3 + b*x**2 + c*x + d

def plot():
    matplotlib.rc('font', family='arial', size = 8)
    fig = plt.figure()

    ax = fig.add_subplot(121)
    R = parameters.R_max
    ax.set_ylim(0, B_max)
    ax.set_xlim(0, R*10000)
    ax.set_ylabel('Buffer (mM)')
    ax.set_xlabel('Distance from the centre ($\mu$m)')

    t_mult = DB/pow(R-r_h, 2)
    t_range = np.linspace(1000, 300000, 10) # time array in milli second for which we want diffusion
    tau_range = [t*t_mult for t in t_range]
    x_range = np.linspace(0, 1, 40)
    h = 1 - (r_h/R)
    roots = check_roots(h)

    colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, len(t_range))]

    for i in range(len(tau_range)):
        arr = [] # buffer concentration with space at tau
        rad = [] # corresponding radius

        for x in x_range:
            term_a = term1(x, h)
            term_b = convergence(x, tau_range[i], roots, h)
            u = term_a - term_b
            r = R + x*(r_h - R)
            arr.append(B_max*u*r_h/r)
            rad.append(r*10000)

        
        ax.plot(rad, arr, color = colors[i])
        ax.text(30, arr[0] + 0.1, f't = {round(t_range[i]/1000)}s', fontsize = 6)

    ax2 = fig.add_subplot(122)
    radius_range = np.linspace(parameters.R_min, parameters.R_max, 100)
    tdiff_arr = []
    for r in radius_range:
        tdiff = cal_tdiff(r)
        tdiff_arr.append(tdiff/1000)
    ax2.plot(radius_range*10000, tdiff_arr, label = 'Analytical solution')

    # get a fit to the diffusion
    popt, _ = optimize.curve_fit(fun, radius_range*10000, np.array(tdiff_arr))

    ax2.plot(radius_range*10000, fun(radius_range*10000, round(popt[0], 4), round(popt[1], 4), \
            round(popt[2], 4), round(popt[3], 4)), ls = '--',\
            label ='y = %5.4f$x^{3}$ + %5.3f$x^{2}$ + %5.3fx + %5.3f' %tuple(popt))

    ax2.set_xlabel('Radius of the cell ($\mu$m)')
    ax2.set_ylabel('$\mathrm{t}_{diff}$ (s)')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.close()

def plot_largest_r():
    matplotlib.rc('font', family='arial', size = 8)
    fig = plt.figure()

    ax = fig.add_subplot(111)
    R = parameters.R_max
    ax.set_ylim(0, B_max)
    ax.set_xlim(0, R*10000)
    ax.set_ylabel('Buffer (mM)')
    ax.set_xlabel('Distance from the centre ($\mu$m)')

    t_mult = DB/pow(R-r_h, 2)
    t_range = np.linspace(1000, 300000, 10) # time array in milli second for which we want diffusion
    tau_range = [t*t_mult for t in t_range]
    x_range = np.linspace(0, 1, 40)
    h = 1 - (r_h/R)
    roots = check_roots(h)

    colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, len(t_range))]

    for i in range(len(tau_range)):
        arr = [] # buffer concentration with space at tau
        rad = [] # corresponding radius

        for x in x_range:
            term_a = term1(x, h)
            term_b = convergence(x, tau_range[i], roots, h)
            u = term_a - term_b
            r = R + x*(r_h - R)
            arr.append(B_max*u*r_h/r)
            rad.append(r*10000)

        
        ax.plot(rad, arr, color = colors[i])
        ax.text(30, arr[0] + 0.1, f't = {round(t_range[i]/1000)}s', fontsize = 6)

    plt.tight_layout()
    plt.show()
    plt.close()

def plot_tau_max():
    matplotlib.rc('font', family='arial', size = 8)
    fig = plt.figure()

    ax2 = fig.add_subplot(111)
    radius_range = np.linspace(parameters.R_min, parameters.R_max, 100)
    tdiff_arr = []
    for r in radius_range:
        tdiff = cal_tdiff(r)
        tdiff_arr.append(tdiff/1000)
    ax2.plot(radius_range*10000, tdiff_arr, label = 'Analytical solution')

    # get a fit to the diffusion
    popt, _ = optimize.curve_fit(fun, radius_range*10000, np.array(tdiff_arr))

    ax2.plot(radius_range*10000, fun(radius_range*10000, round(popt[0], 4), round(popt[1], 4), \
            round(popt[2], 4), round(popt[3], 4)), ls = '--',\
            label ='y = %5.4f$x^{3}$ + %5.3f$x^{2}$ + %5.3fx + %5.3f' %tuple(popt))
    
    print('R min (um)', radius_range[0]*10000, 'tdiff (s)', tdiff_arr[0])
    print('R max (um)', radius_range[-1]*10000, 'tdiff (s)', tdiff_arr[-1])

    ax2.set_xlabel('Radius of the cell ($\mu$m)')
    ax2.set_ylabel('$\mathrm{t}_{diff}$ (s)')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.close()

